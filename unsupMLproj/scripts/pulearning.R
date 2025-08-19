## ===========================================================
## PU-learning with Random Forest for mutation prediction
## 正-未标注（PU）学习 + 随机森林，用于“未报道突变”预测
## Requirements: X_dense (matrix/data.frame, rows=Species, cols=Mutations, values in {0,1})
## ===========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(purrr)
  library(randomForest)
  library(readr)
})

## -------------------------------
## 0) I/O & Global Config 全局参数
## -------------------------------

# Output directory for csv files; will auto-create if not exists
# 输出目录，若不存在则自动创建
if (!exists("pre_dir")) pre_dir <- "predict"
if (!dir.exists(pre_dir)) dir.create(pre_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(42)

# Core hyper-params for PU-RF / 关键超参数
NTREE <- 500           # trees per RF
B_BAGS <- 25           # U-bagging rounds
NEG_RATIO <- 1.0       # unlabeled subsample size relative to #positives (每轮从U中采样与正样本数相当的伪负)
MIN_POS <- 3           # skip mutations with < MIN_POS positives (太少的正样本跳过)
USE_UMAP <- FALSE      # 强烈建议先禁用 UMAP 避免泄漏（若启用需保证对每个目标突变重算UMAP）
DROP_SAME_SITE <- TRUE # 训练某个突变时屏蔽同位点其他等位（避免“捷径特征”/信息泄漏）
ADD_EFFORT <- TRUE     # 添加每个物种的“研究/采样力度”特征：effort = rowSums(X_dense)

# Candidate selection config 候选筛选配置
PTRUE_THRESHOLD <- 0.7  # 概率阈值（未报道且 p_true >= 阈值 → 候选）
TOPK_PER_SPECIES <- 10  # 每个物种导出Top-k候选（未报道依据 p_true 排序）

## -------------------------------------------
## 1) Basic Checks & Coercions 基本检查与整理
## -------------------------------------------

stopifnot(exists("X_dense"))
X_dense <- as.data.frame(X_dense)

# Ensure rownames are present for species / 确保行名是物种名
if (is.null(rownames(X_dense))) stop("X_dense must have rownames as Species / 需要行为物种名")

# Coerce to 0/1 numeric; replace NA with 0 / 转为0/1数值，NA补0
X_dense[is.na(X_dense)] <- 0
X_dense[] <- lapply(X_dense, function(col) {
  v <- suppressWarnings(as.numeric(col))
  v[is.na(v)] <- 0
  pmin(pmax(v, 0), 1)
})
X_dense <- as.data.frame(X_dense)

all_muts <- colnames(X_dense)

## ---------------------------------------------------------
## 2) Name parser for "same-site" masking 位点解析与屏蔽工具
##    Expect names like: rpoB_H526R / gene_REFposALT
## ---------------------------------------------------------

# Parse mutation names into components; if unknown pattern, site_key is NA
# 解析突变命名，若不匹配常见格式则site_key为NA
parse_mut_cols <- function(mut_names) {
  # regex: gene_refAA? pos alt?
  # 例: rpoB_H526R -> gene=rpoB, ref=H, pos=526, alt=R
  mat <- stringr::str_match(mut_names, "^([^_]+)_([A-Z\\*])?(\\d+)([A-Z\\*]+)?$")
  tibble(
    Mutation = mut_names,
    gene = mat[,2],
    refAA = mat[,3],
    pos   = suppressWarnings(as.integer(mat[,4])),
    altAA = mat[,5],
    site_key = ifelse(is.na(mat[,2]) | is.na(mat[,4]),
                      NA_character_,
                      paste0(mat[,2], "_", mat[,4])) # gene_pos
  )
}

col_meta <- parse_mut_cols(all_muts)

# Helper: drop columns that share the same site_key with target_mut
# 工具：训练target时，屏蔽同位点（gene+pos）其它等位突变列
drop_same_site_cols <- function(X_df, col_meta, target_mut) {
  sk <- col_meta$site_key[col_meta$Mutation == target_mut]
  if (length(sk) == 0 || is.na(sk)) return(X_df)  # unknown pattern → do nothing
  drop_set <- col_meta$Mutation[col_meta$site_key == sk & col_meta$Mutation != target_mut]
  keep_cols <- setdiff(colnames(X_df), drop_set)
  X_df[, keep_cols, drop = FALSE]
}

## ---------------------------------------------------------------
## 3) PU-RF for one mutation 单一突变的PU随机森林（含U-bagging）
## ---------------------------------------------------------------

run_pu_rf_for_mut <- function(target_mut,
                              X_dense,
                              col_meta,
                              ntree = NTREE,
                              B = B_BAGS,
                              neg_ratio = NEG_RATIO,
                              use_umap = USE_UMAP,
                              drop_same_site = DROP_SAME_SITE,
                              add_effort = ADD_EFFORT) {
  # Labels and features / 标签与特征
  if (!target_mut %in% colnames(X_dense)) return(NULL)
  y <- as.integer(X_dense[[target_mut]])  # 1=reported (positive), 0=unlabeled
  X_df <- X_dense
  X_df[[target_mut]] <- NULL  # avoid direct leakage
  
  # Optional: drop same-site alleles (prevent "shortcut" features)
  # 可选：屏蔽同位点其它等位，防止“邻位泄漏/捷径”
  if (drop_same_site) {
    X_df <- drop_same_site_cols(X_df, col_meta, target_mut)
  }
  
  # Optional: add sampling/effort feature
  # 可选：加入“研究/采样力度”特征
  if (add_effort) {
    X_df$effort <- rowSums(as.matrix(X_dense))
  }
  
  # Optional: merge UMAP (disabled by default to avoid leakage)
  # 可选：合并UMAP（默认关闭；如需开启，请确保对“去掉target列后的矩阵”重算UMAP）
  if (use_umap) {
    stop("USE_UMAP=TRUE 需先对去掉 target_mut 的矩阵重算 UMAP 后再合并，以避免信息泄漏。")
  }
  
  pos_idx <- which(y == 1)
  unl_idx <- which(y == 0)
  
  # Too few positives → skip
  # 正样本太少则跳过
  if (length(pos_idx) < MIN_POS) {
    return(tibble(
      Species = rownames(X_df),
      p_obs = NA_real_,
      p_true = NA_real_,
      Observed = y,
      Mutation = target_mut
    ))
  }
  
  # U-bagging
  # 多次从未标注中采伪负，训练多个RF取平均概率
  p_obs_mat <- matrix(NA_real_, nrow = nrow(X_df), ncol = B)
  for (b in seq_len(B)) {
    set.seed(100 + b)
    k <- ceiling(length(pos_idx) * neg_ratio)
    neg_sample <- if (length(unl_idx) >= k) sample(unl_idx, k) else unl_idx
    idx <- c(pos_idx, neg_sample)
    
    rf <- randomForest(
      x = X_df[idx, , drop = FALSE],
      y = as.factor(y[idx]),
      ntree = ntree
    )
    p_obs_mat[, b] <- predict(rf, X_df, type = "prob")[, 2]  # p(s=1|x)
  }
  p_obs <- rowMeans(p_obs_mat, na.rm = TRUE)
  
  # Elkan–Noto correction: p_true = p_obs / c
  # 概率校正：把观测概率校正成“真实存在”概率
  c_est <- mean(p_obs[y == 1])
  if (is.na(c_est) || c_est <= 0) {
    p_true <- rep(NA_real_, length(p_obs))
  } else {
    p_true <- p_obs / c_est
    p_true[p_true > 1] <- 1
  }
  
  tibble(
    Species = rownames(X_df),
    p_obs = p_obs,
    p_true = p_true,
    Observed = y,
    Mutation = target_mut
  )
}

## ---------------------------------------------------------
## 4) Run for all mutations 批量对所有突变运行PU-RF
## ---------------------------------------------------------

message("▶ Running PU-RF for all mutations ... / 开始对全部突变进行 PU-RF 预测 …")
res_list <- lapply(all_muts, function(m) {
  tryCatch(
    run_pu_rf_for_mut(m, X_dense, col_meta),
    error = function(e) {
      warning(sprintf("PU-RF failed for %s: %s", m, e$message))
      NULL
    }
  )
})
res <- bind_rows(res_list)

# Save long-format probabilities
# 保存长表（每个 物种-突变 的 p_true 概率）
prob_long_file <- file.path(pre_dir, "PU_RF_prob_long.csv")
write_csv(res, prob_long_file)
message("✔ Saved: ", prob_long_file)

## ---------------------------------------------------------
## 5) Candidate selection 候选突变筛选
## ---------------------------------------------------------

# A) Threshold-based 未报道且 p_true ≥ 阈值
candidates_threshold <- res %>%
  filter(Observed == 0, !is.na(p_true)) %>%
  filter(p_true >= PTRUE_THRESHOLD) %>%
  arrange(desc(p_true))

cand_thr_file <- file.path(pre_dir, "PU_RF_candidates_threshold.csv")
write_csv(candidates_threshold, cand_thr_file)
message("✔ Saved: ", cand_thr_file)

# B) Top-k per species 每个物种选 Top-k 候选（未报道，按 p_true 排序）
candidates_topk <- res %>%
  filter(Observed == 0, !is.na(p_true)) %>%
  group_by(Species) %>%
  arrange(desc(p_true), .by_group = TRUE) %>%
  slice_head(n = TOPK_PER_SPECIES) %>%
  ungroup()

cand_topk_file <- file.path(pre_dir, "PU_RF_candidates_topk.csv")
write_csv(candidates_topk, cand_topk_file)
message("✔ Saved: ", cand_topk_file)

## ---------------------------------------------------------
## 6) (Optional) Wide matrix 生成宽表（可选）
## ---------------------------------------------------------

# 若需要 (物种 × 突变) 的 p_true 宽表，可取消以下注释
# If you need a wide matrix of p_true:
# prob_wide <- res %>%
#   select(Species, Mutation, p_true) %>%
#   pivot_wider(names_from = Mutation, values_from = p_true)
# prob_wide_file <- file.path(pre_dir, "PU_RF_prob_wide.csv")
# write_csv(prob_wide, prob_wide_file)
# message("✔ Saved: ", prob_wide_file)

message("🏁 Done. Use 'PU_RF_candidates_*' to review high-confidence unreported mutations.")
message("完成：请查看候选文件，核对高置信未报道突变。")
