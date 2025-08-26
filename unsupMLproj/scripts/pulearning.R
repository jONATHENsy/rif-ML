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
  library(ggplot2)   # for Recall@K figure / 画 Recall@K 曲线
})

## -------------------------------
## 0) I/O & Global Config 全局参数
## -------------------------------

# Output directory for csv & figures; will auto-create if not exists
# 输出目录，若不存在则自动创建
if (!exists("pre_dir")) pre_dir <- "predict"
if (!dir.exists(pre_dir)) dir.create(pre_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(42)

# Core hyper-params for PU-RF / 关键超参数
NTREE <- 500           # trees per RF
B_BAGS <- 25           # U-bagging rounds
NEG_RATIO <- 1.0       # unlabeled subsample size relative to #positives
MIN_POS <- 3           # skip mutations with < MIN_POS positives
USE_UMAP <- FALSE      # 建议关闭防止泄漏（如启用需按目标突变重算）
DROP_SAME_SITE <- TRUE # 训练目标突变时屏蔽同位点其他等位
ADD_EFFORT <- TRUE     # effort = rowSums(X_dense)

# Candidate selection config 候选筛选配置
PTRUE_THRESHOLD <- 0.7
TOPK_PER_SPECIES <- 10

## -------------------------------------------
## 1) Basic Checks & Coercions 基本检查与整理
## -------------------------------------------

stopifnot(exists("X_dense"))
X_dense <- as.data.frame(X_dense)

if (is.null(rownames(X_dense))) stop("X_dense must have rownames as Species / 需要行为物种名")

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
## ---------------------------------------------------------

parse_mut_cols <- function(mut_names) {
  # e.g., rpoB_H526R -> gene=rpoB, ref=H, pos=526, alt=R
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

drop_same_site_cols <- function(X_df, col_meta, target_mut) {
  sk <- col_meta$site_key[col_meta$Mutation == target_mut]
  if (length(sk) == 0 || is.na(sk)) return(X_df)
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
  if (!target_mut %in% colnames(X_dense)) return(NULL)
  y <- as.integer(X_dense[[target_mut]])  # 1=reported (positive), 0=unlabeled
  X_df <- X_dense
  X_df[[target_mut]] <- NULL  # avoid direct leakage
  
  if (drop_same_site) {
    X_df <- drop_same_site_cols(X_df, col_meta, target_mut)
  }
  if (add_effort) {
    X_df$effort <- rowSums(as.matrix(X_dense))
  }
  if (use_umap) {
    stop("USE_UMAP=TRUE 需先对去掉 target_mut 的矩阵重算 UMAP 后再合并，以避免信息泄漏。")
  }
  
  pos_idx <- which(y == 1)
  unl_idx <- which(y == 0)
  
  if (length(pos_idx) < MIN_POS) {
    return(tibble(
      Species = rownames(X_df),
      p_obs = NA_real_,
      p_true = NA_real_,
      Observed = y,
      Mutation = target_mut
    ))
  }
  
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
  
  # Elkan–Noto correction / 概率校正
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

# B) Top-k per species 每个物种选 Top-k 候选
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

# prob_wide <- res %>%
#   select(Species, Mutation, p_true) %>%
#   pivot_wider(names_from = Mutation, values_from = p_true)
# prob_wide_file <- file.path(pre_dir, "PU_RF_prob_wide.csv")
# write_csv(prob_wide, prob_wide_file)
# message("✔ Saved: ", prob_wide_file)

## ===========================================================
## 7) Simplified PU evaluation: Mask-then-Recover (Recall@K)
##    简化评估：遮蔽再找回，仅报告 Recall@K
## ===========================================================

# Eval hyper-params 评估参数
MASK_FRAC <- 0.20             # fraction of known positives to mask / 遮蔽比例（例如 20%）
N_REPEATS <- 5                # repeats 重复次数
K_SET     <- c(1, 5, 10, 20)  # K for Recall@K

# Helper: run PU-RF for all mutations on a given matrix
run_pu_rf_all <- function(X_in) {
  muts <- colnames(X_in)
  lst <- lapply(muts, function(m) {
    tryCatch(
      run_pu_rf_for_mut(m, X_in, col_meta),
      error = function(e) { warning(sprintf("PU-RF failed for %s: %s", m, e$message)); NULL }
    )
  })
  bind_rows(lst)
}

eval_mask_then_recover <- function(X_in, mask_frac = MASK_FRAC, repeats = N_REPEATS, K_set = K_SET) {
  X_in <- as.data.frame(X_in)
  rs <- list()
  
  for (r in seq_len(repeats)) {
    set.seed(2025 + r)
    
    # 1) Sample positives to mask / 采样需要遮蔽的正样本(物种-突变对)
    P <- which(as.matrix(X_in) == 1, arr.ind = TRUE)
    if (nrow(P) < 2) next
    nmask <- max(1, floor(nrow(P) * mask_frac))
    M <- P[sample(nrow(P), nmask), , drop = FALSE]
    
    # 2) Mask them to 0 (unlabeled) / 遮蔽为0
    X_mask <- X_in
    for (i in seq_len(nrow(M))) {
      X_mask[M[i, 1], M[i, 2]] <- 0
    }
    
    # 3) Train on masked data & predict / 训练并预测
    res_mask <- run_pu_rf_all(X_mask)
    
    # 4) Build masked table / 遮蔽真阳性表
    masked_tbl <- tibble(
      Species = rownames(X_in)[M[, 1]],
      Mutation = colnames(X_in)[M[, 2]]
    )
    
    # 5) Rank among unlabeled per species / 未标注池内按 p_true 排序并排名
    ranks_tbl <- res_mask %>%
      filter(Observed == 0, !is.na(p_true)) %>%
      group_by(Species) %>%
      arrange(desc(p_true), .by_group = TRUE) %>%
      mutate(rank = row_number()) %>%
      ungroup()
    
    hits <- ranks_tbl %>%
      inner_join(masked_tbl, by = c("Species", "Mutation"))
    
    # 6) Compute Recall@K / 计算 Recall@K
    rec_df <- tibble(
      fold = r,
      K = K_set,
      Recall_at_K = sapply(K_set, function(K) mean(hits$rank <= K))
    )
    rs[[r]] <- rec_df
    
    # 保存每折遮蔽与命中（可选，需要的话取消注释）
    # write_csv(masked_tbl, file.path(pre_dir, sprintf("PU_eval_maskedPairs_fold%02d.csv", r)))
    # write_csv(hits,       file.path(pre_dir, sprintf("PU_eval_hits_fold%02d.csv", r)))
  }
  
  out <- bind_rows(rs)
  out_file <- file.path(pre_dir, "PU_RF_eval_recall_at_k.csv")
  write_csv(out, out_file)
  message("✔ Saved: ", out_file)
  
  # Summary with mean & SE / 汇总（均值与标准误）
  summ <- out %>%
    group_by(K) %>%
    summarise(
      mean_recall = mean(Recall_at_K, na.rm = TRUE),
      se_recall   = sd(Recall_at_K, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  summ_file <- file.path(pre_dir, "PU_RF_eval_recall_at_k_summary.csv")
  write_csv(summ, summ_file)
  message("✔ Saved: ", summ_file)
  
  # Figure: Recall@K curve / 画 Recall@K 曲线
  p <- ggplot(summ, aes(x = K, y = mean_recall)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = mean_recall - se_recall, ymax = mean_recall + se_recall), width = 0.2) +
    labs(title = "Mask-then-Recover: Recall@K",
         x = "K (top-K predictions per species)",
         y = "Mean Recall@K ± SE")
  fig_file <- file.path(pre_dir, "PU_RF_recall_at_k.png")
  ggsave(fig_file, p, width = 6.5, height = 4.2, dpi = 300)
  message("✔ Saved figure: ", fig_file)
  
  invisible(list(detail = out, summary = summ))
}
# 1) 每个物种有多少“有效候选池”？
per_species <- res %>%
  group_by(Species) %>%
  summarise(
    n_unlabeled = sum(Observed == 0, na.rm = TRUE),
    n_valid_pred = sum(Observed == 0 & !is.na(p_true)),
    n_above_thr  = sum(Observed == 0 & !is.na(p_true) & p_true >= PTRUE_THRESHOLD)
  ) %>% arrange(n_valid_pred)
print(per_species)

# 2) 哪些突变列被训练为“无效列”（< MIN_POS 或 c_est 失败）？
mut_summary <- res %>%
  group_by(Mutation) %>%
  summarise(
    pos_count = sum(Observed == 1, na.rm = TRUE),
    any_valid = any(!is.na(p_true))
  ) %>% arrange(pos_count)
print(head(mut_summary, 20))           # 看低阳性列
mean(!mut_summary$any_valid)           # 无效列占比

# 3) 哪些物种完全没出现在TopK候选？
missing_in_topk <- setdiff(rownames(X_dense), unique(candidates_topk$Species))
missing_in_topk

message("▶ Running simplified PU evaluation (Mask-then-Recover, Recall@K only) ...")
eval_out <- eval_mask_then_recover(X_dense)
message("🏁 Evaluation done. 查看：PU_RF_eval_recall_at_k_summary.csv / PU_RF_recall_at_k.png")

message("🏁 Done. Use 'PU_RF_candidates_*' and evaluation outputs for figures & text.")
message("完成：候选结果与 Recall@K 评估已生成。")
