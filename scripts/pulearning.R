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
## ===== Force-load master & create working copy =====
# 建议设置 base_dir；没设的话默认当前目录
if (!exists("base_dir")) base_dir <- getwd()

X_DENSE_MASTER_PATH <- file.path(base_dir, "output", "X_dense_master.rds")

# 1) 读取并锁定“母版”（只读；不在脚本里改它）
X_MASTER <- readRDS(X_DENSE_MASTER_PATH)
lockBinding("X_MASTER", .GlobalEnv)     # 以后谁都改不了 X_MASTER

# 2) 生成可修改的“工作副本”，供 PU 学习用
prep_X_for_pu <- function(X) {
  # 稀疏 -> 普通矩阵/数据框
  if (inherits(X, "Matrix")) X <- as.matrix(X)
  X <- as.data.frame(X, check.names = FALSE)
  
  # 行名必须是物种名
  if (is.null(rownames(X))) stop("X_dense must have rownames as Species")
  
  # 按列数值化（更安全地处理字符/因子），并二值化到 {0,1}
  X[] <- lapply(X, function(col) {
    if (is.factor(col)) col <- as.character(col)
    v <- suppressWarnings(as.numeric(col))
    v[is.na(v)] <- 0
    # 如果本来就是 0/1，会保持不变；否则做“>0 视为 1”的硬二值
    as.integer(v > 0)
  })
  X
}

# 工作副本：供脚本后续所有代码使用
X_dense <- prep_X_for_pu(X_MASTER)
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
K_SET     <- c(1,5,10,15,20,25,30,40,50)  # K for Recall@K

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


suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(readr); library(forcats); library(purrr)
})

# 目录与文件
if (!exists("pre_dir")) pre_dir <- "predict"
fig_dir <- file.path(pre_dir, "figs01")


# 读入数据
prob_long   <- read_csv(file.path(pre_dir, "PU_RF_prob_long.csv"), show_col_types = FALSE)
cand_thr    <- read_csv(file.path(pre_dir, "PU_RF_candidates_threshold.csv"), show_col_types = FALSE)
cand_topk   <- read_csv(file.path(pre_dir, "PU_RF_candidates_topk.csv"), show_col_types = FALSE)
eval_detail <- read_csv(file.path(pre_dir, "PU_RF_eval_recall_at_k.csv"), show_col_types = FALSE)
eval_summ   <- read_csv(file.path(pre_dir, "PU_RF_eval_recall_at_k_summary.csv"), show_col_types = FALSE)

# =============== 图1：每物种候选数（阈值法） ===============
p1_df <- cand_thr %>%
  count(Species, name = "n_candidates") %>%
  arrange(desc(n_candidates))

p1 <- ggplot(p1_df, aes(x = reorder(Species, n_candidates), y = n_candidates)) +
  geom_col() + coord_flip() +
  labs(title = "Candidates per species (thresholded)",
       x = "Species", y = "Number of candidates (p_true ≥ τ)") +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir, "fig1_candidates_per_species_threshold.png"), p1,
       width = 7.5, height = 6, dpi = 300)

# =============== 图2：重点物种 Top-K 候选（分面柱状图） ===============
# 请选择你要展示的重点物种（示例）：
focus_species <- c("Mycobacterium tuberculosis",
                   "Bacillus anthracis",
                   "Pseudomonas aeruginosa")
p2_df <- cand_topk %>%
  filter(Species %in% focus_species) %>%
  group_by(Species) %>%
  arrange(desc(p_true), .by_group = TRUE) %>%
  mutate(Mutation = fct_reorder(Mutation, p_true))
p2 <- ggplot(p2_df, aes(x = Mutation, y = p_true)) +
  geom_col() + coord_flip() +
  facet_wrap(~ Species, scales = "free_y") +
  labs(title = "Top-K predicted mutations per species",
       x = "Mutation", y = "p_true") +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir, "fig2_topk_by_species_bar.png"), p2,
       width = 9, height = 6, dpi = 300)

# =============== 图3：全局热图（每物种取 Top-N 未报道） ===============
TOP_N <- 15
p3_df <- prob_long %>%
  filter(Observed == 0, !is.na(p_true)) %>%
  group_by(Species) %>%
  slice_max(p_true, n = TOP_N, with_ties = FALSE) %>%
  ungroup()

# 为了更易读，按每物种的均值排序
species_order <- p3_df %>% group_by(Species) %>% summarise(m = mean(p_true)) %>%
  arrange(desc(m)) %>% pull(Species)
mut_order <- p3_df %>% group_by(Mutation) %>% summarise(m = mean(p_true)) %>%
  arrange(desc(m)) %>% pull(Mutation)

p3 <- ggplot(p3_df, aes(x = factor(Mutation, mut_order),
                        y = factor(Species, species_order),
                        fill = p_true)) +
  geom_tile() +
  scale_fill_gradient(name = "p_true", limits = c(0,1)) +
  labs(title = paste0("Top-", TOP_N, " predicted (unlabeled) per species"),
       x = "Mutation", y = "Species") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(file.path(fig_dir, "fig3_heatmap_topN_per_species.png"), p3,
       width = 12, height = 7.5, dpi = 300)

# =============== 图4：阈值扫描（候选数与物种覆盖） ===============
thr_grid <- seq(0.5, 0.9, by = 0.05)
base_unlabeled <- prob_long %>% filter(Observed == 0, !is.na(p_true))

p4_df <- map_dfr(thr_grid, function(th) {
  d <- base_unlabeled %>% filter(p_true >= th)
  tibble(threshold = th,
         n_candidates = nrow(d),
         n_species = d %>% distinct(Species) %>% nrow())
})

# 两个子图分别画（避免双轴争议）
p4a <- ggplot(p4_df, aes(x = threshold, y = n_candidates)) +
  geom_line() + geom_point() + theme_bw(base_size = 12) +
  labs(title = "Threshold sweep: number of candidates",
       x = "Threshold (p_true ≥ τ)", y = "#Candidates")
p4b <- ggplot(p4_df, aes(x = threshold, y = n_species)) +
  geom_line() + geom_point() + theme_bw(base_size = 12) +
  labs(title = "Threshold sweep: species coverage",
       x = "Threshold (p_true ≥ τ)", y = "#Species with ≥1 candidate")
ggsave(file.path(fig_dir, "fig4a_threshold_sweep_candidates.png"), p4a, width = 6.5, height = 4.2, dpi = 300)
ggsave(file.path(fig_dir, "fig4b_threshold_sweep_species.png"),   p4b, width = 6.5, height = 4.2, dpi = 300)

# =============== 图5：跨物种“热门候选突变”排行（Top-20） ===============
p5_df <- base_unlabeled %>%
  filter(p_true >= 0.6) %>%                   # 可与 PTRUE_THRESHOLD 一致
  distinct(Species, Mutation) %>%
  count(Mutation, name = "n_species") %>%
  arrange(desc(n_species)) %>% slice_head(n = 20) %>%
  mutate(Mutation = fct_reorder(Mutation, n_species))

p5 <- ggplot(p5_df, aes(x = Mutation, y = n_species)) +
  geom_col() + coord_flip() + theme_bw(base_size = 12) +
  labs(title = "Top mutations by cross-species high-confidence predictions",
       x = "Mutation", y = "#Species with p_true ≥ 0.6")
ggsave(file.path(fig_dir, "fig5_hot_mutations_across_species.png"), p5,
       width = 7.5, height = 6, dpi = 300)

# =============== 图6：Recall@K per fold 的箱线图（补充稳定性） ===============
p6 <- ggplot(eval_detail, aes(x = factor(K), y = Recall_at_K)) +
  geom_boxplot(outlier.size = 0.8) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  theme_bw(base_size = 12) +
  labs(title = "Mask-then-Recover: Recall@K per fold",
       x = "K", y = "Recall@K")
ggsave(file.path(fig_dir, "fig6_recall_at_k_box.png"), p6,
       width = 6.5, height = 4.2, dpi = 300)

## ===========================================================
## Focused species compare: Masked vs Unmasked Recall@K
## 指定物种：遮蔽 vs 未遮蔽 的 Recall@K 对比（箱线图）
## 依赖：已定义 X_dense, K_SET, MASK_FRAC, N_REPEATS, run_pu_rf_for_mut(), run_pu_rf_all()
## ===========================================================

suppressPackageStartupMessages({ library(ggplot2); library(forcats) })

# 1) 选你想比较的物种（示例，按需修改）
focus_species <- c(
  "Mycobacterium tuberculosis",
  "Pseudomonas aeruginosa",
  "Bacillus anthracis"
)

# 2) 主函数：对指定物种做遮蔽 vs 未遮蔽对比
eval_focus_species_mask_vs_unmask <- function(X_in, species_focus,
                                              mask_frac = MASK_FRAC,
                                              repeats = N_REPEATS,
                                              K_set = K_SET,
                                              Nout_dir = pre_dir) {
  X_in <- as.data.frame(X_in)
  all_rows <- list()
  
  for (r in seq_len(repeats)) {
    set.seed(3000 + r)
    
    # ---- 采样遮蔽 ----
    P <- which(as.matrix(X_in) == 1, arr.ind = TRUE)
    if (nrow(P) < 2) next
    nmask <- max(1, floor(nrow(P) * mask_frac))
    M <- P[sample(nrow(P), nmask), , drop = FALSE]
    
    X_mask <- X_in
    for (i in seq_len(nrow(M))) X_mask[M[i, 1], M[i, 2]] <- 0
    
    # ---- 训练并预测（遮蔽矩阵）----
    res_mask <- run_pu_rf_all(X_mask)
    
    # 未标注池的排序（每物种）
    ranks_tbl <- res_mask %>%
      dplyr::filter(Observed == 0, !is.na(p_true)) %>%
      dplyr::group_by(Species) %>%
      dplyr::arrange(dplyr::desc(p_true), .by_group = TRUE) %>%
      dplyr::mutate(rank = dplyr::row_number()) %>%
      dplyr::ungroup()
    
    # 遮蔽对（物种, 突变）
    masked_tbl <- tibble::tibble(
      Species = rownames(X_in)[M[, 1]],
      Mutation = colnames(X_in)[M[, 2]]
    )
    
    # 只保留关注物种的正例
    masked_focus <- masked_tbl %>% dplyr::filter(Species %in% species_focus)
    
    # ---- 对每个关注物种计算 Masked 和 Unmasked 的 Recall@K ----
    for (sp in species_focus) {
      # 该物种所有正例（原始）
      pos_all <- which(as.matrix(X_in)[rownames(X_in) == sp, ] == 1)
      mut_all <- colnames(X_in)[pos_all]
      
      # 本折被遮蔽的正例
      masked_sp <- masked_focus %>% dplyr::filter(Species == sp)
      
      # 本折未遮蔽（仍为1）的正例 = 全部正例 - 被遮蔽
      unmasked_mut <- setdiff(mut_all, masked_sp$Mutation)
      unmasked_sp  <- tibble::tibble(Species = sp, Mutation = unmasked_mut)
      
      # 未标注集合在该物种下的 p_true 列表（用于排名参照）
      unl_sp <- ranks_tbl %>% dplyr::filter(Species == sp) %>% dplyr::select(Mutation, p_true)
      
      ## ---- 1) Masked 真实排名（在未标注池里）----
      # 这些在 X_mask 下 Observed=0，所以若有有效 p_true 就会出现在 ranks_tbl
      hits_masked <- masked_sp %>% dplyr::left_join(
        ranks_tbl %>% dplyr::filter(Species == sp) %>% dplyr::select(Mutation, rank),
        by = "Mutation"
      )
      # p_true 为 NA 的样本在 ranks_tbl 匹配不到，视为 rank=Inf（永远不在Top-K）
      hits_masked$rank[is.na(hits_masked$rank)] <- Inf
      
      n_masked <- nrow(hits_masked)
      if (n_masked > 0) {
        for (K in K_set) {
          rec <- sum(hits_masked$rank <= K) / n_masked
          all_rows[[length(all_rows) + 1]] <- tibble::tibble(
            fold = r, Species = sp, Group = "Masked", K = K,
            n_items = n_masked, Recall_at_K = rec
          )
        }
      }
      
      ## ---- 2) Unmasked 反事实排名（若也当作未标注会排到多前）----
      # 取这些未遮蔽正例在 res_mask 里的 p_true 分数
      if (nrow(unmasked_sp) > 0) {
        sc_unmasked <- unmasked_sp %>% dplyr::left_join(
          res_mask %>% dplyr::filter(Species == sp) %>%
            dplyr::select(Mutation, p_true),
          by = "Mutation"
        )
        
        # 用未标注集合 unl_sp 的分布，计算“若加入未标注池”的反事实排名：
        # rank_cf = 1 + # {unlabeled p_true > p_true_of_unmasked}
        if (nrow(unl_sp) == 0) {
          # 该物种无未标注集合，无法计算排名 → 设为 Inf
          sc_unmasked$rank_cf <- Inf
        } else {
          v_unl <- unl_sp$p_true
          sc_unmasked$rank_cf <- purrr::map_dbl(sc_unmasked$p_true, function(s) {
            if (is.na(s)) return(Inf)
            sum(v_unl > s) + 1
          })
        }
        
        n_unmasked <- nrow(sc_unmasked)
        for (K in K_set) {
          rec <- sum(sc_unmasked$rank_cf <= K) / n_unmasked
          all_rows[[length(all_rows) + 1]] <- tibble::tibble(
            fold = r, Species = sp, Group = "Unmasked_cf", K = K,
            n_items = n_unmasked, Recall_at_K = rec
          )
        }
      }
    } # end species loop
  } # end repeats
  
  out <- dplyr::bind_rows(all_rows)
  out_file <- file.path(Nout_dir, "PU_RF_eval_focus_species_mask_vs_unmask.csv")
  readr::write_csv(out, out_file)
  message("✔ Saved: ", out_file)
  
  # 画箱线图：按物种分面
  p <- ggplot(out, aes(x = factor(K), y = Recall_at_K, fill = Group)) +
    geom_boxplot(outlier.size = 0.7, alpha = 0.85, position = position_dodge(width = 0.8)) +
    facet_wrap(~ Species, scales = "free_y") +
    labs(title = "Mask-then-Recover: Recall@K (masked vs unmasked, focus species)",
         x = "K (top-K per species)", y = "Recall@K") +
    theme_bw(base_size = 12)
  fig_file <- file.path(Nout_dir, "PU_RF_focus_species_mask_vs_unmask_box.png")
  ggsave(fig_file, p, width = 9, height = 5.5, dpi = 300)
  message("✔ Saved figure: ", fig_file)
  
  invisible(out)
}

# 运行
focus_out <- eval_focus_species_mask_vs_unmask(X_dense, focus_species)


# ========= 选择一个固定的K来排名与展示 =========
K_SELECT <- 20           # 主文固定K；想看扩展就再跑一次设20
MIN_MASKED_PER_SPECIES <- 6   # 过滤掉样本太少导致不稳定的物种

# 计算“所有物种”的 Masked vs Unmasked_cf（复用之前思路）
eval_mask_unmask_all <- function(X_in, mask_frac = MASK_FRAC, repeats = N_REPEATS, K = K_SELECT) {
  X_in <- as.data.frame(X_in)
  rows <- list()
  
  for (r in seq_len(repeats)) {
    set.seed(4000 + r)
    # ---- 遮蔽（可替换为分层+列保底版本）----
    P <- which(as.matrix(X_in) == 1, arr.ind = TRUE)
    nmask <- max(1, floor(nrow(P) * mask_frac))
    M <- P[sample(nrow(P), nmask), , drop = FALSE]
    X_mask <- X_in; for (i in seq_len(nrow(M))) X_mask[M[i,1], M[i,2]] <- 0
    
    res_mask <- run_pu_rf_all(X_mask)
    
    ranks_tbl <- res_mask %>%
      dplyr::filter(Observed == 0, !is.na(p_true)) %>%
      dplyr::group_by(Species) %>%
      dplyr::arrange(dplyr::desc(p_true), .by_group = TRUE) %>%
      dplyr::mutate(rank = dplyr::row_number()) %>%
      dplyr::ungroup()
    
    masked_tbl <- tibble::tibble(
      Species = rownames(X_in)[M[,1]],
      Mutation = colnames(X_in)[M[,2]]
    )
    
    # ---- 按物种计算 Recall@K（Masked 与 Unmasked_cf）----
    sp_list <- unique(rownames(X_in))
    for (sp in sp_list) {
      # Masked
      masked_sp <- masked_tbl %>% dplyr::filter(Species == sp)
      hits_masked <- masked_sp %>% dplyr::left_join(
        ranks_tbl %>% dplyr::filter(Species == sp) %>% dplyr::select(Mutation, rank),
        by = "Mutation"
      )
      if (nrow(hits_masked) > 0) {
        hits_masked$rank[is.na(hits_masked$rank)] <- Inf
        rec_m <- sum(hits_masked$rank <= K) / nrow(hits_masked)
        rows[[length(rows)+1]] <- tibble::tibble(
          fold=r, Species=sp, Group="Masked", K=K,
          n_items=nrow(hits_masked), Recall_at_K=rec_m
        )
      }
      
      # Unmasked_cf（反事实）
      pos_all <- which(as.matrix(X_in)[rownames(X_in) == sp, ] == 1)
      unmasked_mut <- setdiff(colnames(X_in)[pos_all], masked_sp$Mutation)
      if (length(unmasked_mut) > 0) {
        sc_unmasked <- tibble::tibble(Species=sp, Mutation=unmasked_mut) %>%
          dplyr::left_join(
            res_mask %>% dplyr::filter(Species == sp) %>% dplyr::select(Mutation, p_true),
            by="Mutation"
          )
        unl_sp <- ranks_tbl %>% dplyr::filter(Species==sp) %>% dplyr::select(Mutation, p_true)
        if (nrow(unl_sp) == 0) {
          rank_cf <- rep(Inf, nrow(sc_unmasked))
        } else {
          v_unl <- unl_sp$p_true
          rank_cf <- purrr::map_dbl(sc_unmasked$p_true, function(s) if (is.na(s)) Inf else sum(v_unl > s) + 1)
        }
        rec_u <- sum(rank_cf <= K) / length(rank_cf)
        rows[[length(rows)+1]] <- tibble::tibble(
          fold=r, Species=sp, Group="Unmasked_cf", K=K,
          n_items=length(rank_cf), Recall_at_K=rec_u
        )
      }
    }
  }
  
  out <- dplyr::bind_rows(rows)
  readr::write_csv(out, file.path(pre_dir, sprintf("PU_RF_eval_mask_unmask_all_K%d.csv", K)))
  out
}

all_eval <- eval_mask_unmask_all(X_dense, K = K_SELECT)

# 汇总到物种级别（取均值±SE），并按表现排序
summ_sp <- all_eval %>%
  dplyr::group_by(Species, Group) %>%
  dplyr::summarise(
    mean_recall = mean(Recall_at_K, na.rm = TRUE),
    se_recall   = sd(Recall_at_K, na.rm = TRUE) / sqrt(dplyr::n()),
    n_items_avg = mean(n_items, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(names_from = Group, values_from = c(mean_recall, se_recall, n_items_avg))

# 过滤样本过少的物种
summ_sp <- summ_sp %>% dplyr::filter(n_items_avg_Masked >= MIN_MASKED_PER_SPECIES)

# 定义排序指标：Masked 平均召回高、且与 Unmasked_cf 差距小
summ_sp <- summ_sp %>%
  dplyr::mutate(
    gap = (mean_recall_Unmasked_cf - mean_recall_Masked),      # 越小越好（越不依赖已见监督）
    rank_score = 0.7*mean_recall_Masked - 0.3*pmax(gap, 0)     # 可调权重
  ) %>%
  dplyr::arrange(dplyr::desc(rank_score))

readr::write_csv(summ_sp, file.path(pre_dir, sprintf("PU_RF_species_ranking_K%d.csv", K_SELECT)))

# 选前6个做展示
top_show <- summ_sp %>% dplyr::slice_head(n = 6) %>% dplyr::pull(Species)
top_show
suppressPackageStartupMessages({ library(ggplot2) })
plot_df <- all_eval %>% dplyr::filter(Species %in% top_show)

p <- ggplot(plot_df, aes(x = Group, y = Recall_at_K, fill = Group)) +
  geom_boxplot(outlier.size = 0.7, alpha = 0.85) +
  facet_wrap(~ Species, scales = "free_y") +
  labs(title = sprintf("Recall@%d (Masked vs Unmasked_cf) – Top species", K_SELECT),
       x = "", y = sprintf("Recall@%d", K_SELECT)) +
  theme_bw(base_size = 12)
ggsave(file.path(pre_dir, sprintf("fig_top_species_mask_vs_unmask_K%d.png", K_SELECT)),
       p, width = 9, height = 5.5, dpi = 300)

