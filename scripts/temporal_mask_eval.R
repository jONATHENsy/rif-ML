## ===========================================================
## Temporal Mask → Recover (按年份遮蔽评估)
## 读取 labmuts.csv 提取 (Species, Mutation) 首次发表年份，
## 以年份截点 T_cut 将“未来真阳性”遮成0，比较 Masked vs Unmasked_cf 的 Recall@K。
## 输出：详细结果CSV + 每个K的箱线图（按年份分面）
## Author: you + gpt哥
## ===========================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
  library(purrr); library(stringr); library(ggplot2); library(tibble)
  library(randomForest)
})

## -------------------------
## 0) 路径与全局参数
## -------------------------
# 输入：labmuts.csv 的绝对路径（Windows 建议用 / 斜杠）
labmuts_path <- "C:/Users/user/Desktop/D Drive/2025s1/BIOX7011/rif-ML/unsupMLproj/output/labmuts.csv"

# 输出目录
if (!exists("pre_dir")) pre_dir <- "predict"
out_dir <- file.path(pre_dir, "temporal_eval")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 评估设置
T_GRID  <- c(2012, 2015, 2018, 2020)  # 年份截点，可自行调整
K_SET   <- c(10, 20)                  # 主文K=10，扩展K=20
N_REPEATS <- 5                        # 每个截点重复次数（因列保底会触发采样）
MIN_POS <- 2                          # 遮蔽后每突变列至少保留的正例数（防止整列失效）

# PU-RF 关键超参（与主脚本保持一致）
NTREE <- 800
B_BAGS <- 50
NEG_RATIO <- 2.0
DROP_SAME_SITE <- TRUE
ADD_EFFORT <- TRUE

## -------------------------------------------------
## 1) 基础检查：需要 X_dense 已在内存中
## -------------------------------------------------
stopifnot(exists("X_dense"))
X_dense <- as.data.frame(X_dense)
if (is.null(rownames(X_dense))) stop("X_dense must have rownames as Species.")
X_dense[is.na(X_dense)] <- 0
X_dense[] <- lapply(X_dense, function(col) {
  v <- suppressWarnings(as.numeric(col)); v[is.na(v)] <- 0; pmin(pmax(v,0),1)
})
X_dense <- as.data.frame(X_dense)

## -------------------------------------------------
## 2) 从 labmuts 提取 (Species, Mutation) 的首次年份
## -------------------------------------------------
extract_first_year_tbl <- function(csv_path) {
  d <- readr::read_csv(csv_path, show_col_types = FALSE)
  # 选突变列（按优先顺序）
  mut_cols <- c("AA_mut_name_Ecoli","AA_mut_name","Nt_mut_name_Ecoli","Nt_mut_name")
  mut_col <- mut_cols[mut_cols %in% colnames(d)][1]
  if (is.na(mut_col)) stop("Cannot find a mutation column in labmuts.csv")
  
  # 从 Ref_code 中抽取 4 位年份
  if (!"Ref_code" %in% colnames(d)) stop("labmuts.csv must contain Ref_code")
  d$Year <- stringr::str_extract(d$Ref_code, "(18|19|20)\\d{2}") %>% as.integer()
  
  d %>%
    select(Species, Mutation = all_of(mut_col), Year) %>%
    filter(!is.na(Species), !is.na(Mutation), !is.na(Year)) %>%
    group_by(Species, Mutation) %>%
    summarise(first_year = min(Year), .groups = "drop")
}

first_year_tbl <- extract_first_year_tbl(labmuts_path)

## -------------------------------------------------
## 3) PU-RF（与主代码一致的实现/或直接 source 你的主脚本）
## -------------------------------------------------
parse_mut_cols <- function(mut_names) {
  mat <- stringr::str_match(mut_names, "^([^_]+)_([A-Z\\*])?(\\d+)([A-Z\\*]+)?$")
  tibble(
    Mutation = mut_names,
    gene = mat[,2],
    refAA = mat[,3],
    pos   = suppressWarnings(as.integer(mat[,4])),
    altAA = mat[,5],
    site_key = ifelse(is.na(mat[,2]) | is.na(mat[,4]), NA_character_, paste0(mat[,2], "_", mat[,4]))
  )
}
drop_same_site_cols <- function(X_df, col_meta, target_mut) {
  sk <- col_meta$site_key[col_meta$Mutation == target_mut]
  if (length(sk) == 0 || is.na(sk)) return(X_df)
  drop_set <- col_meta$Mutation[col_meta$site_key == sk & col_meta$Mutation != target_mut]
  keep_cols <- setdiff(colnames(X_df), drop_set)
  X_df[, keep_cols, drop = FALSE]
}
run_pu_rf_for_mut <- function(target_mut, X_dense, col_meta,
                              ntree=NTREE, B=B_BAGS, neg_ratio=NEG_RATIO,
                              drop_same_site=DROP_SAME_SITE, add_effort=ADD_EFFORT) {
  if (!target_mut %in% colnames(X_dense)) return(NULL)
  y <- as.integer(X_dense[[target_mut]])
  X_df <- X_dense; X_df[[target_mut]] <- NULL
  if (drop_same_site) X_df <- drop_same_site_cols(X_df, col_meta, target_mut)
  if (add_effort) X_df$effort <- rowSums(as.matrix(X_dense))
  pos_idx <- which(y==1); unl_idx <- which(y==0)
  if (length(pos_idx) < MIN_POS) {
    return(tibble(Species = rownames(X_df), p_obs = NA_real_, p_true = NA_real_, Observed = y, Mutation = target_mut))
  }
  p_obs_mat <- matrix(NA_real_, nrow=nrow(X_df), ncol=B)
  for (b in seq_len(B)) {
    set.seed(100 + b)
    k <- ceiling(length(pos_idx) * neg_ratio)
    neg_sample <- if (length(unl_idx) >= k) sample(unl_idx, k) else unl_idx
    idx <- c(pos_idx, neg_sample)
    rf <- randomForest(x = X_df[idx, , drop=FALSE], y = as.factor(y[idx]), ntree = ntree)
    p_obs_mat[, b] <- predict(rf, X_df, type="prob")[,2]
  }
  p_obs <- rowMeans(p_obs_mat, na.rm=TRUE)
  # 带收缩的 c 估计
  n_pos <- sum(y==1); c_raw <- mean(p_obs[y==1])
  C_PRIOR <- 0.5; LAMBDA <- if (n_pos <= 3) 10 else 5
  c_est <- ((n_pos * c_raw) + (LAMBDA * C_PRIOR)) / (n_pos + LAMBDA)
  c_est <- max(min(c_est, 0.95), 0.05)
  p_true <- pmin(1, p_obs / c_est)
  
  tibble(Species = rownames(X_df), p_obs=p_obs, p_true=p_true, Observed=y, Mutation=target_mut)
}
run_pu_rf_all <- function(X_in) {
  X_in <- as.data.frame(X_in)
  all_muts <- colnames(X_in)
  col_meta <- parse_mut_cols(all_muts)
  res_list <- lapply(all_muts, function(m) {
    tryCatch(run_pu_rf_for_mut(m, X_in, col_meta), error=function(e) NULL)
  })
  bind_rows(res_list)
}

## -------------------------------------------------
## 4) Temporal Mask → Recover 主流程
## -------------------------------------------------
temporal_mask_vs_unmask <- function(X_in, first_year_tbl, T_grid = T_GRID, K_set = K_SET,
                                    repeats = N_REPEATS, out_dir = out_dir) {
  X_in <- as.data.frame(X_in)
  
  # universe of valid pairs in X_dense
  universe <- tidyr::expand_grid(Species = rownames(X_in), Mutation = colnames(X_in))
  
  all_rows <- list()
  
  for (T_cut in T_grid) {
    message("=== Cut year: ", T_cut, " ===")
    
    # 未来真阳性（将被遮蔽的目标池）
    mask_pool <- first_year_tbl %>%
      filter(first_year > T_cut) %>%
      semi_join(universe, by = c("Species","Mutation"))
    
    if (nrow(mask_pool) == 0) {
      warning("No future positives beyond T_cut = ", T_cut)
      next
    }
    
    pos_per_mut <- colSums(as.matrix(X_in) == 1)                   # 有列名
    cap_per_mut0 <- pmax(0, pos_per_mut - MIN_POS)
    cap_per_mut0 <- setNames(as.numeric(cap_per_mut0), colnames(X_in))  # 关键：恢复列名
    
    
    for (r in seq_len(repeats)) {
      set.seed(7000 + T_cut + r)
      
      # —— 按列容量采样：每个突变列最多遮 cap_per_mut0[col] 个
      cap_per_mut <- cap_per_mut0
      M_list <- list()
      for (mut in unique(mask_pool$Mutation)) {
        cap <- unname(cap_per_mut[mut])
        cap <- ifelse(is.na(cap), 0, cap)
        pool_mut <- mask_pool %>% filter(Mutation == mut)
        if (cap <= 0 || nrow(pool_mut) == 0) next
        take_n <- min(cap, nrow(pool_mut))
        take_idx <- sample(seq_len(nrow(pool_mut)), size = take_n, replace = FALSE)
        take <- pool_mut[take_idx, , drop = FALSE]
        M_list[[mut]] <- take
      }
      M <- bind_rows(M_list)
      
      # 生成遮蔽矩阵
      X_mask <- X_in
      if (nrow(M) > 0) {
        for (i in seq_len(nrow(M))) {
          X_mask[M$Species[i], M$Mutation[i]] <- 0
        }
      }
      
      # 训练并预测
      res_mask <- run_pu_rf_all(X_mask)
      
      # 未标注池的排序
      ranks_tbl <- res_mask %>%
        filter(Observed == 0, !is.na(p_true)) %>%
        group_by(Species) %>%
        arrange(desc(p_true), .by_group = TRUE) %>%
        mutate(rank = row_number()) %>%
        ungroup()
      
      # 逐物种计算 Masked 与 Unmasked_cf 的 Recall@K
      sp_list <- rownames(X_in)
      for (sp in sp_list) {
        # 该物种在 mask_pool 里的被遮蔽候选
        masked_sp <- M %>% filter(Species == sp)
        
        # 全部正例（原始 X_in）
        pos_all <- which(as.matrix(X_in)[rownames(X_in) == sp, ] == 1)
        mut_all <- colnames(X_in)[pos_all]
        unmasked_mut <- setdiff(mut_all, masked_sp$Mutation)
        unl_sp <- ranks_tbl %>% filter(Species == sp) %>% select(Mutation, p_true)
        
        # Masked 真阳性在未标注池里的实际排名
        hits_m <- masked_sp %>% left_join(ranks_tbl %>% filter(Species == sp) %>% select(Mutation, rank),
                                          by = "Mutation")
        hits_m$rank[is.na(hits_m$rank)] <- Inf
        n_masked <- nrow(hits_m)
        
        # Unmasked_cf（反事实）：若也放入未标注池，它们会排在多前
        if (length(unmasked_mut) > 0) {
          sc_u <- tibble(Species = sp, Mutation = unmasked_mut) %>%
            left_join(res_mask %>% filter(Species == sp) %>% select(Mutation, p_true), by="Mutation")
          if (nrow(unl_sp) == 0) {
            rank_cf <- rep(Inf, nrow(sc_u))
          } else {
            v_unl <- unl_sp$p_true
            rank_cf <- map_dbl(sc_u$p_true, function(s) if (is.na(s)) Inf else sum(v_unl > s) + 1)
          }
          n_unmasked <- length(rank_cf)
        } else {
          rank_cf <- integer(0); n_unmasked <- 0
        }
        
        for (K in K_set) {
          if (n_masked > 0) {
            rec_m <- sum(hits_m$rank <= K) / n_masked
            all_rows[[length(all_rows)+1]] <- tibble(
              T_cut=T_cut, fold=r, Species=sp, Group="Masked", K=K,
              n_items=n_masked, Recall_at_K=rec_m
            )
          }
          if (n_unmasked > 0) {
            rec_u <- sum(rank_cf <= K) / n_unmasked
            all_rows[[length(all_rows)+1]] <- tibble(
              T_cut=T_cut, fold=r, Species=sp, Group="Unmasked_cf", K=K,
              n_items=n_unmasked, Recall_at_K=rec_u
            )
          }
        }
      }
    }
  }
  
  out <- bind_rows(all_rows)
  out_file <- file.path(out_dir, "PU_RF_temporal_mask_vs_unmask_detail.csv")
  write_csv(out, out_file); message("✔ Saved: ", out_file)
  
  # 汇总表（均值±SE）
  summ <- out %>%
    group_by(T_cut, K, Group) %>%
    summarise(mean_recall = mean(Recall_at_K, na.rm=TRUE),
              se_recall   = sd(Recall_at_K, na.rm=TRUE)/sqrt(n()),
              n_pairs_avg = mean(n_items, na.rm=TRUE),
              .groups="drop")
  summ_file <- file.path(out_dir, "PU_RF_temporal_mask_vs_unmask_summary.csv")
  write_csv(summ, summ_file); message("✔ Saved: ", summ_file)
  
  # 画图：每个 K 一张箱线图（按 T_cut 分面）
  for (K in K_set) {
    p <- out %>% filter(K==K) %>%
      ggplot(aes(x=Group, y=Recall_at_K, fill=Group)) +
      geom_boxplot(outlier.size=0.7, alpha=0.85) +
      facet_wrap(~ T_cut, nrow = 1) +
      coord_cartesian(ylim = c(0,1)) +
      theme_bw(base_size = 12) +
      labs(title = sprintf("Temporal Mask→Recover: Recall@%d (Masked vs Unmasked_cf)", K),
           x = "", y = sprintf("Recall@%d", K))
    fig_file <- file.path(out_dir, sprintf("fig_temporal_mask_vs_unmask_K%d.png", K))
    ggsave(fig_file, p, width = 10, height = 4, dpi = 300)
    message("✔ Saved figure: ", fig_file)
  }
  
  invisible(list(detail = out, summary = summ))
}

## -------------------------
## 5) 运行
## -------------------------
res_temporal <- temporal_mask_vs_unmask(X_dense, first_year_tbl)
message("🏁 Temporal evaluation done.")
