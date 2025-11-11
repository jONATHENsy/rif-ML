# ===========================================================
# PU-learning with Random Forest for mutation prediction
# - Positive–Unlabeled (PU) learning via U-bagging + RandomForest
# - Input: X_dense (rows = Species, cols = Mutations, values in {0,1})
# - Outputs: probabilities, candidate lists, evaluation (Recall@K) and figures
# ===========================================================

# --------------------------- Packages ---------------------------
suppressPackageStartupMessages({
  pkgs <- c(
    "dplyr","tidyr","tibble","stringr","purrr",
    "randomForest","readr","ggplot2","forcats","Matrix"
  )
  need <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
  lapply(pkgs, library, character.only = TRUE)
})

# --------------------------- Config & I/O ------------------------
# Base directory:
# 1) Prefer option set by main: options(rifml.base_dir = "<path>")
# 2) Otherwise default to "<cwd>/rif-ML/unsupMLproj"
base_dir <- getOption("rifml.base_dir",
                      file.path(getwd(), "rif-ML", "unsupMLproj"))

# Path to master X_dense (read-only). You can override with:
# options(rifml.xdense_master = "<path-to-rds>")
X_DENSE_MASTER_PATH <- getOption(
  "rifml.xdense_master",
  file.path(base_dir, "output", "X_dense_master.rds")
)

# Output directories
pre_dir <- file.path(base_dir, "predict")
fig_dir <- file.path(pre_dir, "figs01")
dir.create(pre_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# RNG seed
set.seed(42)

# --------------------------- Load X_dense ------------------------
# 1) Load the immutable master
X_MASTER <- readRDS(X_DENSE_MASTER_PATH)
lockBinding("X_MASTER", .GlobalEnv)  # prevent accidental overwrite

# 2) Make a working copy for PU learning
prep_X_for_pu <- function(X) {
  if (inherits(X, "Matrix")) X <- as.matrix(X)
  X <- as.data.frame(X, check.names = FALSE)
  
  if (is.null(rownames(X))) stop("X_dense must have rownames as Species.")
  
  X[] <- lapply(X, function(col) {
    if (is.factor(col)) col <- as.character(col)
    v <- suppressWarnings(as.numeric(col))
    v[is.na(v)] <- 0
    as.integer(v > 0)  # hard binarization
  })
  X
}
X_dense <- prep_X_for_pu(X_MASTER)

# --------------------------- Hyper-parameters --------------------
NTREE <- 500        # trees per RF
B_BAGS <- 25        # number of U-bagging rounds
NEG_RATIO <- 1.0    # unlabeled subsample size relative to #positives
MIN_POS <- 3        # skip mutations with < MIN_POS positives
USE_UMAP <- FALSE   # keep FALSE to avoid leakage unless you recompute per target
DROP_SAME_SITE <- TRUE
ADD_EFFORT <- TRUE  # add effort = rowSums(X_dense)

# Candidate selection
PTRUE_THRESHOLD <- 0.7
TOPK_PER_SPECIES <- 10

# Evaluation (Mask-then-Recover)
MASK_FRAC <- 0.20
N_REPEATS <- 5
K_SET     <- c(1,5,10,15,20,25,30,40,50)

# --------------------------- Utilities ---------------------------
stamp <- function(...) message(format(Sys.time(), "%T"), " — ", paste0(...))
safe_write_csv <- function(x, path) { dir.create(dirname(path), FALSE, TRUE); readr::write_csv(x, path) }

# --------------------------- Basic checks ------------------------
X_dense <- as.data.frame(X_dense)
if (is.null(rownames(X_dense))) stop("X_dense must have rownames as Species.")

X_dense[is.na(X_dense)] <- 0
X_dense[] <- lapply(X_dense, function(col) {
  v <- suppressWarnings(as.numeric(col)); v[is.na(v)] <- 0
  pmin(pmax(v, 0), 1)
})
all_muts <- colnames(X_dense)

# --------------------------- Same-site parsing -------------------
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
                      NA_character_, paste0(mat[,2], "_", mat[,4])) # gene_pos
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

# --------------------------- PU-RF core --------------------------
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
  y <- as.integer(X_dense[[target_mut]])  # 1 = observed positive, 0 = unlabeled
  X_df <- X_dense
  X_df[[target_mut]] <- NULL  # prevent leakage
  
  if (drop_same_site) X_df <- drop_same_site_cols(X_df, col_meta, target_mut)
  if (add_effort)     X_df$effort <- rowSums(as.matrix(X_dense))
  if (use_umap) stop("USE_UMAP=TRUE requires recomputing UMAP per-target to avoid leakage.")
  
  pos_idx <- which(y == 1)
  unl_idx <- which(y == 0)
  
  if (length(pos_idx) < MIN_POS) {
    return(tibble(
      Species  = rownames(X_df),
      p_obs    = NA_real_,
      p_true   = NA_real_,
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
  
  # Elkan–Noto correction: p(y=1|x) = p(s=1|x) / c, with c = E[p(s=1|x) | y=1]
  c_est <- mean(p_obs[y == 1])
  if (is.na(c_est) || c_est <= 0) {
    p_true <- rep(NA_real_, length(p_obs))
  } else {
    p_true <- p_obs / c_est
    p_true[p_true > 1] <- 1
  }
  
  tibble(
    Species  = rownames(X_df),
    p_obs    = p_obs,
    p_true   = p_true,
    Observed = y,
    Mutation = target_mut
  )
}

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

# --------------------------- Train & Save -------------------------
stamp("Running PU-RF for all mutations ...")
res_list <- lapply(all_muts, function(m) {
  tryCatch(
    run_pu_rf_for_mut(m, X_dense, col_meta),
    error = function(e) { warning(sprintf("PU-RF failed for %s: %s", m, e$message)); NULL }
  )
})
res <- bind_rows(res_list)

prob_long_file <- file.path(pre_dir, "PU_RF_prob_long.csv")
safe_write_csv(res, prob_long_file)
stamp("Saved: ", prob_long_file)

# --------------------------- Candidates ---------------------------
candidates_threshold <- res %>%
  filter(Observed == 0, !is.na(p_true), p_true >= PTRUE_THRESHOLD) %>%
  arrange(desc(p_true))
safe_write_csv(candidates_threshold, file.path(pre_dir, "PU_RF_candidates_threshold.csv"))

candidates_topk <- res %>%
  filter(Observed == 0, !is.na(p_true)) %>%
  group_by(Species) %>%
  arrange(desc(p_true), .by_group = TRUE) %>%
  slice_head(n = TOPK_PER_SPECIES) %>%
  ungroup()
safe_write_csv(candidates_topk, file.path(pre_dir, "PU_RF_candidates_topk.csv"))

stamp("Candidate lists saved.")

# --------------------------- Evaluation (Mask-then-Recover) -------
eval_mask_then_recover <- function(X_in, mask_frac = MASK_FRAC, repeats = N_REPEATS, K_set = K_SET) {
  X_in <- as.data.frame(X_in)
  rs <- list()
  
  for (r in seq_len(repeats)) {
    set.seed(2025 + r)
    
    # 1) Sample positives to mask
    P <- which(as.matrix(X_in) == 1, arr.ind = TRUE)
    if (nrow(P) < 2) next
    nmask <- max(1, floor(nrow(P) * mask_frac))
    M <- P[sample(nrow(P), nmask), , drop = FALSE]
    
    # 2) Mask to 0 (become unlabeled)
    X_mask <- X_in
    for (i in seq_len(nrow(M))) X_mask[M[i, 1], M[i, 2]] <- 0
    
    # 3) Train on masked data & predict
    res_mask <- run_pu_rf_all(X_mask)
    
    # 4) Build masked positive table
    masked_tbl <- tibble(
      Species  = rownames(X_in)[M[, 1]],
      Mutation = colnames(X_in)[M[, 2]]
    )
    
    # 5) Rank among unlabeled per species
    ranks_tbl <- res_mask %>%
      filter(Observed == 0, !is.na(p_true)) %>%
      group_by(Species) %>%
      arrange(desc(p_true), .by_group = TRUE) %>%
      mutate(rank = row_number()) %>%
      ungroup()
    
    hits <- ranks_tbl %>% inner_join(masked_tbl, by = c("Species", "Mutation"))
    
    # 6) Compute Recall@K
    rec_df <- tibble(
      fold = r,
      K = K_set,
      Recall_at_K = sapply(K_set, function(K) mean(hits$rank <= K))
    )
    rs[[r]] <- rec_df
  }
  
  out <- bind_rows(rs)
  safe_write_csv(out, file.path(pre_dir, "PU_RF_eval_recall_at_k.csv"))
  
  summ <- out %>%
    group_by(K) %>%
    summarise(
      mean_recall = mean(Recall_at_K, na.rm = TRUE),
      se_recall   = sd(Recall_at_K, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  safe_write_csv(summ, file.path(pre_dir, "PU_RF_eval_recall_at_k_summary.csv"))
  
  # Figure
  p <- ggplot(summ, aes(x = K, y = mean_recall)) +
    geom_line() + geom_point() +
    geom_errorbar(aes(ymin = mean_recall - se_recall, ymax = mean_recall + se_recall), width = 0.2) +
    labs(title = "Mask-then-Recover: Recall@K",
         x = "K (top-K predictions per species)",
         y = "Mean Recall@K ± SE") +
    theme_bw(base_size = 12)
  ggsave(file.path(fig_dir, "PU_RF_recall_at_k.png"), p, width = 6.5, height = 4.2, dpi = 300)
  
  invisible(list(detail = out, summary = summ))
}

stamp("Running simplified PU evaluation (Mask-then-Recover, Recall@K) ...")
eval_out <- eval_mask_then_recover(X_dense)
stamp("Evaluation complete.")

# --------------------------- Quick diagnostics --------------------
per_species <- res %>%
  group_by(Species) %>%
  summarise(
    n_unlabeled  = sum(Observed == 0, na.rm = TRUE),
    n_valid_pred = sum(Observed == 0 & !is.na(p_true)),
    n_above_thr  = sum(Observed == 0 & !is.na(p_true) & p_true >= PTRUE_THRESHOLD),
    .groups = "drop"
  ) %>% arrange(n_valid_pred)
print(per_species, n = 10)

mut_summary <- res %>%
  group_by(Mutation) %>%
  summarise(
    pos_count = sum(Observed == 1, na.rm = TRUE),
    any_valid = any(!is.na(p_true)),
    .groups = "drop"
  ) %>% arrange(pos_count)
print(head(mut_summary, 20))
invisible(mean(!mut_summary$any_valid))

missing_in_topk <- setdiff(rownames(X_dense), unique(candidates_topk$Species))
if (length(missing_in_topk)) stamp("Species missing in TopK: ", paste(missing_in_topk, collapse = ", "))

# --------------------------- Figures (global) ---------------------
# Fig 1: candidates per species (threshold)
p1_df <- candidates_threshold %>% count(Species, name = "n_candidates") %>% arrange(desc(n_candidates))
p1 <- ggplot(p1_df, aes(x = reorder(Species, n_candidates), y = n_candidates)) +
  geom_col() + coord_flip() +
  labs(title = "Candidates per species (thresholded)",
       x = "Species", y = "Number of candidates (p_true ≥ τ)") +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir, "fig1_candidates_per_species_threshold.png"), p1, width = 7.5, height = 6, dpi = 300)

# Fig 2: top-K per species (facet)
focus_species <- c("Mycobacterium tuberculosis","Bacillus anthracis","Pseudomonas aeruginosa")
p2_df <- candidates_topk %>%
  filter(Species %in% focus_species) %>%
  group_by(Species) %>%
  arrange(desc(p_true), .by_group = TRUE) %>%
  mutate(Mutation = forcats::fct_reorder(Mutation, p_true))
p2 <- ggplot(p2_df, aes(x = Mutation, y = p_true)) +
  geom_col() + coord_flip() +
  facet_wrap(~ Species, scales = "free_y") +
  labs(title = "Top-K predicted mutations per species",
       x = "Mutation", y = "p_true") +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir, "fig2_topk_by_species_bar.png"), p2, width = 9, height = 6, dpi = 300)

# Fig 3: heatmap-like tile of top-N per species
TOP_N <- 15
p3_df <- res %>%
  filter(Observed == 0, !is.na(p_true)) %>%
  group_by(Species) %>%
  slice_max(p_true, n = TOP_N, with_ties = FALSE) %>%
  ungroup()
species_order <- p3_df %>% group_by(Species) %>% summarise(m = mean(p_true), .groups="drop") %>% arrange(desc(m)) %>% pull(Species)
mut_order <- p3_df %>% group_by(Mutation) %>% summarise(m = mean(p_true), .groups="drop") %>% arrange(desc(m)) %>% pull(Mutation)
p3 <- ggplot(p3_df, aes(x = factor(Mutation, mut_order), y = factor(Species, species_order), fill = p_true)) +
  geom_tile() +
  scale_fill_gradient(name = "p_true", limits = c(0,1)) +
  labs(title = paste0("Top-", TOP_N, " predicted (unlabeled) per species"),
       x = "Mutation", y = "Species") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(file.path(fig_dir, "fig3_heatmap_topN_per_species.png"), p3, width = 12, height = 7.5, dpi = 300)

# Fig 4a/4b: threshold sweeps
thr_grid <- seq(0.5, 0.9, by = 0.05)
base_unlabeled <- res %>% filter(Observed == 0, !is.na(p_true))
p4_df <- purrr::map_dfr(thr_grid, function(th) {
  d <- base_unlabeled %>% filter(p_true >= th)
  tibble(threshold = th, n_candidates = nrow(d), n_species = d %>% distinct(Species) %>% nrow())
})
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

# Fig 5: cross-species “hot” mutations (Top-20 by coverage)
p5_df <- base_unlabeled %>%
  filter(p_true >= 0.6) %>%
  distinct(Species, Mutation) %>%
  count(Mutation, name = "n_species") %>%
  arrange(desc(n_species)) %>%
  slice_head(n = 20) %>%
  mutate(Mutation = forcats::fct_reorder(Mutation, n_species))
p5 <- ggplot(p5_df, aes(x = Mutation, y = n_species)) +
  geom_col() + coord_flip() + theme_bw(base_size = 12) +
  labs(title = "Top mutations by cross-species high-confidence predictions",
       x = "Mutation", y = "#Species with p_true ≥ 0.6")
ggsave(file.path(fig_dir, "fig5_hot_mutations_across_species.png"), p5, width = 7.5, height = 6, dpi = 300)

# Fig 6: Recall@K per fold (boxplot)
eval_detail <- readr::read_csv(file.path(pre_dir, "PU_RF_eval_recall_at_k.csv"), show_col_types = FALSE)
p6 <- ggplot(eval_detail, aes(x = factor(K), y = Recall_at_K)) +
  geom_boxplot(outlier.size = 0.8) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  theme_bw(base_size = 12) +
  labs(title = "Mask-then-Recover: Recall@K per fold",
       x = "K", y = "Recall@K")
ggsave(file.path(fig_dir, "fig6_recall_at_k_box.png"), p6, width = 6.5, height = 4.2, dpi = 300)

# --------------------------- Focus species (optional) ------------
focus_species <- c("Mycobacterium tuberculosis","Pseudomonas aeruginosa","Bacillus anthracis")

eval_focus_species_mask_vs_unmask <- function(X_in, species_focus,
                                              mask_frac = MASK_FRAC,
                                              repeats = N_REPEATS,
                                              K_set = K_SET,
                                              Nout_dir = pre_dir) {
  X_in <- as.data.frame(X_in); all_rows <- list()
  
  for (r in seq_len(repeats)) {
    set.seed(3000 + r)
    P <- which(as.matrix(X_in) == 1, arr.ind = TRUE)
    if (nrow(P) < 2) next
    nmask <- max(1, floor(nrow(P) * mask_frac))
    M <- P[sample(nrow(P), nmask), , drop = FALSE]
    X_mask <- X_in; for (i in seq_len(nrow(M))) X_mask[M[i, 1], M[i, 2]] <- 0
    
    res_mask <- run_pu_rf_all(X_mask)
    
    ranks_tbl <- res_mask %>%
      dplyr::filter(Observed == 0, !is.na(p_true)) %>%
      dplyr::group_by(Species) %>%
      dplyr::arrange(dplyr::desc(p_true), .by_group = TRUE) %>%
      dplyr::mutate(rank = dplyr::row_number()) %>%
      dplyr::ungroup()
    
    masked_tbl <- tibble::tibble(
      Species = rownames(X_in)[M[, 1]],
      Mutation = colnames(X_in)[M[, 2]]
    )
    masked_focus <- masked_tbl %>% dplyr::filter(Species %in% species_focus)
    
    for (sp in species_focus) {
      pos_all <- which(as.matrix(X_in)[rownames(X_in) == sp, ] == 1)
      mut_all <- colnames(X_in)[pos_all]
      masked_sp <- masked_focus %>% dplyr::filter(Species == sp)
      
      hits_masked <- masked_sp %>% dplyr::left_join(
        ranks_tbl %>% dplyr::filter(Species == sp) %>% dplyr::select(Mutation, rank),
        by = "Mutation"
      )
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
      
      unmasked_mut <- setdiff(mut_all, masked_sp$Mutation)
      if (length(unmasked_mut) > 0) {
        sc_unmasked <- tibble::tibble(Species = sp, Mutation = unmasked_mut) %>%
          dplyr::left_join(
            res_mask %>% dplyr::filter(Species == sp) %>% dplyr::select(Mutation, p_true),
            by = "Mutation"
          )
        unl_sp <- ranks_tbl %>% dplyr::filter(Species == sp) %>% dplyr::select(Mutation, p_true)
        if (nrow(unl_sp) == 0) {
          rank_cf <- rep(Inf, nrow(sc_unmasked))
        } else {
          v_unl <- unl_sp$p_true
          rank_cf <- purrr::map_dbl(sc_unmasked$p_true, function(s) if (is.na(s)) Inf else sum(v_unl > s) + 1)
        }
        for (K in K_set) {
          rec <- sum(rank_cf <= K) / length(rank_cf)
          all_rows[[length(all_rows) + 1]] <- tibble::tibble(
            fold = r, Species = sp, Group = "Unmasked_cf", K = K,
            n_items = length(rank_cf), Recall_at_K = rec
          )
        }
      }
    }
  }
  
  out <- dplyr::bind_rows(all_rows)
  readr::write_csv(out, file.path(Nout_dir, "PU_RF_eval_focus_species_mask_vs_unmask.csv"))
  
  p <- ggplot(out, aes(x = factor(K), y = Recall_at_K, fill = Group)) +
    geom_boxplot(outlier.size = 0.7, alpha = 0.85, position = position_dodge(width = 0.8)) +
    facet_wrap(~ Species, scales = "free_y") +
    labs(title = "Mask-then-Recover: Recall@K (masked vs unmasked, focus species)",
         x = "K (top-K per species)", y = "Recall@K") +
    theme_bw(base_size = 12)
  ggsave(file.path(Nout_dir, "PU_RF_focus_species_mask_vs_unmask_box.png"),
         p, width = 9, height = 5.5, dpi = 300)
  
  invisible(out)
}
focus_out <- eval_focus_species_mask_vs_unmask(X_dense, focus_species)

# --------------------------- All-species at fixed K --------------
K_SELECT <- 20
MIN_MASKED_PER_SPECIES <- 6

eval_mask_unmask_all <- function(X_in, mask_frac = MASK_FRAC, repeats = N_REPEATS, K = K_SELECT) {
  X_in <- as.data.frame(X_in)
  rows <- list()
  
  for (r in seq_len(repeats)) {
    set.seed(4000 + r)
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
    
    for (sp in unique(rownames(X_in))) {
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
      
      # Unmasked_cf
      pos_all <- which(as.matrix(X_in)[rownames(X_in) == sp, ] == 1)
      unmasked_mut <- setdiff(colnames(X_in)[pos_all], masked_sp$Mutation)
      if (length(unmasked_mut) > 0) {
        sc_unmasked <- tibble::tibble(Species=sp, Mutation=unmasked_mut) %>%
          dplyr::left_join(
            res_mask %>% dplyr::filter(Species == sp) %>% dplyr::select(Mutation, p_true),
            by="Mutation"
          )
        unl_sp <- ranks_tbl %>% dplyr::filter(Species==sp) %>% dplyr::select(Mutation, p_true)
        rank_cf <- if (nrow(unl_sp) == 0) rep(Inf, nrow(sc_unmasked)) else {
          v_unl <- unl_sp$p_true
          purrr::map_dbl(sc_unmasked$p_true, function(s) if (is.na(s)) Inf else sum(v_unl > s) + 1)
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

summ_sp <- all_eval %>%
  dplyr::group_by(Species, Group) %>%
  dplyr::summarise(
    mean_recall = mean(Recall_at_K, na.rm = TRUE),
    se_recall   = sd(Recall_at_K, na.rm = TRUE) / sqrt(dplyr::n()),
    n_items_avg = mean(n_items, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(names_from = Group, values_from = c(mean_recall, se_recall, n_items_avg)) %>%
  dplyr::filter(n_items_avg_Masked >= MIN_MASKED_PER_SPECIES) %>%
  dplyr::mutate(
    gap = (mean_recall_Unmasked_cf - mean_recall_Masked),
    rank_score = 0.7*mean_recall_Masked - 0.3*pmax(gap, 0)
  ) %>%
  dplyr::arrange(dplyr::desc(rank_score))

readr::write_csv(summ_sp, file.path(pre_dir, sprintf("PU_RF_species_ranking_K%d.csv", K_SELECT)))

top_show <- summ_sp %>% dplyr::slice_head(n = 6) %>% dplyr::pull(Species)
plot_df <- all_eval %>% dplyr::filter(Species %in% top_show)

p <- ggplot(plot_df, aes(x = Group, y = Recall_at_K, fill = Group)) +
  geom_boxplot(outlier.size = 0.7, alpha = 0.85) +
  facet_wrap(~ Species, scales = "free_y") +
  labs(title = sprintf("Recall@%d (Masked vs Unmasked_cf) – Top species", K_SELECT),
       x = "", y = sprintf("Recall@%d", K_SELECT)) +
  theme_bw(base_size = 12)
ggsave(file.path(pre_dir, sprintf("fig_top_species_mask_vs_unmask_K%d.png", K_SELECT)),
       p, width = 9, height = 5.5, dpi = 300)

stamp("PU-learning finished. Artifacts in: 'predict/' and 'predict/figs01/'.")
