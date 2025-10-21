# =============================================================
# compareclustering_plus_all.R
# -------------------------------------------------------------
# 同时对原始 X_dense、残差 X_resid、筛选子集 X_dense_filtered 做聚类分析
# 输出 cluster label + silhouette score 表格，用于性能比较
# 动态调整 n_neighbors 以适配小样本数据（如 core species）
# =============================================================

message("\n>>> compareclustering_plus_all.R started")

suppressPackageStartupMessages({
  library(tidyverse)
  library(uwot)
  library(dbscan)
  library(cluster)
  library(proxy)
  library(mclust)
})

if (!exists("base_dir"))
  base_dir <- "C:/Users/user/Desktop/D Drive/2025s1/BIOX7011/rif-ML/unsupMLproj"

out_dir <- file.path(base_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------
# 1. Helper function
# -------------------------------------------------------------
run_pipeline <- function(mat, metric = "euclidean",
                         cluster_method = "hdbscan",
                         minPts = 5, seed = 123) {
  set.seed(seed)
  n_neighbors <- min(15, max(2, nrow(mat) - 1))  # 动态调整
  
  if (metric == "jaccard") {
    dmat <- proxy::dist(mat, method = "Jaccard")
    um <- uwot::umap(as.matrix(dmat), input = "dist", metric = "precomputed",
                     n_neighbors = n_neighbors, min_dist = 0.3, verbose = FALSE)
  } else {
    um <- uwot::umap(mat, metric = metric, n_neighbors = n_neighbors,
                     min_dist = 0.3, verbose = FALSE)
  }
  clust <- switch(cluster_method,
                  "hdbscan" = dbscan::hdbscan(as.data.frame(um), minPts = minPts)$cluster,
                  "kmeans"  = kmeans(um, centers = 4)$cluster,
                  "dbscan"  = dbscan::dbscan(um, eps = 1.2, minPts = minPts)$cluster,
                  "gmm"     = Mclust(um)$classification)
  list(umap = um, cluster = clust)
}

# -------------------------------------------------------------
# 2. Clustering schemes
# -------------------------------------------------------------
schemes <- expand.grid(
  metric = c("euclidean", "manhattan", "cosine"),
  method = c("hdbscan", "kmeans", "dbscan", "gmm"),
  stringsAsFactors = FALSE
) %>% mutate(name = paste(toupper(metric), toupper(method), sep = "_"))

# -------------------------------------------------------------
# 3. Run clustering for each matrix variant
# -------------------------------------------------------------
variant_list <- list(
  original   = if (exists("X_dense")) X_dense else stop("X_dense not found"),
  residual   = if (exists("X_resid")) X_resid else stop("X_resid not found"),
  filtered   = if (exists("X_dense_filtered")) X_dense_filtered else stop("X_dense_filtered not found")
)

for (variant in names(variant_list)) {
  cat("\n>>> Running clustering on:", variant, "\n")
  mat <- variant_list[[variant]]
  
  if (nrow(mat) < 4) {
    warning(paste("Skipped", variant, "because it has < 4 rows (too small for clustering)"))
    next
  }
  
  results <- pmap(schemes, ~run_pipeline(mat, metric = ..1, cluster_method = ..2))
  unit_names <- rownames(mat)
  
  cluster_tbl <- purrr::imap_dfc(results, function(res, i) {
    tibble(!!schemes$name[i] := res$cluster)
  }) %>% mutate(Unit = unit_names, .before = 1)
  
  mean_sil <- map2_dbl(results, schemes$name, function(res, nm) {
    cl <- res$cluster
    um <- res$umap
    valid_idx <- which(!is.na(cl) & cl > 0)
    if (length(valid_idx) < 2 || length(unique(cl[valid_idx])) < 2) {
      return(NA_real_)
    }
    sil <- silhouette(cl[valid_idx], dist(um[valid_idx, ]))
    mean(sil[, 3])
  })
  
  sil_df <- tibble(
    Scheme = schemes$name,
    Mean_Silhouette = round(mean_sil, 3),
    Valid = !is.na(mean_sil),
    Source = variant
  )
  
  write_csv(cluster_tbl, file.path(out_dir, paste0("cluster_labels_", variant, ".csv")))
  write_csv(sil_df,      file.path(out_dir, paste0("silhouette_scores_", variant, ".csv")))
}

message("\n✓ compareclustering_plus_all.R finished: all variant results saved")
