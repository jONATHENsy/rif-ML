# =============================================================
# compareclustering.R
# -------------------------------------------------------------
# Run UMAP + multiple clustering algorithms (HDBSCAN, KMeans, DBSCAN, GMM)
# with several distance metrics and compare the resulting cluster assignments.
# Executed automatically by main.R (after X_dense is created).
# =============================================================

message("\n>>> compareclustering.R started")

if (!exists("X_dense"))
  stop("X_dense not found – run clusteringplot.R (or main.R) first")

if (!exists("base_dir")) {
  base_dir <- "C:/Users/user/Desktop/D Drive/2025s1/BIOX7011/rif-ML/unsupMLproj"
}

# -------------------------------------------------------------
# 0. Libraries
# -------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(uwot)
  library(dbscan)
  library(cluster)
  library(proxy)
  library(mclust)   # For GMM
})

# -------------------------------------------------------------
# 1. helper – run one UMAP + clustering pipeline
# -------------------------------------------------------------
run_pipeline <- function(mat, metric = "euclidean",
                         cluster_method = "hdbscan",
                         minPts = 5, seed = 123) {
  
  set.seed(seed)
  
  if (metric == "jaccard") {
    dmat <- proxy::dist(mat, method = "Jaccard")
    um   <- uwot::umap(
      as.matrix(dmat),
      input       = "dist",
      metric      = "precomputed",
      n_neighbors = 15,
      min_dist    = 0.3,
      verbose     = FALSE
    )
  } else {
    um   <- uwot::umap(
      mat,
      metric      = metric,
      n_neighbors = 15,
      min_dist    = 0.3,
      verbose     = FALSE
    )
  }
  
  clust <- switch(cluster_method,
                  "hdbscan" = dbscan::hdbscan(as.data.frame(um), minPts = minPts)$cluster,
                  "kmeans"  = kmeans(um, centers = 4)$cluster,
                  "dbscan"  = dbscan::dbscan(um, eps = 1.2, minPts = minPts)$cluster,
                  "gmm"     = Mclust(um)$classification)
  
  list(umap = um, cluster = clust)
}

# -------------------------------------------------------------
# 2. metric × algorithm schemes to evaluate
# -------------------------------------------------------------
schemes <- expand.grid(
  metric = c("euclidean", "manhattan", "cosine"),
  method = c("hdbscan", "kmeans", "dbscan", "gmm"),
  stringsAsFactors = FALSE
)
schemes <- schemes %>% mutate(name = paste(toupper(metric), toupper(method), sep = "_"))

# -------------------------------------------------------------
# 3. run all combinations
# -------------------------------------------------------------
results <- pmap(schemes, ~run_pipeline(X_dense, metric = ..1, cluster_method = ..2))

# combine cluster labels into one data-frame
unit_names <- rownames(X_dense)

cluster_tbl <- purrr::imap_dfc(results, function(res, i) {
  tibble(!!schemes$name[i] := res$cluster)
}) %>% mutate(Unit = unit_names, .before = 1)

# -------------------------------------------------------------
# 4. print basic summaries
# -------------------------------------------------------------
walk(schemes$name, \(nm) {
  cat("\n###", nm, "\n")
  print(table(cluster_tbl[[nm]]))
})

# -------------------------------------------------------------
# 5. silhouette scores (robust version)
# -------------------------------------------------------------
mean_sil <- map2_dbl(results, schemes$name, function(res, nm) {
  cl <- res$cluster
  um <- res$umap
  valid_idx <- which(!is.na(cl) & cl > 0)
  
  if (length(valid_idx) < 2 || length(unique(cl[valid_idx])) < 2) {
    return(NA_real_)  # Skip invalid clustering
  }
  
  sil <- silhouette(cl[valid_idx], dist(um[valid_idx, ]))
  mean(sil[, 3])
})

sil_df <- tibble(
  Scheme = schemes$name,
  Mean_Silhouette = round(mean_sil, 3),
  Valid = !is.na(mean_sil)
)
print(sil_df)

# -------------------------------------------------------------
# 6. save outputs
# -------------------------------------------------------------
out_dir <- file.path(base_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write_csv(cluster_tbl, file.path(out_dir, "cluster_labels_all_methods.csv"))
write_csv(sil_df,      file.path(out_dir, "silhouette_all_methods.csv"))

message("✓ compareclustering.R finished – results written to ",
        normalizePath(out_dir))
