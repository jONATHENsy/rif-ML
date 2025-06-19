# =============================================================
# compareclustering.R
# -------------------------------------------------------------
# Run UMAP + HDBSCAN with several distance metrics and
# compare the resulting cluster assignments.
# Executed automatically by main.R (after X_dense is created).
# =============================================================

message("\n>>> compareclustering.R started")

if (!exists("X_dense"))
  stop("X_dense not found – run clusteringplot.R (or main.R) first")

library(tidyverse)
library(uwot)
library(dbscan)
library(cluster)        # silhouette()
library(proxy)          # extra distance functions

# -------------------------------------------------------------
# 1. helper – run one UMAP + HDBSCAN pipeline
# -------------------------------------------------------------
run_pipeline <- function(mat, metric = "euclidean",
                         minPts = 5, seed = 123) {
  
  set.seed(seed)
  
  if (metric == "jaccard") {
    # Jaccard needs a pre-computed distance matrix
    dmat <- proxy::dist(mat, method = "Jaccard")
    um   <- uwot::umap(
      dmat,
      input       = "dist",
      n_neighbors = 15,
      min_dist    = 0.3,
      verbose     = FALSE)
  } else {
    um   <- uwot::umap(
      mat,
      metric      = metric,
      n_neighbors = 15,
      min_dist    = 0.3,
      verbose     = FALSE)
  }
  
  hdb  <- dbscan::hdbscan(as.data.frame(um), minPts = minPts)
  
  list(
    umap    = um,                 # 2-D coordinates
    cluster = hdb$cluster         # integer labels (0 = noise)
  )
}

# -------------------------------------------------------------
# 2. distance schemes to evaluate
# -------------------------------------------------------------
schemes <- tribble(
  ~name,        ~metric,
  "EUCLIDEAN",  "euclidean",
  "MANHATTAN",  "manhattan",
  "COSINE",     "cosine",
  "JACCARD",    "jaccard"
)

# -------------------------------------------------------------
# 3. run them all
# -------------------------------------------------------------
results <- map(schemes$metric, ~run_pipeline(X_dense, .x))

# combine cluster labels into one data-frame
unit_names <- rownames(X_dense)

cluster_tbl <- purrr::imap_dfc(results, \(res, i) {
  tibble(!!schemes$name[i] := res$cluster)
}) |>
  mutate(Unit = unit_names, .before = 1)

# -------------------------------------------------------------
# 4. print basic summaries
# -------------------------------------------------------------
walk(schemes$name, \(nm) {
  cat("\n###", nm, "\n")
  print(table(cluster_tbl[[nm]]))
})

cat("\nContingency (EUCLIDEAN vs MANHATTAN):\n")
print(table(cluster_tbl$EUCLIDEAN, cluster_tbl$MANHATTAN))

# -------------------------------------------------------------
# 5. mean silhouette for each scheme (higher = better)
#    silhouette is computed on 2-D UMAP Euclidean distances
# -------------------------------------------------------------
mean_sil <- map2_dbl(results, schemes$name, \(res, nm) {
  sil <- silhouette(res$cluster, dist(res$umap))
  mean(sil[, 3])
})

sil_df <- tibble(Scheme = schemes$name,
                 Mean_Silhouette = round(mean_sil, 3))
print(sil_df)

# -------------------------------------------------------------
# 6. save outputs (cluster table + silhouette table)
# -------------------------------------------------------------
out_dir <- file.path(base_dir, "output")  # change if you like
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write_csv(cluster_tbl,
          file.path(out_dir, "cluster_labels_all_metrics.csv"))
write_csv(sil_df,
          file.path(out_dir, "silhouette_summary.csv"))

message("✓ compareclustering.R finished – results written to ",
        normalizePath(out_dir))

