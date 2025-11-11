# =============================================================
# compareandplot.R — Auto-select top clustering schemes & plot
# - Run multiple (metric × algo) on X_dense → mean silhouette
# - Pick top_k schemes (optionally restricted to one metric)
# - Plot UMAP & Heatmap with a shared cluster color map
# - Save tables (silhouettes, labels) and figures
# =============================================================

suppressPackageStartupMessages({
  pkgs <- c(
    "tidyverse","Matrix","uwot","ComplexHeatmap","circlize",
    "RColorBrewer","grid","dbscan","cluster","proxy","mclust",
    "patchwork","ggplotify","UpSetR"
  )
  need <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
  lapply(pkgs, library, character.only = TRUE)
})

message("\n🔎 compareandplot.R started")

# ----------------------------- Config -----------------------------
# Read paths from main if available; otherwise use a safe default.
base_dir <- getOption("rifml.base_dir",
                      file.path(getwd(), "rif-ML", "unsupMLproj"))
fig_dir  <- getOption("rifml.fig_dir", file.path(base_dir, "figures"))
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# X_dense: prefer option set by main (checked version), else try local default.
xdense_path_opt <- getOption("rifml.xdense", file.path(base_dir, "data", "X_dense_checked.RDS"))
if (file.exists(xdense_path_opt)) {
  X_dense <- readRDS(xdense_path_opt)
} else if (exists("X_dense")) {
  # allow running interactively where X_dense is already in memory
  X_dense <- get("X_dense")
} else {
  stop("Cannot find X_dense. Expected RDS at: ", xdense_path_opt)
}

# Ensure numeric 0/1 matrix with species rownames
X_dense <- as.matrix(X_dense)
mode(X_dense) <- "numeric"
if (is.null(rownames(X_dense))) stop("X_dense must have rownames = species.")
X_dense[is.na(X_dense)] <- 0
unit_names <- rownames(X_dense)

# Analysis knobs
restrict_metric  <- NULL       # e.g., "manhattan"; or NULL = across all metrics
top_k            <- 3          # how many schemes to show
min_cluster_pts  <- 5          # HDBSCAN/DBSCAN minPts
dbscan_eps       <- 1.2        # DBSCAN eps
umap_neighbors   <- 15
umap_min_dist    <- 0.3
set.seed(123)

# Heatmap padding tweak (mm): nudge right to avoid left labels being clipped
shift_right_mm <- 28

# Output folder for csv tables from this script
out_dir <- file.path(base_dir, "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------- Clustering pipeline ----------------------
run_pipeline <- function(mat, metric = "euclidean",
                         cluster_method = "hdbscan",
                         minPts = 5, eps = 1.2, seed = 123,
                         n_neighbors = 15, min_dist = 0.3) {
  set.seed(seed)
  # UMAP on distances if using Jaccard; otherwise on data
  if (metric == "jaccard") {
    dmat <- proxy::dist(mat, method = "Jaccard")
    um   <- uwot::umap(as.matrix(dmat), input = "dist", metric = "precomputed",
                       n_neighbors = n_neighbors, min_dist = min_dist, verbose = FALSE)
  } else {
    um   <- uwot::umap(mat, metric = metric, n_neighbors = n_neighbors,
                       min_dist = min_dist, verbose = FALSE)
  }
  
  clust <- switch(tolower(cluster_method),
                  "hdbscan" = dbscan::hdbscan(as.data.frame(um), minPts = minPts)$cluster,
                  "kmeans"  = kmeans(um, centers = 4)$cluster,
                  "dbscan"  = dbscan::dbscan(um, eps = eps, minPts = minPts)$cluster,
                  "gmm"     = mclust::Mclust(um)$classification,
                  stop("Unknown cluster_method: ", cluster_method)
  )
  
  list(umap = um, cluster = clust)
}

schemes <- expand.grid(
  metric = c("euclidean","manhattan","cosine"),
  method = c("hdbscan","kmeans","dbscan","gmm"),
  stringsAsFactors = FALSE
) |>
  dplyr::mutate(name = paste(toupper(metric), toupper(method), sep = "_"))

results <- purrr::pmap(schemes, ~run_pipeline(
  X_dense, metric = ..1, cluster_method = ..2,
  minPts = min_cluster_pts, eps = dbscan_eps, seed = 123,
  n_neighbors = umap_neighbors, min_dist = umap_min_dist
))

# ----------------------- Silhouette selection ---------------------
mean_sil <- purrr::map2_dbl(results, schemes$name, function(res, nm) {
  cl <- res$cluster
  um <- res$umap
  valid_idx <- which(!is.na(cl) & cl > 0)
  if (length(valid_idx) < 2 || length(unique(cl[valid_idx])) < 2) return(NA_real_)
  sil <- cluster::silhouette(cl[valid_idx], dist(um[valid_idx, ]))
  mean(sil[, 3])
})

sil_df <- tibble::tibble(
  Scheme = schemes$name,
  Metric = toupper(schemes$metric),
  Method = toupper(schemes$method),
  Mean_Silhouette = round(mean_sil, 4),
  Valid = !is.na(mean_sil)
)

# Select top methods (optionally restricted to a metric)
sel_df <- dplyr::filter(sil_df, Valid)
if (!is.null(restrict_metric)) {
  sel_df <- dplyr::filter(sel_df, Metric == toupper(restrict_metric))
}
if (nrow(sel_df) == 0) {
  warning("No valid schemes under restrict_metric = ", restrict_metric,
          ". Falling back to all metrics.")
  sel_df <- dplyr::filter(sil_df, Valid)
}
n_take <- min(as.integer(top_k), nrow(sel_df))
top_methods <- sel_df |>
  dplyr::slice_max(order_by = Mean_Silhouette, n = n_take, with_ties = FALSE) |>
  dplyr::pull(Scheme)

message("⭐ Top methods: ", paste(top_methods, collapse = ", "))

# Save tables (full silhouettes + all cluster labels)
write_csv(sil_df, file.path(out_dir, "silhouette_all_methods_auto.csv"))

cluster_tbl <- tibble::tibble(Unit = unit_names)
for (i in seq_len(nrow(schemes))) {
  cluster_tbl[[schemes$name[i]]] <- results[[i]]$cluster
}
write_csv(cluster_tbl, file.path(out_dir, "cluster_labels_all_methods_auto.csv"))

# ------------------- Unified cluster color mapping ----------------
all_clusters <- unique(unlist(lapply(results, function(r) r$cluster)))
all_clusters <- sort(unique(na.omit(as.character(all_clusters))))
base_pal <- RColorBrewer::brewer.pal(max(3, min(9, length(all_clusters))), "Set1")
cluster_colors <- setNames(rep(base_pal, length.out = length(all_clusters)), all_clusters)
if ("0" %in% names(cluster_colors)) cluster_colors[["0"]] <- "grey70"  # HDBSCAN noise
message("🎨 Cluster→color: ", paste(names(cluster_colors), collapse = ","))

# --------------------------- Plot helpers -------------------------
build_umap_plot <- function(method_name) {
  i <- which(schemes$name == method_name)
  res <- results[[i]]
  umap_df <- as_tibble(res$umap, .name_repair = "unique") |>
    setNames(c("UMAP1","UMAP2")) |>
    mutate(Species = unit_names,
           Cluster = factor(as.character(res$cluster), levels = names(cluster_colors)))
  ggplot(umap_df, aes(UMAP1, UMAP2, colour = Cluster)) +
    geom_point(size = 3, alpha = .85) +
    scale_colour_manual(values = cluster_colors, na.value = "grey80") +
    theme_classic(base_size = 14) +
    theme(plot.title = element_blank(), legend.position = "none")
}

build_heatmap_topN <- function(method_name, topN = 30) {
  # Use cluster blocks on rows and plot topN most frequent mutations overall
  i <- which(schemes$name == method_name)
  res <- results[[i]]
  
  umap_df <- as_tibble(res$umap, .name_repair = "unique") |>
    setNames(c("UMAP1","UMAP2")) |>
    mutate(Species = unit_names,
           Cluster = factor(as.character(res$cluster), levels = names(cluster_colors)))
  
  # Align matrix rows to umap_df order
  idx <- match(umap_df$Species, rownames(X_dense))
  keep <- which(!is.na(idx))
  umap_df <- umap_df[keep, , drop = FALSE]
  sub_mat <- X_dense[idx[keep], , drop = FALSE]
  
  # Pick topN columns by global frequency
  mut_freq <- colSums(sub_mat, na.rm = TRUE)
  topN <- min(topN, ncol(sub_mat))
  top_cols <- names(sort(mut_freq, decreasing = TRUE))[seq_len(topN)]
  sub_mat <- sub_mat[, top_cols, drop = FALSE]
  
  # Discrete "off/on_cluster" matrix so the 1's inherit the cluster color
  state_mat <- matrix("off",
                      nrow = nrow(sub_mat), ncol = ncol(sub_mat),
                      dimnames = dimnames(sub_mat))
  row_clusters <- as.character(umap_df$Cluster)
  for (r in seq_len(nrow(state_mat))) {
    on_idx <- which(sub_mat[r, ] == 1)
    if (length(on_idx)) state_mat[r, on_idx] <- paste0("on_", row_clusters[r])
  }
  
  col_map <- c(off = "white")
  for (cl in levels(umap_df$Cluster)) {
    if (is.na(cl)) next
    col_map[paste0("on_", cl)] <- if (cl == "0") "grey70" else cluster_colors[[cl]]
  }
  
  row_anno <- ComplexHeatmap::rowAnnotation(
    Cluster = umap_df$Cluster,
    col = list(Cluster = cluster_colors),
    show_annotation_name = FALSE,
    width = unit(4, "mm")
  )
  
  g <- grid::grid.grabExpr(
    ComplexHeatmap::draw(
      ComplexHeatmap::Heatmap(
        state_mat,
        name = "Mut",
        col = col_map,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_split = umap_df$Cluster,
        show_row_names = TRUE,  row_names_gp = grid::gpar(fontsize = 9),
        show_column_names = TRUE, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7),
        right_annotation = row_anno,
        na_col = "white"
      ),
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      padding = unit(c(6, 6, 6, shift_right_mm), "mm")
    )
  )
  ggplotify::as.ggplot(g)
}

# ----------------------- Render patchwork figs --------------------
methods_to_show <- top_methods[1:min(3, length(top_methods))]

# 3×UMAP panel
umap_list <- lapply(methods_to_show, build_umap_plot)
fig_umap <- wrap_plots(umap_list, nrow = 1, guides = "collect") +
  plot_annotation(tag_levels = "A")
ggsave(file.path(fig_dir, "FIG_UMAP_patchwork_abc.png"),
       fig_umap, width = 18, height = 5.5, dpi = 300)
ggsave(file.path(fig_dir, "FIG_UMAP_patchwork_abc.pdf"),
       fig_umap, width = 18, height = 5.5)

# 3×Top-30 heatmaps panel
heat_list <- lapply(methods_to_show, build_heatmap_topN)
fig_heat <- wrap_plots(heat_list, nrow = 1, guides = "collect") +
  plot_annotation(tag_levels = "A")
ggsave(file.path(fig_dir, "FIG_HeatmapTop30_patchwork_abc.png"),
       fig_heat, width = 24, height = 8, dpi = 300)
ggsave(file.path(fig_dir, "FIG_HeatmapTop30_patchwork_abc.pdf"),
       fig_heat, width = 24, height = 8)

# ----------------------- Optional: single full heatmaps -----------
# If you still want one full heatmap per top method (all columns),
# flip the flag below to TRUE. Otherwise skip to UpSet.
make_full_heatmaps <- FALSE
if (make_full_heatmaps) {
  for (method in methods_to_show) {
    i <- which(schemes$name == method)
    res <- results[[i]]
    umap_df <- as_tibble(res$umap, .name_repair = "unique") |>
      setNames(c("UMAP1", "UMAP2")) |>
      mutate(Species = unit_names,
             Cluster = factor(as.character(res$cluster), levels = names(cluster_colors)))
    idx <- match(umap_df$Species, rownames(X_dense))
    keep <- which(!is.na(idx))
    umap_df <- umap_df[keep, , drop = FALSE]
    sub_mat <- X_dense[idx[keep], , drop = FALSE]
    
    row_anno <- ComplexHeatmap::rowAnnotation(
      Cluster = umap_df$Cluster,
      col = list(Cluster = cluster_colors),
      show_annotation_name = FALSE,
      width = unit(4, "mm")
    )
    
    png_w <- max(3600, ncol(sub_mat) * 28)
    png_h <- max(1600, nrow(sub_mat) * 22)
    
    png(file.path(fig_dir, paste0("heatmap_full_", method, ".png")),
        width = png_w, height = png_h, res = 300, type = "cairo-png")
    ComplexHeatmap::draw(
      ComplexHeatmap::Heatmap(
        sub_mat, name = "Mut",
        col = c("0" = "white", "1" = "steelblue"),
        cluster_rows = FALSE, cluster_columns = TRUE,
        row_split = umap_df$Cluster,
        show_row_names = TRUE, row_names_side = "left", row_names_gp = gpar(fontsize = 10),
        show_column_names = FALSE,
        right_annotation = row_anno,
        use_raster = TRUE, raster_device = "png"
      ),
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      padding = unit(c(6, 6, 6, shift_right_mm), "mm")
    )
    dev.off()
  }
}

# ----------------------- Compact UpSet overview -------------------
# A single, quick UpSet of all mutations across species.
X_t <- t(X_dense)
mutation_df <- as.data.frame(X_t)
mutation_df$Mutation <- rownames(mutation_df)
mutation_df <- dplyr::relocate(mutation_df, Mutation)
png(file.path(fig_dir, "mutation_species_upset.png"), width = 1600, height = 1000, res = 120)
UpSetR::upset(mutation_df[, -1], nsets = min(25, ncol(mutation_df) - 1),
              nintersects = 30, keep.order = TRUE,
              sets.bar.color = "steelblue", order.by = "freq")
dev.off()

message("✅ Done — auto-selected top methods: ", paste(top_methods, collapse = ", "))
