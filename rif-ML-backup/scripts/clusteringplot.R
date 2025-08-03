# =============================================================
# clusteringplot.R – UMAP + HDBSCAN + heat-map visualisations
# -------------------------------------------------------------
# Can be sourced by main.R or run independently
# =============================================================

message("\n>>> clusteringplot.R started")

## -----------------------------------------------------------
## 1. Dependencies (already loaded in main.R, but keep for safety)
## -----------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(uwot)
  library(dbscan)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(ggrepel)
})

## -----------------------------------------------------------
## 2. Locate input objects / files
## -----------------------------------------------------------
if (!exists("base_dir"))
  base_dir <- "D:/2025s1/BIOX7011/rif-ML/unsupMLproj"   # fallback

mut_file <- file.path(base_dir, "output", "labmuts.csv")

# if X_dense & muts_completed already exist (main.R) we reuse them
if (!exists("muts_completed"))
  muts_completed <- read_csv(mut_file, show_col_types = FALSE)

if (!exists("X_dense")) {
  # ---------- rebuild 0/1 matrix ----------
  unit_tbl <- muts_completed %>%
    mutate(Unit = paste(Species, RefSeq_ID, sep = "__")) %>%
    select(Unit, Orthologous_AA = AA_mut_name_Ecoli) %>%
    distinct()
  
  mut_index <- unit_tbl %>%
    distinct(Orthologous_AA) %>%
    arrange(Orthologous_AA) %>%
    mutate(idx = row_number())
  
  enc <- inner_join(unit_tbl, mut_index, by = "Orthologous_AA")
  
  X <- sparseMatrix(
    i = as.integer(factor(enc$Unit, levels = unique(enc$Unit))),
    j = enc$idx,
    x = 1,
    dims = c(n_distinct(enc$Unit), nrow(mut_index)),
    dimnames = list(unique(enc$Unit), mut_index$Orthologous_AA)
  )
  
  X_dense <- as.matrix(X)
  X_dense[X_dense > 0] <- 1
}

## -----------------------------------------------------------
## 3.  UMAP  +  HDBSCAN
## -----------------------------------------------------------
set.seed(123)
umap_res <- uwot::umap(
  X_dense,
  metric       = "manhattan",   # change if you prefer
  n_neighbors  = 15,
  min_dist     = 0.3,
  verbose      = TRUE
)

umap_df <- as_tibble(umap_res, .name_repair = "unique") %>%
  setNames(c("UMAP1", "UMAP2")) %>%
  mutate(Unit = rownames(X_dense))

hdb <- hdbscan(umap_df[, c("UMAP1", "UMAP2")], minPts = 5)
umap_df$Cluster <- factor(hdb$cluster)

## -----------------------------------------------------------
## 4.  UMAP scatter plot (saved)
## -----------------------------------------------------------
p_umap <- ggplot(umap_df, aes(UMAP1, UMAP2, colour = Cluster)) +
  geom_point(size = 3, alpha = .85) +
  scale_colour_brewer(palette = "Set1", na.translate = FALSE) +
  theme_classic(base_size = 14) +
  labs(title = "UMAP + HDBSCAN clusters")

fig_dir <- file.path(base_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(fig_dir, "umap_clusters.png"),
       p_umap, width = 6, height = 5, dpi = 300)

## -----------------------------------------------------------
## 5.  Heat-map of 30 most frequent mutations
## -----------------------------------------------------------
mut_freq    <- colSums(X_dense)
top30       <- names(sort(mut_freq, decreasing = TRUE))[1:30]
sub_mat     <- X_dense[umap_df$Unit, top30, drop = FALSE]

cluster_ids <- sort(unique(umap_df$Cluster))
cluster_cols <- setNames(brewer.pal(length(cluster_ids), "Set1")[seq_along(cluster_ids)],
                         cluster_ids)

row_anno <- rowAnnotation(
  Cluster = umap_df$Cluster,
  col     = list(Cluster = cluster_cols),
  show_annotation_name = FALSE
)

highlight_species <- c("Escherichia coli")
row_labels <- gsub("__.*$", "", umap_df$Unit)
label_gp   <- gpar(fontface = ifelse(row_labels %in% highlight_species,
                                     "bold", "plain"),
                   fontsize = 6)

ht <- Heatmap(
  sub_mat,
  name             = "Mut",
  col              = c("0" = "white", "1" = "steelblue"),
  cluster_rows     = FALSE,
  cluster_columns  = TRUE,
  show_row_names   = TRUE,
  row_names_gp     = label_gp,
  left_annotation  = row_anno,
  row_split        = umap_df$Cluster,
  width            = unit(8, "cm"),
  height           = unit(12, "cm")
)

png(file.path(fig_dir, "cluster_heatmap.png"), width = 1000, height = 1400, res = 120)
draw(ht)
dev.off()

## -----------------------------------------------------------
## 6.  Write data tables for supervised phase
## -----------------------------------------------------------
write_csv(umap_df,
          file.path(base_dir, "output", "umap_coords.csv"))

cluster_summary <- umap_df %>%
  count(Cluster, name = "Num_Units") %>%
  arrange(Cluster)

write_csv(cluster_summary,
          file.path(base_dir, "output", "cluster_summary.csv"))

message("✓ clusteringplot.R finished – artefacts saved to ",
        normalizePath(fig_dir))
