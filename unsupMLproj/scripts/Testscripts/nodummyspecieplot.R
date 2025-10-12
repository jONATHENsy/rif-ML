# =============================================================
# clusteringplot_filtered_full.R – 使用 Manhattan 最佳三法，绘制所有突变的完整图
# =============================================================

message("\n🔍 clusteringplot_filtered_full.R started")

suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(uwot)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(UpSetR)
})

# 设置路径
base_dir <- "C:/Users/user/Desktop/D Drive/2025s1/BIOX7011/rif-ML/unsupMLproj"
fig_dir <- file.path(base_dir, "figures", "filtered_nodummy")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# 数据路径
mat_path <- file.path(base_dir, "output", "X_dense_midhigh.RDS")
label_file <- file.path(base_dir, "output", "cluster_labels_filtered.csv")

# 👇 读取并过滤物种
remove_species <- c("Vibrio parahaemolyticus", "Vibrio vulnificus", "Streptomyces lividans")

X_dense <- readRDS(mat_path)
mode(X_dense) <- "numeric"
X_dense[is.na(X_dense)] <- 0
X_dense <- X_dense[!(rownames(X_dense) %in% remove_species), , drop = FALSE]

cluster_labels <- read_csv(label_file, show_col_types = FALSE) %>%
  mutate(Unit = as.character(Unit)) %>%
  filter(!(Unit %in% remove_species))

# 执行 UMAP
set.seed(123)
umap_res <- umap(
  X_dense,
  metric = "manhattan",
  n_neighbors = min(15, nrow(X_dense)-1),
  min_dist = 0.3,
  verbose = TRUE
)

umap_df_base <- as_tibble(umap_res, .name_repair = "unique") %>%
  setNames(c("UMAP1", "UMAP2")) %>%
  mutate(Species = rownames(X_dense))

# 使用 Manhattan 距离下评分最高的三种方法
top_methods <- c("MANHATTAN_HDBSCAN", "MANHATTAN_GMM", "MANHATTAN_KMEANS")

for (method in top_methods) {
  umap_df <- umap_df_base %>%
    left_join(cluster_labels %>% select(Species = Unit, Cluster = !!sym(method)), by = "Species") %>%
    mutate(Cluster = factor(Cluster))
  
  # UMAP plot
  p <- ggplot(umap_df, aes(UMAP1, UMAP2, colour = Cluster)) +
    geom_point(size = 3, alpha = .85) +
    scale_colour_brewer(palette = "Set1", na.translate = FALSE) +
    theme_classic(base_size = 14) +
    labs(title = paste("UMAP +", method, "clusters (filtered, all mutations)"))
  
  ggsave(file.path(fig_dir, paste0("umap_", method, ".png")),
         p, width = 6, height = 5, dpi = 300)
  
  # Heatmap：全部突变
  sub_mat <- X_dense[umap_df$Species, , drop = FALSE]
  
  cluster_ids <- sort(unique(umap_df$Cluster))
  palette_colors <- brewer.pal(max(3, length(cluster_ids)), "Set1")
  cluster_cols <- setNames(palette_colors[seq_along(cluster_ids)], cluster_ids)
  
  row_anno <- rowAnnotation(
    Cluster = umap_df$Cluster,
    col = list(Cluster = cluster_cols),
    show_annotation_name = FALSE
  )
  
  ht <- Heatmap(
    sub_mat,
    name = "Mut",
    col = c("0" = "white", "1" = "steelblue"),
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    left_annotation = row_anno,
    row_split = umap_df$Cluster,
    width = unit(33, "cm"),
    height = unit(12, "cm")
  )
  
  png(file.path(fig_dir, paste0("heatmap_", method, "_allmut.png")), width = 1600, height = 1400, res = 120)
  draw(ht)
  dev.off()
}

# Upset Plot 1：物种集合交集
X_transpose <- t(X_dense)
mutation_df <- as.data.frame(X_transpose)
mutation_df$Mutation <- rownames(mutation_df)
upset_df <- mutation_df %>% relocate(Mutation, .after = last_col())

png(file.path(fig_dir, "mutation_species_upset.png"), width = 1600, height = 1000, res = 120)
upset(upset_df[, -ncol(upset_df)], nsets = 10, nintersects = 30,
      keep.order = TRUE, sets.bar.color = "steelblue", order.by = "freq")
dev.off()

# Upset Plot 2：突变在物种中的分布（前10）
mut_freq <- sort(colSums(X_dense), decreasing = TRUE)
top_mutations <- names(mut_freq)[1:min(10, length(mut_freq))]
submat <- X_dense[, top_mutations, drop = FALSE]
inverted_df <- as.data.frame(submat)
inverted_df$Species <- rownames(submat)
inverted_df <- inverted_df %>% relocate(Species, .after = last_col())

png(file.path(fig_dir, "top10mut_species_upset.png"), width = 1600, height = 1000, res = 120)
upset(inverted_df[, -ncol(inverted_df)], nsets = 10, nintersects = 30,
      keep.order = TRUE, sets.bar.color = "tomato", order.by = "freq",
      mainbar.y.label = "Number of species", sets.x.label = "Top mutations")
dev.off()

message("✅ Done plotting with top Manhattan clustering methods.")
