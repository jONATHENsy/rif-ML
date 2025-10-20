# =============================================================
# clusteringplot_plus_all.R – 针对不同数据矩阵 (原始/残差/筛选)
# 分别绘制 UMAP + 热图 + upset plot
# =============================================================

message("\n>>> clusteringplot_plus_all.R started")

suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(uwot)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(UpSetR)
})

# 设置基本路径
if (!exists("base_dir"))
  # =============================================================
# clusteringplot_plus_all.R – 针对不同数据矩阵 (原始/残差/筛选)
# 分别绘制 UMAP + 热图 + upset plot
# =============================================================

message("\n>>> clusteringplot_plus_all.R started")

suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(uwot)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(UpSetR)
})

# 设置基本路径
if (!exists("base_dir"))
  base_dir <- "D:/2025s1/BIOX7011/rif-ML/unsupMLproj"

label_variants <- c("original", "residual", "filtered")

# ========== ⬇️ 新增：重建 clean X_dense 并保存为 RDS ==========

message("🧬 Rebuilding clean X_dense from muts_completed ...")
mut_file <- file.path(base_dir, "output", "labmuts.csv")
muts_completed <- read_csv(mut_file, show_col_types = FALSE)

unit_tbl <- muts_completed %>%
  mutate(Species_clean = str_replace(Strain_ID, "_[^_]+$", ""),
         Species_clean = str_replace_all(Species_clean, "_", " ")) %>%
  select(Species = Species_clean, Orthologous_AA = AA_mut_name_Ecoli) %>%
  filter(!is.na(Species), !is.na(Orthologous_AA)) %>%
  distinct()

mut_index <- unit_tbl %>%
  distinct(Orthologous_AA) %>%
  arrange(Orthologous_AA) %>%
  mutate(idx = row_number())

enc <- inner_join(unit_tbl, mut_index, by = "Orthologous_AA")

X <- sparseMatrix(
  i = as.integer(factor(enc$Species, levels = unique(enc$Species))),
  j = enc$idx,
  x = 1,
  dims = c(n_distinct(enc$Species), nrow(mut_index)),
  dimnames = list(unique(enc$Species), mut_index$Orthologous_AA)
)

X_dense_clean <- as.matrix(X)
mode(X_dense_clean) <- "numeric"
X_dense_clean[X_dense_clean > 0] <- 1

# 覆盖保存
saveRDS(X_dense_clean, file.path(base_dir, "output", "X_dense.RDS"))
message("✅ Saved clean X_dense to output/X_dense.RDS")

# ========== ⬇️ 原有分析流程（不需改动） ==========

# 自定义函数：读取并清洗矩阵
sanitize_matrix <- function(mat_path) {
  mat <- readRDS(mat_path)
  if (!is.matrix(mat)) {
    mat <- as.matrix(mat)
  }
  mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  rownames(mat) <- rownames(readRDS(mat_path))
  return(mat)
}

for (variant in label_variants) {
  message("\n>>> Processing:", variant)
  fig_dir <- file.path(base_dir, "figures", variant)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  mat_path <- switch(variant,
                     original = file.path(base_dir, "output", "X_dense.RDS"),
                     residual = file.path(base_dir, "output", "X_resid.RDS"),
                     filtered = file.path(base_dir, "output", "X_dense_filtered.RDS"))

  label_file <- file.path(base_dir, "output", paste0("cluster_labels_", variant, ".csv"))
  X_dense <- sanitize_matrix(mat_path)
  cluster_labels <- read_csv(label_file, show_col_types = FALSE) %>%
    mutate(Unit = as.character(Unit))

  set.seed(123)
  umap_res <- uwot::umap(
    X_dense,
    metric = "manhattan",
    n_neighbors = min(15, nrow(X_dense)-1),
    min_dist = 0.3,
    verbose = TRUE
  )

  umap_df_base <- as_tibble(umap_res, .name_repair = "unique") %>%
    setNames(c("UMAP1", "UMAP2")) %>%
    mutate(Species = rownames(X_dense))

  top_methods <- c("EUCLIDEAN_KMEANS", "EUCLIDEAN_GMM", "COSINE_KMEANS")

  for (method in top_methods) {
    umap_df <- umap_df_base %>%
      left_join(cluster_labels %>% select(Species = Unit, Cluster = !!sym(method)), by = "Species") %>%
      mutate(Cluster = factor(Cluster))

    # 如果 cluster 数量少于 3，手动调整颜色数
    n_cluster <- length(unique(umap_df$Cluster))
    n_colors <- max(3, n_cluster)

    # UMAP plot
    p <- ggplot(umap_df, aes(UMAP1, UMAP2, colour = Cluster)) +
      geom_point(size = 3, alpha = .85) +
      scale_colour_brewer(palette = "Set1", na.translate = FALSE) +
      theme_classic(base_size = 14) +
      labs(title = paste("UMAP +", method, "clusters (", variant, ")"))

    ggsave(file.path(fig_dir, paste0("umap_", method, ".png")),
           p, width = 6, height = 5, dpi = 300)

    # Heatmap
    mut_freq <- colSums(X_dense)
    top30 <- names(sort(mut_freq, decreasing = TRUE))[1:min(30, length(mut_freq))]
    sub_mat <- X_dense[umap_df$Species, top30, drop = FALSE]

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
      width = unit(8, "cm"),
      height = unit(12, "cm")
    )

    png(file.path(fig_dir, paste0("heatmap_", method, ".png")), width = 1000, height = 1400, res = 120)
    draw(ht)
    dev.off()
  }

  # Upset plot 1：物种集合交集
  X_transpose <- t(X_dense)
  mutation_df <- as.data.frame(X_transpose)
  mutation_df$Mutation <- rownames(mutation_df)
  upset_df <- mutation_df %>% relocate(Mutation, .after = last_col())

  png(file.path(fig_dir, "mutation_species_upset.png"), width = 1600, height = 1000, res = 120)
  upset(upset_df[, -ncol(upset_df)], nsets = 10, nintersects = 30,
        keep.order = TRUE, sets.bar.color = "steelblue", order.by = "freq")
  dev.off()

  # Upset plot 2：突变在物种中的分布
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
}

message("✓ clusteringplot_plus_all.R completed")


label_variants <- c("original", "residual", "filtered")

# ========== ⬇️ 新增：重建 clean X_dense 并保存为 RDS ==========

message("🧬 Rebuilding clean X_dense from muts_completed ...")
mut_file <- file.path(base_dir, "output", "labmuts.csv")
muts_completed <- read_csv(mut_file, show_col_types = FALSE)

unit_tbl <- muts_completed %>%
  mutate(Species_clean = str_replace(Strain_ID, "_[^_]+$", ""),
         Species_clean = str_replace_all(Species_clean, "_", " ")) %>%
  select(Species = Species_clean, Orthologous_AA = AA_mut_name_Ecoli) %>%
  filter(!is.na(Species), !is.na(Orthologous_AA)) %>%
  distinct()

mut_index <- unit_tbl %>%
  distinct(Orthologous_AA) %>%
  arrange(Orthologous_AA) %>%
  mutate(idx = row_number())

enc <- inner_join(unit_tbl, mut_index, by = "Orthologous_AA")

X <- sparseMatrix(
  i = as.integer(factor(enc$Species, levels = unique(enc$Species))),
  j = enc$idx,
  x = 1,
  dims = c(n_distinct(enc$Species), nrow(mut_index)),
  dimnames = list(unique(enc$Species), mut_index$Orthologous_AA)
)

X_dense_clean <- as.matrix(X)
mode(X_dense_clean) <- "numeric"
X_dense_clean[X_dense_clean > 0] <- 1

# 覆盖保存
saveRDS(X_dense_clean, file.path(base_dir, "output", "X_dense.RDS"))
message("✅ Saved clean X_dense to output/X_dense.RDS")

# ========== ⬇️ 原有分析流程（不需改动） ==========

# 自定义函数：读取并清洗矩阵
sanitize_matrix <- function(mat_path) {
  mat <- readRDS(mat_path)
  if (!is.matrix(mat)) {
    mat <- as.matrix(mat)
  }
  mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  rownames(mat) <- rownames(readRDS(mat_path))
  return(mat)
}

for (variant in label_variants) {
  message("\n>>> Processing:", variant)
  fig_dir <- file.path(base_dir, "figures", variant)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  
  mat_path <- switch(variant,
                     original = file.path(base_dir, "output", "X_dense.RDS"),
                     residual = file.path(base_dir, "output", "X_resid.RDS"),
                     filtered = file.path(base_dir, "output", "X_dense_filtered.RDS"))
  
  label_file <- file.path(base_dir, "output", paste0("cluster_labels_", variant, ".csv"))
  X_dense <- sanitize_matrix(mat_path)
  cluster_labels <- read_csv(label_file, show_col_types = FALSE) %>%
    mutate(Unit = as.character(Unit))
  
  set.seed(123)
  umap_res <- uwot::umap(
    X_dense,
    metric = "manhattan",
    n_neighbors = min(15, nrow(X_dense)-1),
    min_dist = 0.3,
    verbose = TRUE
  )
  
  umap_df_base <- as_tibble(umap_res, .name_repair = "unique") %>%
    setNames(c("UMAP1", "UMAP2")) %>%
    mutate(Species = rownames(X_dense))
  
  top_methods <- c("EUCLIDEAN_KMEANS", "EUCLIDEAN_GMM", "COSINE_KMEANS")
  
  for (method in top_methods) {
    umap_df <- umap_df_base %>%
      left_join(cluster_labels %>% select(Species = Unit, Cluster = !!sym(method)), by = "Species") %>%
      mutate(Cluster = factor(Cluster))
    
    # 如果 cluster 数量少于 3，手动调整颜色数
    n_cluster <- length(unique(umap_df$Cluster))
    n_colors <- max(3, n_cluster)
    
    # UMAP plot
    p <- ggplot(umap_df, aes(UMAP1, UMAP2, colour = Cluster)) +
      geom_point(size = 3, alpha = .85) +
      scale_colour_brewer(palette = "Set1", na.translate = FALSE) +
      theme_classic(base_size = 14) +
      labs(title = paste("UMAP +", method, "clusters (", variant, ")"))
    
    ggsave(file.path(fig_dir, paste0("umap_", method, ".png")),
           p, width = 6, height = 5, dpi = 300)
    
    # Heatmap
    mut_freq <- colSums(X_dense)
    top30 <- names(sort(mut_freq, decreasing = TRUE))[1:min(30, length(mut_freq))]
    sub_mat <- X_dense[umap_df$Species, top30, drop = FALSE]
    
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
      width = unit(8, "cm"),
      height = unit(12, "cm")
    )
    
    png(file.path(fig_dir, paste0("heatmap_", method, ".png")), width = 1000, height = 1400, res = 120)
    draw(ht)
    dev.off()
  }
  
  # Upset plot 1：物种集合交集
  X_transpose <- t(X_dense)
  mutation_df <- as.data.frame(X_transpose)
  mutation_df$Mutation <- rownames(mutation_df)
  upset_df <- mutation_df %>% relocate(Mutation, .after = last_col())
  
  png(file.path(fig_dir, "mutation_species_upset.png"), width = 1600, height = 1000, res = 120)
  upset(upset_df[, -ncol(upset_df)], nsets = 10, nintersects = 30,
        keep.order = TRUE, sets.bar.color = "steelblue", order.by = "freq")
  dev.off()
  
  # Upset plot 2：突变在物种中的分布
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
}

message("✓ clusteringplot_plus_all.R completed")


