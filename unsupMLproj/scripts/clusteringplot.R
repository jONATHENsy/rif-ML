# =============================================================
# clusteringplot.R – UMAP + multi-method cluster + heatmap + upset plot
# -------------------------------------------------------------
# Can be sourced by main.R or run independently
# =============================================================

message("\n>>> clusteringplot.R started")

## -----------------------------------------------------------
## 1. Dependencies
## -----------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(uwot)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  
})

## -----------------------------------------------------------
## 2. Locate input objects / files
## -----------------------------------------------------------
if (!exists("base_dir"))
  base_dir <- "C:/Users/user/Desktop/D Drive/2025s1/BIOX7011/rif-ML/unsupMLproj"

mut_file <- file.path(base_dir, "output", "labmuts.csv")
label_file <- file.path(base_dir, "output", "cluster_labels_all_methods.csv")
fig_dir <- file.path(base_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (!exists("muts_completed"))
  muts_completed <- read_csv(mut_file, show_col_types = FALSE)

## -----------------------------------------------------------
## 3. Create binary mutation matrix based on Species only
## -----------------------------------------------------------
if (!exists("X_dense")) {
  unit_tbl <- muts_completed %>%
    mutate(Species_clean = str_replace(Strain_ID, "_[^_]+$", ""),         # 去掉菌株
           Species_clean = str_replace_all(Species_clean, "_", " ")) %>%  # 恢复空格
    select(Species = Species_clean, Orthologous_AA = AA_mut_name_Ecoli) %>%
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
  
  X_dense <- as.matrix(X)
  X_dense[X_dense > 0] <- 1
}

## -----------------------------------------------------------
## 4. UMAP projection (shared for all methods)
## -----------------------------------------------------------
set.seed(123)
umap_res <- uwot::umap(
  X_dense,
  metric = "manhattan",
  n_neighbors = 15,
  min_dist = 0.3,
  verbose = TRUE
)

umap_df_base <- as_tibble(umap_res, .name_repair = "unique") %>%
  setNames(c("UMAP1", "UMAP2")) %>%
  mutate(Species = rownames(X_dense))

## -----------------------------------------------------------
## 5. Generate UMAP and heatmaps for top 3 methods
## -----------------------------------------------------------
top_methods <- c("EUCLIDEAN_KMEANS", "EUCLIDEAN_GMM", "COSINE_KMEANS")
cluster_labels <- read_csv(label_file, show_col_types = FALSE)

for (method in top_methods) {
  umap_df <- umap_df_base %>%
    left_join(cluster_labels %>% select(Species = Unit, Cluster = !!sym(method)), by = "Species") %>%
    mutate(Cluster = factor(Cluster))
  
  ## Save UMAP plot
  p <- ggplot(umap_df, aes(UMAP1, UMAP2, colour = Cluster)) +
    geom_point(size = 3, alpha = .85) +
    scale_colour_brewer(palette = "Set1", na.translate = FALSE) +
    theme_classic(base_size = 14) +
    labs(title = paste("UMAP +", method, "clusters"))
  
  ggsave(file.path(fig_dir, paste0("umap_", method, ".png")),
         p, width = 6, height = 5, dpi = 300)
  
  ## Save mutation heatmap
  mut_freq <- colSums(X_dense)
  top30 <- names(sort(mut_freq, decreasing = TRUE))[1:30]
  sub_mat <- X_dense[umap_df$Species, top30, drop = FALSE]
  
  cluster_ids <- sort(unique(umap_df$Cluster))
  cluster_cols <- setNames(brewer.pal(length(cluster_ids), "Set1")[seq_along(cluster_ids)], cluster_ids)
  
  row_anno <- rowAnnotation(
    Cluster = umap_df$Cluster,
    col = list(Cluster = cluster_cols),
    show_annotation_name = FALSE
  )
  
  label_gp <- gpar(fontface = ifelse(umap_df$Species == "Escherichia coli", "bold", "plain"),
                   fontsize = 6)
  
  ht <- Heatmap(
    sub_mat,
    name = "Mut",
    col = c("0" = "white", "1" = "steelblue"),
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    row_names_gp = label_gp,
    left_annotation = row_anno,
    row_split = umap_df$Cluster,
    width = unit(8, "cm"),
    height = unit(12, "cm")
  )
  
  png(file.path(fig_dir, paste0("heatmap_", method, ".png")), width = 1000, height = 1400, res = 120)
  draw(ht)
  dev.off()
}

## -----------------------------------------------------------
## 6. Upset Plot of mutation presence across species
##     Using classic and reliable UpSetR package
## -----------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(UpSetR)  # ✅ 轻量经典方案
})

# 转置：行 = mutation，列 = species
X_transpose <- t(X_dense)  # 原始行是物种，现在每行为一个突变
mutation_df <- as.data.frame(X_transpose)
mutation_df$Mutation <- rownames(mutation_df)

# 将 Mutation 列移到最后（UpSetR 要求 binary 列在前）
upset_df <- mutation_df %>%
  relocate(Mutation, .after = last_col())

# 画图（注意 UpSetR 使用 pdf/png 需在 dev.off() 之后才能看到）
png(file.path(fig_dir, "mutation_species_upset.png"), width = 1600, height = 1000, res = 120)

upset(
  upset_df[, -ncol(upset_df)],   # 去掉最后的 Mutation 列，只保留 0/1 列
  nsets = 10,                    # 显示前10个物种集合
  nintersects = 30,             # 显示前30个交集
  keep.order = TRUE,
  sets.bar.color = "steelblue",
  order.by = "freq"
)

dev.off()

## -----------------------------------------------------------
## 7. Upset Plot of top mutations across species
##     (inverted: sets = mutations, elements = species)
## -----------------------------------------------------------

# 统计每个突变在多少物种中出现
mut_freq <- sort(colSums(X_dense), decreasing = TRUE)

# 选取前10个突变
top_mutations <- names(mut_freq)[1:10]

# 创建一个新矩阵：行 = 物种，列 = 这10个突变
submat <- X_dense[, top_mutations]

# 转置：列为突变（sets），行是物种（元素）
inverted_df <- as.data.frame(submat)
inverted_df$Species <- rownames(submat)

# 将物种放到最后一列
inverted_df <- inverted_df %>%
  relocate(Species, .after = last_col())

# 保存图像
png(file.path(fig_dir, "top10mut_species_upset.png"), width = 1600, height = 1000, res = 120)

upset(
  inverted_df[, -ncol(inverted_df)],   # 只保留突变列作为 sets
  nsets = 10,                          # 前10个突变
  nintersects = 30,                   # 最多显示30个交集组合
  keep.order = TRUE,
  sets.bar.color = "tomato",
  order.by = "freq",
  mainbar.y.label = "Number of species",
  sets.x.label = "Top mutations"
)

dev.off()

