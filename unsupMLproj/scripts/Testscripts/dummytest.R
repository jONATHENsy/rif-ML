# =============================================================
# validate_dummy_clustering.R – 在人工构建数据上测试聚类方法
# =============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(uwot)
  library(cluster)
  library(RColorBrewer)
  library(circlize)
})

message("\n🧪 Dummy clustering validation started")

# -------------------------------------------------------------
# 1. 构建 dummy binary matrix
# -------------------------------------------------------------
dummy_data <- tibble::tibble(
  Species = c("Species_1", "Species_2", "Species_3", "Species_4", "Species_5"),
  Mutations = c("A,B,C,D,E,Z", "M,N,O,P,Q,Z", "U,V,W,X,Y,Z", "M,N,O,R,S,Z", "A,B,C,F,G,Z")
)

# 提取所有突变（合并后唯一）
all_muts <- dummy_data %>%
  pull(Mutations) %>%
  strsplit(",") %>%
  unlist() %>%
  unique() %>%
  sort()

# 构建 binary 矩阵
X_dummy <- matrix(0, nrow = nrow(dummy_data), ncol = length(all_muts),
                  dimnames = list(dummy_data$Species, all_muts))

for (i in 1:nrow(dummy_data)) {
  muts <- strsplit(dummy_data$Mutations[i], ",")[[1]]
  X_dummy[i, muts] <- 1
}

X_dummy <- as.matrix(X_dummy)

# -------------------------------------------------------------
# 2. 运行 KMeans 聚类（聚成 2 类）
# -------------------------------------------------------------
set.seed(123)
k <- 2
clustering <- kmeans(X_dummy, centers = k)
cluster_labels <- as.factor(clustering$cluster)

# -------------------------------------------------------------
# 3. UMAP 可视化
# -------------------------------------------------------------
umap_res <- umap(X_dummy, n_neighbors = 3, min_dist = 0.1, metric = "manhattan")

umap_df <- as_tibble(umap_res, .name_repair = "unique") %>%
  setNames(c("UMAP1", "UMAP2")) %>%
  mutate(Species = rownames(X_dummy),
         Cluster = cluster_labels)

ggplot(umap_df, aes(UMAP1, UMAP2, color = Cluster, label = Species)) +
  geom_point(size = 4) +
  geom_text(nudge_y = 0.1, size = 4) +
  theme_minimal(base_size = 14) +
  scale_color_brewer(palette = "Set1") +
  ggtitle("UMAP Projection of Dummy Dataset with KMeans Clustering")

# -------------------------------------------------------------
# 4. 绘制 heatmap
# -------------------------------------------------------------
n_cluster <- length(unique(cluster_labels))
palette_colors <- brewer.pal(max(3, n_cluster), "Set1")[1:n_cluster]
cluster_levels <- levels(cluster_labels)
names(palette_colors) <- cluster_levels

row_anno <- rowAnnotation(
  Cluster = cluster_labels,
  col = list(Cluster = palette_colors),
  show_annotation_name = FALSE
)

ht <- Heatmap(
  X_dummy,
  name = "Presence",
  col = c("0" = "white", "1" = "steelblue"),
  cluster_columns = TRUE,
  cluster_rows = FALSE,
  left_annotation = row_anno,
  row_names_side = "left"
)

draw(ht)

message("✅ Done.")

