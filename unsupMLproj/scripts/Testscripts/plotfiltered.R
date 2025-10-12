# =============================================================
# clusteringplot_filtered_full.R – 使用 Manhattan 最佳三法，绘制所有突变的完整图
# 关键改动（已全局应用）：
# 1) 行注释移到右侧（right_annotation），释放左侧空间
# 2) 动态测量行名宽度 + 追加裕量，保证 y 轴物种名称不被截断
# 3) 画布尺寸/留白按行列动态扩展；导出高分辨率 PNG 与矢量 PDF
# =============================================================

message("\n🔍 clusteringplot_filtered_full.R started")

suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(uwot)
  library(ComplexHeatmap)
  library(circlize)
  library(magick)
  library(RColorBrewer)
  library(UpSetR)
  library(grid)
})

# ---------------- 路径 ----------------
base_dir <- "C:/Users/user/Desktop/D Drive/2025s1/BIOX7011/rif-ML/unsupMLproj"
fig_dir  <- file.path(base_dir, "figures", "filtered_mhtest02")  # 或改成 filtered_full
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------- 数据 ----------------
mat_path   <- file.path(base_dir, "output", "X_dense_midhigh.RDS")
label_file <- file.path(base_dir, "output", "cluster_labels_filtered.csv")

X_dense <- readRDS(mat_path)
mode(X_dense) <- "numeric"
X_dense[is.na(X_dense)] <- 0

cluster_labels <- read_csv(label_file, show_col_types = FALSE) %>%
  mutate(Unit = as.character(Unit))

# ---------------- UMAP ----------------
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

# 使用 Manhattan 下评分较高的三法
top_methods <- c("MANHATTAN_HDBSCAN", "MANHATTAN_GMM", "MANHATTAN_KMEANS")

for (method in top_methods) {
  umap_df <- umap_df_base %>%
    left_join(cluster_labels %>% select(Species = Unit, Cluster = !!sym(method)), by = "Species") %>%
    mutate(Cluster = factor(Cluster))
  
  # ---------- UMAP 图 ----------
  p <- ggplot(umap_df, aes(UMAP1, UMAP2, colour = Cluster)) +
    geom_point(size = 3, alpha = .85) +
    scale_colour_brewer(palette = "Set1", na.translate = FALSE) +
    theme_classic(base_size = 14) +
    labs(title = paste("UMAP +", method, "clusters (filtered, all mutations)"))
  
  ggsave(file.path(fig_dir, paste0("umap_", method, ".png")),
         p, width = 6, height = 5, dpi = 300)
  
  # ---------- Heatmap：全部突变 ----------
  sub_mat <- X_dense[umap_df$Species, , drop = FALSE]
  
  cluster_ids <- sort(unique(umap_df$Cluster))
  palette_colors <- RColorBrewer::brewer.pal(max(3, length(cluster_ids)), "Set1")
  cluster_cols <- stats::setNames(palette_colors[seq_along(cluster_ids)], cluster_ids)
  
  # 行注释放到右侧，避免压缩左侧行名空间
  row_anno <- ComplexHeatmap::rowAnnotation(
    Cluster = umap_df$Cluster,
    col = list(Cluster = cluster_cols),
    show_annotation_name = FALSE,
    width = unit(3, "mm")
  )
  
  # ---- 动态设备尺寸与留白 ----
  n_mut <- ncol(sub_mat)
  n_sp  <- nrow(sub_mat)
  
  # 画布尺寸（像素）
  png_w <- max(3200, n_mut * 18)      # 每列 ~18 px；列多时可调到 20~24
  png_h <- max(2400, n_sp * 50)       # 每物种 ~50 px，确保 y 轴可读
  dpi   <- 300
  
  # 行名/列名所需最大宽度
  rowlab_w <- ComplexHeatmap::max_text_width(rownames(sub_mat), gp = gpar(fontsize = 10))
  collab_h <- ComplexHeatmap::max_text_width(colnames(sub_mat), gp = gpar(fontsize = 6))
  
  # 给行名宽度加足裕量（6mm，可按需调大）
  row_names_max <- rowlab_w + unit(6, "mm")
  
  # Heatmap 对象
  ht <- ComplexHeatmap::Heatmap(
    sub_mat,
    name = "Mut",
    col = c("0" = "white", "1" = "steelblue"),
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    row_split = umap_df$Cluster,
    
    # ✅ 行名显示
    show_row_names = TRUE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 10),
    row_names_max_width = row_names_max,
    
    # ✅ 列名竖排
    show_column_names = TRUE,
    column_names_rot  = 90,
    column_names_gp   = gpar(fontsize = 6),
    column_names_centered = TRUE,
    column_names_max_height = collab_h + unit(2, "mm"),
    
    # ✅ 注释改到右侧
    right_annotation = row_anno,
    
    width  = unit(1, "npc"),
    height = unit(1, "npc"),
    
    use_raster = TRUE,
    raster_device = "png"
  )
  
  # 计算四边 padding（mm）——左侧至少 32mm，保证行名不截断
  left_pad_mm   <- max(convertWidth(row_names_max, "mm", valueOnly = TRUE) + 10, 32)
  bottom_pad_mm <- convertWidth(collab_h, "mm", valueOnly = TRUE) + 8
  shift_right_mm <- 28
  
  # --- 导出 PNG ---
  png(file.path(fig_dir, paste0("heatmap_", method, "_allmut.png")),
      width = png_w, height = png_h, res = dpi, type = "cairo-png")
  
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    padding = unit(c(6, 6, bottom_pad_mm, left_pad_mm), "mm") # top,right,bottom,left
  )
  dev.off()
  
  # --- 导出 PDF（矢量，便于无限放大） ---
  pdf_w_in <- max(8, n_mut * 0.18)
  pdf_h_in <- max(8, n_sp  * 0.45)
  pdf(file.path(fig_dir, paste0("heatmap_", method, "_allmut.pdf")),
      width = pdf_w_in, height = pdf_h_in)
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    padding = unit(c(6, 6, bottom_pad_mm, left_pad_mm), "mm")
  )
  dev.off()
}

# ---------------- UpSet 1：物种集合交集 ----------------
X_transpose <- t(X_dense)
mutation_df <- as.data.frame(X_transpose)
mutation_df$Mutation <- rownames(mutation_df)
upset_df <- mutation_df %>% relocate(Mutation, .after = last_col())

png(file.path(fig_dir, "mutation_species_upset.png"), width = 1600, height = 1000, res = 120)
upset(upset_df[, -ncol(upset_df)], nsets = 10, nintersects = 30,
      keep.order = TRUE, sets.bar.color = "steelblue", order.by = "freq")
dev.off()

# ---------------- UpSet 2：top10 突变的物种分布 ----------------
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

    
  