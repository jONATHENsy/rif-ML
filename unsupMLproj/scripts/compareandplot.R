# =============================================================
# compareandplot.R – 自动选最优聚类并作图（全突变）
# - 先在当前 X_dense 上对多种 metric×算法 做聚类 → 计算平均轮廓系数
# - 选出 top_methods（默认限制为 metric="manhattan"，更贴合你之前习惯）
# - 用各自的 UMAP 与聚类标签画 UMAP & Heatmap
# - 已修复：整体向右“平移”保证左侧物种名不被遮挡
# =============================================================

message("\n🔍 clusteringplot_filtered_full.R started")

suppressPackageStartupMessages({
  library(tidyverse)   # includes purrr/dplyr/readr/ggplot2
  library(Matrix)
  library(uwot)
  library(ComplexHeatmap)
  library(circlize)
  library(magick)
  library(RColorBrewer)
  library(UpSetR)
  library(grid)
  # for clustering + evaluation
  library(dbscan)
  library(cluster)
  library(proxy)
  library(mclust)
})

# ---------------- 参数区（可改） ----------------
# 仅从这个 metric 里挑前K名；设为 NULL 则在所有 metric 中挑
restrict_metric <- NULL   # 常用："manhattan" / NULL
top_k          <- 3              # 选前K个方案
min_cluster_pts <- 5             # HDBSCAN/DBSCAN 的 minPts
dbscan_eps      <- 1.2           # DBSCAN 的 eps
set.seed(123)

# 额外向右平移（mm），避免左侧行名被遮挡
shift_right_mm <- 28  # 24~40 之间调；不改其它配置

# ---------------- 路径 ----------------
base_dir <- "C:/Users/user/Desktop/D Drive/2025s1/BIOX7011/rif-ML/unsupMLproj"
fig_dir  <- file.path(base_dir, "figures", "filtered_h02")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------- 数据 ----------------
mat_path   <- file.path(base_dir, "output", "X_dense_high.RDS")
# 可选：如果有外部聚类标签文件，这里给出路径；没有也没关系
label_file <- file.path(base_dir, "output", "cluster_labels_filtered.csv")

# 👇 读取并过滤物种
remove_species <- c(
  "Vibrio parahaemolyticus",
  "Vibrio vulnificus",
  "Streptomyces lividans"
)

X_dense <- readRDS(mat_path)
mode(X_dense) <- "numeric"
X_dense[is.na(X_dense)] <- 0

# 过滤前记录哪些物种会被移除（仅用于日志，可删）
to_remove_now <- intersect(remove_species, rownames(X_dense))

# 过滤矩阵行
X_dense <- X_dense[!(rownames(X_dense) %in% remove_species), , drop = FALSE]

# 过滤后更新物种名（供后续使用）
unit_names <- rownames(X_dense)

# 若存在外部聚类标签文件，则同步过滤（可选）
if (file.exists(label_file)) {
  cluster_labels <- readr::read_csv(label_file, show_col_types = FALSE) %>%
    dplyr::mutate(Unit = as.character(Unit)) %>%
    dplyr::filter(!(Unit %in% remove_species))
}

if (length(to_remove_now) > 0) {
  message("Filtered species: ", paste(to_remove_now, collapse = ", "))
}

# =============================================================
# 一、compareclustering（内联版）
# =============================================================

run_pipeline <- function(mat, metric = "euclidean",
                         cluster_method = "hdbscan",
                         minPts = 5, eps = 1.2, seed = 123) {
  set.seed(seed)
  
  if (metric == "jaccard") {
    dmat <- proxy::dist(mat, method = "Jaccard")
    um   <- uwot::umap(as.matrix(dmat), input = "dist", metric = "precomputed",
                       n_neighbors = 15, min_dist = 0.3, verbose = FALSE)
  } else {
    um   <- uwot::umap(mat, metric = metric, n_neighbors = 15,
                       min_dist = 0.3, verbose = FALSE)
  }
  
  clust <- switch(cluster_method,
                  "hdbscan" = dbscan::hdbscan(as.data.frame(um), minPts = minPts)$cluster,
                  "kmeans"  = kmeans(um, centers = 4)$cluster,
                  "dbscan"  = dbscan::dbscan(um, eps = eps, minPts = minPts)$cluster,
                  "gmm"     = mclust::Mclust(um)$classification,
                  stop("Unknown cluster_method: ", cluster_method)
  )
  
  list(umap = um, cluster = clust)
}

schemes <- expand.grid(
  metric = c("euclidean", "manhattan", "cosine"),
  method = c("hdbscan", "kmeans", "dbscan", "gmm"),
  stringsAsFactors = FALSE
) %>% mutate(name = paste(toupper(metric), toupper(method), sep = "_"))

# 跑全部方案
results <- purrr::pmap(schemes, ~run_pipeline(
  X_dense, metric = ..1, cluster_method = ..2,
  minPts = min_cluster_pts, eps = dbscan_eps, seed = 123
))

# 平均轮廓系数（只在有效聚类上算）
mean_sil <- purrr::map2_dbl(results, schemes$name, function(res, nm) {
  cl <- res$cluster
  um <- res$umap
  valid_idx <- which(!is.na(cl) & cl > 0)
  if (length(valid_idx) < 2 || length(unique(cl[valid_idx])) < 2) return(NA_real_)
  sil <- cluster::silhouette(cl[valid_idx], dist(um[valid_idx, ]))
  mean(sil[, 3])
})

sil_df <- tibble(
  Scheme = schemes$name,
  Metric = toupper(schemes$metric),
  Method = toupper(schemes$method),
  Mean_Silhouette = round(mean_sil, 4),
  Valid = !is.na(mean_sil)
)

# 选择 top_methods
sel_df <- sil_df %>% filter(Valid)
if (!is.null(restrict_metric)) {
  sel_df <- sel_df %>% filter(Metric == toupper(restrict_metric))
}

# 若过滤后没有有效方案，回退到全体有效方案
if (nrow(sel_df) == 0) {
  warning("No valid schemes under restrict_metric = ", restrict_metric,
          ". Falling back to all metrics.")
  sel_df <- sil_df %>% filter(Valid)
}

# slice_head() 的 n 必须是常量
n_take <- min(as.integer(top_k), nrow(sel_df))

top_methods <- sel_df %>%
  slice_max(order_by = Mean_Silhouette, n = n_take, with_ties = FALSE) %>%
  pull(Scheme)

message("⭐ Top methods: ", paste(top_methods, collapse = ", "))


# 把所有方案的聚类标签放成一张表（备用/记录）
cluster_tbl <- tibble(Unit = unit_names)
for (i in seq_len(nrow(schemes))) {
  cluster_tbl[[schemes$name[i]]] <- results[[i]]$cluster
}

# 存档（可选）
out_dir <- file.path(base_dir, "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(sil_df,      file.path(out_dir, "silhouette_all_methods_auto.csv"))
write_csv(cluster_tbl, file.path(out_dir, "cluster_labels_all_methods_auto.csv"))

# =============================================================
# 二、按 top_methods 作图
# =============================================================

for (method in top_methods) {
  i <- which(schemes$name == method)
  res <- results[[i]]
  
  umap_df <- as_tibble(res$umap, .name_repair = "unique") %>%
    setNames(c("UMAP1", "UMAP2")) %>%
    mutate(Species = unit_names,
           Cluster = factor(res$cluster))
  
  # ---------- UMAP ----------
  p <- ggplot(umap_df, aes(UMAP1, UMAP2, colour = Cluster)) +
    geom_point(size = 3, alpha = .85) +
    scale_colour_brewer(palette = "Set1", na.translate = FALSE) +
    theme_classic(base_size = 14) +
    labs(title = paste("UMAP +", method, "clusters (filtered, all mutations)"))
  
  ggsave(file.path(fig_dir, paste0("umap_", method, ".png")),
         p, width = 6, height = 5, dpi = 300)
  
  # ---------- Heatmap ----------
  sub_mat <- X_dense[umap_df$Species, , drop = FALSE]
  
  # 颜色：若包含“0”（噪声），给一个中性灰
  cluster_ids <- sort(unique(umap_df$Cluster))
  pal <- RColorBrewer::brewer.pal(max(3, length(cluster_ids)), "Set1")
  cluster_cols <- setNames(pal[seq_along(cluster_ids)], cluster_ids)
  if ("0" %in% names(cluster_cols)) cluster_cols[["0"]] <- "grey70"
  
  row_anno <- ComplexHeatmap::rowAnnotation(
    Cluster = umap_df$Cluster,
    col = list(Cluster = cluster_cols),
    show_annotation_name = FALSE,
    width = unit(4, "mm")    # 小幅加粗注释带
  )
  
  n_mut <- ncol(sub_mat)
  n_sp  <- nrow(sub_mat)
  
  # ====== 画布尺寸（更宽）======
  col_px_full <- 32  # 👈 全量热图：每列像素（原来是 18，太瘦）
  row_px_full <- 26  # 行高像素（原来 50 很高）
  png_w <- max(3600, n_mut * col_px_full)
  png_h <- max(1600, n_sp  * row_px_full)
  dpi   <- 300
  # ============================
  
  rowlab_w <- ComplexHeatmap::max_text_width(rownames(sub_mat), gp = gpar(fontsize = 10))
  collab_h <- ComplexHeatmap::max_text_width(colnames(sub_mat), gp = gpar(fontsize = 8))
  row_names_max <- rowlab_w + unit(12, "mm")
  
  ht <- ComplexHeatmap::Heatmap(
    sub_mat,
    name = "Mut",
    col = c("0" = "white", "1" = "steelblue"),
    cluster_rows = FALSE,
    cluster_columns = TRUE,        # 保持列聚类（得到列树）
    row_split = umap_df$Cluster,   # 行分面（用聚类标签）
    show_row_names = TRUE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 10),
    row_names_max_width = row_names_max,
    show_column_names = TRUE,
    column_names_rot  = 90,
    column_names_gp   = gpar(fontsize = 8),
    column_names_centered = TRUE,
    column_names_max_height = collab_h + unit(2, "mm"),
    right_annotation = row_anno,
    width  = unit(1, "npc"),
    height = unit(1, "npc"),
    use_raster = TRUE,
    raster_device = "png"
  )
  
  left_pad_mm   <- max(convertWidth(row_names_max, "mm", valueOnly = TRUE) + 10, 34)
  bottom_pad_mm <- convertWidth(collab_h, "mm", valueOnly = TRUE) + 10
  
  # --- PNG（全量，更宽）---
  png(file.path(fig_dir, paste0("heatmap_", method, "_allmut.png")),
      width = png_w, height = png_h, res = dpi, type = "cairo-png")
  ht_drawn <- ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    padding = unit(c(6, 6, bottom_pad_mm, left_pad_mm + shift_right_mm), "mm")
  )
  dev.off()
  
  # --- PDF（矢量）---
  # 按列/行像素估算英寸：1 英寸≈90px 的经验换算
  pdf_w_in <- max(8,  n_mut * (col_px_full/90))
  pdf_h_in <- max(6,  n_sp  * (row_px_full/90))
  out_pdf <- file.path(fig_dir, sprintf("heatmap_%s_allmut.pdf", method))
  if (file.exists(out_pdf)) {
    out_pdf <- file.path(fig_dir, sprintf("heatmap_%s_allmut_%s.pdf",
                                          method, format(Sys.time(), "%Y%m%d_%H%M%S")))
  }
  pdf(out_pdf, width = pdf_w_in, height = pdf_h_in)
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    padding = unit(c(6, 6, bottom_pad_mm, left_pad_mm + shift_right_mm), "mm")
  )
  dev.off()
  message("PDF saved: ", out_pdf)
  
  # ---------- Top-30-only clear heatmap (fixed width & padding) ----------
  ord_try <- try(ComplexHeatmap::column_order(ht_drawn), silent = TRUE)
  if (!inherits(ord_try, "try-error")) {
    ord_vec <- if (is.list(ord_try)) ord_try[[1]] else ord_try
    ordered_names <- colnames(sub_mat)[ord_vec]
    keep_cols <- ordered_names[ordered_names %in% top30]
  } else {
    keep_cols <- intersect(colnames(sub_mat), top30)
  }
  
  if (length(keep_cols) == 0) {
    message("No overlap between matrix columns and top30 for ", method, "; skip top30 plot.")
  } else {
    sub_mat_top <- sub_mat[, keep_cols, drop = FALSE]
    
    # —— 关键：用“像素”控制列宽/行高，并把左/右留白精确转为像素 —— #
    n_top <- ncol(sub_mat_top)
    n_sp  <- nrow(sub_mat_top)
    dpi   <- 300
    
    col_px <- 70     # 每一列的目标像素宽度（更宽更清爽；可调 60–90）
    row_px <- 26     # 每一行的目标像素高度（可调 22–30）
    
    left_pad_mm  <- 6    # top30 不再额外右移，避免重复留白
    right_pad_mm <- 45   # 为右侧图例和注释条留足空间（很重要）
    top_pad_mm   <- 6
    
    # 把 mm 换算为像素
    mm2px <- function(mm, dpi) mm / 25.4 * dpi
    left_pad_px   <- mm2px(left_pad_mm,  dpi)
    right_pad_px  <- mm2px(right_pad_mm, dpi)
    top_pad_px    <- mm2px(top_pad_mm,   dpi)
    
    # 画布宽高（像素）：= 左留白 + 列宽*n_top + 右留白
    png_w_t <- ceiling(left_pad_px + n_top * col_px + right_pad_px)
    # 高度：给行名 + 热图 + (旋转列名占的“高度”)
    # 列名高度由 ComplexHeatmap 处理，这里主要保证行方向有足够空间
    png_h_t <- max(1200, ceiling(n_sp * row_px + top_pad_px + 80))
    
    # 文本样式与最大宽/高（不过度放大）
    rowlab_w_t <- ComplexHeatmap::max_text_width(rownames(sub_mat_top),
                                                 gp = gpar(fontsize = 10))
    collab_h_t <- unit(12, "mm")  # 给列名留 12mm 高度（旋转 90°）
    
    ht_top <- ComplexHeatmap::Heatmap(
      sub_mat_top,
      name = "Mut",
      col = c("0" = "white", "1" = "steelblue"),
      cluster_rows = FALSE,
      cluster_columns = FALSE,      # 固定使用完整图里的列顺序
      row_split = umap_df$Cluster,
      
      show_row_names = TRUE,
      row_names_side = "left",
      row_names_gp   = gpar(fontsize = 10),
      row_names_max_width = rowlab_w_t,
      
      show_column_names = TRUE,
      column_names_rot  = 90,
      column_names_gp   = gpar(fontsize = 9),
      column_names_centered = TRUE,
      column_names_max_height = collab_h_t,
      
      right_annotation = row_anno,   # 复用完整图的行注释
      width  = unit(1, "npc"),
      height = unit(1, "npc"),
      use_raster = TRUE,
      raster_device = "png"
    )
    
    # —— 真正绘图：右侧给足 padding，避免图例挤到热图上 —— #
    png(file.path(fig_dir, sprintf("heatmap_%s_top30.png", method)),
        width = png_w_t, height = png_h_t, res = dpi, type = "cairo-png")
    ComplexHeatmap::draw(
      ht_top,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      # padding 顺序：top, right, bottom, left
      padding = unit(c(top_pad_mm, right_pad_mm, 8, left_pad_mm), "mm"),
      merge_legend = TRUE
    )
    dev.off()
    
    # 同步生成 PDF（宽高与 PNG 对齐换算）
    pdf(file.path(fig_dir, sprintf("heatmap_%s_top30.pdf", method)),
        width  = png_w_t / dpi,
        height = png_h_t / dpi)
    ComplexHeatmap::draw(
      ht_top,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      padding = unit(c(top_pad_mm, right_pad_mm, 8, left_pad_mm), "mm"),
      merge_legend = TRUE
    )
    dev.off()
    
    message("Top-30 clear heatmap saved for ", method,
            "  [", n_top, " cols, ", round(png_w_t/dpi,2), "×", round(png_h_t/dpi,2), " in].")
  }
  # ---------- end of Top-30-only clear heatmap ----------
  
  
  
}

# =============================================================
# 三、UpSet（与聚类无关，保持原逻辑）
# =============================================================
X_transpose <- t(X_dense)
mutation_df <- as.data.frame(X_transpose)
mutation_df$Mutation <- rownames(mutation_df)
upset_df <- mutation_df %>% relocate(Mutation, .after = dplyr::last_col())

png(file.path(fig_dir, "mutation_species_upset.png"), width = 1600, height = 1000, res = 120)
upset(upset_df[, -ncol(upset_df)], nsets = 10, nintersects = 30,
      keep.order = TRUE, sets.bar.color = "steelblue", order.by = "freq")
dev.off()

mut_freq <- sort(colSums(X_dense), decreasing = TRUE)
top_mutations <- names(mut_freq)[1:min(10, length(mut_freq))]
submat <- X_dense[, top_mutations, drop = FALSE]
inverted_df <- as.data.frame(submat)
inverted_df$Species <- rownames(submat)
inverted_df <- inverted_df %>% relocate(Species, .after = dplyr::last_col())

png(file.path(fig_dir, "top10mut_species_upset.png"), width = 1600, height = 1000, res = 120)
upset(inverted_df[, -ncol(inverted_df)], nsets = 10, nintersects = 30,
      keep.order = TRUE, sets.bar.color = "tomato", order.by = "freq",
      mainbar.y.label = "Number of species", sets.x.label = "Top mutations")
dev.off()

message("✅ Done – auto-selected top methods: ", paste(top_methods, collapse = ", "))
