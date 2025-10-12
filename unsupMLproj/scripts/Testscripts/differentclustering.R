# ============================================================
#  EXTRA: run UMAP + HDBSCAN with multiple distance measures
# ============================================================

library(tidyverse)
library(uwot)
library(dbscan)
library(proxy)        # for Jaccard distance
library(cluster)      # silhouette

## ---------- 函数：单次运行 ----------------------------------
run_pipeline <- function(X_dense, metric = "euclidean",
                         minPts = 5, seed = 123) {
  set.seed(seed)
  
  if (metric == "jaccard") {
    jac <- proxy::dist(X_dense, method = "Jaccard")
    um  <- uwot::umap(jac,          # ← 指明 uwot
                      input = "dist",
                      n_neighbors = 15,
                      min_dist    = 0.3,
                      verbose     = FALSE)
  } else {
    um  <- uwot::umap(X_dense,
                      metric       = metric,
                      n_neighbors  = 15,
                      min_dist     = 0.3,
                      verbose      = FALSE)
  }
  
  hdb <- dbscan::hdbscan(as.data.frame(um), minPts = minPts)
  list(umap = um, cluster = hdb$cluster)
}

## ---------- 1. 需要评估的方案列表 ---------------------------
schemes <- tribble(
  ~name,       ~metric,
  "EUCLIDEAN", "euclidean",
  "MANHATTAN", "manhattan",
  "COSINE",    "cosine",
  "JACCARD",   "jaccard"      # 预计算
)

## ---------- 2. 逐个跑并收集结果 -----------------------------
results <- map(schemes$metric, ~run_pipeline(X_dense, .x))

# 将结果整理到同一 data.frame
unit_names <- rownames(X_dense)
cluster_tbl <- purrr::imap_dfc(results, function(res, i){
  tibble(!!schemes$name[i] := res$cluster)
}) %>% mutate(Unit = unit_names, .before = 1)

## ---------- 3. 查看簇规模 & 交叉对比 -------------------------
walk(schemes$name, ~{
  cat("\n###", .x, "\n")
  print(table(cluster_tbl[[.x]]))
})

# 交叉表示例：EUCLIDEAN vs MANHATTAN
cat("\nContingency (EUCLIDEAN vs MANHATTAN):\n")
print(table(cluster_tbl$EUCLIDEAN, cluster_tbl$MANHATTAN))

## ---------- 4. (可选) 平均轮廓系数 --------------------------
sil_scores <- map2_dbl(results, schemes$name, function(res, nm){
  sil <- silhouette(res$cluster, dist(res$umap))
  mean(sil[, 3])
})
sil_df <- tibble(Scheme = schemes$name, Mean_Silhouette = round(sil_scores, 3))
print(sil_df)

