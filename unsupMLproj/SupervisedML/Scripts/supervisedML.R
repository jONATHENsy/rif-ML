## -----------------------------------------------------------
## Supervised ML to predict presence of one mutation
## -----------------------------------------------------------

library(randomForest)
library(caret)   # for confusion matrix
library(dplyr)

# 1. 选择目标突变（这里举例第一个突变）
mut_freq <- sort(colSums(X_dense), decreasing = TRUE)
target_mut <- names(mut_freq)[1]   # 例如 rpoB_H526Y

# 2. 准备数据集
X_df <- as.data.frame(X_dense)
y <- as.factor(X_df[[target_mut]])   # 0/1 label
X_df[[target_mut]] <- NULL           # 避免信息泄漏

# 3. 可选加入其他特征：UMAP, cluster
if (exists("umap_df_base")) {
  umap_coords <- umap_df_base %>%
    select(Species, UMAP1, UMAP2)
  X_df <- X_df %>%
    rownames_to_column("Species") %>%
    left_join(umap_coords, by = "Species") %>%
    column_to_rownames("Species")
}

# 4. 构建随机森林模型
set.seed(42)
rf_model <- randomForest(x = X_df, y = y, ntree = 500, importance = TRUE)

# 5. 输出重要特征（有哪些突变与目标突变最强相关）
importance_df <- as.data.frame(importance(rf_model)) %>%
  arrange(desc(MeanDecreaseGini)) %>%
  rownames_to_column("Feature")

write_csv(importance_df, file.path(fig_dir, paste0("RF_importance_", target_mut, ".csv")))

# 6. 评估模型表现
pred <- predict(rf_model, X_df)
conf <- confusionMatrix(pred, y)
print(conf)

# 7. 保存预测标签和真实值
predict_tbl <- tibble(
  Species = rownames(X_df),
  Observed = y,
  Predicted = pred
)

write_csv(predict_tbl, file.path(fig_dir, paste0("RF_prediction_", target_mut, ".csv")))

