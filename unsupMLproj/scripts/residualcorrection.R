# -------------------------------------
# ✅ 加载依赖包
# -------------------------------------
library(readr)
library(dplyr)
library(vegan)     # for adonis2
library(umap)      # for visualization (optional)
library(tibble)    # for column_to_rownames

# -------------------------------------
# ✅ Step 1: 读取 labmuts 数据和构建突变矩阵
# -------------------------------------
base_dir <- "C:/Users/user/Desktop/D Drive/2025s1/BIOX7011/rif-ML/unsupMLproj"
labmuts_path <- file.path(base_dir, "output", "labmuts.csv")
labmuts <- read_csv(labmuts_path)

# 清理物种名
labmuts <- labmuts %>%
  filter(!is.na(Species)) %>%
  mutate(Species = trimws(Species))

# 构建突变 presence 矩阵
X_dense <- labmuts %>%
  mutate(presence = 1) %>%
  select(Species, Mutation = AA_mut_name_Ecoli, presence) %>%
  distinct() %>%
  tidyr::pivot_wider(names_from = Mutation, values_from = presence, values_fill = 0)

X_mat <- X_dense %>% column_to_rownames("Species")

# -------------------------------------
# ✅ Step 2: 构建 study_score，划分 study_group
# -------------------------------------
confounder_df <- labmuts %>%
  group_by(Species) %>%
  summarise(
    num_mutations = n(),
    num_studies = n_distinct(Ref_code),
    study_score = log1p(num_mutations) + log1p(num_studies)
  )

# 对 Species 排序以对齐矩阵
confounder_df <- confounder_df %>%
  filter(Species %in% rownames(X_mat)) %>%
  arrange(match(Species, rownames(X_mat)))

# 添加 study_score 和 group
study_score <- confounder_df$study_score
names(study_score) <- confounder_df$Species

study_group <- ifelse(study_score >= quantile(study_score, 0.75), "High",
                      ifelse(study_score <= quantile(study_score, 0.25), "Low", NA))
study_group <- factor(study_group)

# 只保留 High 和 Low 的物种
valid_species <- names(study_group[!is.na(study_group)])
X_mat_sub <- X_mat[valid_species, ]
study_group_sub <- study_group[valid_species]
study_score_sub <- study_score[valid_species]

# -------------------------------------
# ✅ Step 3: Residualization（去除 study_score 影响）
# -------------------------------------
X_resid <- X_mat_sub
for (i in 1:ncol(X_mat_sub)) {
  fit <- lm(X_mat_sub[, i] ~ study_score_sub)
  X_resid[, i] <- residuals(fit)
}

# -------------------------------------
# ✅ Step 4: PERMANOVA 检验是否还存在 confounding bias
# -------------------------------------
adonis_resid <- adonis2(X_resid ~ study_group_sub, permutations = 999)
cat("\n✅ Residualized PERMANOVA result:\n")
print(adonis_resid)

write.csv(X_resid, file = file.path(base_dir, "output", "X_residualized.csv"), row.names = TRUE)

