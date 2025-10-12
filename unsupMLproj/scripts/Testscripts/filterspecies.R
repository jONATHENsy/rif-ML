# =============================================================
# filter_species.R – 筛选研究充分度适中、结构有代表性的物种
# 输出核心物种 + 子矩阵
# =============================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# 设置基本路径
if (!exists("base_dir"))
  base_dir <- "C:/Users/user/Desktop/D Drive/2025s1/BIOX7011/rif-ML/unsupMLproj"

labmuts_path <- file.path(base_dir, "output", "labmuts.csv")
labmuts <- read_csv(labmuts_path, show_col_types = FALSE)

# 清理物种名
labmuts <- labmuts %>%
  filter(!is.na(Species)) %>%
  mutate(Species = trimws(Species))

# 计算研究充分度指标（confounder_score）
confounder_df <- labmuts %>%
  group_by(Species) %>%
  summarise(
    num_mutations = n(),
    num_studies = n_distinct(Ref_code),
    study_score = log1p(num_mutations) + log1p(num_studies),
    .groups = "drop"
  ) %>%
  arrange(desc(study_score))

# 修改筛选标准（放宽要求）
q_low <- quantile(confounder_df$study_score, 0.7, na.rm = TRUE)
q_high <- quantile(confounder_df$study_score, 1, na.rm = TRUE)

core_species <- confounder_df %>%
  filter(num_mutations >= 3, num_studies >= 1, study_score >= q_low, study_score <= q_high) %>%
  pull(Species)

# 保存核心物种列表
write_csv(tibble(Species = core_species), file.path(base_dir, "output", "core_species.csv"))

# 构建 X_dense_filtered 矩阵
unit_tbl <- labmuts %>%
  mutate(Species_clean = str_replace(Strain_ID, "_[^_]+$", ""),
         Species_clean = str_replace_all(Species_clean, "_", " ")) %>%
  select(Species = Species_clean, Orthologous_AA = AA_mut_name_Ecoli) %>%
  distinct()

unit_tbl <- unit_tbl %>% filter(Species %in% core_species)

mut_index <- unit_tbl %>%
  distinct(Orthologous_AA) %>%
  arrange(Orthologous_AA) %>%
  mutate(idx = row_number())

enc <- inner_join(unit_tbl, mut_index, by = "Orthologous_AA")

X <- Matrix::sparseMatrix(
  i = as.integer(factor(enc$Species, levels = unique(enc$Species))),
  j = enc$idx,
  x = 1,
  dims = c(n_distinct(enc$Species), nrow(mut_index)),
  dimnames = list(unique(enc$Species), mut_index$Orthologous_AA)
)

X_dense_filtered <- as.matrix(X)
X_dense_filtered[X_dense_filtered > 0] <- 1

saveRDS(X_dense_filtered, file = file.path(base_dir, "output", "X_dense_high.RDS"))
message("✓ filter_species.R completed. Core species count:", nrow(X_dense_filtered))
unit_tbl_filtered <- unit_tbl %>% filter(Species %in% core_species)
message("📈 对应突变数：", unit_tbl_filtered %>% distinct(Orthologous_AA) %>% nrow())

