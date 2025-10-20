# ---- 脚本说明 ----
# 输入: labmuts.csv（包含物种与文献数据）
# 输出: confounder_score.csv（每个物种的研究充分度评分）

library(dplyr)
library(readr)

# Step 1: 读取数据
labmuts <- read_csv("labmuts.csv")  # 请根据实际路径调整

# Step 2: 标准化物种名列（避免 NA 或拼写误差）
labmuts <- labmuts %>%
  filter(!is.na(Species)) %>%
  mutate(Species = trimws(Species))

# Step 3: 聚合物种研究充分度指标
confounder_df <- labmuts %>%
  group_by(Species) %>%
  summarise(
    num_mutations = n(),  # 突变条目总数
    num_studies = n_distinct(Ref_code),  # 涉及文献数
    year_range = max(Year, na.rm = TRUE) - min(Year, na.rm = TRUE),  # 研究时间跨度
    study_score = log1p(num_mutations) + log1p(num_studies) + log1p(1 + year_range)
  ) %>%
  arrange(desc(study_score))

# Step 4: 导出评分表
write_csv(confounder_df, "confounder_score.csv")
