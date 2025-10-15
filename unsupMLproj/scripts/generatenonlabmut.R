suppressPackageStartupMessages({
  library(dplyr); library(readr); library(googlesheets4)
  library(ALJEbinf); library(Biostrings)
})

# 1) 基础路径
if (!exists("base_dir")) base_dir <- getwd()
data_dir <- file.path(base_dir, "data")
ref_dir  <- file.path(base_dir, "reference")
out_dir  <- file.path(base_dir, "output")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)

# 2) 复用坐标（不存在就报错）
coords_file <- file.path(out_dir, "coordinates.RData")
if (!file.exists(coords_file)) {
  stop("未找到 ", coords_file, "。请先用你原先脚本生成一次 coordinates.RData。")
}
load(coords_file)  # 提供对象：coordinates

# 3) 读取与坐标匹配的参考（必须与生成 coordinates 用的是同一 FASTA）
refs_tbl <- read_csv(file.path(data_dir, "rpoB_reference.csv"), show_col_types = FALSE)
seqs     <- readDNAStringSet(file.path(ref_dir, "rpoB_references.fasta"))

# 4) 下载“非 Lab mutant”记录并最小清洗
gs4_deauth()
norm_origin <- function(x) trimws(tolower(x))

nonlab_df <- read_sheet(
  "https://docs.google.com/spreadsheets/d/1-OihBy1WlfU46_npueE-ZJhmo1eH4eJf7BslsAVh6jc",
  sheet = "Data", col_types = "c"
) |>
  filter(
    !startsWith(Ref_code, "ALJE"),
    Gene == "rpoB",
    is.na(Origin) | norm_origin(Origin) != "lab mutant"
  ) |>
  mutate(across(where(is.character), ~ na_if(., ""))) |>
  mutate(
    AA_pos        = suppressWarnings(as.integer(AA_pos)),
    AA_pos_Ecoli  = suppressWarnings(as.integer(AA_pos_Ecoli)),
    Nt_pos        = suppressWarnings(as.integer(Nt_pos)),
    # 核心：补 Nt_pos，避免 fillMutationsTable 在 NA 上报错
    Nt_pos = ifelse(is.na(Nt_pos) & !is.na(AA_pos),       3L * AA_pos       - 2L, Nt_pos),
    Nt_pos = ifelse(is.na(Nt_pos) & !is.na(AA_pos_Ecoli), 3L * AA_pos_Ecoli - 2L, Nt_pos)
  ) |>
  filter(!is.na(Nt_pos)) |>
  arrange(Species, Ref_code, AA_pos, AA_pos_Ecoli)

write_csv(nonlab_df, file.path(data_dir, "nonlab_reported_mutations.csv"))

# 5) 坐标化（用“旧坐标”）
nonlab_done <- nonlab_df |>
  ALJEbinf::fillMutationsTable(refs_tbl, seqs, coordinates) |>
  filter(!is.na(AA_mut_name_Ecoli))

# 6) 输出
write_csv(nonlab_done, file.path(out_dir, "nonlabmuts.csv"))
write_lines(unique(na.omit(nonlab_done$AA_mut_name_Ecoli)),
            file.path(out_dir, "seen_ecoli_mutations.txt"))

# 保险：确认没混入 Lab mutant
stopifnot(!any(norm_origin(nonlab_done$Origin) == "lab mutant", na.rm = TRUE))


