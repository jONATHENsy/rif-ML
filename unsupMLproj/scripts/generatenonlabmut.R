suppressPackageStartupMessages({
  library(googlesheets4)
  library(dplyr)
  library(readr)
  library(ALJEbinf)
  library(Biostrings)
  library(openssl)
})

# ---------- 可选：基准目录 ----------
if (!exists("base_dir")) base_dir <- getwd()
data_dir  <- file.path(base_dir, "data")
out_dir   <- file.path(base_dir, "output")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)

gs4_deauth()   # 公开表，无需 OAuth

norm_origin <- function(x) trimws(tolower(x))

# -------------------------------------------------
# 1) 下载并保存「非 Lab Mutant」原始行
# -------------------------------------------------
muts <- read_sheet(
  "https://docs.google.com/spreadsheets/d/1-OihBy1WlfU46_npueE-ZJhmo1eH4eJf7BslsAVh6jc",
  sheet = "Data", col_types = "c"
) %>%
  filter(
    !startsWith(Ref_code, "ALJE"),
    Gene == "rpoB",
    is.na(Origin) | norm_origin(Origin) != "lab mutant"  # 反选（鲁棒）
  ) %>%
  mutate(
    across(everything(), ~na_if(., "")),
    AA_pos       = suppressWarnings(as.integer(AA_pos)),
    AA_pos_Ecoli = suppressWarnings(as.integer(AA_pos_Ecoli))
  ) %>%
  arrange(Species, Ref_code, AA_pos, AA_pos_Ecoli)

write_csv(muts, file.path(data_dir, "allother_mutations.csv"))
message(sprintf("非 Lab Mutant 原始行数：%d", nrow(muts)))

# -------------------------------------------------
# 2) 下载并保存 References（去掉末尾句号）
# -------------------------------------------------
refs <- read_sheet(
  "https://docs.google.com/spreadsheets/d/1-OihBy1WlfU46_npueE-ZJhmo1eH4eJf7BslsAVh6jc",
  sheet = "References", col_types = "c"
) %>%
  select(-Comments) %>%
  mutate(DOI = sub("\\.$", "", Link)) %>%
  select(-Link) %>%
  arrange(Ref_code)

if (length(refs$Ref_code) != length(unique(refs$Ref_code))) {
  warning("Duplicated Ref_code entries detected in References!")
}
write_csv(refs, file.path(data_dir, "reports_references.csv"))

# -------------------------------------------------
# 3) 读取本地 reference FASTA / CSV
# -------------------------------------------------
refs_tbl <- read_csv(file.path(data_dir, "rpoB_reference.csv"),
                     show_col_types = FALSE)
seqs  <- readDNAStringSet(file.path(data_dir, "rpoB_references.fasta"))

# -------------------------------------------------
# 4) 生成或读取坐标（按 fasta 的 sha1 缓存）
# -------------------------------------------------
fastahash_file <- file.path(data_dir, "fastahash.Rds")
coords_file    <- file.path(out_dir,  "coordinates.RData")
fasta_sha1     <- as.character(sha1(file(file.path(data_dir, "rpoB_references.fasta"))))

if (file.exists(fastahash_file) && identical(readRDS(fastahash_file), fasta_sha1)) {
  message("Sequences unchanged, loading saved coordinates.")
  load(coords_file)  # 加载对象：coordinates
} else {
  message("Sequences updated, rebuilding coordinates.")
  coordinates <- ALJEbinf::getAllCoordinates(
    seqs, "rpoB_Escherichia_coli_MG1655"
  )
  saveRDS(fasta_sha1, fastahash_file)
  save(coordinates, file = coords_file)
}

# -------------------------------------------------
# 5) 用非 Lab Mutant 数据生成完整突变表
#    关键修复：删除核苷酸层面的列，避免 NA 触发 if(NA) 报错
# -------------------------------------------------
mutation_list_reports <- muts %>%
  select(-any_of(c("Nt_pos","Nt_from","Nt_to","Nt_ref","Nt_mut",
                   "Nucleotide_change","Codon_pos")))

muts_completed <- mutation_list_reports %>%
  ALJEbinf::fillMutationsTable(refs_tbl, seqs, coordinates) %>%
  filter(!is.na(AA_mut_name_Ecoli))

# 保险起见：再次确认未混入 Lab Mutant
if (any(norm_origin(muts_completed$Origin) == "lab mutant", na.rm = TRUE)) {
  stop("检测到 Lab Mutant 混入，请检查上游筛选。")
}

# 输出
out_file <- file.path(out_dir, "nonlabmutant_results.csv")
write_csv(muts_completed, out_file)
message(sprintf("已写出：%s（行数：%d）", out_file, nrow(muts_completed)))
