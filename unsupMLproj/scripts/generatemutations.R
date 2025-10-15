# -------------------------------------------------
# 0. 基础目录（只需改这一行）
# -------------------------------------------------
base_dir <- "D:/2025s1/BIOX7011/processed"
dir.create(file.path(base_dir, "data"),   showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "output"), showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------
# 1. 加载包
# -------------------------------------------------
library(googlesheets4)
library(tidyverse)
library(ALJEbinf)
library(Biostrings)

gs4_deauth()   # 无需 OAuth

# -------------------------------------------------
# 2. 下载突变表
# -------------------------------------------------
muts <- read_sheet(
  "https://docs.google.com/spreadsheets/d/1-OihBy1WlfU46_npueE-ZJhmo1eH4eJf7BslsAVh6jc",
  sheet = "Data", col_types = "c") |>
  filter(substr(Ref_code, 1, 4) != "ALJE",
         Gene == "rpoB",
         Origin %in% "Isolate") |>
  arrange(Species, Ref_code, AA_pos, AA_pos_Ecoli)

write_csv(muts, file.path(base_dir, "data", "allother_mutations.csv"))

# -------------------------------------------------
# 3. 下载 References
# -------------------------------------------------
refs <- read_sheet(
  "https://docs.google.com/spreadsheets/d/1-OihBy1WlfU46_npueE-ZJhmo1eH4eJf7BslsAVh6jc",
  sheet = "References", col_types = "c") |>
  select(-Comments) |>
  mutate(DOI = sub("\\.$", "", Link)) |>
  select(-Link) |>
  arrange(Ref_code)

if (length(refs$Ref_code) != length(unique(refs$Ref_code)))
  warning("Duplicated Ref_code entries detected!")

write_csv(refs, file.path(base_dir, "data", "reports_references.csv"))

# -------------------------------------------------
# 4. 读取 reference FASTA / CSV
# -------------------------------------------------
refs_tbl <- read_csv(file.path(base_dir, "data", "rpoB_reference.csv"),
                     show_col_types = FALSE)
seqs  <- readDNAStringSet(file.path(base_dir, "data", "rpoB_references.fasta"))
mutation_list_reports <- read_csv(
  file.path(base_dir, "data", "lab_reported_mutations.csv"),
  show_col_types = FALSE)

# -------------------------------------------------
# 5. 生成或读取坐标
# -------------------------------------------------
fastahash_file <- file.path(base_dir, "data", "fastahash.Rds")
coords_file    <- file.path(base_dir, "output", "coordinates.RData")

if (file.exists(fastahash_file) &&
    as.character(openssl::sha1(file(file.path(base_dir, "data", "rpoB_references.fasta")))) ==
    readRDS(fastahash_file)) {
  
  message("Sequences unchanged, loading saved coordinates")
  load(coords_file)
  
} else {
  saveRDS(as.character(openssl::sha1(
    file(file.path(base_dir, "data", "rpoB_references.fasta")))),
    fastahash_file)
  
  message("Sequences updated, rebuilding coordinates")
  coordinates <- ALJEbinf::getAllCoordinates(
    seqs, "rpoB_Escherichia_coli_MG1655")
  save(coordinates, file = coords_file)
}

# -------------------------------------------------
# 6. 填充突变表并输出
# -------------------------------------------------
muts_completed <- mutation_list_reports |>
  fillMutationsTable(refs_tbl, seqs, coordinates) |>
  filter(!is.na(AA_mut_name_Ecoli))

write_csv(muts_completed, file.path(base_dir, "output", "othermuts.csv"))
