# ==============================================================
#  main.R  –  Master workflow (download → clean → unsupervised)
# ==============================================================

## --------------------------------------------------------------
## 0. Base directory (edit if you move the project)
## --------------------------------------------------------------
base_dir <- "D:/2025s1/BIOX7011/rif-ML/unsupMLproj"

dir.create(file.path(base_dir, "data"),   showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "output"), showWarnings = FALSE, recursive = TRUE)

setwd(base_dir)

## --------------------------------------------------------------
## 1. Packages
## --------------------------------------------------------------
pkgs <- c(
  # data wrangling
  "googlesheets4", "tidyverse", "Biostrings", "openssl",
  # ALJE - mutation utilities
  "ALJEbinf",
  # ML & visualisation
  "Matrix", "uwot", "dbscan", "ComplexHeatmap",
  "patchwork", "ggrepel", "cluster","RColorBrewer"
)
to_ins <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_ins)) install.packages(to_ins)
lapply(pkgs, library, character.only = TRUE)

gs4_deauth()

## --------------------------------------------------------------
## 2. Download mutation sheet (lab mutants only)
## --------------------------------------------------------------
muts_path <- file.path("data", "lab_reported_mutations.csv")

muts_raw <- read_sheet(
  "https://docs.google.com/spreadsheets/d/1-OihBy1WlfU46_npueE-ZJhmo1eH4eJf7BslsAVh6jc",
  sheet     = "Data",
  col_types = "c"
) |>
  filter(substr(Ref_code, 1, 4) != "ALJE",
         Gene == "rpoB",
         Origin %in% "Lab mutant") |>
  arrange(Species, Ref_code, AA_pos, AA_pos_Ecoli)

write_csv(muts_raw, muts_path)
message(format(Sys.time(), "%T"), " – mutations saved.")

## --------------------------------------------------------------
## 3. Download reference list
## --------------------------------------------------------------
refs_path <- file.path("data", "reports_references.csv")

refs <- read_sheet(
  "https://docs.google.com/spreadsheets/d/1-OihBy1WlfU46_npueE-ZJhmo1eH4eJf7BslsAVh6jc",
  sheet     = "References",
  col_types = "c"
) |>
  select(-Comments) |>
  mutate(DOI = sub("\\.$", "", Link)) |>
  select(-Link) |>
  arrange(Ref_code)

write_csv(refs, refs_path)
message(format(Sys.time(), "%T"), " – references saved.")

