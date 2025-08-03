# ==============================================================
#  main.R  –  Master workflow (download → clean → unsupervised)
# ==============================================================

## --------------------------------------------------------------
## 0. Base directory (edit if you move the project)
## --------------------------------------------------------------
base_dir <- "D:/2025s1/BIOX7011/rif-ML/unsupMLproj"

dir.create(file.path(base_dir, "data"),   FALSE, TRUE)
dir.create(file.path(base_dir, "output"), FALSE, TRUE)
dir.create(file.path(base_dir, "figures"), FALSE, TRUE)
setwd(base_dir)

## --------------------------------------------------------------
## 1. Packages
## --------------------------------------------------------------
pkgs <- c("googlesheets4", "tidyverse", "Biostrings", "openssl",
          "ALJEbinf", "Matrix", "uwot", "dbscan", "ComplexHeatmap",
          "patchwork", "ggrepel", "cluster", "RColorBrewer")
new <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(new)) install.packages(new)
invisible(lapply(pkgs, library, character.only = TRUE))

# No authentication needed
gs4_deauth()

## --------------------------------------------------------------
## 2. Download & clean mutation sheet
## --------------------------------------------------------------
muts_path <- file.path("data", "lab_reported_mutations.csv")

muts_raw <- read_sheet(
  "https://docs.google.com/spreadsheets/d/1-OihBy1WlfU46_npueE-ZJhmo1eH4eJf7BslsAVh6jc",
  sheet = "Data", col_types = "c") |>
  filter(substr(Ref_code, 1, 4) != "ALJE",
         Gene == "rpoB",
         Origin == "Lab mutant") |>
  arrange(Species, Ref_code, AA_pos, AA_pos_Ecoli)

write_csv(muts_raw, muts_path)
message(format(Sys.time(), "%T"), " – mutations saved to ", muts_path)

# Clean Nt_pos to avoid downstream NA issues
muts_clean <- muts_raw %>%
  filter(!is.na(Nt_pos), grepl("^[0-9]+$", Nt_pos)) %>%
  mutate(Nt_pos = as.integer(Nt_pos))

## --------------------------------------------------------------
## 3. Download reference metadata
## --------------------------------------------------------------
refs_path <- file.path("data", "reports_references.csv")

refs <- read_sheet(
  "https://docs.google.com/spreadsheets/d/1-OihBy1WlfU46_npueE-ZJhmo1eH4eJf7BslsAVh6jc",
  sheet = "References", col_types = "c") |>
  select(-Comments) |>
  mutate(DOI = sub("\\.$", "", Link)) |>
  select(-Link) |>
  arrange(Ref_code)

write_csv(refs, refs_path)
message(format(Sys.time(), "%T"), " – references saved to ", refs_path)

## --------------------------------------------------------------
## 4. Load reference FASTA and coordinates
## --------------------------------------------------------------
refs_tbl <- read_csv(file.path(base_dir, "reference", "rpoB_reference.csv"),
                     show_col_types = FALSE)
seqs  <- readDNAStringSet(file.path(base_dir, "reference", "rpoB_references.fasta"))
mutation_list_reports <- read_csv(
  file.path(base_dir, "data", "lab_reported_mutations.csv"),
  show_col_types = FALSE)

fastahash_file <- file.path(base_dir, "data", "fastahash.Rds")
coords_file    <- file.path(base_dir, "output", "coordinates.RData")

if (file.exists(fastahash_file) &&
    as.character(openssl::sha1(file(file.path(base_dir, "reference", "rpoB_references.fasta")))) ==
    readRDS(fastahash_file)) {
  
  message("Sequences unchanged, loading saved coordinates")
  load(coords_file)
  
} else {
  saveRDS(as.character(openssl::sha1(
    file(file.path(base_dir, "reference", "rpoB_references.fasta")))),
    fastahash_file)
  
  message("Sequences updated, rebuilding coordinates")
  coordinates <- ALJEbinf::getAllCoordinates(
    seqs, "rpoB_Escherichia_coli_MG1655")
  save(coordinates, file = coords_file)
}
## --------------------------------------------------------------
## 5. Run fillMutationsTable and write output
## --------------------------------------------------------------
muts_completed <- muts_clean %>%
  fillMutationsTable(refs_tbl, seqs, coordinates) %>%
  filter(!is.na(AA_mut_name_Ecoli))

write_csv(muts_completed, "output/labmuts.csv")
message(format(Sys.time(), "%T"), " – labmuts.csv written")

## --------------------------------------------------------------
## 6. Run downstream analysis scripts
## --------------------------------------------------------------
script_dir <- file.path(base_dir, "scripts")

for (script in c("compareclustering.R", "clusteringplot.R")) {
  scr_path <- file.path(script_dir, script)
  if (file.exists(scr_path)) {
    message("---- Running ", script, " ----")
    source(scr_path, echo = TRUE)
  } else {
    warning("Script not found: ", scr_path)
  }
}

message("\n✅  Pipeline complete – all scripts executed.")
