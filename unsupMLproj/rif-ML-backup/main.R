# ==============================================================
#  main.R  –  Master workflow (download → clean → unsupervised)
# ==============================================================

## --------------------------------------------------------------
## 0. Base directory (edit if you move the project)
## --------------------------------------------------------------
base_dir <- "D:/2025s1/BIOX7011/rif-ML/unsupMLproj"

dir.create(file.path(base_dir, "data"),   FALSE, TRUE)
dir.create(file.path(base_dir, "output"), FALSE, TRUE)
dir.create(file.path(base_dir, "figures"), FALSE, TRUE)   # for later plots
setwd(base_dir)

## --------------------------------------------------------------
## 1. Packages  (install once, then load)
## --------------------------------------------------------------
pkgs <- c(
  # data / I/O
  "googlesheets4", "tidyverse", "Biostrings", "openssl",
  # mutation helpers
  "ALJEbinf",
  # ML & viz
  "Matrix", "uwot", "dbscan", "ComplexHeatmap",
  "patchwork", "ggrepel", "cluster", "RColorBrewer"
)
new <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(new)) install.packages(new)
invisible(lapply(pkgs, library, character.only = TRUE))

gs4_deauth()                 # public Sheet → no OAuth popup

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
         Origin == "Lab mutant") |>
  arrange(Species, Ref_code, AA_pos, AA_pos_Ecoli)

write_csv(muts_raw, muts_path)
message(format(Sys.time(), "%T"), " – mutations saved to ", muts_path)

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
message(format(Sys.time(), "%T"), " – references saved to ", refs_path)

## --------------------------------------------------------------
## 4. Reference FASTA  & coordinate generation
## --------------------------------------------------------------
refs_tbl <- read_csv("data/rpoB_reference.csv", show_col_types = FALSE)
seqs     <- readDNAStringSet("data/rpoB_references.fasta")
mutation_list_reports <- read_csv(muts_path, show_col_types = FALSE)

fastahash_file <- "data/fastahash.Rds"
coords_file    <- "output/coordinates.RData"

if (file.exists(fastahash_file) &&
    as.character(openssl::sha1(file("data/rpoB_references.fasta"))) ==
    readRDS(fastahash_file)) {
  
  load(coords_file)
  message("→ Coordinates loaded from cache")
  
} else {
  saveRDS(as.character(openssl::sha1(file("data/rpoB_references.fasta"))),
          fastahash_file)
  message("→ Re-building coordinates …")
  coordinates <- ALJEbinf::getAllCoordinates(
    seqs, "rpoB_Escherichia_coli_MG1655")
  save(coordinates, file = coords_file)
}

## --------------------------------------------------------------
## 5. Complete mutation table & export
## --------------------------------------------------------------
muts_completed <- mutation_list_reports |>
  fillMutationsTable(refs_tbl, seqs, coordinates) |>
  filter(!is.na(AA_mut_name_Ecoli))

write_csv(muts_completed, "output/labmuts.csv")
message(format(Sys.time(), "%T"), " – labmuts.csv written")

## ==============================================================
## 6.  RUN downstream analysis scripts
##    (each script assumes setwd(base_dir) already)
## ==============================================================

script_dir <- file.path(base_dir, "scripts")

for (scr in c("generatemutations.R",
              "clusteringplot.R",
              "compareclustering.R")) {
  scr_path <- file.path(script_dir, scr)
  if (file.exists(scr_path)) {
    message("---- Running ", scr, " ----")
    source(scr_path, echo = TRUE)   # echo=TRUE prints lines to console
  } else {
    warning("Script not found: ", scr_path)
  }
}

message("\n✅  Pipeline complete – all scripts executed.")


