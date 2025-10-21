#install.packages("devtools")
#devtools::install_github("JanEngelstaedter/ALJEbinf", build_vignettes = TRUE)

library(googlesheets4)
library(tidyverse)
library(ALJEbinf)

gs4_deauth() #function to inactive Google authorization
muts <- read_sheet("https://docs.google.com/spreadsheets/d/1-OihBy1WlfU46_npueE-ZJhmo1eH4eJf7BslsAVh6jc/edit?usp=sharing",
                   sheet = "Data",
                   col_types = "c") |>
  filter(substr(Ref_code, 1, 4) != "ALJE") |>
  filter(Gene == "rpoB") |>
  filter(Origin %in% c("Isolate")) |>
  arrange(Species, Ref_code, AA_pos, AA_pos_Ecoli)
write_csv(muts, "./data/reported_mutations.csv")

refs <- read_sheet("https://docs.google.com/spreadsheets/d/1-OihBy1WlfU46_npueE-ZJhmo1eH4eJf7BslsAVh6jc/edit?usp=sharing",
                   sheet = "References",
                   col_types = "c") |>
  select(-Comments) |>
  mutate(DOI = sub("\\.$", "", Link)) |>
  select(-Link) |>
  arrange(Ref_code)

if (length(refs$Ref_code) != (length(unique(refs$Ref_code)))) {
  warning("Duplicated Ref_code entries detected!")
}

write_csv(refs, "./data/reports_references.csv")

#load reference sequences and their information 
refs <- read_csv("./data/rpoB_references.csv", show_col_types = FALSE)
seqs <- readDNAStringSet("./data/rpoB_references.fasta")
mutation_list_reports <- read_csv("./data/reported_mutations.csv", show_col_types = FALSE)

# Check whether the above files have been changed and hence the coordinates need to be updated 
if(file.exists("./data/fastahash.Rds") && as.character(openssl::sha1(file("./data/rpoB_references.fasta"))) == readRDS("./data/fastahash.Rds")){
  print("Sequences file has not changed, loading original coordinates")
  load(file = "./output/coordinates.RData")
} else {
  old_fastahash <- as.character(openssl::sha1(file("./data/rpoB_references.fasta")))
  saveRDS(old_fastahash, 
          file = "./data/fastahash.Rds")
  print("Sequences file has changed, regenerating coordinates")
  #get coordinates
  coordinates <- ALJEbinf::getAllCoordinates(seqs, "rpoB_Escherichia_coli_MG1655")
  save(coordinates, file = "./output/coordinates.RData")
}

#2 load and complete table of reported mutations:
muts <- mutation_list_reports |>
  fillMutationsTable(refs, seqs, coordinates) |>
  filter(!is.na(AA_mut_name_Ecoli))  # filter out "bad" entries that couldn't be mapped to E. coli

write_csv(muts, "./output/muts.csv")

