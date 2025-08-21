## -----------------------------------------------------------
## Packages
## -----------------------------------------------------------
library(tidyverse)   # dplyr / tidyr / ggplot2
library(Matrix)      # sparse matrices
library(FactoMineR)  # PCA (optional)
library(umap)        # legacy UMAP (not used here)
# install.packages("uwot")   # install once if missing
library(uwot)        # current UMAP implementation
library(dbscan)      # HDBSCAN clustering
library(patchwork)   # multi-plot layout
library(ggrepel)     # nicer text labels

## -----------------------------------------------------------
## File paths
## -----------------------------------------------------------
base_dir <- "D:/2025s1/BIOX7011/processed"
mut_file <- file.path(base_dir, "output", "labmuts.csv")

## -----------------------------------------------------------
## 1.  Read mutation table and build 0/1 sparse matrix
## -----------------------------------------------------------
muts <- read_csv(mut_file, show_col_types = FALSE)

unit_id <- muts %>%                # define analysis unit = Species__RefSeq_ID
  mutate(Unit = paste(Species, RefSeq_ID, sep = "__"))

analysis_tbl <- unit_id %>%        # keep only Unit + orthologous AA mutation
  select(Unit, Orthologous_AA = AA_mut_name_Ecoli) %>% 
  distinct()

mut_index <- analysis_tbl %>%      # give each unique mutation an integer idx
  distinct(Orthologous_AA) %>% 
  arrange(Orthologous_AA) %>% 
  mutate(idx = row_number())

enc <- analysis_tbl %>%            # merge indices back
  inner_join(mut_index, by = "Orthologous_AA")

X <- sparseMatrix(                 # rows = Units, cols = mutations
  i = as.integer(factor(enc$Unit, levels = unique(enc$Unit))),
  j = enc$idx,
  x = 1,
  dims = c(n_distinct(enc$Unit), nrow(mut_index)),
  dimnames = list(unique(enc$Unit), mut_index$Orthologous_AA)
)

## -----------------------------------------------------------
## 2.  UMAP on dense matrix (works for ≤ ~30 M cells)
##     – convert sparse → dense and keep 0/1
## -----------------------------------------------------------
X_dense <- as.matrix(X)
X_dense[X_dense > 0] <- 1           # ensure binary

set.seed(123)
umap_res <- uwot::umap(
  X_dense,
  n_neighbors = 15,
  min_dist    = 0.3,
  metric      = "euclidean",        # "cosine" also fine
  verbose     = TRUE
)

umap_df <- as_tibble(umap_res, .name_repair = "unique") %>%  # add col names
  setNames(c("UMAP1", "UMAP2")) %>%
  mutate(Unit = rownames(X_dense))

## -----------------------------------------------------------
## 3.  Quick UMAP scatter plot
## -----------------------------------------------------------
ggplot(umap_df, aes(UMAP1, UMAP2)) +
  geom_point(size = 3, alpha = .8) +
  theme_classic()

## -----------------------------------------------------------
## 4.  HDBSCAN clustering in UMAP space
## -----------------------------------------------------------
hdb <- hdbscan(umap_df[, c("UMAP1", "UMAP2")], minPts = 5)
umap_df$Cluster <- factor(hdb$cluster)

ggplot(umap_df, aes(UMAP1, UMAP2, colour = Cluster)) +
  geom_point(size = 3) +
  scale_colour_brewer(palette = "Set1") +
  theme_classic()

## -----------------------------------------------------------
## 5.  List top mutations per cluster
## -----------------------------------------------------------
top_muts <- function(cl){
  rows <- umap_df$Unit[umap_df$Cluster == cl]
  colSums(X[rows, , drop = FALSE]) %>% 
    sort(decreasing = TRUE) %>% head(8)
}
lapply(levels(factor(umap_df$Cluster))[-1], top_muts)   # skip noise (0)

## -----------------------------------------------------------
## 6.  Save figure to repo
## -----------------------------------------------------------
p_umap <- ggplot(umap_df, aes(UMAP1, UMAP2, colour = Cluster)) +
  geom_point(size = 3, alpha = .85) +
  scale_colour_brewer(palette = "Set1", na.translate = FALSE) +
  theme_classic(base_size = 14) +
  labs(title = "UMAP + HDBSCAN clusters")

out_dir  <- "C:/Users/stupidZHUO/rif-ML/processeddata"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(out_dir, "umap_clusters.png"),
  plot     = p_umap,
  width    = 6,     # inches
  height   = 5,
  dpi      = 300
)

cat("✅ Figure saved to", file.path(out_dir, "umap_clusters.png"), "\n")

## -----------------------------------------------------------
## 7.  Heat-map of 30 most frequent mutations
## -----------------------------------------------------------
library(ComplexHeatmap)

mut_freq     <- colSums(X)                              # column frequencies
top_30_cols  <- names(sort(mut_freq, decreasing = TRUE))[1:30]

sub_mat <- as.matrix(X[umap_df$Unit, top_30_cols])      # keep same row order
sub_mat[sub_mat > 0] <- 1

Heatmap(
  sub_mat,
  name            = "Mut",
  col             = c("0" = "white", "1" = "steelblue"),
  cluster_rows    = FALSE,
  cluster_columns = TRUE,
  show_row_names  = FALSE,
  row_split       = umap_df$Cluster
)

