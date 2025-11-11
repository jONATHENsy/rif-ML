============================================================
rif-ML — Unsupervised and PU-Learning Pipeline
============================================================

Overview
------------------------------------------------------------
This project analyzes rpoB mutations related to rifampicin resistance
across many bacterial species. The workflow integrates both unsupervised
clustering and positive–unlabeled (PU) supervised learning, coordinated
through a single master script: main.R

------------------------------------------------------------
Directory Structure
------------------------------------------------------------

unsupMLproj/
│
├─ main.R                  – Master workflow (downloads, cleaning, analysis)
│
├─ scripts/
│   ├─ compareandplot.R     – Unsupervised clustering and visualization (UMAP + Heatmap)
│   ├─ compareclustering.R  – test clustering comparison across metrics and methods
│   ├─ confounderscore.R    – Computes per-species research coverage score
│   ├─ dummytest.R          – Dummy dataset test for clustering functions
│   ├─ generatemutations.R  – Generates lab and non-lab mutation tables
│   ├─ packages.R           – Ensure all required packages can be installed
│   ├─ pulearning.R         – PU-learning (supervised) analysis and plotting
│
├─ data/                   – Raw and intermediate CSV / RDS files(include specie/strain rference data)
├─ output/                 – Generated and saved data/ results (labmuts.csv, nonlabmuts.csv, coordinates, etc.)
├─ figures/                – AUMAP, heatmaps, UpSet, etc.
└─ reference/              – Reference FASTA 

------------------------------------------------------------
Main Script: main.R
------------------------------------------------------------

The entry point for the entire pipeline.

Workflow steps:

0) Setup
1) Download and clean mutation records (lab + nonlab)
2) Load reference FASTA and coordinate table (with FASTA hash cache)
3) Run fillMutationsTable() → produce labmuts.csv & nonlabmuts.csv
4) Optionally normalize X_dense column names to “rpoB_” prefix
5) Run downstream analysis: unsupervised clustering + PU-learning

Each downstream script automatically inherits file paths and directories
from main.R.



------------------------------------------------------------
Credits
------------------------------------------------------------

If you use this pipeline or any derivative work, please cite the ALJEbinf
R package for coordinate mapping and the rpoB mutation dataset obtained
from the shared Google Sheet.

------------------------------------------------------------
End of README.txt