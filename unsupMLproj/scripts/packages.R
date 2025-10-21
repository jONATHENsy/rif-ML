## —— 0) 可选：固定一个本地库路径，避免权限问题（自行修改路径）——
# dir.create("C:/Rlibs", showWarnings = FALSE, recursive = TRUE)
# .libPaths(c("C:/Rlibs", .libPaths()))

## —— 1) 安装/加载 BiocManager（管理 Bioconductor 包）——
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
BiocManager::version()  # 看看Bioconductor版本是否与R 4.5.x匹配

## —— 2) 你项目常用的包清单 —— 
cran_pkgs <- c(
  # tidyverse 及常用扩展
  "tidyverse","ggplot2","scales","data.table","Matrix","stringr","readr","readxl","writexl",
  "purrr","tibble","forcats",
  # 可视化与排版
  "cowplot","ggrepel","patchwork","gridExtra","RColorBrewer","viridis",
  # 你脚本中出现过
  "UpSetR","uwot","circlize"
)

bioc_pkgs <- c(
  # 你现在就报缺的：
  "ComplexHeatmap","Biostrings",
  # ComplexHeatmap/生信常用依赖（Bioc会自动装依赖，但列出来更直观）
  "IRanges","S4Vectors","XVector","BiocGenerics",
  # 你的作业/项目里提到过的
  "DECIPHER"            # 序列/功能注释里可能用到
)

## —— 3) 安装 CRAN 包（跳过已安装）——
to_install_cran <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(to_install_cran)) {
  install.packages(to_install_cran, repos = "https://cloud.r-project.org")
}

## —— 4) 安装 Bioconductor 包（跳过已安装）——
to_install_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(to_install_bioc)) {
  BiocManager::install(to_install_bioc, update = TRUE, ask = FALSE)
}

## —— 5) 验证是否都装好了 —— 
ok <- sapply(c(cran_pkgs, bioc_pkgs), requireNamespace, quietly = TRUE)
sort(ok)
