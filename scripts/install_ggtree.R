#!/usr/bin/env Rscript
# Install ggtree from Bioconductor

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos="https://cloud.r-project.org")
}

BiocManager::install("ggtree", update=FALSE, ask=FALSE)

cat("ggtree installation complete\n")
