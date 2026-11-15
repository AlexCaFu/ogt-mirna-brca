#!/usr/bin/env Rscript
# =====================================================================
# make_repro_metadata.R
# Captures the computational environment (session information and package
# versions) of the R run.
#
# RUN THIS from the repository root, ideally in the SAME R environment
# that ran the analysis (so the recorded package versions are the real
# ones). Uses base R only; no packages are installed.
#
#   setwd("/path/to/ogt-mirna-brca")
#   source("make_repro_metadata.R")
#
# The checksum list of the deposit (MANIFEST.sha256) is generated separately.
# =====================================================================

root <- getwd()
if (!dir.exists("results")) dir.create("results")
message("Repository root: ", root)

## 1. Full session information -----------------------------------------
si_path <- file.path("results", "sessionInfo.txt")
writeLines(capture.output(sessionInfo()), si_path)
message("[1/2] Wrote ", si_path)

## 2. Exact versions of the analysis packages --------------------------
key_pkgs <- c(
  "multiMiR", "TCGAbiolinks", "SummarizedExperiment", "edgeR", "limma",
  "GEOquery", "survival", "readr", "dplyr", "tidyr", "ggplot2",
  "Biobase", "BiocManager", "matrixStats"
)
ip <- as.data.frame(installed.packages()[, c("Package", "Version")],
                    stringsAsFactors = FALSE)
sel <- ip[ip$Package %in% key_pkgs, , drop = FALSE]
sel <- sel[order(sel$Package), ]
pv_path <- file.path("results", "package_versions.csv")
utils::write.csv(sel, pv_path, row.names = FALSE)
message("[2/2] Wrote ", pv_path, " (", nrow(sel), " packages recorded)")
missing <- setdiff(key_pkgs, sel$Package)
if (length(missing))
  message("      note: not installed in this session: ",
          paste(missing, collapse = ", "))

message("\nDone. Environment captured.")

## Optional, advanced: a full lockfile with renv
## (leaves an renv/ folder and renv.lock; only if you want dependency pinning)
##   if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
##   renv::init(bare = TRUE); renv::snapshot(prompt = FALSE)
