#!/usr/bin/env Rscript
# =====================================================================
# make_repro_metadata.R
# Captures the computational environment and generates an integrity
# manifest so the deposited version is fully reproducible and verifiable.
#
# RUN THIS from the repository root, ideally in the SAME R environment
# that ran the analysis (so the recorded package versions are the real
# ones). Uses base R only; no packages are installed.
#
#   setwd("/path/to/ogt-mirna-brca")
#   source("make_repro_metadata.R")
#
# After it finishes you can close R. The remaining steps (git, Zenodo)
# are in DEPOSIT_STEPS.md and do not require R.
# =====================================================================

root <- getwd()
if (!dir.exists("results")) dir.create("results")
message("Repository root: ", root)

## 1. Full session information -----------------------------------------
si_path <- file.path("results", "sessionInfo.txt")
writeLines(capture.output(sessionInfo()), si_path)
message("[1/3] Wrote ", si_path)

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
message("[2/3] Wrote ", pv_path, " (", nrow(sel), " packages recorded)")
missing <- setdiff(key_pkgs, sel$Package)
if (length(missing))
  message("      note: not installed in this session: ",
          paste(missing, collapse = ", "))

## 3. Checksum manifest (MD5, base R) ----------------------------------
search_dirs <- c("R", "tests", "results", "figures", "data")
search_dirs <- search_dirs[dir.exists(search_dirs)]
files <- list.files(search_dirs, recursive = TRUE, full.names = TRUE)
files <- files[file.exists(files) & !dir.exists(files)]
# never checksum ignored/large content or the manifest itself
drop <- grepl("GDCdata|/geo/|/raw/|/cache/|\\.DS_Store|MANIFEST\\.md5", files)
files <- files[!drop]

if (length(files)) {
  md5 <- tools::md5sum(files)
  man <- data.frame(md5 = unname(md5), file = names(md5),
                    stringsAsFactors = FALSE)
  man <- man[order(man$file), ]
  writeLines(sprintf("%s  %s", man$md5, man$file), "MANIFEST.md5")
  utils::write.csv(man, "MANIFEST.md5.csv", row.names = FALSE)
  message("[3/3] Wrote MANIFEST.md5 (", nrow(man), " files hashed)")
} else {
  message("[3/3] No files found to hash yet. Add your R/, tests/, and ",
          "results/ files, then re-run this script.")
}

message("\nDone. Environment captured and manifest written.")
message("You can now close R. Continue with DEPOSIT_STEPS.md (git + Zenodo).")

## Optional, advanced: a full lockfile with renv
## (leaves an renv/ folder and renv.lock; only if you want dependency pinning)
##   if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
##   renv::init(bare = TRUE); renv::snapshot(prompt = FALSE)
