# =====================================================================
# R/00_setup.R    packages and shared helpers
# =====================================================================
cran_pkgs <- c("dplyr", "tibble", "readr", "tidyr", "stringr", "ggplot2",
               "patchwork", "survival", "survminer")
bioc_pkgs <- c("TCGAbiolinks", "SummarizedExperiment", "GEOquery", "limma", "multiMiR")

if (isTRUE(getOption("ogt.install", FALSE))) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  miss <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
  miss <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) BiocManager::install(miss, update = FALSE, ask = FALSE)
}

suppressMessages({
  for (p in c(cran_pkgs, bioc_pkgs)) {
    if (!requireNamespace(p, quietly = TRUE))
      stop(sprintf("Package '%s' not installed. Run once with options(ogt.install=TRUE).", p))
    library(p, character.only = TRUE)
  }
})

## Reduce a full TCGA barcode to the 12-character patient submitter ID.
patient_id <- function(x) substr(x, 1, 12)
## Primary-tumor test (sample type "01").
is_primary <- function(x) substr(x, 14, 15) == "01"
log2p <- function(x, pseudo = 1) log2(x + pseudo)

## Map a mature miRNA name to its stem-loop base for TCGA lookup.
## NOTE: TCGA "miRNA Expression Quantification" is stem-loop level; the mature
## guide arm (-3p/-5p) is used as the label but the measured quantity is the
## stem-loop reads-per-million (summed across loci). This is stated in the methods.
stem_of <- function(x) tolower(sub("-(3p|5p)$", "", sub("hsa-miR", "hsa-mir", x)))
