# =====================================================================
# R/00_setup.R
# Libraries and helper functions used throughout the pipeline.
#
# The statistical and string helpers below are the ones used to produce
# the reported numbers (recovered from the analysis session); the small
# barcode / transform helpers are defined explicitly here.
# =====================================================================

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(multiMiR)
  library(limma)
  library(GEOquery)
  library(Biobase)
  library(survival)
  library(readr)
  library(dplyr)
})

## ---------- barcode / transform helpers ----------
# TCGA patient barcode = first 12 characters (TCGA-XX-XXXX).
patient_id <- function(x) substr(as.character(x), 1, 12)

# log2 transform with a pseudo-count.
log2p <- function(x, pseudo = 1) log2(x + pseudo)

# Map a mature miRNA name to its precursor stem, to match the miRNA_ID
# rows of the TCGA miRNA matrix (e.g. hsa-miR-200c-3p -> hsa-mir-200c).
stem_of <- function(x) {
  x <- tolower(x)
  x <- gsub("-[35]p$", "", x)
  x
}

## ---------- microRNA / gene name canonicalization ----------
canon_mir <- function(x) {
  x <- tolower(x)
  x <- gsub("hsa-", "", x)
  x <- gsub("^mir-?", "", x)
  x <- gsub("-[35]p$", "", x)
  x <- gsub("\\*$", "", x)
  x <- gsub("-[0-9]+$", "", x)
  trimws(x)
}

canon_tokens <- function(s) {
  toks <- unlist(strsplit(as.character(s), "[ ,;/|]+")); toks <- toks[nzchar(toks)]
  unique(canon_mir(toks))
}

gene_tokens <- function(s) {
  toks <- unlist(strsplit(as.character(s), "\\s*///\\s*|[;,]"))
  tolower(trimws(toks[nzchar(toks)]))
}

# Map probe/feature rows of a GEO ExpressionSet to gene or miRNA names.
feature_map <- function(eset, kind = c("gene", "mirna")) {
  kind <- match.arg(kind); fd <- fData(eset); cn <- colnames(fd)
  pick <- function(patterns) {
    for (p in patterns) { h <- grep(p, cn, ignore.case = TRUE, value = TRUE); if (length(h)) return(h[1]) }
    NA_character_ }
  col <- if (kind == "gene") pick(c("^gene.?symbol$", "gene.?symbol", "^symbol$"))
         else pick(c("mirna_id", "mirna", "mir_id", "transcript.?id\\(", "transcript.?id", "^name$"))
  if (is.na(col)) return(setNames(rownames(fd), rownames(fd)))
  setNames(as.character(fd[[col]]), rownames(fd))
}

# Log2-transform a GEO expression matrix only if it appears linear-scale.
to_log2 <- function(ex) if (max(ex, na.rm = TRUE) > 100) log2(ex + 1) else ex

## ---------- statistical helpers (unit-tested in tests/test_layer6.R) ----------
# Partial Spearman via the residual-rank method: rank-transform x and y,
# regress each on the covariates (factors allowed), correlate the residuals.
partial_spearman <- function(x, y, cov) {
  d <- data.frame(x = x, y = y, cov, check.names = FALSE)
  d <- d[stats::complete.cases(d), , drop = FALSE]
  if (nrow(d) < 20) return(list(rho = NA, p = NA, n = nrow(d)))
  d$rx <- rank(d$x); d$ry <- rank(d$y)
  cn <- setdiff(colnames(cov), c("x", "y"))
  rhs <- paste(sprintf("`%s`", cn), collapse = " + ")
  ex <- residuals(lm(as.formula(paste("rx ~", rhs)), data = d))
  ey <- residuals(lm(as.formula(paste("ry ~", rhs)), data = d))
  ct <- suppressWarnings(cor.test(ex, ey, method = "pearson"))
  list(rho = unname(ct$estimate), p = ct$p.value, n = nrow(d))
}

# Fisher-z CI for Spearman with the Bonett-Wright standard error.
spearman_ci_fisher <- function(rho, n, conf = 0.95) {
  se <- sqrt((1 + rho^2 / 2) / (n - 3)); z <- atanh(rho); zc <- qnorm(1 - (1 - conf) / 2)
  c(lower = tanh(z - zc * se), upper = tanh(z + zc * se))
}

# Bootstrap CI for Spearman.
boot_spearman_ci <- function(x, y, B = 2000, conf = 0.95) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  rhos <- replicate(B, { i <- sample.int(n, n, replace = TRUE)
                         suppressWarnings(cor(x[i], y[i], method = "spearman")) })
  stats::quantile(rhos, c((1 - conf) / 2, 1 - (1 - conf) / 2), na.rm = TRUE)
}

# Benjamini-Hochberg helper.
bh <- function(p) p.adjust(p, method = "BH")

if (!dir.exists("results")) dir.create("results")
if (!dir.exists("figures")) dir.create("figures")
if (!dir.exists("tables"))  dir.create("tables")
message("Setup complete: libraries loaded, helpers defined.")
