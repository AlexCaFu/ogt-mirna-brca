# =====================================================================
# R/config.R
# Central configuration for the four-layer integrative analysis.
# Nothing here is a result; these are the inputs and the documented,
# auditable definitions of every data access.
# =====================================================================

## ---- Candidate microRNA -> enzyme pairs (the hypotheses being tested) ----
## target_claim uses ";" for a microRNA claimed to target more than one enzyme.
## Each microRNA-to-enzyme pair becomes one tested row.
suppressMessages({library(dplyr); library(tibble); library(readr)})

CANDIDATES <- tribble(
  ~mirna,            ~target_claim, ~drug,         ~model,
  "hsa-miR-200c-3p", "OGA;EOGT",    "Tamoxifen",   "MCF-7",
  "hsa-miR-148a-3p", "OGT",         "Tamoxifen",   "MCF-7",
  "hsa-miR-152-3p",  "EOGT",        "Tamoxifen",   "MCF-7",
  "hsa-miR-24-3p",   "OGT",         "Tamoxifen",   "MCF-7",
  "hsa-miR-222-3p",  "OGA;EOGT",    "Doxorubicin", "MCF-7",
  "hsa-miR-7-5p",    "OGT;OGA",     "Doxorubicin", "MCF-7",
  "hsa-miR-181b-5p", "OGA",         "Doxorubicin", "MDA-MB-231",
  "hsa-miR-181d-5p", "OGT",         "Doxorubicin", "(ovarian)",
  "hsa-miR-195-5p",  "OGT",         "Doxorubicin", "MCF-7"
)

GENES <- c("OGT", "OGA", "EOGT")
# OGA is the current HGNC symbol for the gene formerly named MGEA5; both are
# accepted when resolving against a platform annotation.
GENE_ALIASES <- c(OGA = "MGEA5", MGEA5 = "OGA")

## ---- TCGA access (Layers 2 and 3) ----
TCGA_PROJECT <- "TCGA-BRCA"
# Gene expression: STAR - Counts workflow, TPM assay (tpm_unstrand).
# microRNA expression: reads_per_million per stem-loop.
# Samples: primary tumor only (barcode sample-type code "01").

## ---- GEO access (Layer 4), documented from the series records ----
## IMPORTANT: the "resistant" and "sensitive" keyword hints below are only a
## first guess. The pipeline prints the full sample table for each series so
## the assignment is confirmed against the real metadata before any fold
## change is computed. Fold change convention throughout: log2(resistant) minus
## log2(sensitive); a positive value means higher in the resistant state.
GEO_SERIES <- list(
  GSE26459 = list(
    accession   = "GSE26459",
    measures    = "mRNA",          # enzyme fold changes (OGT, OGA, EOGT)
    channels    = "single",        # single-channel array
    description = "Subcloned MCF-7, tamoxifen-sensitive vs resistant (Gonzalez-Malerva 2011)",
    resistant_hint = "resist|/tr|tamr|tam-r|tam r",
    sensitive_hint = "sensitive|parent|wild|wt|mcf-?7$|mcf-?7 "
  ),
  GSE28267 = list(
    accession   = "GSE28267",
    measures    = "miRNA",         # microRNA fold changes
    channels    = "two",           # TWO-COLOR miRCURY LNA (Hy3/Hy5 log-ratios)
    description = "MCF-7 sensitive vs LY2 resistant, two-color miRCURY LNA (miRbase 11.0)",
    resistant_hint = "ly2|resist",
    sensitive_hint = "mcf-?7|sensitive|parent"
  ),
  GSE66607 = list(
    accession   = "GSE66607",
    measures    = "miRNA",         # microRNA fold changes (second model)
    channels    = "single",        # Affymetrix GeneChip miRNA 3.0
    description = "MCF7 parent (MCF7-Pa) vs tamoxifen-resistant (MCF7-Re), Affymetrix miRNA 3.0",
    resistant_hint = "re$|resist|mcf7-re|mcf-7-re",
    sensitive_hint = "pa$|parent|mcf7-pa|mcf-7-pa|sensitive"
  )
)

## ---- Parameters ----
PARAMS <- list(
  cor_method    = "spearman",
  fdr_method    = "BH",
  fdr_threshold = 0.05,
  pseudo        = 1,          # log2(x + 1)
  seed          = 1234
)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
set.seed(PARAMS$seed)
