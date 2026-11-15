# =====================================================================
# R/config.R
# Central configuration for the pipeline: parameters, data sources, and
# the candidate microRNA -> O-GlcNAc enzyme interactions to be tested.
# Edit values here rather than inside the analysis scripts.
# =====================================================================

## Reproducibility parameters -----------------------------------------
PARAMS <- list(
  seed  = 1,        # controls the bootstrap resampling in 06_robustness.R
  pseudo = 1        # pseudo-count for log2(x + pseudo) transforms
)

## Project / cohort ---------------------------------------------------
TCGA_PROJECT <- "TCGA-BRCA"

## O-GlcNAc cycling enzymes (gene symbols). OGA is annotated as MGEA5 in
## many resources; both symbols are handled via GENE_ALIASES below.
ENZYMES <- c("OGT", "MGEA5", "EOGT")

## Symbol aliases so probe-to-gene matching is robust across platforms.
GENE_ALIASES <- c(
  "OGA"   = "MGEA5",
  "MGEA5" = "MGEA5",
  "NCOAT" = "MGEA5",
  "OGT"   = "OGT",
  "EOGT"  = "EOGT",
  "AER61" = "EOGT"
)

## Candidate miRNA -> enzyme interactions evaluated in the study.
## These are the resistance-associated microRNAs compiled from the
## reviewed literature and screened against the enzyme 3' UTRs; the set
## corresponds to the pairs graded in Table 2.
CANDIDATES <- rbind(
  data.frame(mirna = "hsa-miR-24-3p",   enzyme = "OGT",   drug = "Tamoxifen"),
  data.frame(mirna = "hsa-miR-7-5p",    enzyme = "OGT",   drug = "Doxorubicin"),
  data.frame(mirna = "hsa-miR-200c-3p", enzyme = "EOGT",  drug = "Tamoxifen"),
  data.frame(mirna = "hsa-miR-222-3p",  enzyme = "MGEA5", drug = "Doxorubicin"),
  data.frame(mirna = "hsa-miR-181b-5p", enzyme = "MGEA5", drug = "Doxorubicin"),
  data.frame(mirna = "hsa-miR-181d-5p", enzyme = "OGT",   drug = "Doxorubicin"),
  data.frame(mirna = "hsa-miR-148a-3p", enzyme = "OGT",   drug = "Tamoxifen"),
  data.frame(mirna = "hsa-miR-200c-3p", enzyme = "MGEA5", drug = "Tamoxifen"),
  data.frame(mirna = "hsa-miR-7-5p",    enzyme = "MGEA5", drug = "Doxorubicin"),
  data.frame(mirna = "hsa-miR-195-5p",  enzyme = "OGT",   drug = "Doxorubicin"),
  data.frame(mirna = "hsa-miR-222-3p",  enzyme = "EOGT",  drug = "Doxorubicin"),
  data.frame(mirna = "hsa-miR-152-3p",  enzyme = "EOGT",  drug = "Tamoxifen"),
  stringsAsFactors = FALSE
)

## GEO resistant-model series.
## `orient` records which sample group is the resistant state; the actual
## per-sample assignment is read from results/geo_<name>_samples.csv
## (a two-column file: geo_accession, group_guess in {sensitive,resistant}).
GEO_SERIES <- list(
  GSE26459 = list(measures = "mRNA",
                  description = "Tamoxifen-sensitive vs resistant MCF7 (GPL570)"),
  GSE28267 = list(measures = "miRNA",
                  description = "Tamoxifen-sensitive MCF7 vs resistant LY2"),
  GSE66607 = list(measures = "miRNA",
                  description = "Isogenic MCF7 vs tamoxifen-resistant MCF7 (single replicate)")
)

## False-discovery-rate threshold for calling an interaction supported.
FDR_THRESHOLD <- 0.05
