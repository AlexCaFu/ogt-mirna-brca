#!/usr/bin/env Rscript
# =====================================================================
# run_all.R  -- master reproduction script
# Runs the full pipeline end to end, in order, from the repository root.
#
#   setwd("/path/to/ogt-mirna-brca")
#   source("run_all.R")
#
# The first run downloads several GB from GDC and GEO and may take a while;
# later runs reuse the local cache. After it completes, compare the numbers
# with the result tables in results/. Scripts 01/02/03/04 are reconstructions
# of the documented method; 06 and the tests are recovered from the original
# analysis session.
# =====================================================================

t0 <- Sys.time()

pipeline <- c(
  "R/00_setup.R",                     # libraries, parameters, helper functions
  "R/01_multimir.R",                  # seed-based miRNA -> enzyme target compilation
  "R/02_tcga_correlation_survival.R", # TCGA-BRCA correlation + survival (feeds Table 2 + Fig 2)
  "R/03_geo_foldchange.R",            # GEO resistant-model differential expression
  "R/04_assemble_table2.R",           # assemble Table 2 from results
  "R/06_robustness.R"                 # partial corr, stratification, CI, tumor/normal, Cox
)

run_one <- function(path) {
  if (!file.exists(path)) { message(">> SKIP (not found): ", path); return(invisible(FALSE)) }
  message("\n==================== RUN: ", path, " ====================")
  source(path, echo = TRUE, max.deparse.length = Inf)
  invisible(TRUE)
}
invisible(lapply(pipeline, run_one))

## ---- validation tests (synthetic data; should all pass) ----
tests <- c("tests/smoke_test_layer4.R", "tests/test_layer4b.R", "tests/test_layer6.R")
message("\n-------------------- VALIDATION TESTS --------------------")
invisible(lapply(tests, run_one))

message("\nPipeline finished in ",
        round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
        " min. Outputs in results/, tables/, and figures/.")
message("Confirm values against VERIFICATION.md, then run make_repro_metadata.R and close R.")
