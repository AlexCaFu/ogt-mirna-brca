# =====================================================================
# R/01_multimir.R    LAYER 1: target-database provenance
# Records validated and predicted support for each microRNA -> enzyme pair.
# This does not change the tested pairs; it documents the evidence behind them.
# Network: multiMiR web service.
# =====================================================================
source("R/config.R")

l1 <- bind_rows(lapply(CANDIDATES$mirna, function(m) {
  pr <- tryCatch(multiMiR::get_multimir(org = "hsa", mirna = m, target = GENES,
                                        table = "predicted", summary = TRUE)@data,
                 error = function(e) NULL)
  va <- tryCatch(multiMiR::get_multimir(org = "hsa", mirna = m, target = GENES,
                                        table = "validated", summary = TRUE)@data,
                 error = function(e) NULL)
  tibble(
    mirna = m,
    n_pred_db       = if (!is.null(pr) && nrow(pr)) length(unique(pr$database)) else 0L,
    validated_in    = if (!is.null(va) && nrow(va)) paste(sort(unique(va$database)), collapse = ";") else "",
    validated_pmids = if (!is.null(va) && nrow(va)) paste(unique(va$pubmed_id), collapse = ";") else ""
  )
}))
l1$call <- ifelse(l1$validated_in != "", "VALIDATED",
                  ifelse(l1$n_pred_db >= 3, "PREDICTED>=3DB", "WEAK<3DB"))
write_csv(l1, "results/01_target_consensus.csv")
message("Layer 1 complete -> results/01_target_consensus.csv")
