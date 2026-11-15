# =====================================================================
# R/01_multimir.R
# Compile predicted and experimentally validated microRNA -> enzyme
# interactions for the candidate set, using the multiMiR package
# (which aggregates TargetScan and other prediction and validation
# databases). Records, for each candidate pair, whether the interaction
# is validated and/or predicted. Output: results/01_targets.csv
#
# NOTE: reconstructed from the documented method. multiMiR queries live
# databases, so counts may drift with database updates; the archived
# results/ table is the record of the version used.
# =====================================================================
source("R/config.R")
source("R/00_setup.R")

enzyme_of_interest <- unique(CANDIDATES$enzyme)
mirs <- unique(CANDIDATES$mirna)

# Pull validated and predicted target tables for the candidate miRNAs.
val <- tryCatch(
  multiMiR::get_multimir(org = "hsa", mirna = mirs, table = "validated",
                         summary = FALSE)@data,
  error = function(e) { message("validated query failed: ", conditionMessage(e)); NULL })
prd <- tryCatch(
  multiMiR::get_multimir(org = "hsa", mirna = mirs, table = "predicted",
                         summary = FALSE)@data,
  error = function(e) { message("predicted query failed: ", conditionMessage(e)); NULL })

hit <- function(tbl, mir, gene) {
  if (is.null(tbl) || !nrow(tbl)) return(FALSE)
  g <- toupper(tbl$target_symbol)
  m <- tolower(tbl$mature_mirna_id)
  aliases <- toupper(unique(c(gene, names(GENE_ALIASES)[GENE_ALIASES == gene], unname(GENE_ALIASES[gene]))))
  any(tolower(m) == tolower(mir) & g %in% aliases)
}

targets <- do.call(rbind, lapply(seq_len(nrow(CANDIDATES)), function(i) {
  mir <- CANDIDATES$mirna[i]; enz <- CANDIDATES$enzyme[i]
  data.frame(
    mirna       = mir,
    enzyme      = enz,
    drug        = CANDIDATES$drug[i],
    validated   = hit(val, mir, enz),
    predicted   = hit(prd, mir, enz),
    stringsAsFactors = FALSE)
}))
targets$interaction_status <- ifelse(targets$validated, "Validated",
                              ifelse(targets$predicted, "Predicted", "Not recovered"))

readr::write_csv(targets, "results/01_targets.csv")
message("Wrote results/01_targets.csv (", nrow(targets), " candidate pairs)")
print(targets)
