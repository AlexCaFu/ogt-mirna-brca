# =====================================================================
# R/03_geo_inspect.R    LAYER 4a: INSPECT (download and audit, no fold changes)
# For each GEO series, download it and write out the platform, the full sample
# phenotype table, and the probe-annotation column names, plus the keyword-based
# GUESS of resistant vs sensitive. You then confirm or correct the assignment
# in R/04_geo_foldchange.R before any fold change is computed.
# Network: NCBI GEO.
# =====================================================================
source("R/config.R")
suppressMessages(library(GEOquery))

geo_objects <- list()
for (nm in names(GEO_SERIES)) {
  s <- GEO_SERIES[[nm]]
  message("\n=====================================================")
  message(sprintf("Series %s  |  %s  |  %s-channel", s$accession, s$measures, s$channels))
  message(s$description)
  g <- tryCatch(getGEO(s$accession, GSEMatrix = TRUE, getGPL = TRUE),
                error = function(e) { message("  download failed: ", conditionMessage(e)); NULL })
  if (is.null(g)) next
  eset <- if (length(g) > 1) g[[1]] else g[[1]]
  geo_objects[[nm]] <- eset

  ph <- pData(eset)
  # descriptive columns most likely to carry the sensitive/resistant label
  desc_cols <- grep("title|characteristic|source|description|label|group",
                    colnames(ph), ignore.case = TRUE, value = TRUE)
  samp <- data.frame(geo_accession = rownames(ph), ph[, desc_cols, drop = FALSE],
                     check.names = FALSE, row.names = NULL)

  # keyword guess (to be verified, never trusted blindly)
  txt <- tolower(apply(samp[, -1, drop = FALSE], 1, paste, collapse = " | "))
  guess <- ifelse(grepl(s$resistant_hint, txt), "resistant",
                  ifelse(grepl(s$sensitive_hint, txt), "sensitive", "UNASSIGNED"))
  samp$group_guess <- guess

  write.csv(samp, sprintf("results/geo_%s_samples.csv", nm), row.names = FALSE)
  message(sprintf("  %d samples; platform %s; %d probes",
                  ncol(eset), annotation(eset), nrow(eset)))
  message("  sample table written -> results/geo_", nm, "_samples.csv")
  message("  keyword group guess:"); print(table(guess))
  message("  probe-annotation columns available for mapping:")
  print(utils::head(colnames(fData(eset)), 30))
}
saveRDS(geo_objects, "results/geo_objects.rds")
message("\nLayer 4a (inspect) complete. Open the results/geo_*_samples.csv files and\n",
        "confirm the group_guess column before running R/04_geo_foldchange.R.")
