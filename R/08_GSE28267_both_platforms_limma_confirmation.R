# GSE28267, both array batches, named-arm probes selected by miRBase accession. limma in R.
# Run on 15 November 2026; output saved as results/08b_GSE28267_limma_R_confirmation.csv.
suppressMessages({library(GEOquery); library(limma)})
arm <- c("miR-200c-3p"="MIMAT0000617","miR-200c-5p"="MIMAT0004657","miR-148a-3p"="MIMAT0000243",
         "miR-152-3p"="MIMAT0000438","miR-24-3p"="MIMAT0000080","miR-222-3p"="MIMAT0000279",
         "miR-7-5p"="MIMAT0000252","miR-181b-5p"="MIMAT0000257","miR-181d-5p"="MIMAT0002821",
         "miR-195-5p"="MIMAT0000461")
gl <- getGEO("GSE28267", GSEMatrix = TRUE, getGPL = TRUE)
out <- list()
for (k in seq_along(gl)) {
  es <- gl[[k]]; plat <- annotation(es); ph <- pData(es)
  cc <- grep("characteristics", colnames(ph), value = TRUE)
  txt <- apply(ph[, cc, drop = FALSE], 1, paste, collapse = " | ")
  keep <- !grepl("agent: E2", txt)
  grp <- factor(ifelse(grepl("resistant", txt[keep]), "resistant", "sensitive"),
                levels = c("sensitive", "resistant"))
  cat("\n", plat, ": arrays used =", sum(keep), "\n"); print(table(grp))
  ex  <- exprs(es)[, keep, drop = FALSE]
  fit <- eBayes(lmFit(ex, model.matrix(~ grp)))
  tt  <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  acc <- as.character(fData(es)[rownames(tt), "Accession"])
  for (nm in names(arm)) {
    hit <- which(vapply(strsplit(acc, "/"), function(a) arm[[nm]] %in% a, logical(1)))
    for (i in hit) out[[length(out) + 1]] <- data.frame(
      platform = plat, mirna = nm, probe = rownames(tt)[i],
      log2FC = round(tt$logFC[i], 3), p = signif(tt$P.Value[i], 3),
      adj_p = signif(tt$adj.P.Val[i], 3))
  }
}
res <- do.call(rbind, out)
print(res, row.names = FALSE)
write.csv(res, "results/08b_GSE28267_limma_R_confirmation.csv", row.names = FALSE)
