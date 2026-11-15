# Deaths per univariable Cox model and Schoenfeld residual tests, TCGA-BRCA.
# Run on 15 November 2026 after R/07 (uses the matrix gx); output saved as
# results/09_survival_events_schoenfeld.csv.
suppressMessages({library(TCGAbiolinks); library(survival)})
qm <- GDCquery(project = "TCGA-BRCA", data.category = "Transcriptome Profiling",
               data.type = "miRNA Expression Quantification")
GDCdownload(qm); mi <- GDCprepare(qm)
rc   <- grep("reads_per_million", colnames(mi), value = TRUE)
mraw <- as.matrix(mi[, rc]); rownames(mraw) <- mi$miRNA_ID
full <- sub("reads_per_million_miRNA_mapped_", "", rc)
keep <- substr(full, 14, 15) == "01"
mraw <- mraw[, keep, drop = FALSE]; colnames(mraw) <- substr(full[keep], 1, 12)
mraw <- mraw[, !duplicated(colnames(mraw)), drop = FALSE]
mir_vec <- function(stem) {
  h <- grep(paste0("^", stem, "(-[0-9]+)?$"), rownames(mraw), value = TRUE)
  v <- if (length(h) == 1) mraw[h, ] else colSums(mraw[h, , drop = FALSE])
  log2(v + 1)
}
cl0 <- GDCquery_clinic("TCGA-BRCA", "clinical")
idc <- if ("bcr_patient_barcode" %in% names(cl0)) "bcr_patient_barcode" else "submitter_id"
cl  <- data.frame(id = cl0[[idc]],
                  time  = ifelse(is.na(cl0$days_to_death), cl0$days_to_last_follow_up, cl0$days_to_death),
                  event = ifelse(cl0$vital_status == "Dead", 1, 0))
cl  <- cl[!is.na(cl$time) & cl$time > 0, ]
one <- function(vec, label) {
  df <- merge(data.frame(id = names(vec), expr = as.numeric(vec)), cl, by = "id")
  df <- df[!is.na(df$expr), ]
  fit <- coxph(Surv(time, event) ~ expr, data = df)
  s <- summary(fit)$coefficients
  data.frame(feature = label, n = nrow(df), deaths = sum(df$event),
             HR = round(exp(s[1, "coef"]), 3), p = signif(s[1, "Pr(>|z|)"], 3),
             schoenfeld_p = signif(cox.zph(fit)$table[1, "p"], 3))
}
stems <- c("miR-200c-3p"="hsa-mir-200c","miR-148a-3p"="hsa-mir-148a","miR-152-3p"="hsa-mir-152",
           "miR-24-3p"="hsa-mir-24","miR-222-3p"="hsa-mir-222","miR-7-5p"="hsa-mir-7",
           "miR-181b-5p"="hsa-mir-181b","miR-181d-5p"="hsa-mir-181d","miR-195-5p"="hsa-mir-195")
out <- do.call(rbind, c(
  lapply(c("OGT","OGA","EOGT"), function(g) one(log2(gx[g, ] + 1), g)),
  lapply(names(stems), function(m) one(mir_vec(stems[[m]]), m))))
print(out, row.names = FALSE)
write.csv(out, "results/09_survival_events_schoenfeld.csv", row.names = FALSE)
