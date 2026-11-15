# =====================================================================
# R/06_robustness.R
# Robustness / sensitivity analyses for the Perspective, run on the
# TCGA-BRCA data (reads from the local GDC cache; re-downloads if absent).
#
#   A. miR-200c-3p vs EOGT: unadjusted, partial (adjusting for tumor purity,
#      ER status, PAM50 subtype), ER-positive-restricted, and subtype-stratified.
#   B. 95% confidence intervals for the priority correlation (bootstrap + Fisher z).
#   C. Rank of miR-200c-3p among ALL predicted EOGT-targeting microRNAs.
#   D. EOGT expression, tumor vs adjacent normal, in TCGA-BRCA.
#   E. Multivariable Cox for EOGT and miR-200c-3p, adjusting for age and stage.
#
# Every number is computed here; nothing is hard-coded. Outputs go to results/.
# =====================================================================
source("R/config.R")
source("R/00_setup.R")
set.seed(PARAMS$seed)

## ---------- load TCGA-BRCA gene expression ----------
message("Preparing TCGA-BRCA gene expression ...")
qg <- GDCquery(project = TCGA_PROJECT, data.category = "Transcriptome Profiling",
               data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")
GDCdownload(qg); se <- GDCprepare(qg)
tpm <- assay(se, "tpm_unstrand"); rownames(tpm) <- rowData(se)$gene_name
tpm <- tpm[!duplicated(rownames(tpm)), , drop = FALSE]
stopifnot("EOGT" %in% rownames(tpm))
type_all <- substr(colnames(tpm), 14, 15)

## ---------- D. EOGT tumor vs adjacent normal (sample types 01 vs 11) ----------
eogt_all <- log2p(tpm["EOGT", ], PARAMS$pseudo)
tum <- eogt_all[type_all == "01"]; nor <- eogt_all[type_all == "11"]
D <- if (length(nor) >= 3) {
  wt <- suppressWarnings(wilcox.test(tum, nor))
  data.frame(analysis = "EOGT tumor vs adjacent normal (log2 TPM)",
             n_tumor = length(tum), n_normal = length(nor),
             median_tumor = round(median(tum), 3), median_normal = round(median(nor), 3),
             log2FC = round(median(tum) - median(nor), 3), wilcox_p = signif(wt$p.value, 3))
} else data.frame(analysis = "EOGT tumor vs normal", note = "insufficient normal samples")
readr::write_csv(D, "results/06D_EOGT_tumor_vs_normal.csv"); print(D)

## primary-tumor gene matrix, one sample per patient, log2 scale
gx <- tpm[, type_all == "01", drop = FALSE]
colnames(gx) <- patient_id(colnames(gx)); gx <- gx[, !duplicated(colnames(gx)), drop = FALSE]
gx <- log2p(gx, PARAMS$pseudo)

## ---------- load TCGA-BRCA microRNA expression ----------
message("Preparing TCGA-BRCA microRNA ...")
qm <- GDCquery(project = TCGA_PROJECT, data.category = "Transcriptome Profiling",
               data.type = "miRNA Expression Quantification")
GDCdownload(qm); mi <- GDCprepare(qm)
rc <- grep("reads_per_million", colnames(mi), value = TRUE)
mraw <- as.matrix(mi[, rc]); rownames(mraw) <- mi$miRNA_ID
full <- sub("reads_per_million_miRNA_mapped_", "", rc); keep <- substr(full, 14, 15) == "01"
mraw <- mraw[, keep, drop = FALSE]; colnames(mraw) <- patient_id(full[keep])
mraw <- mraw[, !duplicated(colnames(mraw)), drop = FALSE]

# Sum precursors of a mature miRNA, then log2-transform.
mir_vec2 <- function(mature) {
  h <- grep(paste0("^", stem_of(mature), "(-[0-9]+)?$"), rownames(mraw), value = TRUE)
  if (!length(h)) return(NULL)
  log2p(if (length(h) == 1) mraw[h, ] else colSums(mraw[h, , drop = FALSE]), PARAMS$pseudo)
}

## align miR-200c-3p and EOGT
m200 <- mir_vec2("hsa-miR-200c-3p")
com <- intersect(names(m200), colnames(gx)); x <- m200[com]; y <- gx["EOGT", com]
message(sprintf("Paired tumors for miR-200c-3p vs EOGT: %d", length(com)))

## ---------- annotations: purity, ER status, PAM50 subtype ----------
purity <- setNames(rep(NA_real_, length(com)), com)
pu <- tryCatch({ e <- new.env(); utils::data("Tumor.purity", package = "TCGAbiolinks", envir = e)
                 as.data.frame(e$Tumor.purity) }, error = function(e) NULL)
if (!is.null(pu)) {
  cpe_col <- grep("CPE|ESTIMATE", colnames(pu), ignore.case = TRUE, value = TRUE)[1]
  pid <- patient_id(pu$Sample.ID); val <- as.numeric(gsub(",", ".", pu[[cpe_col]]))
  val <- tapply(val, pid, function(v) mean(v, na.rm = TRUE))
  purity[com] <- val[com]
  message(sprintf("Tumor purity matched for %d/%d tumors (%s)", sum(!is.na(purity)), length(com), cpe_col))
}
pam50 <- setNames(rep(NA_character_, length(com)), com); er <- pam50
sb <- tryCatch(as.data.frame(TCGAquery_subtype("brca")), error = function(e) NULL)
if (!is.null(sb)) {
  pcol  <- grep("patient|sample", colnames(sb), ignore.case = TRUE, value = TRUE)[1]
  pam_c <- grep("PAM50|mRNA.*subtype|Subtype_PAM50", colnames(sb), ignore.case = TRUE, value = TRUE)[1]
  er_c  <- grep("^ER|ER.?Status|ER_Status", colnames(sb), ignore.case = TRUE, value = TRUE)[1]
  key <- patient_id(as.character(sb[[pcol]]))
  if (!is.na(pam_c)) pam50[com] <- setNames(as.character(sb[[pam_c]]), key)[com]
  if (!is.na(er_c))  er[com]    <- setNames(as.character(sb[[er_c]]),  key)[com]
  message(sprintf("PAM50 matched %d, ER matched %d", sum(!is.na(pam50)), sum(!is.na(er))))
}
er_bin <- ifelse(grepl("pos", er, ignore.case = TRUE), "Positive",
          ifelse(grepl("neg", er, ignore.case = TRUE), "Negative", NA))

## ---------- A. unadjusted, partial, ER-positive, subtype-stratified ----------
ct0 <- suppressWarnings(cor.test(x, y, method = "spearman"))
cov <- data.frame(purity = purity[com], ER = factor(er_bin[com]), PAM50 = factor(pam50[com]))
cov <- cov[, vapply(cov, function(z) sum(!is.na(z)) > 20 && length(unique(na.omit(z))) > 1, logical(1)), drop = FALSE]
part <- if (ncol(cov)) partial_spearman(x, y, cov) else list(rho = NA, p = NA, n = NA)
er_ok <- er_bin[com] == "Positive" & !is.na(er_bin[com])
ctER <- if (sum(er_ok) >= 20) suppressWarnings(cor.test(x[er_ok], y[er_ok], method = "spearman")) else NULL

A <- data.frame(
  model = c("Unadjusted", paste0("Partial (", paste(colnames(cov), collapse = "+"), ")"), "ER-positive only"),
  rho = round(c(ct0$estimate, part$rho, if (!is.null(ctER)) ctER$estimate else NA), 3),
  p   = signif(c(ct0$p.value, part$p, if (!is.null(ctER)) ctER$p.value else NA), 3),
  n   = c(length(x), part$n, if (!is.null(ctER)) sum(er_ok) else NA))
readr::write_csv(A, "results/06A_partial_correlation.csv"); print(A)

strat <- do.call(rbind, lapply(split(seq_along(com), pam50[com]), function(idx) {
  if (length(idx) < 15) return(NULL)
  ct <- suppressWarnings(cor.test(x[idx], y[idx], method = "spearman"))
  data.frame(subtype = pam50[com][idx][1], n = length(idx),
             rho = round(unname(ct$estimate), 3), p = signif(ct$p.value, 3))
}))
readr::write_csv(strat, "results/06A_subtype_stratified.csv"); print(strat)

## ---------- B. confidence intervals for the priority correlation ----------
fis <- spearman_ci_fisher(unname(ct0$estimate), length(x))
boo <- boot_spearman_ci(x, y, B = 2000)
B <- data.frame(rho = round(unname(ct0$estimate), 3),
                fisher_lower = round(fis["lower"], 3), fisher_upper = round(fis["upper"], 3),
                boot_lower = round(boo[[1]], 3), boot_upper = round(boo[[2]], 3), n = length(x))
readr::write_csv(B, "results/06B_correlation_CI.csv"); print(B)

## ---------- C. rank of miR-200c-3p among predicted EOGT-targeting microRNAs ----------
pred <- tryCatch(multiMiR::get_multimir(org = "hsa", target = "EOGT", table = "predicted",
                                        summary = FALSE)@data, error = function(e) NULL)
if (!is.null(pred) && nrow(pred)) {
  predmirs <- unique(pred$mature_mirna_id)
  rk <- do.call(rbind, lapply(predmirs, function(mm) {
    v <- mir_vec2(mm); if (is.null(v)) return(NULL)
    cc <- intersect(names(v), colnames(gx)); if (length(cc) < 50) return(NULL)
    data.frame(mirna = mm, rho = suppressWarnings(cor(v[cc], gx["EOGT", cc], method = "spearman")))
  }))
  rk <- rk[is.finite(rk$rho), ]; rk <- rk[order(rk$rho), ]; rk$rank <- seq_len(nrow(rk))
  pos <- grep("200c-3p", rk$mirna, ignore.case = TRUE)
  C <- data.frame(n_predicted_EOGT_miRNAs_tested = nrow(rk),
                  miR_200c_3p_rank = if (length(pos)) rk$rank[pos][1] else NA,
                  miR_200c_3p_rho = if (length(pos)) round(rk$rho[pos][1], 3) else NA,
                  percentile = if (length(pos)) round(100 * rk$rank[pos][1] / nrow(rk), 1) else NA)
  readr::write_csv(rk, "results/06C_EOGT_miRNA_ranking_full.csv")
  readr::write_csv(C, "results/06C_miR200c_rank_summary.csv"); print(C)
} else message("multiMiR EOGT prediction unavailable; ranking skipped.")

## ---------- E. multivariable Cox ----------
E <- tryCatch({
  cl <- GDCquery_clinic(TCGA_PROJECT, "clinical")
  cl <- data.frame(patient = cl$bcr_patient_barcode,
                   time = ifelse(is.na(cl$days_to_death), cl$days_to_last_follow_up, cl$days_to_death),
                   event = ifelse(cl$vital_status == "Dead", 1, 0),
                   age = suppressWarnings(as.numeric(cl$age_at_index)),
                   stage = as.character(cl$ajcc_pathologic_stage))
  cl <- cl[is.finite(cl$time) & cl$time > 0, ]
  cl$stage2 <- gsub("[ABC]$", "", gsub("Stage ", "", cl$stage))
  mk <- function(vec, name) {
    df <- merge(data.frame(patient = names(vec), v = as.numeric(vec)), cl, by = "patient")
    df$PAM50 <- pam50[df$patient]
    df <- df[stats::complete.cases(df[, c("time","event","v","age","stage2")]), ]
    fit <- summary(survival::coxph(survival::Surv(time, event) ~ v + age + factor(stage2), data = df))
    data.frame(feature = name, HR_adjusted = round(exp(fit$coefficients["v","coef"]), 3),
               p_adjusted = signif(fit$coefficients["v","Pr(>|z|)"], 3), n = nrow(df))
  }
  rbind(mk(y, "EOGT"), mk(x, "miR-200c-3p"))
}, error = function(e) { message("Cox step skipped: ", conditionMessage(e)); NULL })
if (!is.null(E)) { readr::write_csv(E, "results/06E_multivariable_cox.csv"); print(E) }

writeLines(capture.output(sessionInfo()), "results/sessionInfo_layer6.txt")
message("\nLayer 6 complete. Results in results/06*.csv")
