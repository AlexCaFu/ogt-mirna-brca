# =====================================================================
# R/02_tcga_correlation_survival.R
# Patient-level evaluation of the candidate interactions in TCGA-BRCA:
#   - Spearman correlation between each microRNA and its target enzyme,
#     with Benjamini-Hochberg correction across all tested pairs.
#   - Cox proportional-hazards regression for each microRNA and enzyme
#     (overall survival), with BH correction.
# Outputs: results/02_tcga_correlation.csv, results/02_tcga_survival.csv
#
# NOTE: reconstructed from the documented method. The data loading mirrors
# 06_robustness.R. Run and confirm the values match VERIFICATION.md; the
# authoritative values are the archived results/ tables.
# =====================================================================
source("R/config.R")
source("R/00_setup.R")
set.seed(PARAMS$seed)

## ---------- load TCGA-BRCA gene + miRNA matrices ----------
qg <- GDCquery(project = TCGA_PROJECT, data.category = "Transcriptome Profiling",
               data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")
GDCdownload(qg); se <- GDCprepare(qg)
tpm <- assay(se, "tpm_unstrand"); rownames(tpm) <- rowData(se)$gene_name
tpm <- tpm[!duplicated(rownames(tpm)), , drop = FALSE]
type_all <- substr(colnames(tpm), 14, 15)
gx <- tpm[, type_all == "01", drop = FALSE]
colnames(gx) <- patient_id(colnames(gx)); gx <- gx[, !duplicated(colnames(gx)), drop = FALSE]
gx <- log2p(gx, PARAMS$pseudo)

qm <- GDCquery(project = TCGA_PROJECT, data.category = "Transcriptome Profiling",
               data.type = "miRNA Expression Quantification")
GDCdownload(qm); mi <- GDCprepare(qm)
rc <- grep("reads_per_million", colnames(mi), value = TRUE)
mraw <- as.matrix(mi[, rc]); rownames(mraw) <- mi$miRNA_ID
full <- sub("reads_per_million_miRNA_mapped_", "", rc); keep <- substr(full, 14, 15) == "01"
mraw <- mraw[, keep, drop = FALSE]; colnames(mraw) <- patient_id(full[keep])
mraw <- mraw[, !duplicated(colnames(mraw)), drop = FALSE]

mir_vec2 <- function(mature) {
  h <- grep(paste0("^", stem_of(mature), "(-[0-9]+)?$"), rownames(mraw), value = TRUE)
  if (!length(h)) return(NULL)
  log2p(if (length(h) == 1) mraw[h, ] else colSums(mraw[h, , drop = FALSE]), PARAMS$pseudo)
}

# Resolve an enzyme's row in the gene matrix, allowing symbol aliases.
resolve_gene <- function(sym) {
  cands <- toupper(unique(c(sym, names(GENE_ALIASES)[GENE_ALIASES == sym], unname(GENE_ALIASES[sym]))))
  hit <- rownames(gx)[toupper(rownames(gx)) %in% cands]
  if (length(hit)) hit[1] else NA_character_
}

## ---------- A. per-pair Spearman correlation ----------
corr <- do.call(rbind, lapply(seq_len(nrow(CANDIDATES)), function(i) {
  mir <- CANDIDATES$mirna[i]; enz <- CANDIDATES$enzyme[i]
  gsym <- resolve_gene(enz); v <- mir_vec2(mir)
  if (is.na(gsym) || is.null(v)) return(data.frame(mirna = mir, enzyme = enz,
        rho = NA, p = NA, n = NA))
  com <- intersect(names(v), colnames(gx))
  ct <- suppressWarnings(cor.test(v[com], gx[gsym, com], method = "spearman"))
  data.frame(mirna = mir, enzyme = enz, rho = unname(ct$estimate),
             p = ct$p.value, n = length(com))
}))
corr$FDR <- bh(corr$p)
corr$supported <- with(corr, is.finite(rho) & rho < 0 & FDR < FDR_THRESHOLD)
corr$rho <- round(corr$rho, 3); corr$p <- signif(corr$p, 3); corr$FDR <- signif(corr$FDR, 3)
readr::write_csv(corr, "results/02_tcga_correlation.csv")
message("Wrote results/02_tcga_correlation.csv"); print(corr)

## ---------- B. Cox proportional hazards (overall survival) ----------
cl <- GDCquery_clinic(TCGA_PROJECT, "clinical")
cl <- data.frame(patient = cl$bcr_patient_barcode,
                 time = ifelse(is.na(cl$days_to_death), cl$days_to_last_follow_up, cl$days_to_death),
                 event = ifelse(cl$vital_status == "Dead", 1, 0))
cl <- cl[is.finite(cl$time) & cl$time > 0, ]

cox_feature <- function(vec, name) {
  df <- merge(data.frame(patient = names(vec), v = as.numeric(vec)), cl, by = "patient")
  df <- df[stats::complete.cases(df), ]
  if (nrow(df) < 30) return(data.frame(feature = name, HR = NA, p = NA, n = nrow(df)))
  fit <- summary(survival::coxph(survival::Surv(time, event) ~ v, data = df))
  data.frame(feature = name, HR = unname(exp(fit$coefficients["v", "coef"])),
             p = unname(fit$coefficients["v", "Pr(>|z|)"]), n = nrow(df))
}

feat_list <- c(
  setNames(lapply(unique(CANDIDATES$mirna), mir_vec2), unique(CANDIDATES$mirna)),
  setNames(lapply(unique(CANDIDATES$enzyme), function(e) {
    g <- resolve_gene(e); if (is.na(g)) return(NULL)
    setNames(gx[g, ], colnames(gx)) }), unique(CANDIDATES$enzyme))
)
surv <- do.call(rbind, lapply(names(feat_list), function(nm) {
  v <- feat_list[[nm]]; if (is.null(v)) return(NULL); cox_feature(v, nm) }))
surv$FDR <- bh(surv$p)
surv$HR <- round(surv$HR, 3); surv$p <- signif(surv$p, 3); surv$FDR <- signif(surv$FDR, 3)
readr::write_csv(surv, "results/02_tcga_survival.csv")
message("Wrote results/02_tcga_survival.csv"); print(surv)

message("\nTCGA correlation and survival complete.")
