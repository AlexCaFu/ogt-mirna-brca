# =====================================================================
# R/02_tcga_correlation_survival.R
# LAYER 2: Spearman correlation of each candidate microRNA against its target
#          enzyme across paired TCGA-BRCA primary tumors, Benjamini-Hochberg FDR.
# LAYER 3: Cox overall-survival hazard ratio per enzyme and per microRNA
#          (continuous predictor), with a median-split log-rank p for reference.
# Network: NIH Genomic Data Commons (large download; cached after first run).
# =====================================================================
source("R/config.R")

## ---------- Download gene expression (STAR - Counts, TPM) ----------
message("TCGA-BRCA gene expression download (large)...")
qg <- GDCquery(project = TCGA_PROJECT, data.category = "Transcriptome Profiling",
               data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")
GDCdownload(qg); se <- GDCprepare(qg)
gx <- assay(se, "tpm_unstrand"); rownames(gx) <- rowData(se)$gene_name
gx <- gx[, is_primary(colnames(gx)), drop = FALSE]
colnames(gx) <- patient_id(colnames(gx))
gx <- gx[, !duplicated(colnames(gx)), drop = FALSE]      # one aliquot per patient
gx <- gx[!duplicated(rownames(gx)), , drop = FALSE]      # one probe per gene symbol
gx <- log2p(gx, PARAMS$pseudo)                           # (irrelevant to Spearman; kept for clarity)

resolve <- function(s) {
  if (s %in% rownames(gx)) return(s)
  alt <- GENE_ALIASES[s]
  if (!is.na(alt) && alt %in% rownames(gx)) return(unname(alt))
  NA_character_
}

## ---------- Download microRNA expression (reads per million) ----------
message("TCGA-BRCA microRNA download...")
qm <- GDCquery(project = TCGA_PROJECT, data.category = "Transcriptome Profiling",
               data.type = "miRNA Expression Quantification")
GDCdownload(qm); mi <- GDCprepare(qm)
rc <- grep("reads_per_million", colnames(mi), value = TRUE)
mraw <- as.matrix(mi[, rc]); rownames(mraw) <- mi$miRNA_ID
full <- sub("reads_per_million_miRNA_mapped_", "", rc)
keep <- is_primary(full)
mraw <- mraw[, keep, drop = FALSE]; colnames(mraw) <- patient_id(full[keep])
mraw <- mraw[, !duplicated(colnames(mraw)), drop = FALSE]

## Expression vector for one mature microRNA: sum reads-per-million across its
## stem-loop loci, then log2. Regex is anchored so hsa-mir-24 matches
## hsa-mir-24-1 and hsa-mir-24-2 but not hsa-mir-240.
mir_vec <- function(mature) {
  s  <- stem_of(mature)
  rx <- paste0("^", s, "(-[0-9]+)?$")
  h  <- grep(rx, rownames(mraw), value = TRUE)
  if (!length(h)) return(NULL)
  v <- if (length(h) == 1) mraw[h, ] else colSums(mraw[h, , drop = FALSE])
  log2p(v, PARAMS$pseudo)
}

## ---------- LAYER 2: correlation ----------
l2 <- list()
for (i in seq_len(nrow(CANDIDATES))) {
  mv <- mir_vec(CANDIDATES$mirna[i]); if (is.null(mv)) next
  for (g in unique(strsplit(CANDIDATES$target_claim[i], ";")[[1]])) {
    rg <- resolve(g); if (is.na(rg)) next
    com <- intersect(names(mv), colnames(gx))
    ct  <- suppressWarnings(cor.test(mv[com], gx[rg, com], method = PARAMS$cor_method))
    l2[[paste(CANDIDATES$mirna[i], g)]] <- tibble(
      mirna = CANDIDATES$mirna[i], target = g, symbol_used = rg,
      rho = unname(ct$estimate), p = ct$p.value, n = length(com))
  }
}
l2 <- bind_rows(l2)
l2$fdr <- p.adjust(l2$p, PARAMS$fdr_method)
l2$inverse_sig <- l2$rho < 0 & l2$fdr < PARAMS$fdr_threshold
write_csv(l2, "results/02_tcga_correlation.csv")
message(sprintf("Layer 2 complete: %d pairs, n paired = %d, %d significant inverse -> results/02_tcga_correlation.csv",
                nrow(l2), l2$n[1], sum(l2$inverse_sig)))

## ---------- LAYER 3: survival ----------
cl <- GDCquery_clinic(TCGA_PROJECT, "clinical")
cl <- transmute(cl, id = bcr_patient_barcode,
                time  = ifelse(is.na(days_to_death), days_to_last_follow_up, days_to_death),
                event = ifelse(vital_status == "Dead", 1, 0))
cl <- filter(cl, !is.na(time), time > 0)

surv_one <- function(vec, label) {
  df <- inner_join(tibble(id = names(vec), expr = as.numeric(vec)), cl, by = "id")
  df <- filter(df, !is.na(expr)); if (nrow(df) < 20) return(NULL)
  cox <- summary(coxph(Surv(time, event) ~ expr, data = df))$coefficients
  df$grp <- ifelse(df$expr >= median(df$expr, na.rm = TRUE), "high", "low")
  lr <- tryCatch(survminer::surv_pvalue(survfit(Surv(time, event) ~ grp, data = df), df)$pval,
                 error = function(e) NA)
  tibble(feature = label, HR = exp(cox[1, "coef"]), HR_p = cox[1, "Pr(>|z|)"],
         logrank_p = lr, n = nrow(df))
}
s3 <- list()
for (g in GENES) { rg <- resolve(g); if (!is.na(rg)) s3[[g]] <- surv_one(gx[rg, ], paste0(g, "(", rg, ")")) }
for (i in seq_len(nrow(CANDIDATES))) {
  mv <- mir_vec(CANDIDATES$mirna[i]); if (!is.null(mv)) s3[[CANDIDATES$mirna[i]]] <- surv_one(mv, CANDIDATES$mirna[i])
}
s3 <- bind_rows(s3); write_csv(s3, "results/03_survival.csv")
message("Layer 3 complete -> results/03_survival.csv")

saveRDS(list(gx = gx, mir_vec = mir_vec, resolve = resolve), "results/tcga_objects.rds")
