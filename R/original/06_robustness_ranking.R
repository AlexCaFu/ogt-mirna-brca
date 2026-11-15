# =====================================================================
# R/06_robustness_ranking.R    LAYER 5: Robustness, CIs, and ranking
#
# Analyses:
#  1. Bootstrap 95% CIs on all 12 candidate Spearman correlations
#  2. Partial Spearman: miR-200c-3p to EOGT controlling for PAM50,
#     ER status, and tumor purity (confirms the correlation is not
#     a subtype or purity artifact)
#  3. Subtype-stratified correlations (ER-positive only, Luminal A,
#     Luminal B) - the clinically relevant subset for tamoxifen
#  4. EOGT expression: primary tumors vs solid-tissue-normal (TCGA-BRCA)
#  5. Rank of miR-200c-3p among ALL miRNAs measured in TCGA-BRCA
#     for inverse correlation with EOGT
#
# Reads:  results/tcga_objects.rds   (from R/02_tcga_correlation_survival.R)
#         results/02_tcga_correlation.csv
# Net:    GDCquery_clinic for subtypes (small); GDC for normal tissue
# Writes: results/06_bootstrap_ci.csv
#         results/06_partial_correlation.csv
#         results/06_stratified_correlation.csv
#         results/06_eogt_normal_tumor.csv
#         results/06_eogt_ranking.csv
#         results/06_summary.txt   (plain-text numbers for the manuscript)
#
# ABSOLUTE RULE: every number written here came from the actual data.
# No number is asserted without computation.
# =====================================================================
source("R/00_setup.R")
source("R/config.R")
set.seed(PARAMS$seed)

## ---- Load TCGA objects ----
message("Loading TCGA objects...")
objs  <- readRDS("results/tcga_objects.rds")
gx    <- objs$gx          # log2(TPM+1), gene x patient
l2    <- readr::read_csv("results/02_tcga_correlation.csv", show_col_types = FALSE)

## Recover full miRNA RPM matrix from the mir_vec closure environment.
## This is preserved by R's serialization of closures.
mraw <- environment(objs$mir_vec)$mraw   # stem-loop x patient, raw RPM

common   <- intersect(colnames(gx), colnames(mraw))
gx_com   <- gx[, common, drop = FALSE]
mraw_com <- mraw[, common, drop = FALSE]
n_paired <- length(common)
EOGT_v   <- as.numeric(gx_com["EOGT", ])
message(sprintf("Paired cohort: %d patients, %d miRNA stem-loops, %d genes",
                n_paired, nrow(mraw_com), nrow(gx_com)))

## ---- Helper: mature miRNA vector (sum loci, log2+1) ----
get_mir <- function(mature) {
  s  <- stem_of(mature)
  rx <- paste0("^", s, "(-[0-9]+)?$")
  h  <- grep(rx, rownames(mraw_com), value = TRUE)
  if (!length(h)) return(NULL)
  v  <- if (length(h) == 1) mraw_com[h, ] else colSums(mraw_com[h, , drop = FALSE])
  log2p(v, PARAMS$pseudo)
}

## ---- Helper: partial Spearman via rank residuals ----
## Method: regress rank(x) and rank(y) on covariates Z; correlate
## residuals. Equivalent to the partial Spearman of Hollander & Wolfe.
## p-value uses t distribution with n - ncol(Z) - 2 df.
partial_spearman <- function(x, y, Z, label = "") {
  ok  <- complete.cases(x, y, Z)
  x   <- x[ok]; y <- y[ok]; Z <- Z[ok, , drop = FALSE]
  rx  <- rank(x); ry <- rank(y)
  ex  <- residuals(lm(rx ~ Z)); ey <- residuals(lm(ry ~ Z))
  r   <- cor(ex, ey)
  n   <- length(ex); k <- ncol(Z)
  tst <- r * sqrt((n - k - 2) / (1 - r^2))
  pv  <- 2 * pt(-abs(tst), df = n - k - 2)
  dplyr::tibble(label = label, partial_rho = round(r, 4),
                p = signif(pv, 3), n = n, n_covariates = k)
}

## ---- Helper: bootstrap 95% CI for one Spearman rho ----
boot_spearman_ci <- function(x, y, R = 2000) {
  ok <- complete.cases(x, y); x <- x[ok]; y <- y[ok]; n <- length(x)
  rhos <- replicate(R, {
    i <- sample(n, replace = TRUE)
    suppressWarnings(cor(x[i], y[i], method = "spearman"))
  })
  quantile(rhos, c(0.025, 0.975), na.rm = TRUE)
}

# ===================================================================
# ANALYSIS 1: Bootstrap 95% CIs for all 12 candidate pairs
# ===================================================================
message("\n=== Analysis 1: Bootstrap 95% CIs (R = 2000) ===")
ci_rows <- lapply(seq_len(nrow(l2)), function(i) {
  mv <- get_mir(l2$mirna[i]); if (is.null(mv)) return(NULL)
  rg <- l2$symbol_used[i]
  if (!rg %in% rownames(gx_com)) return(NULL)
  ci <- boot_spearman_ci(mv, gx_com[rg, ])
  dplyr::tibble(mirna = l2$mirna[i], target = l2$target[i],
                rho = l2$rho[i], ci_lo = round(ci[1], 4),
                ci_hi = round(ci[2], 4), fdr = l2$fdr[i],
                inverse_sig = l2$inverse_sig[i])
})
ci_tab <- dplyr::bind_rows(ci_rows)
readr::write_csv(ci_tab, "results/06_bootstrap_ci.csv")
pri <- dplyr::filter(ci_tab, mirna == "hsa-miR-200c-3p", target == "EOGT")
message(sprintf("  Priority node rho = %.4f [%.4f, %.4f]",
                pri$rho, pri$ci_lo, pri$ci_hi))
print(as.data.frame(ci_tab))

# ===================================================================
# ANALYSIS 2 + 3: Partial Spearman and subtype-stratified correlations
# ===================================================================
message("\n=== Analysis 2/3: Partial Spearman and stratification ===")

## Get PAM50 and ER status via TCGAquery_subtype (most reliable source)
subtype_df <- tryCatch({
  s  <- TCGAbiolinks::TCGAquery_subtype(tumor = "brca")
  dplyr::tibble(
    patient = as.character(s$patient),
    PAM50   = as.character(s$BRCA_Subtype_PAM50),
    ER      = as.character(s[["ER.Status"]])
  )
}, error = function(e) {
  message("  TCGAquery_subtype unavailable: ", conditionMessage(e))
  tryCatch({
    cl <- GDCquery_clinic(TCGA_PROJECT, "clinical")
    dplyr::tibble(
      patient = cl$bcr_patient_barcode,
      PAM50   = if ("paper_BRCA_Subtype_PAM50" %in% colnames(cl))
                  as.character(cl$paper_BRCA_Subtype_PAM50) else NA_character_,
      ER      = if ("er_status_by_ihc" %in% colnames(cl))
                  as.character(cl$er_status_by_ihc) else NA_character_
    )
  }, error = function(e2) {
    message("  Clinical data unavailable; skipping partial/stratified analyses")
    NULL
  })
})

## Purity: from PanCancerAtlas_subtypes (has CPE purity column)
purity_df <- tryCatch({
  pc   <- TCGAbiolinks::PanCancerAtlas_subtypes()
  brca <- pc[tolower(pc$cancer.type) == "brca", ]
  dplyr::tibble(
    patient = substr(as.character(brca$pan.samplesID), 1, 12),
    purity  = suppressWarnings(as.numeric(brca$purity))
  )
}, error = function(e) {
  message("  Purity data unavailable; partial correlation will not include purity")
  NULL
})

if (!is.null(subtype_df)) {
  m200c <- get_mir("hsa-miR-200c-3p")
  df_p  <- dplyr::tibble(patient = common,
                          mir     = as.numeric(m200c),
                          eogt    = EOGT_v)
  df_p  <- dplyr::left_join(df_p, subtype_df, by = "patient")
  if (!is.null(purity_df))
    df_p <- dplyr::left_join(df_p, purity_df, by = "patient")

  df_p$er_bin  <- as.integer(tolower(df_p$ER) %in% c("positive","pos"))
  df_p$pam_lum <- as.integer(df_p$PAM50 %in% c("Luminal A","Luminal B"))

  message("  PAM50 distribution:");  print(table(df_p$PAM50, useNA = "ifany"))
  message("  ER distribution:");     print(table(df_p$ER,    useNA = "ifany"))

  partial_rows <- list()
  Z1 <- matrix(df_p$er_bin, ncol = 1, dimnames = list(NULL, "er_bin"))
  partial_rows[[1]] <- partial_spearman(df_p$mir, df_p$eogt, Z1,
                                        "Adjusted for ER status")
  Z2 <- cbind(er_bin = df_p$er_bin, pam_lum = df_p$pam_lum)
  partial_rows[[2]] <- partial_spearman(df_p$mir, df_p$eogt, Z2,
                                        "Adjusted for ER + Luminal PAM50")
  if (!is.null(purity_df) && "purity" %in% colnames(df_p)) {
    Z3 <- cbind(er_bin = df_p$er_bin, pam_lum = df_p$pam_lum,
                purity = df_p$purity)
    partial_rows[[3]] <- partial_spearman(df_p$mir, df_p$eogt, Z3,
                                          "Adjusted for ER + PAM50 + purity")
  }
  partial_tab <- dplyr::bind_rows(partial_rows)
  readr::write_csv(partial_tab, "results/06_partial_correlation.csv")
  message("  Partial correlations:"); print(as.data.frame(partial_tab))

  ## Subtype-stratified
  strat_rows <- list()
  for (sg in c("Positive", "Luminal A", "Luminal B")) {
    idx <- if (sg == "Positive") tolower(df_p$ER) %in% c("positive", "pos")
           else df_p$PAM50 == sg
    idx[is.na(idx)] <- FALSE
    if (sum(idx) < 30) { message("  ", sg, ": n=", sum(idx), " (too few; skipped)"); next }
    ct <- suppressWarnings(
      cor.test(df_p$mir[idx], df_p$eogt[idx], method = "spearman"))
    ci <- boot_spearman_ci(df_p$mir[idx], df_p$eogt[idx])
    strat_rows[[sg]] <- dplyr::tibble(
      subgroup = sg, n = sum(idx),
      rho = round(unname(ct$estimate), 4),
      ci_lo = round(ci[1], 4), ci_hi = round(ci[2], 4),
      p = ct$p.value)
  }
  if (length(strat_rows)) {
    strat_tab <- dplyr::bind_rows(strat_rows)
    strat_tab$fdr_strat <- p.adjust(strat_tab$p, "BH")
    readr::write_csv(strat_tab, "results/06_stratified_correlation.csv")
    message("  Stratified correlations:"); print(as.data.frame(strat_tab))
  }
} else {
  message("  Skipped: no subtype/ER data available.")
}

# ===================================================================
# ANALYSIS 4: EOGT tumor vs solid-tissue-normal
# ===================================================================
message("\n=== Analysis 4: EOGT tumor vs solid-tissue-normal ===")
cache_norm <- "results/eogt_normal_cache.rds"
if (!file.exists(cache_norm)) {
  message("  Downloading normal tissue samples from GDC (first run only)...")
  qn <- tryCatch(
    GDCquery(project = TCGA_PROJECT,
             data.category = "Transcriptome Profiling",
             data.type = "Gene Expression Quantification",
             workflow.type = "STAR - Counts",
             sample.type = "Solid Tissue Normal"),
    error = function(e) { message("  GDCquery failed: ", conditionMessage(e)); NULL })
  if (!is.null(qn)) {
    GDCdownload(qn)
    se_n <- GDCprepare(qn)
    eogt_n <- log2p(assay(se_n, "tpm_unstrand")["EOGT", ], PARAMS$pseudo)
    saveRDS(eogt_n, cache_norm)
    rm(se_n); gc()
  } else {
    saveRDS(NULL, cache_norm)
  }
}
eogt_norm <- readRDS(cache_norm)
if (!is.null(eogt_norm)) {
  eogt_tum <- EOGT_v
  wt  <- wilcox.test(eogt_tum, eogt_norm, conf.int = TRUE)
  fc  <- mean(eogt_tum, na.rm = TRUE) - mean(eogt_norm, na.rm = TRUE)
  norm_tab <- dplyr::tibble(
    n_tumor  = sum(!is.na(eogt_tum)),
    n_normal = sum(!is.na(eogt_norm)),
    mean_tumor_log2tpm  = round(mean(eogt_tum,  na.rm = TRUE), 3),
    mean_normal_log2tpm = round(mean(eogt_norm, na.rm = TRUE), 3),
    log2FC_tumor_vs_normal = round(fc, 3),
    wilcox_p = signif(wt$p.value, 3),
    direction = if (fc > 0) "EOGT higher in tumor" else "EOGT lower in tumor"
  )
  readr::write_csv(norm_tab, "results/06_eogt_normal_tumor.csv")
  message("  EOGT tumor vs normal:"); print(as.data.frame(norm_tab))
} else {
  message("  Skipped: normal tissue data unavailable.")
}

# ===================================================================
# ANALYSIS 5: Rank miR-200c-3p among ALL miRNAs correlated with EOGT
# ===================================================================
message("\n=== Analysis 5: Pan-miRNome ranking for EOGT correlation ===")
message("  Computing Spearman rho for all ", nrow(mraw_com),
        " stem-loops vs EOGT (this takes 2-5 min)...")

rho_all <- apply(mraw_com, 1, function(x) {
  ok <- is.finite(x) & is.finite(EOGT_v) & (x > 0 | EOGT_v > 0)
  if (sum(ok) < 100) return(NA_real_)
  # require detectable expression in at least 10% of samples
  if (sum(x[ok] > 0.5) < sum(ok) * 0.10) return(NA_real_)
  suppressWarnings(cor(x[ok], EOGT_v[ok], method = "spearman"))
})
rho_all <- rho_all[!is.na(rho_all)]
n_total  <- length(rho_all)
n_neg    <- sum(rho_all < 0)
message("  Tested ", n_total, " stem-loops; ", n_neg, " have negative correlation with EOGT")

## miR-200c stems
mir200c_stems <- grep("^hsa-mir-200c$|^hsa-mir-200c-", names(rho_all), value = TRUE)
if (!length(mir200c_stems)) mir200c_stems <- grep("hsa-mir-200c", names(rho_all), value = TRUE)
if (!length(mir200c_stems)) stop("hsa-mir-200c not found in mraw")
message("  miR-200c stems found: ", paste(mir200c_stems, collapse = ", "))
mir200c_rho <- min(rho_all[mir200c_stems], na.rm = TRUE)   # most inverse locus
rank_most_inv <- sum(rho_all <= mir200c_rho, na.rm = TRUE) # rank among all negatives
pct_inv       <- round(100 * rank_most_inv / n_total, 1)

rank_tab <- dplyr::tibble(
  mirna                  = "hsa-miR-200c-3p",
  rho_eogt               = round(mir200c_rho, 4),
  rank_most_inverse      = rank_most_inv,
  n_stem_loops_tested    = n_total,
  n_negatively_correlated = n_neg,
  percentile_most_inverse = pct_inv,
  interpretation = sprintf("miR-200c-3p is among the %d most inversely correlated miRNAs
    out of %d tested (%.1f percentile)", rank_most_inv, n_total, pct_inv)
)
readr::write_csv(rank_tab, "results/06_eogt_ranking.csv")
message("  Ranking:"); print(as.data.frame(dplyr::select(rank_tab, -interpretation)))
message("  ", rank_tab$interpretation)

# ===================================================================
# SUMMARY: plain-text numbers ready to insert into the manuscript
# ===================================================================
summary_lines <- c(
  "=== Layer 5 summary for manuscript ===",
  "",
  "1. Bootstrap 95% CIs (all 12 pairs):",
  apply(ci_tab, 1, function(r)
    sprintf("   %s to %s: rho=%.4f [%.4f, %.4f]  FDR=%.3g",
            r["mirna"], r["target"], as.numeric(r["rho"]),
            as.numeric(r["ci_lo"]), as.numeric(r["ci_hi"]), as.numeric(r["fdr"]))),
  "",
  "2. Partial Spearman (miR-200c-3p to EOGT):"
)
if (file.exists("results/06_partial_correlation.csv")) {
  pt <- readr::read_csv("results/06_partial_correlation.csv", show_col_types = FALSE)
  summary_lines <- c(summary_lines,
    apply(pt, 1, function(r)
      sprintf("   %s: partial_rho=%.4f  p=%.3g", r["label"],
              as.numeric(r["partial_rho"]), as.numeric(r["p"]))))
}
if (file.exists("results/06_stratified_correlation.csv")) {
  st <- readr::read_csv("results/06_stratified_correlation.csv", show_col_types = FALSE)
  summary_lines <- c(summary_lines, "",
    "3. Subtype-stratified (miR-200c-3p to EOGT):",
    apply(st, 1, function(r)
      sprintf("   %s (n=%s): rho=%.4f [%.4f, %.4f]  p=%.3g",
              r["subgroup"], r["n"], as.numeric(r["rho"]),
              as.numeric(r["ci_lo"]), as.numeric(r["ci_hi"]),
              as.numeric(r["p"]))))
}
if (file.exists("results/06_eogt_normal_tumor.csv")) {
  nt <- readr::read_csv("results/06_eogt_normal_tumor.csv", show_col_types = FALSE)
  summary_lines <- c(summary_lines, "",
    "4. EOGT tumor vs normal:",
    sprintf("   n_tumor=%d  n_normal=%d  log2FC_tumor_vs_normal=%.3f  Wilcoxon_p=%.3g",
            nt$n_tumor, nt$n_normal,
            nt$log2FC_tumor_vs_normal, nt$wilcox_p))
}
summary_lines <- c(summary_lines, "",
  "5. miR-200c-3p ranking among all miRNAs vs EOGT:",
  sprintf("   rho=%.4f  rank=%d of %d tested  (%.1f percentile most inverse)",
          rank_tab$rho_eogt, rank_tab$rank_most_inverse,
          rank_tab$n_stem_loops_tested, rank_tab$percentile_most_inverse))
writeLines(summary_lines, "results/06_summary.txt")
message("\nLayer 5 complete. Verify results/06_summary.txt before inserting numbers.")
message("ABSOLUTE CHECK: open each CSV, confirm values, then add to manuscript.")
