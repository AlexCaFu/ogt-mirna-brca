# =====================================================================
# R/04_assemble_table2.R
# Assemble the quantitative columns of Table 2 by joining the target
# compilation (01), the TCGA correlation and survival (02), and the GEO
# fold changes (03), then apply the tier logic. Writes a regenerated
# table for verification against the curated tables/Table_2.csv.
#
# Tier logic (as defined in the manuscript):
#   Tier 1: luciferase-validated AND inverse TCGA correlation (FDR < 0.05)
#   Tier 2: inverse TCGA correlation (FDR < 0.05), no gold-standard validation
#   Tier 3: null, positive, or directionally inconsistent in breast tissue
# =====================================================================
source("R/config.R")
source("R/00_setup.R")

rd <- function(f) if (file.exists(f)) readr::read_csv(f, show_col_types = FALSE) else NULL
tg   <- rd("results/01_targets.csv")
corr <- rd("results/02_tcga_correlation.csv")
surv <- rd("results/02_tcga_survival.csv")
geo  <- rd("results/03_geo_foldchanges.csv")

stopifnot(!is.null(tg), !is.null(corr))

key <- CANDIDATES
key$rho <- NA; key$FDR <- NA; key$validated <- FALSE; key$predicted <- FALSE
key$HR <- NA;  key$HR_p <- NA; key$mir_dFC <- NA; key$enz_dFC <- NA

for (i in seq_len(nrow(key))) {
  m <- key$mirna[i]; e <- key$enzyme[i]
  ct <- corr[corr$mirna == m & corr$enzyme == e, ]
  if (nrow(ct)) { key$rho[i] <- ct$rho[1]; key$FDR[i] <- ct$FDR[1] }
  tt <- tg[tg$mirna == m & tg$enzyme == e, ]
  if (nrow(tt)) { key$validated[i] <- isTRUE(tt$validated[1]); key$predicted[i] <- isTRUE(tt$predicted[1]) }
  if (!is.null(surv)) {
    sm <- surv[surv$feature == m, ]; if (nrow(sm)) { key$HR[i] <- sm$HR[1]; key$HR_p[i] <- sm$p[1] }
  }
  if (!is.null(geo)) {
    # miRNA fold change from a miRNA series; enzyme fold change from an mRNA series
    gm <- geo[geo$measure == "miRNA" & tolower(geo$target) == tolower(m), ]
    if (nrow(gm)) key$mir_dFC[i] <- gm$log2FC[which.max(abs(gm$log2FC))]
    ge <- geo[geo$measure == "mRNA" & toupper(geo$target) == toupper(e), ]
    if (nrow(ge)) key$enz_dFC[i] <- ge$log2FC[which.max(abs(ge$log2FC))]
  }
}

inverse_sig <- with(key, is.finite(rho) & rho < 0 & is.finite(FDR) & FDR < FDR_THRESHOLD)
key$tier <- ifelse(key$validated & inverse_sig, "Tier 1",
            ifelse(inverse_sig, "Tier 2", "Tier 3"))

# order like the manuscript: tiers, then by |rho|
key <- key[order(key$tier, -abs(replace(key$rho, is.na(key$rho), 0))), ]
readr::write_csv(key, "tables/Table_2_regenerated.csv")
message("Wrote tables/Table_2_regenerated.csv")
message("Compare against the curated tables/Table_2.csv and VERIFICATION.md.")
print(key[, c("tier","mirna","enzyme","drug","rho","FDR","HR","HR_p","mir_dFC","enz_dFC")])
