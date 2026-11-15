# Arm-specific correlation of miR-200c-3p and miR-200c-5p with EOGT in TCGA-BRCA
# (Isoform Expression Quantification). Run on 15 November 2026; output saved as
# results/07_miR200c_arm_specific.csv. Run from the repository root.
suppressMessages({library(TCGAbiolinks); library(SummarizedExperiment)})
qg <- GDCquery(project = "TCGA-BRCA", data.category = "Transcriptome Profiling",
               data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")
GDCdownload(qg); se <- GDCprepare(qg)
gx <- assay(se, "tpm_unstrand"); rownames(gx) <- rowData(se)$gene_name
gx <- gx[, substr(colnames(gx), 14, 15) == "01", drop = FALSE]
colnames(gx) <- substr(colnames(gx), 1, 12)
gx <- gx[, !duplicated(colnames(gx)), drop = FALSE]
eogt <- log2(gx["EOGT", ] + 1)

qi <- GDCquery(project = "TCGA-BRCA", data.category = "Transcriptome Profiling",
               data.type = "Isoform Expression Quantification")
GDCdownload(qi); iso <- GDCprepare(qi)
arms <- c("miR-200c-3p" = "MIMAT0000617", "miR-200c-5p" = "MIMAT0004657")
iso$mimat <- sub("^mature,", "", iso$miRNA_region)
iso <- iso[iso$mimat %in% arms & substr(iso$barcode, 14, 15) == "01", ]
res <- lapply(names(arms), function(a) {
  d  <- iso[iso$mimat == arms[[a]], ]
  v  <- tapply(d$reads_per_million_miRNA_mapped, d$barcode, sum)
  pt <- substr(names(v), 1, 12); v <- v[!duplicated(pt)]; names(v) <- pt[!duplicated(pt)]
  com <- intersect(names(v), names(eogt))
  ct  <- suppressWarnings(cor.test(log2(v[com] + 1), eogt[com], method = "spearman"))
  data.frame(arm = a, mimat = arms[[a]], n = length(com), median_RPM = round(median(v[com]), 2),
             rho = round(unname(ct$estimate), 4), p = signif(ct$p.value, 3))
})
res <- do.call(rbind, res)
print(res, row.names = FALSE)
write.csv(res, "results/07_miR200c_arm_specific.csv", row.names = FALSE)
