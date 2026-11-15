# results/

| File | Content |
|---|---|
| `01_target_consensus.csv` | multiMiR prediction and validation support per candidate microRNA |
| `02_tcga_correlation.csv` | Spearman correlation for the 12 tested pairs in TCGA-BRCA (n = 1076), BH across 12 pairs. microRNA abundance is precursor-level |
| `03_survival.csv`, `03_survival_with_ci.csv` | univariable Cox models; confidence intervals derived from the hazard ratio and Wald p value |
| `04_geo_foldchange.csv` | first analysis of the resistant models (July 2026). Enzyme rows (GSE26459) are current. microRNA rows are SUPERSEDED: only the GPL7722 batch of GSE28267 was used and six values came from the opposite-arm probe |
| `08_geo_named_arm_foldchange.csv` | current resistant-model table: named-arm probes, both batches of GSE28267, GSE66607, and the enzyme rows of GSE26459 |
| `07_miR200c_arm_specific.csv` | correlation of miR-200c-3p and miR-200c-5p (isoform-level quantification) with EOGT |
| `08b_GSE28267_limma_R_confirmation.csv` | limma in R for both batches of GSE28267, named-arm probes; source of the manuscript values |
| `09_survival_events_schoenfeld.csv` | rerun of the univariable Cox models with deaths per model and Schoenfeld tests |
| `06A` to `06E` | sensitivity analyses for miR-200c/EOGT (partial correlation, subtype strata, confidence intervals, ranking among predicted microRNAs, tumor versus normal, multivariable Cox) |
| `sessionInfo.txt`, `package_versions.csv` | R environment of the July 2026 run |
