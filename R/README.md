# R/

- `config.R`, `00_setup.R`, `01_multimir.R`, `02_tcga_correlation_survival.R`, `03_geo_foldchange.R`, `04_assemble_table2.R`: reconstructions of the analysis that reproduce the deposited result tables.
- `06_robustness.R`: the authors' own sensitivity-analysis script.
- `original/`: the authors' original scripts as run in July 2026 (`03_geo_inspect.R`, `04_geo_foldchange.R`, `05_assemble.R`, `06_robustness_ranking.R` and others). `03_geo_inspect.R` keeps only the first platform batch of GSE28267, and `04_geo_foldchange.R` merges both arms of a microRNA and keeps the probe with the largest absolute change. Both choices are superseded for microRNAs by `python/08_geo_named_arm.py`.
- `07_arm_specific_miR200c_EOGT.R`: arm-specific (isoform-level) TCGA correlation for miR-200c. Run on 15 November 2026 (`results/07_miR200c_arm_specific.csv`).
- `08_GSE28267_both_platforms_limma_confirmation.R`: limma in R for both batches of GSE28267, run on 15 November 2026 (`results/08b_GSE28267_limma_R_confirmation.csv`). Source of the GSE28267 values in the manuscript.
- `09_survival_events_and_schoenfeld.R`: event counts and proportional-hazards tests. Run on 15 November 2026 (`results/09_survival_events_schoenfeld.csv`).
