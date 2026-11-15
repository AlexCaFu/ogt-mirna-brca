# Candidate microRNA regulation of O-GlcNAc-related enzymes in breast cancer drug resistance: analysis code and processed results

Code, processed result tables, figure scripts and supplementary tables for the manuscript

> Hernández-Huerta MT, Vasquez Martínez IP, Cabrera-Fuentes HA. Candidate microRNA regulation of O-GlcNAc-related enzymes in breast cancer drug resistance: an integrative analysis prioritizing the miR-200c-3p/EOGT node. Submitted to Drug Discovery Today, 2026.

This release (2.1.0, 14 November 2026) corresponds to manuscript version 7.

## What the analysis does

Twelve candidate microRNA-enzyme pairs (nine microRNAs; OGT, OGA/MGEA5, EOGT) are evaluated in three layers:

1. Target prediction (multiMiR; TargetScanHuman 8.0 site tables by microRNA arm).
2. Co-expression and overall survival in TCGA-BRCA primary tumors (Spearman correlation with Benjamini-Hochberg correction across the 12 pairs, n = 1076; Cox regression).
3. Expression change in public resistant cell models from GEO.

The analysis is observational. It prioritizes miR-200c/EOGT as a hypothesis for functional testing and does not demonstrate direct targeting.

## Points a reader must know

- TCGA microRNA values in the main analysis are stem-loop (precursor) counts, which do not separate the 5p and 3p arms. For miR-200c the arms were quantified separately from the isoform-level files (`R/07`, `results/07_miR200c_arm_specific.csv`): the correlation with EOGT is carried by miR-200c-3p (rho -0.255), and miR-200c-5p is about 300-fold less abundant and correlates weakly (rho -0.111).
- TargetScanHuman 8.0 predicts no site for miR-200c-3p in the EOGT 3' untranslated region. One non-conserved 7mer-m8 site exists for miR-200c-5p (Supplementary Table S6).
- GSE28267 was deposited as two array batches, GPL7722 (2 MCF-7 and 4 LY2 arrays) and GPL7723 (10 MCF-7 and 4 LY2 arrays). The first analysis (July 2026, `results/04_geo_foldchange.csv`) used only GPL7722 and, for some microRNAs, reported the opposite-arm probe. The current analysis (`python/08_geo_named_arm.py`, `results/08_geo_named_arm_foldchange.csv`) analyzes both batches separately and reads each microRNA from the probe of its named arm. The fold changes were confirmed with limma in R (`R/08`, `results/08b_GSE28267_limma_R_confirmation.csv`), which is the source of the values in the manuscript. The batches differ in magnitude and, for miR-24-3p, miR-222-3p, miR-148a-3p and miR-181b-5p, in direction.
- EOGT (GSE26459) and the microRNAs (GSE28267) were measured in different cell models.
- R scripts 00 to 04 in `R/` are reconstructions that reproduce the result tables. The authors' original scripts are in `R/original/`. Scripts `R/07`, `R/08` and `R/09` were run on 15 November 2026; their outputs are `results/07`, `results/08b` and `results/09`. Hazard ratios in `results/09` differ from `results/03_survival.csv` (July 2026) by at most 0.009; the manuscript reports the July values together with the death counts and Schoenfeld tests of the rerun.

## Layout

```
README.md, LICENSE, CITATION.cff, .zenodo.json, MANIFEST.sha256, VERIFICATION.md
run_all.R, make_repro_metadata.R
R/            analysis scripts (see R/README.md); R/original/ holds the authors' original scripts
python/       08_geo_named_arm.py (named-arm analysis of the resistant models), requirements.txt
results/      processed result tables (see results/README.md)
figures/      Figure 1, Figure 2, Graphical Abstract (PNG, PDF), LEGENDS.md, src/ (matplotlib scripts)
tables/       Table_2.csv, Supplementary_Tables_S1-S6.xlsx
tests/        tests of the original pipeline
data/         README only; no third-party data are redistributed
```

## Reproducing

R (version and packages in `results/sessionInfo.txt`): open R in the repository folder and run `source("run_all.R")`. The first run downloads TCGA-BRCA and GEO data.

Resistant models, named arm (Python 3.10 or later, packages in `python/requirements.txt`):

```
python3 python/08_geo_named_arm.py
```

Figures, from the result tables:

```
python3 figures/src/fig1.py results figures
python3 figures/src/fig2.py results figures
python3 figures/src/ga.py results figures
```

`VERIFICATION.md` lists the values printed in the manuscript next to the deposited file and value they come from.

## Data sources

TCGA-BRCA (NIH Genomic Data Commons, via TCGAbiolinks); GEO series GSE26459 (GPL570), GSE28267 (GPL7722 and GPL7723), GSE66607 (GPL16384); TargetScanHuman 8.0 site tables; multiMiR. No third-party data are stored here.

## License

Code: MIT (see LICENSE). Tables and figures: CC BY 4.0.

## Contact

Héctor A. Cabrera-Fuentes, ORCID 0000-0002-8502-7098.
