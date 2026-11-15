# Data

No raw data are stored in this repository. Every dataset is public and is
re-downloaded programmatically by the pipeline, so results are reproducible
from source. This keeps the repository small and respects the data providers'
terms of use.

## The Cancer Genome Atlas (TCGA-BRCA)

- Project: `TCGA-BRCA`
- Assays used: gene expression (STAR counts), miRNA quantification, clinical
  follow-up, and the PAM50 molecular-subtype and tumor-purity annotations.
- Retrieved with the `TCGAbiolinks` package, which downloads into a local
  `GDCdata/` cache, which is not stored in the repository.
- Access and use are governed by the NIH Genomic Data Sharing Policy. The
  expression and clinical data used here are open-access.

## Gene Expression Omnibus (GEO)

| Accession | Contrast | Platform |
|---|---|---|
| `GSE26459` | Tamoxifen-sensitive vs resistant MCF7 (mRNA) | GPL570 |
| `GSE28267` | Tamoxifen-sensitive MCF-7 vs endocrine-resistant LY2 (miRNA), two array batches | GPL7722 (6 arrays) and GPL7723 (14 arrays), miRCURY LNA, two-color |
| `GSE66607` | Isogenic MCF7 vs tamoxifen-resistant MCF7 (miRNA) | GPL16384 |

Retrieved with `GEOquery` (R scripts) or directly from the GEO FTP site (`python/08_geo_named_arm.py`). `GSE66607` has a single replicate per state and is
treated as descriptive.

## Reproducing from source

The download and processing steps are performed by `run_all.R` and by the scripts listed in the main README.
The local caches (`GDCdata/`, `data/geo/`, `data/cache/`) are not part of the repository.
