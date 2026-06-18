# Microarray Analysis CTEPH vs PAH

This repository contains the R script used to assess associations between
ribosomal protein gene expression (microarray) and clinical/hemodynamic
parameters in a validation cohort of CTEPH, PAH, and healthy control subjects.

---

## Overview

Gene expression data were obtained from Affymetrix microarrays and normalized
using the Robust Multi-array Average (RMA) algorithm (log2 summarization).
Associations between candidate gene expression and clinical parameters were
assessed using **partial Spearman correlation**, with sex included as a
controlling variable to account for its confounding effect on gene expression.

Two levels of analysis are performed:

1. **Gene-level** — individual ribosomal protein genes (RPL and RPS family members)
2. **Family-level** — mean expression across all RPL or RPS candidate genes

Given the exploratory nature of the analysis and the limited sample size (n = 14),
no correction for multiple comparisons was applied. Results are considered
hypothesis-generating.

---

## Cohort

| Group   | n |
|---------|---|
| CTEPH   | 5 |
| PAH     | 5 |
| Control | 4 |
| **Total** | **14** |

Sex distribution: CTEPH 1F/4M · PAH 2F/3M · Control 2F/2M

---

## Requirements

R (>= 4.1.0) and the following packages:

```r
install.packages(c(
  "ggplot2", "dplyr", "tidyr", "patchwork",
  "openxlsx", "pheatmap", "ppcor", "ggh4x"
))
```

---

## Input data

The raw microarray and clinical data used in this study contain patient-level
information and are not publicly available due to ethical and privacy restrictions.

**Data availability:** The datasets are available from the corresponding author
upon reasonable request, subject to ethical review and data sharing agreement.
Please contact mgratacos@idibgi.org.

To reproduce the analysis, place the input files in `01_data/correlations/`:

| File | Description |
|------|-------------|
| `genes.xlsx` | Expression matrix: rows = genes, columns = samples (RMA log2 values) |
| `clinica.xlsx` | Clinical metadata: rows = samples, columns = clinical/hemodynamic parameters |

Sample naming convention:
- CTEPH samples: `CTEPH1`, `CTEPH2`, ...
- PAH samples: `PAH1`, `PAH2`, ...
- Control samples: `control1`, `control2`, ...

---

## Usage

```r
source("correlation_analysis_CTEPH_github.R")
```

The script runs sequentially. Key parameters that can be adjusted at the top:

```r
candidate_genes       <- c("RPS", "RPL23A", "RPL4", "RPL21", "RPSA", "RPL32")
MISSING_THRESHOLD_PCT <- 30    # max % missing values allowed per clinical variable
n_per_family          <- 4     # max scatter plots per family (RPL / RPS)
```

---

## Output files

```
02_results/analisi_final/sex_sensitivity/
├── correlation_results_partial_spearman.csv   # gene-level results
└── correlation_results_grouped_families.csv   # family-level results

03_plots/sex_sensitivity/
├── Fig_S1_missing_values.png
├── Fig_heatmap_partial_spearman.pdf           # gene-level heatmap (|rho| >= 0.5)
├── Fig_heatmap_partial_spearman_filtered.pdf  # filtered: params with >= 1 gene at |rho| >= 0.5
└── Fig_scatter_partial_spearman.pdf           # individual gene scatter plots

Fig_heatmap_grouped_families.png               # family-level heatmap
Fig_heatmap_grouped_families_filtered.png      # family-level heatmap (filtered)
Fig_scatter_grouped_families.png               # family-level scatter plots
```

---

## Statistical methods

- **Partial Spearman correlation** (`pcor.test`, ppcor package) with sex as covariate
- Minimum 6 complete observations required per gene–parameter pair
- Significance thresholds: `*` p < 0.05, `**` p < 0.01, `***` p < 0.001
- Variables with > 30% missing values excluded from analysis
- Family-level scores computed as the arithmetic mean across all candidate genes of each family

---

## Session info

```r
sessionInfo()
```

Results were generated with R 4.4.1. Full session information is available
upon request.

---

## Contact

1) Olga Tura-Ceide - olgaturac@gmail.com (correspondance author and PI)
1) Míriam Peracaula Domínguez - mperacaula@idibgi.org (1st author)
1,2) Miquel Gratacós i Aurich — mgratacos@idibgi.org (bionformatitian)

1) Translational Research on Cardiovascular Respiratory Diseases (CAREs), IDIBGI
2) Grup de Recerca en Estadística i Anàlisi de Dades Composicionals (GR-EADC), UdG
