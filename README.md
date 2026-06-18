# Microarray Expression Analysis — CTEPH vs PAH

This repository contains the R scripts used to perform differential expression
analysis and correlation analysis between ribosomal protein gene expression
(microarray) and clinical/hemodynamic parameters in a cohort of CTEPH, PAH,
and healthy control subjects.

---

## Repository structure

```
microarray_analysis_CTEPH_vs_PAH/
├── README.md
├── limma_analysis_sex_subanalysis_v2.R        # Script 1: differential expression
├── correlation_analysis_CTEPH_github.R        # Script 2: correlation analysis
├── 01_data/
│   └── correlations/                          # input data (not publicly available)
├── 02_results/                                # generated automatically
└── 03_plots/                                  # generated automatically
```

---

## Scripts

### Script 1 — `limma_analysis_sex_subanalysis_v2.R`

Differential expression analysis (DEA) comparing CTEPH vs Control, PAH vs Control,
and CTEPH vs PAH using the **limma** framework on RMA-normalized Affymetrix microarray data.

Includes:
- ComBat batch correction (preserving biological variability)
- Sex-stratified subanalyses (males / females separately)
- Disease × Sex interaction model to test whether differential expression
  of ribosomal genes is sex-dependent
- Boxplots of candidate ribosomal genes stratified by sex and disease group

**Main outputs:**
```
02_results/analisi_final/sex_sensitivity/
├── interaction_Disease_x_Sex.xlsx
└── logFC_ribosomal_sex_stratified.xlsx

03_plots/sex_sensitivity/
├── Fig_ribosomal_CTEPH_vs_CTL_by_sex.pdf
├── Fig_ribosomal_all_groups_by_sex.pdf
└── Fig_heatmap_logFC_sex_stratified.pdf
```

---

### Script 2 — `correlation_analysis_CTEPH_github.R`

Correlation analysis between candidate ribosomal gene expression and
clinical/hemodynamic parameters in the validation cohort (n = 14).

Uses **partial Spearman correlation** (ppcor) with sex as a controlling variable.
Two levels of analysis are performed:
1. **Gene-level** — individual RPL and RPS genes
2. **Family-level** — mean expression across all RPL or RPS candidate genes

**Main outputs:**
```
02_results/analisi_final/sex_sensitivity/
├── correlation_results_partial_spearman.csv   # gene-level results
└── correlation_results_grouped_families.csv   # family-level results

03_plots/sex_sensitivity/
├── Fig_S1_missing_values.png
├── Fig_heatmap_partial_spearman.pdf
├── Fig_heatmap_partial_spearman_filtered.pdf
└── Fig_scatter_partial_spearman.pdf

Fig_heatmap_grouped_families.png
Fig_heatmap_grouped_families_filtered.png
Fig_scatter_grouped_families.png
```

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
  # Differential expression
  "limma", "sva",
  # Correlation analysis
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
Please contact olgaturac@gmail.com or mperacaula@igibgi.org

To reproduce the correlation analysis, place the input files in `01_data/correlations/`:

| File | Description |
|------|-------------|
| `genes.xlsx` | Expression matrix: rows = genes, columns = samples (RMA log2 values) |
| `clinica.xlsx` | Clinical metadata: rows = samples, columns = clinical/hemodynamic parameters |

Sample naming convention:
- CTEPH samples: `CTEPH1`, `CTEPH2`, ...
- PAH samples: `PAH1`, `PAH2`, ...
- Control samples: `control1`, `control2`, ...

---

## Statistical methods

**Differential expression (Script 1)**
- RMA normalization + ComBat batch correction
- limma-voom with empirical Bayes moderation (eBayes)
- Significance: adjusted p-value < 0.05 (Benjamini-Hochberg FDR), |log2FC| > 1
- Disease × Sex interaction term to assess sex-dependent effects

**Correlation analysis (Script 2)**
- Partial Spearman correlation (`pcor.test`, ppcor) with sex as covariate
- Minimum 6 complete observations required per gene–parameter pair
- Significance thresholds: `*` p < 0.05, `**` p < 0.01, `***` p < 0.001
- Variables with > 30% missing values excluded
- Family-level scores: arithmetic mean across all candidate genes per family
- No correction for multiple comparisons (exploratory analysis)

---

## Session info

Results were generated with R 4.4.1. Full session information is available
upon request.

---

## Contact

1) Olga Tura-Ceide - olgaturac@gmail.com (correspondance author and PI)
1) Míriam Peracaula Domínguez - mperacaula@idibgi.org (1st author)
1,2) Miquel Gratacós i Aurich — mgratacos@idibgi.org (bionformatitian)

1) Translational Research on Cardiovascular Respiratory Diseases (CAREs), IDIBGI
2) Grup de Recerca en Estadística i Anàlisi de Dades Composicionals (GR-EADC), UdG
