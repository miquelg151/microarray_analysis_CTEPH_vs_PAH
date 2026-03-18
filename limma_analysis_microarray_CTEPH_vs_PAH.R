# =============================================================================
# Microarray Gene Expression Analysis Pipeline
# =============================================================================
# Description:
#   End-to-end pipeline for processing and analysing Affymetrix microarray data.
#   Steps include:
#     1. Raw data loading and RMA normalisation
#     2. Batch effect correction using ComBat (SVA)
#     3. Differential expression analysis via limma (three contrasts)
#     4. Volcano plot visualisation
#     5. PCA before and after batch correction
#
# Contrasts analysed:
#   - CTEPH vs CTL  (Chronic Thromboembolic Pulmonary Hypertension vs Control / PIE samples)
#   - PAH vs CTL    (Pulmonary Arterial Hypertension vs Control / PIE samples)
#   - CTEPH vs PAH
#
# Input files (place in 01_data/):
#   - Sample_Info_sex+age_2.xlsx   (sample metadata)
#   - geneArray_TAC.txt            (pre-processed expression matrix from TAC)
#   - *.CEL files                  (raw Affymetrix CEL files)
#
# Output directories:
#   - 02_results/   expression matrices and DEG tables (.xlsx)
#   - 03_plots/     volcano plots and PCA figures (.pdf)
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Environment setup
# -----------------------------------------------------------------------------

renv::restore()   # Restore project-specific package versions via renv

library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)        # Combine multiple ggplot2 plots
library(showtext)         # Google Fonts support in plots

library(affy)             # Reading Affymetrix CEL files
library(oligo)            # Oligonucleotide array processing + RMA normalisation

library(limma)            # Linear models for microarray differential expression
library(sva)              # Surrogate variable analysis / ComBat batch correction

library(ggrepel)          # Non-overlapping text labels in ggplot2


# --- Output directories -------------------------------------------------------

dir.create("02_results/analisi_final/age_sex", recursive = TRUE, showWarnings = FALSE)
dir.create("03_plots/pca",                      recursive = TRUE, showWarnings = FALSE)


# -----------------------------------------------------------------------------
# 1. Custom theme and colour palette
# -----------------------------------------------------------------------------

font_add_google(name = "Lato", family = "Lato")
showtext_auto()

#' Minimal ggplot2 theme with Lato font
theme_microbiome <- function() {
	theme_minimal() +
		theme(
			plot.title      = element_text(family = "Lato", face = "bold",
																		 hjust = 0.5, margin = margin(b = 20)),
			plot.subtitle   = element_text(family = "Lato", hjust = 0.5,
																		 vjust = 1.5, face = "italic", color = "grey40"),
			axis.title.x    = element_text(family = "Lato", vjust = -5,
																		 color = "grey25", face = "italic"),
			axis.title.y    = element_text(family = "Lato", vjust = 6,
																		 color = "grey25", face = "italic"),
			axis.text       = element_text(family = "Lato"),
			legend.title    = element_text(family = "Lato", face = "bold"),
			legend.text     = element_text(family = "Lato"),
			panel.grid      = element_blank(),
			strip.text      = element_text(family = "Lato", face = "bold",
																		 margin = margin(t = 6, b = 6)),
			strip.background = element_rect(fill = "grey95", color = NA),
			plot.margin     = margin(20, 20, 20, 20),
			axis.line       = element_line(color = "grey50", linewidth = 0.3),
			axis.ticks      = element_line(color = "grey60", linewidth = 0.2),
			plot.background = element_rect(fill = "white", color = NA),
			panel.background = element_rect(fill = "white", color = NA),
			legend.key.size = unit(0.7, "lines"),
			legend.spacing.y = unit(0.4, "lines"),
			legend.background = element_blank()
		)
}

# Six-colour palette used throughout all figures
colors_custom <- c(
	"#018dfa",   # blue
	"#ff3230",   # red
	"#37d275",   # green
	"#a1085c",   # magenta
	"#c0e15c",   # lime
	"#7d49e2"    # purple
)


# -----------------------------------------------------------------------------
# 2. Load sample metadata
# -----------------------------------------------------------------------------

phenodata <- readxl::read_xlsx("01_data/Sample_Info_sex+age_2.xlsx")

# Rename and type-cast the batch variable (Array Date)
phenodata$ArrayDate       <- as.factor(phenodata$`Array Date`)
phenodata$`Array Date`    <- NULL

# Factor-encode grouping variables for model matrix construction
phenodata$Disease <- factor(phenodata$Disease)
phenodata$Sex     <- factor(phenodata$Sex)


# -----------------------------------------------------------------------------
# 3. Raw CEL file loading and RMA normalisation
# -----------------------------------------------------------------------------

cel_files <- list.celfiles("01_data/", full.names = TRUE)
raw       <- read.celfiles(cel_files)

# Visual QC: log2 intensity distributions before normalisation
boxplot(
	log2(exprs(raw) + 1),
	las = 2, outline = FALSE,
	main = "Before normalisation (raw, log2)",
	ylab = "log2 expression"
)

# Apply RMA (Robust Multi-array Average) normalisation:
# background correction + quantile normalisation + log2 summarisation
norm      <- rma(raw, normalize = TRUE)
exprs_mat <- exprs(norm)

# Remove the bottom 10% lowest-variance probesets (uninformative features)
cv  <- apply(exprs_mat, 1, sd) / rowMeans(exprs_mat)
cut <- quantile(cv, 0.10, na.rm = TRUE)
exprs_mat <- exprs_mat[cv > cut, ]

# Visual QC: distributions after normalisation
boxplot(
	exprs_mat,
	las = 2, outline = FALSE,
	main = "After RMA normalisation"
)

# Save the filtered, normalised matrix for later use in ComBat
exprs_mat_rma <- exprs_mat

# -----------------------------------------------------------------------------
# 4. Load TAC expression matrix and align annotations
# -----------------------------------------------------------------------------
#
# The TAC software (Thermo Fisher) exports a tab-separated expression file
# with probe IDs, gene symbols, and one column per sample.

exprs_tac <- read.table(
	"01_data/geneArray_TAC.txt",
	header          = TRUE,
	sep             = "\t",
	dec             = ",",
	quote           = "",
	stringsAsFactors = FALSE,
	check.names     = FALSE
)

# Split annotation (columns 1-2) from expression values (columns 3+)
annot     <- exprs_tac[, 1:2]          # columns: ID, Gene Symbol
exprs_mat <- as.matrix(exprs_tac[, -(1:2)])
rownames(exprs_mat) <- annot$ID
storage.mode(exprs_mat) <- "numeric"

# Strip the TAC-appended suffix from sample column names
colnames(exprs_mat) <- sub("\\.rma-gene-full\\.chp Signal", "",
													 colnames(exprs_mat))

# Keep only probes with a valid gene symbol annotation
gene_symbol    <- annot$`Gene Symbol`[match(rownames(exprs_mat), annot$ID)]
keep           <- !is.na(gene_symbol) & gene_symbol != ""
exprs_mat_filt <- exprs_mat[keep, ]

cat(sprintf("Probes retained after annotation filter: %d / %d\n",
						sum(keep), length(keep)))


# -----------------------------------------------------------------------------
# 5. ComBat batch correction
# -----------------------------------------------------------------------------
#
# ComBat adjusts for known batch effects (here: array processing date) while
# preserving biological variance specified in the model (Disease, Sex, Age).

mod_combat <- model.matrix(
	~ Disease + Sex + Age,
	data = phenodata
)

# --- Batch-correct the RMA-filtered matrix (used for DEG analysis) -----------
exprs_combat <- ComBat(
	dat        = exprs_mat_filt,   # probes x samples matrix (annotated probes only)
	batch      = phenodata$ArrayDate,
	mod        = mod_combat,
	par.prior  = TRUE,             # use parametric empirical Bayes priors
	prior.plots = FALSE
)

# Visual QC: distributions after ComBat
boxplot(exprs_combat,
				outline = FALSE, las = 2,
				main = "After ComBat batch correction")

plotMDS(exprs_combat, labels = phenodata$Disease)

# --- Also batch-correct and export the full RMA matrix (all probes) ----------
exprs_combat_full <- ComBat(
	dat        = exprs_mat_rma,
	batch      = phenodata$ArrayDate,
	mod        = mod_combat,
	par.prior  = TRUE,
	prior.plots = FALSE
)

# Annotate and export the full batch-corrected matrix for downstream analysis and interpretation of the results
exprs_combat_full_df <- as.data.frame(exprs_combat_full)
exprs_combat_full_df$ID <- rownames(exprs_combat_full_df)
exprs_combat_full_df <- exprs_combat_full_df[, c("ID", setdiff(colnames(exprs_combat_full_df), "ID"))]
rownames(exprs_combat_full_df) <- NULL

exprs_combat_full_df$GeneSymbol <-
	annot$`Gene Symbol`[match(exprs_combat_full_df$ID, annot$ID)]

exprs_combat_full_df <- exprs_combat_full_df |>
	mutate(GeneSymbol = gsub('"', '', GeneSymbol)) |>
	separate_rows(GeneSymbol, sep = ";\\s*")

openxlsx::write.xlsx(
	exprs_combat_full_df,
	"02_results/exprs_mat_raw_normalized_batch_corr.xlsx"
)


# -----------------------------------------------------------------------------
# 6. Differential expression analysis (limma)
# -----------------------------------------------------------------------------
#
# Model: ~ Disease + Sex + Age
#   - Disease reference level: CTL (first alphabetical factor level)
#   - Two direct contrasts from model coefficients: CTEPH vs CTL, PAH vs CTL
#   - One indirect contrast built with makeContrasts: CTEPH vs PAH

design <- model.matrix(~ Disease + Sex + Age, data = phenodata)
colnames(design) <- make.names(colnames(design))   # ensure valid R names

fit <- lmFit(exprs_combat, design)   # fit linear model per probeset
fit <- eBayes(fit)                   # empirical Bayes moderation of variances


# --- Extract results for each contrast ---------------------------------------

res_CTEPH_vs_CTL <- topTable(fit, coef = "DiseaseCTEPH", n = Inf)
res_PAH_vs_CTL   <- topTable(fit, coef = "DiseasePAH",   n = Inf)

# For the CTEPH-vs-PAH contrast, we need an explicit contrast matrix
contrast_matrix <- makeContrasts(
	CTEPH_vs_PAH = DiseaseCTEPH - DiseasePAH,
	levels = design
)
fit_contrast      <- contrasts.fit(fit, contrast_matrix)
fit_contrast      <- eBayes(fit_contrast)
res_CTEPH_vs_PAH  <- topTable(fit_contrast, coef = "CTEPH_vs_PAH", n = Inf)


# --- Annotate results with gene symbols --------------------------------------
#
# Each results data frame gets the probe ID added as a column, then we map
# to gene symbols, expand multi-gene probes (";"-separated), and keep only
# the probe with the largest absolute t-statistic per gene (most informative).

annotate_results <- function(res_df, annot) {
	res_df$ID <- rownames(res_df)
	res_df$GeneSymbol <- annot$`Gene Symbol`[match(res_df$ID, annot$ID)]
	res_df <- res_df |>
		mutate(GeneSymbol = gsub('"', '', GeneSymbol)) |>
		separate_rows(GeneSymbol, sep = ";\\s*") |>
		group_by(GeneSymbol) |>
		slice_max(order_by = abs(t), n = 1, with_ties = FALSE) |>
		ungroup()
	rownames(res_df) <- res_df$GeneSymbol
	res_df
}

res_CTEPH_vs_CTL <- annotate_results(res_CTEPH_vs_CTL, annot)
res_PAH_vs_CTL   <- annotate_results(res_PAH_vs_CTL,   annot)
res_CTEPH_vs_PAH <- annotate_results(res_CTEPH_vs_PAH, annot)


# --- Export DEG tables -------------------------------------------------------

out_dir <- "02_results/analisi_final/age_sex/"

openxlsx::write.xlsx(res_CTEPH_vs_PAH,
										 file.path(out_dir, "DEGs_CTEPH_vs_PAH_limma_batch.xlsx"),
										 rowNames = TRUE)
openxlsx::write.xlsx(res_CTEPH_vs_CTL,
										 file.path(out_dir, "DEGs_CTL_vs_CTEPH_limma_batch.xlsx"),
										 rowNames = TRUE)
openxlsx::write.xlsx(res_PAH_vs_CTL,
										 file.path(out_dir, "DEGs_CTL_vs_PAH_limma_batch.xlsx"),
										 rowNames = TRUE)

# -----------------------------------------------------------------------------
# 7. Volcano plots
# -----------------------------------------------------------------------------

# Re-import from disk to ensure reproducibility from the saved files
res_CTEPH_vs_PAH <- openxlsx::read.xlsx(
	file.path(out_dir, "DEGs_CTEPH_vs_PAH_limma_batch.xlsx"), rowNames = TRUE)
res_CTEPH_vs_CTL <- openxlsx::read.xlsx(
	file.path(out_dir, "DEGs_CTL_vs_CTEPH_limma_batch.xlsx"), rowNames = TRUE)
res_PAH_vs_CTL   <- openxlsx::read.xlsx(
	file.path(out_dir, "DEGs_CTL_vs_PAH_limma_batch.xlsx"),   rowNames = TRUE)


#' Draw a volcano plot for a limma topTable result
#'
#' @param df          Data frame from topTable (rows = genes).
#' @param logfc_cutoff Absolute log2FC threshold for significance (default 0.585 ≈ 1.5-fold).
#' @param padj_cutoff  Adjusted p-value threshold (default 0.05).
#' @param title        Plot title string.
#' @param label_top    Number of top genes (by adj.P.Val) to label.
#' @return A ggplot2 object.

plot_volcano <- function(df,
												 logfc_cutoff = 0.585,
												 padj_cutoff  = 0.05,
												 title        = "",
												 label_top    = 10) {
	
	df <- df |>
		mutate(
			gene       = rownames(df),
			log10padj  = -log10(adj.P.Val),
			status     = case_when(
				adj.P.Val < padj_cutoff & logFC >  logfc_cutoff ~ "Up",
				adj.P.Val < padj_cutoff & logFC < -logfc_cutoff ~ "Down",
				TRUE ~ "NS"
			)
		)
	
	# Label only the most significant genes
	top_genes <- df |> arrange(adj.P.Val) |> slice(1:label_top)
	
	ggplot(df, aes(logFC, log10padj)) +
		geom_point(aes(color = status), alpha = 0.6, size = 1.5) +
		ggrepel::geom_text_repel(data  = top_genes,
														 aes(label = gene),
														 size  = 3,
														 max.overlaps = Inf) +
		scale_color_manual(
			values = c("Down" = "#3B4CC0", "NS" = "grey85", "Up" = "#B40426")
		) +
		geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff),
							 linetype = "dashed", linewidth = 0.4) +
		geom_hline(yintercept = -log10(padj_cutoff),
							 linetype = "dashed", linewidth = 0.4) +
		theme_microbiome() +
		labs(
			title = title,
			x     = expression(log[2]~Fold~Change),
			y     = expression(-log[10]~adj.~italic(p))
		)
}

volcano_cteph_vs_pah <- plot_volcano(res_CTEPH_vs_PAH, title = "CTEPH vs PAH")
volcano_cteph_vs_ctl <- plot_volcano(res_CTEPH_vs_CTL, title = "CTEPH vs CTL")
volcano_pah_vs_ctl   <- plot_volcano(res_PAH_vs_CTL,   title = "PAH vs CTL")

# Save each volcano plot to a PDF
ggsave("03_plots/volcano_cteph_vs_pah.pdf",
			 plot = volcano_cteph_vs_pah, width = 10, height = 10)
ggsave("03_plots/volcano_cteph_vs_ctl.pdf",
			 plot = volcano_cteph_vs_ctl, width = 10, height = 10)
ggsave("03_plots/volcano_pah_vs_ctl.pdf",
			 plot = volcano_pah_vs_ctl,   width = 10, height = 10)


# -----------------------------------------------------------------------------
# 8. PCA: before vs after ComBat batch correction
# -----------------------------------------------------------------------------
#
# We run PCA on the transposed expression matrix (samples as rows) to visualise
# whether samples cluster by disease group or by batch before and after correction.

# Helper to extract variance explained for axis labels
pct_var <- function(pca_obj, pc) {
	round(100 * summary(pca_obj)$importance[2, pc], 1)
}

# --- PCA on annotated probes BEFORE batch correction -------------------------

pca_before    <- prcomp(t(exprs_mat_filt), center = TRUE, scale. = FALSE)
pca_before_df <- data.frame(
	PC1        = pca_before$x[, 1],
	PC2        = pca_before$x[, 2],
	Disease    = phenodata$Disease,
	Array.Date = phenodata$ArrayDate,
	Sex        = phenodata$Sex
)

pca_b <- ggplot(pca_before_df, aes(PC1, PC2, color = Disease)) +
	geom_point(size = 4, alpha = 0.7) +
	stat_ellipse(type = "norm", level = 0.75, linetype = 1, linewidth = 1, alpha = 0.7) +
	scale_color_manual(values = colors_custom[1:3]) +
	labs(
		title = "PCA before ComBat batch correction",
		x = paste0("PC1 (", pct_var(pca_before, 1), "%)"),
		y = paste0("PC2 (", pct_var(pca_before, 2), "%)")
	) +
	theme_microbiome()

pca_b_batch <- ggplot(pca_before_df, aes(PC1, PC2, color = Disease, shape = Array.Date)) +
	geom_point(alpha = 0.7, size = 4) +
	stat_ellipse(type = "norm", level = 0.75, linewidth = 1, linetype = 1, alpha = 0.7) +
	scale_color_manual(values = colors_custom[1:3]) +
	labs(
		title    = "PCA before ComBat batch correction",
		subtitle = "Shape = array batch",
		color    = "Disease",
		shape    = "Batch",
		x = paste0("PC1 (", pct_var(pca_before, 1), "%)"),
		y = paste0("PC2 (", pct_var(pca_before, 2), "%)")
	) +
	theme_microbiome()


# --- PCA on annotated probes AFTER ComBat correction -------------------------

pca_after    <- prcomp(t(exprs_combat), center = TRUE, scale. = FALSE)
pca_after_df <- data.frame(
	PC1        = pca_after$x[, 1],
	PC2        = pca_after$x[, 2],
	Disease    = phenodata$Disease,
	Array.Date = phenodata$ArrayDate,
	Sex        = phenodata$Sex
)

pca_a <- ggplot(pca_after_df, aes(PC1, PC2, color = Disease)) +
	geom_point(size = 4, alpha = 0.7) +
	stat_ellipse(type = "norm", level = 0.75, linetype = 1, linewidth = 1, alpha = 0.7) +
	scale_color_manual(values = colors_custom[1:3]) +
	labs(
		title = "PCA after ComBat batch correction",
		x = paste0("PC1 (", pct_var(pca_after, 1), "%)"),
		y = paste0("PC2 (", pct_var(pca_after, 2), "%)")
	) +
	theme_microbiome()

pca_a_b <- ggplot(pca_after_df, aes(PC1, PC2, color = Disease, shape = Array.Date)) +
	geom_point(size = 4, alpha = 0.7) +
	stat_ellipse(type = "norm", level = 0.75, linewidth = 1, linetype = 1, alpha = 0.7) +
	scale_color_manual(values = colors_custom[1:3]) +
	labs(
		title    = "PCA after ComBat batch correction",
		subtitle = "Shape = array batch",
		color    = "Disease",
		shape    = "Batch",
		x = paste0("PC1 (", pct_var(pca_after, 1), "%)"),
		y = paste0("PC2 (", pct_var(pca_after, 2), "%)")
	) +
	theme_microbiome()


# --- Combine and export PCA figures ------------------------------------------

tag_theme <- theme(plot.tag = element_text(size = 16, face = "bold"))

combined_pca_before <- (pca_b / pca_b_batch) +
	plot_annotation(tag_levels = "A", theme = tag_theme)

combined_pca_after <- (pca_a / pca_a_b) +
	plot_annotation(tag_levels = "A", theme = tag_theme)

ggsave("03_plots/pca/PCA_before_batch_corr.pdf",
			 plot = combined_pca_before, height = 14, width = 10)
ggsave("03_plots/pca/PCA_after_batch_corr.pdf",
			 plot = combined_pca_after,  height = 14, width = 10)

message("Pipeline complete. Results written to 02_results/ and 03_plots/.")