# =============================================================================
# MICROARRAY EXPRESSION vs CLINICAL PARAMETERS — PARTIAL SPEARMAN CORRELATION
# =============================================================================
# METHOD:
#   Partial Spearman correlation (pcor.test) with sex as a controlling variable.
#   Measures the gene~parameter association while mathematically removing the
#   effect of sex, without eliminating biological variability across groups.
#
# STUDY CONTEXT:
#   Validation cohort: n=14 samples (CTEPH, PAH, Control)
#   Expression data: Affymetrix microarray, RMA-normalized log2 values
#   Clinical data:   hemodynamic and laboratory parameters
#
# STATISTICAL NOTE:
#   Given the exploratory nature of the analysis and the limited sample size,
#   no correction for multiple comparisons was applied. Results are considered
#   hypothesis-generating. Partial Spearman correlation was chosen to account
#   for the confounding effect of sex while preserving the non-parametric
#   nature of the analysis.
#
# OUTPUT FILES:
#   - 02_results/analisi_final/sex_sensitivity/correlation_results_partial_spearman.csv
#   - 02_results/analisi_final/sex_sensitivity/correlation_results_grouped_families.csv
#   - 03_plots/sex_sensitivity/Fig_heatmap_partial_spearman.pdf
#   - 03_plots/sex_sensitivity/Fig_heatmap_partial_spearman_filtered.pdf
#   - 03_plots/sex_sensitivity/Fig_scatter_partial_spearman.pdf
#   - Fig_heatmap_grouped_families.png
#   - Fig_heatmap_grouped_families_filtered.png
#   - Fig_scatter_grouped_families.png
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(openxlsx)
library(pheatmap)
library(ppcor)      # partial correlation — install.packages("ppcor") if needed
library(ggh4x)      # nested axis labels — install.packages("ggh4x") if needed


# =============================================================================
# 1. DATA LOADING
# =============================================================================

expr_matrix   <- read.xlsx("01_data/correlations/genes.xlsx",   check.names = FALSE, rowNames = TRUE)
clinical_data <- read.xlsx("01_data/correlations/clinica.xlsx", check.names = FALSE, rowNames = TRUE)

# Save Sex BEFORE numeric conversion (factor must be preserved)
sex_vector <- as.factor(clinical_data$Sex)

# Rename duplicate column names (adds _1 suffix to second occurrence)
# Duplicates arise from pre/post-intervention hemodynamic measurements
names(clinical_data) <- make.unique(names(clinical_data), sep = "_")

# Convert all columns to numeric (Sex will become NA, recovered below)
clinical_data <- clinical_data %>%
  dplyr::mutate(dplyr::across(dplyr::everything(),
                               ~ suppressWarnings(as.numeric(as.character(.)))))

# Restore Sex as factor
clinical_data$Sex <- sex_vector

cat("Sex distribution:\n")
print(table(clinical_data$Sex, useNA = "always"))
cat("\n")

# Remove post-intervention hemodynamic duplicate columns (keep pre-intervention only)
clinical_data$`PAPd.(mmHg)_1` <- NULL
clinical_data$`PAPm.(mmHg)_1` <- NULL
clinical_data$`PAPs.(mmHg)_1` <- NULL
clinical_data$`Qt.(l/min)_1`  <- NULL
clinical_data$`CI.(l/min)_1`  <- NULL
clinical_data$`RAP.(mmHg)_1`  <- NULL
clinical_data$`PCWP.(mmHg)_1` <- NULL


# =============================================================================
# 2. CANDIDATE GENES  ← modify with your gene list
# =============================================================================

candidate_genes <- c("RPS", "RPL23A", "RPL4", "RPL21", "RPSA", "RPL32")


# =============================================================================
# 3. SAMPLE GROUPS AND COLORS
# =============================================================================

sample_names <- colnames(expr_matrix)

group_df <- data.frame(
  sample = sample_names,
  group  = dplyr::case_when(
    grepl("^CTEPH",   sample_names, ignore.case = TRUE) ~ "CTEPH",
    grepl("^PAH",     sample_names, ignore.case = TRUE) ~ "PAH",
    grepl("^control", sample_names, ignore.case = TRUE) ~ "Control",
    TRUE ~ "Other"
  ),
  stringsAsFactors = FALSE
)

group_colors <- c("CTEPH" = "#B2182B", "PAH" = "#2166AC", "Control" = "#4DAC26")

cat("Samples per group:\n")
print(table(group_df$group))

cat("\nSex distribution by group:\n")
print(table(
  group_df$group[match(rownames(clinical_data), group_df$sample)],
  clinical_data$Sex, useNA = "always"
))
cat("\n")


# =============================================================================
# 4. MISSING VALUE ASSESSMENT AND VARIABLE SELECTION
# =============================================================================

common_samples <- intersect(rownames(clinical_data), colnames(expr_matrix))
cat("Common samples between expression matrix and clinical data:", length(common_samples), "\n\n")

clinical_sub <- clinical_data[common_samples, , drop = FALSE]
clinical_sub$group  <- group_df$group[match(rownames(clinical_sub), group_df$sample)]
clinical_sub$sample <- rownames(clinical_sub)

# Exclude Sex, group, and sample from missing value calculation
vars_for_missing <- setdiff(names(clinical_sub), c("Sex", "group", "sample"))

missing_count <- colSums(is.na(clinical_sub[, vars_for_missing]))
missing_pct   <- round(100 * missing_count / nrow(clinical_sub), 1)

missing_summary <- data.frame(
  variable    = names(missing_count),
  n_missing   = missing_count,
  pct_missing = missing_pct,
  n_complete  = nrow(clinical_sub) - missing_count
) %>% dplyr::arrange(n_missing, dplyr::desc(n_complete))

print(missing_summary)
write.csv(missing_summary, "missing_values_summary.csv", row.names = FALSE)

# Selection threshold: keep variables with <= 30% missing values
# Controls lack many hemodynamic parameters (not measured in healthy subjects)
MISSING_THRESHOLD_PCT <- 30

params_selected <- missing_summary %>%
  dplyr::filter(pct_missing <= MISSING_THRESHOLD_PCT) %>%
  dplyr::pull(variable)

# Exclude covariables from correlation targets
params_selected <- setdiff(params_selected, c("Sex", "group", "sample"))

cat(sprintf("\nSelected variables (<=%.0f%% missing): %d\n",
            MISSING_THRESHOLD_PCT, length(params_selected)))
cat(params_selected, sep = "\n")
cat("\n")

# Missing value heatmap (Supplementary Figure)
png("03_plots/sex_sensitivity/Fig_S1_missing_values.png", width = 2400, height = 1800, res = 300)
pheatmap(
  t(is.na(clinical_sub[, vars_for_missing]) * 1),
  color         = c("white", "#B2182B"),
  cluster_rows  = TRUE, cluster_cols = FALSE,
  fontsize_row  = 7, fontsize_col = 8,
  legend_breaks = c(0, 1), legend_labels = c("Present", "Missing"),
  main          = "Missing values per variable and sample", angle_col = 45
)
dev.off()


# =============================================================================
# 5. GENE AVAILABILITY CHECK
# =============================================================================

genes_available <- candidate_genes[candidate_genes %in% rownames(expr_matrix)]
genes_missing_g <- candidate_genes[!candidate_genes %in% rownames(expr_matrix)]

if (length(genes_missing_g) > 0) {
  cat("WARNING: Candidate genes not found in expression matrix:", genes_missing_g, "\n\n")
}

expr_sub <- expr_matrix[genes_available, common_samples, drop = FALSE]

cat(sprintf("Samples with Sex information: %d / %d\n\n",
            sum(!is.na(clinical_sub[common_samples, "Sex"])),
            length(common_samples)))

# =============================================================================
# 6. PARTIAL SPEARMAN CORRELATION (sex-controlled)
#
#    pcor.test(x = expression, y = parameter, z = sex, method = "spearman")
#    Returns partial rho: gene~parameter association with sex effect removed.
#
#    Minimum 6 complete observations required (n >= p + 3, with p = 1 covariate)
# =============================================================================

cor_results <- data.frame()

for (gene in genes_available) {

  gene_expr <- as.numeric(expr_sub[gene, common_samples])
  sex_num   <- as.numeric(clinical_sub[common_samples, "Sex"])  # pcor.test requires numeric

  for (param in params_selected) {

    param_vals   <- as.numeric(clinical_sub[common_samples, param])
    complete_idx <- complete.cases(gene_expr, param_vals, sex_num)
    n_complete   <- sum(complete_idx)

    if (n_complete >= 6) {
      tryCatch({
        pct <- pcor.test(
          x      = gene_expr[complete_idx],
          y      = param_vals[complete_idx],
          z      = sex_num[complete_idx],
          method = "spearman"
        )

        cor_results <- rbind(cor_results, data.frame(
          Gene      = gene,
          Parameter = param,
          rho       = round(pct$estimate, 3),
          p.value   = round(pct$p.value,  6),
          n         = n_complete
        ))
      }, error = function(e) {
        cat(sprintf("  Error: %s ~ %s — %s\n", gene, param, e$message))
      })
    }
  }
}

# Sort by descending |rho|
cor_results <- cor_results %>% dplyr::arrange(dplyr::desc(abs(rho)))

cat("=== TOP CORRELATIONS (|rho| > 0.6) ===\n")
print(cor_results %>% dplyr::filter(abs(rho) > 0.6))

write.csv(
  cor_results,
  "02_results/analisi_final/sex_sensitivity/correlation_results_partial_spearman.csv",
  row.names = FALSE
)


# =============================================================================
# 7. HEATMAP — ALL GENES (full, unfiltered)
#    Genes grouped by family (RPL / RPS / snoRNA)
#    Cells: rho value + significance star below
#    Ordered by mean rho (most positive to most negative)
# =============================================================================

if (nrow(cor_results) > 0) {

  # Order parameters and genes by mean rho (descending)
  param_order <- cor_results %>%
    dplyr::group_by(Parameter) %>%
    dplyr::summarise(mean_rho = mean(rho, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_rho)) %>%
    dplyr::pull(Parameter)

  gene_order <- cor_results %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(mean_rho = mean(rho, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_rho)) %>%
    dplyr::pull(Gene)

  # Assign gene family labels
  heatmap_data <- cor_results %>%
    dplyr::mutate(
      Parameter = factor(Parameter, levels = param_order),
      Gene      = factor(Gene,      levels = gene_order),
      family    = dplyr::case_when(
        grepl("^RPL", Gene)  ~ "RPL",
        grepl("^RPS", Gene)  ~ "RPS",
        grepl("^SNOR", Gene) ~ "snoRNA",
        TRUE                 ~ "Other"
      ),
      sig_label = dplyr::case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE            ~ ""
      ),
      cell_label = ifelse(
        sig_label != "",
        paste0(rho, "\n \n", sig_label),
        as.character(rho)
      )
    )

  # Family display order (top to bottom on Y axis)
  family_order  <- c("RPL", "RPS", "snoRNA", "Other")
  family_colors <- c("RPL" = "#B2182B", "RPS" = "#2166AC",
                     "snoRNA" = "#4DAC26", "Other" = "grey40")

  # Filter: keep only rows where |rho| >= 0.5
  heatmap_data <- heatmap_data %>%
    dplyr::filter(abs(rho) >= 0.5)

  # Re-order genes grouping by family, then by rho within each family
  gene_order_family <- heatmap_data %>%
    dplyr::distinct(Gene, family) %>%
    dplyr::mutate(family = factor(family, levels = family_order)) %>%
    dplyr::arrange(family, Gene) %>%
    dplyr::pull(Gene) %>%
    as.character()

  heatmap_data <- heatmap_data %>%
    dplyr::mutate(
      Gene   = factor(Gene,   levels = gene_order_family),
      family = factor(family, levels = family_order)
    )

  p_heatmap <- ggplot(heatmap_data, aes(x = Parameter, y = Gene, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(
      aes(label = cell_label),
      size = 2, lineheight = 0.85,
      color = ifelse(abs(heatmap_data$rho) > 0.6, "white", "black")
    ) +

    # Dashed separator line between gene families
    geom_hline(
      data = heatmap_data %>%
        dplyr::distinct(Gene, family) %>%
        dplyr::arrange(Gene) %>%
        dplyr::mutate(gene_num = as.numeric(Gene)) %>%
        dplyr::group_by(family) %>%
        dplyr::summarise(max_num = max(gene_num), .groups = "drop") %>%
        dplyr::filter(max_num < max(as.numeric(heatmap_data$Gene))),
      aes(yintercept = max_num + 0.5),
      color = "grey40", linewidth = 0.6, linetype = "dashed",
      inherit.aes = FALSE
    ) +

    scale_fill_gradientn(
      colors   = c("#2166AC", "white", "#B2182B"),
      limits   = c(-1, 1),
      name     = expression("Partial Spearman " * rho),
      na.value = "grey90"
    ) +

    # Nested Y axis (ggh4x)
    scale_y_discrete(guide = guide_axis_nested(delim = "_")) +

    # Family labels on the left, colour-coded by family
    annotate(
      "text",
      x     = rep(-Inf, length(unique(heatmap_data$family))),
      y     = heatmap_data %>%
        dplyr::group_by(family) %>%
        dplyr::summarise(mid = mean(as.numeric(Gene)), .groups = "drop") %>%
        dplyr::arrange(match(family, family_order)) %>%
        dplyr::pull(mid),
      label = heatmap_data %>%
        dplyr::group_by(family) %>%
        dplyr::summarise(mid = mean(as.numeric(Gene)), .groups = "drop") %>%
        dplyr::arrange(match(family, family_order)) %>%
        dplyr::pull(family) %>% as.character(),
      hjust    = 1.3,
      size     = 3,
      fontface = "bold",
      color    = family_colors[
        heatmap_data %>%
          dplyr::group_by(family) %>%
          dplyr::summarise(.groups = "drop") %>%
          dplyr::pull(family) %>%
          as.character()
      ]
    ) +

    labs(
      title    = "Partial Spearman correlations — controlling for sex",
      subtitle = "* p<0.05  ** p<0.01  *** p<0.001",
      caption  = paste0("All samples (n=", length(common_samples), "). |rho| >= 0.5 shown.")
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y     = element_text(size = 9, face = "italic"),
      axis.title      = element_blank(),
      plot.title      = element_text(face = "bold", size = 11),
      plot.subtitle   = element_text(size = 9, color = "grey40"),
      plot.caption    = element_text(size = 8, color = "grey50"),
      legend.position = "right",
      legend.title    = element_text(family = "sans", hjust = 0.5,
                                     margin = margin(b = 5)),
      panel.grid      = element_blank(),
      plot.margin     = margin(10, 10, 10, 40)
    ) +
    coord_cartesian(clip = "off")
}

pdf(
  "03_plots/sex_sensitivity/Fig_heatmap_partial_spearman.pdf",
  width  = 14,
  height = max(4, length(unique(heatmap_data$Gene)) * 0.5 + 2)
)
p_heatmap
dev.off()


# =============================================================================
# 7B. HEATMAP — FILTERED (parameters where at least one gene has |rho| >= 0.5)
#     All genes shown for retained parameters (no blank cells within columns)
# =============================================================================

if (nrow(cor_results) > 0) {

  param_order <- cor_results %>%
    dplyr::group_by(Parameter) %>%
    dplyr::summarise(mean_rho = mean(rho, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_rho)) %>%
    dplyr::pull(Parameter)

  gene_order <- cor_results %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(mean_rho = mean(rho, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_rho)) %>%
    dplyr::pull(Gene)

  heatmap_data <- cor_results %>%
    dplyr::mutate(
      Parameter = factor(Parameter, levels = param_order),
      Gene      = factor(Gene,      levels = gene_order),
      family    = dplyr::case_when(
        grepl("^RPL", Gene)  ~ "RPL",
        grepl("^RPS", Gene)  ~ "RPS",
        grepl("^SNOR", Gene) ~ "snoRNA",
        TRUE                 ~ "Other"
      ),
      sig_label = dplyr::case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE            ~ ""
      ),
      cell_label = ifelse(
        sig_label != "",
        paste0(rho, "\n \n", sig_label),
        as.character(rho)
      )
    )

  family_order  <- c("RPL", "RPS", "snoRNA", "Other")
  family_colors <- c("RPL" = "#B2182B", "RPS" = "#2166AC",
                     "snoRNA" = "#4DAC26", "Other" = "grey40")

  # Keep parameters where at least one gene reaches |rho| >= 0.5
  # All genes are shown for the retained parameters (avoids blank cells)
  params_to_show <- heatmap_data %>%
    dplyr::group_by(Parameter) %>%
    dplyr::summarise(max_rho = max(abs(rho), na.rm = TRUE), .groups = "drop") %>%
    dplyr::filter(max_rho >= 0.5) %>%
    dplyr::pull(Parameter)

  heatmap_data <- heatmap_data %>%
    dplyr::filter(Parameter %in% params_to_show)

  gene_order_family <- heatmap_data %>%
    dplyr::distinct(Gene, family) %>%
    dplyr::mutate(family = factor(family, levels = family_order)) %>%
    dplyr::arrange(family, Gene) %>%
    dplyr::pull(Gene) %>%
    as.character()

  heatmap_data <- heatmap_data %>%
    dplyr::mutate(
      Gene      = factor(Gene,      levels = gene_order_family),
      family    = factor(family,    levels = family_order),
      Parameter = factor(Parameter, levels = param_order[param_order %in% params_to_show])
    )

  p_heatmap <- ggplot(heatmap_data, aes(x = Parameter, y = Gene, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(
      aes(label = cell_label),
      size = 2, lineheight = 0.85,
      color = ifelse(abs(heatmap_data$rho) > 0.6, "white", "black")
    ) +
    geom_hline(
      data = heatmap_data %>%
        dplyr::distinct(Gene, family) %>%
        dplyr::arrange(Gene) %>%
        dplyr::mutate(gene_num = as.numeric(Gene)) %>%
        dplyr::group_by(family) %>%
        dplyr::summarise(max_num = max(gene_num), .groups = "drop") %>%
        dplyr::filter(max_num < max(as.numeric(heatmap_data$Gene))),
      aes(yintercept = max_num + 0.5),
      color = "grey40", linewidth = 0.6, linetype = "dashed",
      inherit.aes = FALSE
    ) +
    scale_fill_gradientn(
      colors   = c("#2166AC", "white", "#B2182B"),
      limits   = c(-1, 1),
      name     = expression("Partial Spearman " * rho),
      na.value = "grey90"
    ) +
    scale_y_discrete(guide = guide_axis_nested(delim = "_")) +
    annotate(
      "text",
      x     = rep(-Inf, length(unique(heatmap_data$family))),
      y     = heatmap_data %>%
        dplyr::group_by(family) %>%
        dplyr::summarise(mid = mean(as.numeric(Gene)), .groups = "drop") %>%
        dplyr::arrange(match(family, family_order)) %>%
        dplyr::pull(mid),
      label = heatmap_data %>%
        dplyr::group_by(family) %>%
        dplyr::summarise(mid = mean(as.numeric(Gene)), .groups = "drop") %>%
        dplyr::arrange(match(family, family_order)) %>%
        dplyr::pull(family) %>% as.character(),
      hjust    = 1.3,
      size     = 3,
      fontface = "bold",
      color    = family_colors[
        heatmap_data %>%
          dplyr::group_by(family) %>%
          dplyr::summarise(.groups = "drop") %>%
          dplyr::pull(family) %>%
          as.character()
      ]
    ) +
    labs(
      title    = "Partial Spearman correlations — controlling for sex (filtered)",
      subtitle = "* p<0.05  ** p<0.01  *** p<0.001",
      caption  = paste0(
        "All samples (n=", length(common_samples), "). ",
        "Parameters shown where at least one gene reaches |rho| >= 0.5. ",
        "All genes displayed for retained parameters."
      )
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y     = element_text(size = 9, face = "italic"),
      axis.title      = element_blank(),
      plot.title      = element_text(face = "bold", size = 11),
      plot.subtitle   = element_text(size = 9, color = "grey40"),
      plot.caption    = element_text(size = 8, color = "grey50"),
      legend.position = "right",
      legend.title    = element_text(family = "sans", hjust = 0.5,
                                     margin = margin(b = 5)),
      panel.grid      = element_blank(),
      plot.margin     = margin(10, 10, 10, 40)
    ) +
    coord_cartesian(clip = "off")
}

pdf(
  "03_plots/sex_sensitivity/Fig_heatmap_partial_spearman_filtered.pdf",
  width  = 14,
  height = max(4, length(unique(heatmap_data$Gene)) * 0.5 + 2)
)
p_heatmap
dev.off()


# =============================================================================
# 8. SCATTER PLOTS — TOP HITS (individual genes)
#    Y axis: raw RMA-normalized log2 expression
#    Sex adjustment is internal to pcor.test (not visualised as residuals)
#    Point shape indicates sex for transparency
# =============================================================================

top_hits <- cor_results %>%
  dplyr::filter(abs(rho) >= 0.5) %>%
  dplyr::arrange(dplyr::desc(abs(rho))) %>%
  head(9)

if (nrow(top_hits) == 0) {
  cat("WARNING: No correlations with |rho| > 0.5. Showing top 9 by absolute rho.\n")
  top_hits <- cor_results %>% dplyr::arrange(dplyr::desc(abs(rho))) %>% head(9)
}

scatter_list <- list()

for (i in seq_len(nrow(top_hits))) {
  gene  <- top_hits$Gene[i]
  param <- top_hits$Parameter[i]
  rho   <- top_hits$rho[i]
  pval  <- top_hits$p.value[i]
  n_obs <- top_hits$n[i]

  plot_data <- data.frame(
    x      = as.numeric(clinical_sub[common_samples, param]),
    y      = as.numeric(expr_sub[gene, common_samples]),
    sample = common_samples,
    group  = clinical_sub[common_samples, "group"],
    sex    = clinical_sub[common_samples, "Sex"]
  ) %>%
    dplyr::filter(complete.cases(.))

  p <- ggplot(plot_data, aes(x = x, y = y, color = group)) +
    geom_smooth(method = "lm", se = TRUE,
                color = "grey50", linetype = "dashed", linewidth = 0.8,
                inherit.aes = FALSE, aes(x = x, y = y)) +
    geom_point(aes(shape = sex), size = 3.5, alpha = 0.9) +
    geom_text(aes(label = sample),
              vjust = -0.9, size = 2.2, show.legend = FALSE) +
    scale_color_manual(values = group_colors) +
    scale_shape_manual(
      values = c("0" = 16, "1" = 17),  # adjust if Sex levels differ
      name   = "Sex",
      labels = c("0" = "Female", "1" = "Male")
    ) +
    annotate("text",
             x = -Inf, y = Inf, hjust = -0.1, vjust = 1.6,
             label = paste0("\u03c1 = ", rho, "  p = ", pval, "\nn = ", n_obs),
             size = 3, fontface = "italic", color = "grey30") +
    labs(
      x     = param,
      y     = paste0(gene, " (RMA-normalized log2 expression)"),
      title = paste0(gene, " ~ ", param),
      color = "Group"
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title      = element_text(face = "bold", size = 9),
      legend.position = "bottom",
      legend.key.size = unit(0.4, "cm"),
      legend.text     = element_text(size = 7)
    )

  scatter_list[[i]] <- p
}

if (length(scatter_list) > 0) {
  combined_plot <- wrap_plots(scatter_list, ncol = 3) +
    plot_annotation(
      title   = "Partial Spearman correlations — gene expression vs clinical parameters",
      caption = paste0(
        "\u03c1: partial Spearman's rho (ppcor). ",
        "Dashed line: linear fit (all samples). n=", length(common_samples), "."
      ),
      theme = theme(
        plot.title   = element_text(face = "bold", size = 12),
        plot.caption = element_text(size = 8, color = "grey50")
      )
    )
}

pdf("03_plots/sex_sensitivity/Fig_scatter_partial_spearman.pdf",
    width = 14, height = 12)
combined_plot
dev.off()


# =============================================================================
# 9. FAMILY-LEVEL ANALYSIS (RPL_mean and RPS_mean)
# =============================================================================
# Expression values of all RPL and RPS candidate genes are averaged to generate
# family-level scores. This reduces probe-level noise and increases robustness.
# Partial Spearman correlation is then applied to these family scores.
# =============================================================================

# ── 9.1 Compute mean expression per family ──────────────────────────────────

# Identify probes for each family from available genes
genes_RPL <- genes_available[grepl("^RPL|^RPSA", genes_available)]
genes_RPS <- genes_available[grepl("^RPS",        genes_available)]

cat("RPL genes included in family mean:", genes_RPL, "\n")
cat("RPS genes included in family mean:", genes_RPS, "\n\n")

# Compute per-sample mean across all family members
RPL_mean <- colMeans(expr_sub[genes_RPL, common_samples, drop = FALSE], na.rm = TRUE)
RPS_mean <- colMeans(expr_sub[genes_RPS, common_samples, drop = FALSE], na.rm = TRUE)

# Combine into a 2 x n matrix
expr_grouped <- rbind(RPL = RPL_mean, RPS = RPS_mean)


# ── 9.2 Partial Spearman correlation for family scores ──────────────────────

sex_num <- as.numeric(clinical_sub[common_samples, "Sex"])

cor_grouped <- data.frame()

for (family in rownames(expr_grouped)) {

  family_expr <- as.numeric(expr_grouped[family, common_samples])

  for (param in params_selected) {

    param_vals   <- as.numeric(clinical_sub[common_samples, param])
    complete_idx <- complete.cases(family_expr, param_vals, sex_num)
    n_complete   <- sum(complete_idx)

    if (n_complete >= 6) {
      tryCatch({
        pct <- pcor.test(
          x      = family_expr[complete_idx],
          y      = param_vals[complete_idx],
          z      = sex_num[complete_idx],
          method = "spearman"
        )
        cor_grouped <- rbind(cor_grouped, data.frame(
          Gene      = family,   # "RPL" or "RPS"
          Parameter = param,
          rho       = round(pct$estimate, 3),
          p.value   = signif(pct$p.value, 3),
          n         = n_complete
        ))
      }, error = function(e) {
        cat(sprintf("Error: %s ~ %s\n", family, param))
      })
    }
  }
}

cor_grouped <- cor_grouped %>% dplyr::arrange(dplyr::desc(abs(rho)))

cat("=== TOP FAMILY-LEVEL CORRELATIONS (|rho| > 0.6) ===\n")
print(cor_grouped %>% dplyr::filter(abs(rho) > 0.6))

write.csv(cor_grouped, "correlation_results_grouped_families.csv", row.names = FALSE)


# ── 9.3 Heatmap — family scores, unfiltered ─────────────────────────────────

param_order_g <- cor_grouped %>%
  dplyr::group_by(Parameter) %>%
  dplyr::summarise(mean_rho = mean(rho, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(mean_rho)) %>%
  dplyr::pull(Parameter)

heatmap_grouped <- cor_grouped %>%
  dplyr::mutate(
    Parameter = factor(Parameter, levels = param_order_g),
    Gene      = factor(Gene, levels = c("RPL", "RPS")),
    sig_label = dplyr::case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    cell_label = ifelse(
      sig_label != "",
      paste0(rho, "\n \n", sig_label),
      as.character(rho)
    )
  ) %>%
  dplyr::filter(abs(rho) >= 0.5)

p_heatmap_grouped <- ggplot(heatmap_grouped,
                             aes(x = Parameter, y = Gene, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = cell_label),
    size = 2.8, lineheight = 0.85,
    color = ifelse(abs(heatmap_grouped$rho) > 0.6, "white", "black")
  ) +
  scale_fill_gradientn(
    colors   = c("#2166AC", "white", "#B2182B"),
    limits   = c(-1, 1),
    name     = expression("Partial Spearman " * rho),
    na.value = "grey90"
  ) +
  labs(
    title    = "Partial Spearman correlations — RPL and RPS family means vs clinical parameters",
    subtitle = "* p<0.05  ** p<0.01  *** p<0.001  |  Values = mean expression across all family members",
    caption  = paste0(
      "RPL family (n=", length(genes_RPL), " genes): ", paste(genes_RPL, collapse = ", "), "\n",
      "RPS family (n=", length(genes_RPS), " genes): ", paste(genes_RPS, collapse = ", "), "\n",
      "All samples (n=", length(common_samples), "). Partial correlation controlling for sex. |rho| >= 0.5 shown."
    )
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y     = element_text(size = 11, face = "bold.italic"),
    axis.title      = element_blank(),
    plot.title      = element_text(face = "bold", size = 11),
    plot.subtitle   = element_text(size = 9, color = "grey40"),
    plot.caption    = element_text(size = 7, color = "grey50", hjust = 0),
    legend.position = "right",
    legend.title    = element_text(family = "sans", hjust = 0, margin = margin(b = 5)),
    panel.grid      = element_blank()
  )

ggsave("Fig_heatmap_grouped_families.png",
       p_heatmap_grouped, width = 12, height = 3, dpi = 300)


# ── 9.4 Heatmap — family scores, filtered ───────────────────────────────────
# Parameters shown where at least one family reaches |rho| >= 0.5
# Both families always shown for retained parameters (no blank cells)

params_to_show <- cor_grouped %>%
  dplyr::group_by(Parameter) %>%
  dplyr::summarise(max_rho = max(abs(rho), na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(max_rho >= 0.5) %>%
  dplyr::pull(Parameter)

heatmap_grouped_filt <- cor_grouped %>%
  dplyr::filter(Parameter %in% params_to_show) %>%
  dplyr::mutate(
    Parameter = factor(Parameter,
                       levels = param_order_g[param_order_g %in% params_to_show]),
    Gene      = factor(Gene, levels = c("RPL", "RPS")),
    sig_label = dplyr::case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    cell_label = ifelse(
      sig_label != "",
      paste0(rho, "\n \n", sig_label),
      as.character(rho)
    )
  )

p_heatmap_grouped_filt <- ggplot(heatmap_grouped_filt,
                                  aes(x = Parameter, y = Gene, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = cell_label),
    size = 2.8, lineheight = 0.85,
    color = ifelse(abs(heatmap_grouped_filt$rho) > 0.6, "white", "black")
  ) +
  scale_fill_gradientn(
    colors   = c("#2166AC", "white", "#B2182B"),
    limits   = c(-1, 1),
    name     = expression("Partial Spearman " * rho),
    na.value = "grey90"
  ) +
  labs(
    title    = "Partial Spearman correlations — RPL and RPS family means vs clinical parameters (filtered)",
    subtitle = "* p<0.05  ** p<0.01  *** p<0.001  |  Parameters shown where at least one family has |rho| >= 0.5",
    caption  = paste0(
      "RPL family (n=", length(genes_RPL), " genes): ", paste(genes_RPL, collapse = ", "), "\n",
      "RPS family (n=", length(genes_RPS), " genes): ", paste(genes_RPS, collapse = ", "), "\n",
      "All samples (n=", length(common_samples), "). Partial correlation controlling for sex."
    )
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y     = element_text(size = 11, face = "bold.italic"),
    axis.title      = element_blank(),
    plot.title      = element_text(face = "bold", size = 11),
    plot.subtitle   = element_text(size = 9, color = "grey40"),
    plot.caption    = element_text(size = 7, color = "grey50", hjust = 0),
    legend.position = "right",
    legend.title    = element_text(family = "sans", hjust = 0, margin = margin(b = 5)),
    panel.grid      = element_blank()
  )

ggsave("Fig_heatmap_grouped_families_filtered.png",
       p_heatmap_grouped_filt, width = 12, height = 3, dpi = 300)


# ── 9.5 Scatter plots — family scores vs top parameters ─────────────────────
# Top 4 parameters per family (RPL first, then RPS), ordered by |rho|
# Family title colour: RPL = red (#B2182B), RPS = blue (#2166AC)

n_per_family <- 4

top_grouped <- cor_grouped %>%
  dplyr::filter(abs(rho) > 0.5) %>%
  dplyr::mutate(family_order = ifelse(Gene == "RPL", 1, 2)) %>%
  dplyr::group_by(Gene) %>%
  dplyr::slice_max(order_by = abs(rho), n = n_per_family, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(family_order, dplyr::desc(abs(rho))) %>%
  dplyr::select(-family_order)

# Fallback if no correlations exceed threshold
if (nrow(top_grouped) == 0) {
  top_grouped <- cor_grouped %>%
    dplyr::mutate(family_order = ifelse(Gene == "RPL", 1, 2)) %>%
    dplyr::arrange(family_order, dplyr::desc(abs(rho))) %>%
    dplyr::select(-family_order) %>%
    head(6)
}

cat(sprintf("Family scatter panels: RPL = %d, RPS = %d, Total = %d\n",
            sum(top_grouped$Gene == "RPL"),
            sum(top_grouped$Gene == "RPS"),
            nrow(top_grouped)))

scatter_grouped_list <- list()

for (i in seq_len(nrow(top_grouped))) {
  family <- as.character(top_grouped$Gene[i])
  param  <- top_grouped$Parameter[i]
  rho    <- top_grouped$rho[i]
  pval   <- top_grouped$p.value[i]
  n_obs  <- top_grouped$n[i]

  y_vals <- as.numeric(expr_grouped[family, common_samples])

  plot_data <- data.frame(
    x      = as.numeric(clinical_sub[common_samples, param]),
    y      = y_vals,
    sample = common_samples,
    group  = clinical_sub[common_samples, "group"],
    sex    = clinical_sub[common_samples, "Sex"]
  ) %>% dplyr::filter(complete.cases(.))

  # Title colour indicates family: RPL = red, RPS = blue
  family_color <- if (family == "RPL") "#B2182B" else "#2166AC"

  p <- ggplot(plot_data, aes(x = x, y = y)) +
    geom_smooth(method = "lm", se = TRUE,
                color = "grey50", linetype = "dashed", linewidth = 0.8,
                inherit.aes = FALSE, aes(x = x, y = y)) +
    geom_point(aes(shape = group, color = group), size = 3.5, alpha = 0.9) +
    geom_text(aes(label = sample, color = group),
              vjust = -0.9, size = 2.2, show.legend = FALSE) +
    scale_color_manual(values = group_colors) +
    scale_shape_manual(
      values = c("CTEPH" = 16, "PAH" = 17, "Control" = 15),
      name   = "Group"
    ) +
    annotate("text",
             x = -Inf, y = Inf, hjust = -0.1, vjust = 1.6,
             label = paste0("\u03c1 = ", rho, "  p = ", pval, "\nn = ", n_obs),
             size = 3, fontface = "italic", color = "grey30") +
    labs(
      x     = param,
      y     = paste0(family, " family\n(mean log2 expression)"),
      title = paste0(family, " family ~ ", param),
      color = "Group"
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title      = element_text(face = "bold", size = 9, color = family_color),
      legend.position = "bottom",
      legend.key.size = unit(0.4, "cm"),
      legend.text     = element_text(size = 7)
    )

  scatter_grouped_list[[i]] <- p
}

if (length(scatter_grouped_list) > 0) {
  combined_grouped <- wrap_plots(scatter_grouped_list, ncol = 4) +
    plot_annotation(
      title    = "Ribosomal family mean expression vs clinical parameters",
      subtitle = paste0(
        "RPL family (red titles): mean of ", length(genes_RPL), " genes. ",
        "RPS family (blue titles): mean of ", length(genes_RPS), " genes.\n",
        "Partial Spearman \u03c1 controlling for sex. Shape = disease group."
      ),
      caption  = "\u03c1: partial Spearman's rho (ppcor). Dashed line: linear fit.",
      theme    = theme(
        plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        plot.caption  = element_text(size = 8, color = "grey50")
      )
    )

  ggsave("Fig_scatter_grouped_families.png",
         combined_grouped,
         width  = 18,
         height = ceiling(nrow(top_grouped) / 4) * 4.5,
         dpi    = 300)
  cat("Grouped family scatter plots saved.\n\n")
}
