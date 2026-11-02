# ==============================================================================
# Script 01: Methylation-Expression Correlation Analysis
#
# Author:   Getachew Derbew Belay
# Date:     October 2026
# Platform: R 4.5.2
# Contact:  derbew1953@gmail.com / d610113006@tmu.edu.tw
#
# Purpose:  Compute Spearman correlations between promoter methylation and
#           gene expression overall and by subgroup; produce Figure 2 (2x4
#           scatter), Figure S1 (all genes), Figure 3 (stratified heatmap).
# Produces: Table_Methylation_Expression_Correlation.csv, TableS2_Stratified_Correlations.csv, and,
#.          Figure2_Correlation_Scatter.pdf, FigureS1_Correlation_All_Genes.pdf, Figure3_Stratified_Correlation_Heatmap.pdf
#
# Save as:  01_correlation.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/
#
# Depends:  nsclc_base.rds  (from Script 00)
#
# Seeds:    none
# ==============================================================================

setwd("~/TCGA_Local_Storage/TCGA NSCLC Project")
base <- readRDS("nsclc_base.rds")
list2env(base, envir = environment())
library(ggplot2); library(patchwork); library(pheatmap); library(dplyr)

## ---- Align expression and methylation matrices -------------------------------
common_genes    <- intersect(rownames(gene_meth), rownames(master_nsclc_matrix))
common_patients <- intersect(colnames(gene_meth), colnames(master_nsclc_matrix))
meth_sub <- gene_meth[common_genes, common_patients]
expr_sub <- master_nsclc_matrix[common_genes, common_patients]
cat("Aligned:", length(common_genes), "genes x", length(common_patients), "patients\n")

## ---- Overall Spearman correlation --------------------------------------------
cor_results <- data.frame(gene = common_genes, rho = NA_real_, p_value = NA_real_,
                          n = NA_integer_, stringsAsFactors = FALSE)
for (i in seq_along(common_genes)) {
  g <- common_genes[i]
  m <- as.numeric(meth_sub[g, ]); e <- as.numeric(expr_sub[g, ])
  ok <- is.finite(m) & is.finite(e)
  if (sum(ok) < 10) next
  if (sd(m[ok]) == 0 || sd(e[ok]) == 0) next
  tt <- tryCatch(cor.test(m[ok], e[ok], method = "spearman", exact = FALSE),
                 error = function(e) NULL)
  if (is.null(tt)) next
  cor_results$rho[i]     <- tt$estimate
  cor_results$p_value[i] <- tt$p.value
  cor_results$n[i]       <- sum(ok)
}
cor_results$fdr <- p.adjust(cor_results$p_value, method = "BH")
cor_results <- cor_results[order(cor_results$fdr), ]
write.csv(cor_results, "Table_Methylation_Expression_Correlation.csv", row.names = FALSE)
print(cor_results)

## ---- Figure 2: 2x4 scatter grid for the 8 highest-|rho| genes ---------------
sel8 <- c("CHRNB4","DAPK1","KDM1A","SPARC","WIF1","SHOX2","KDM2A","PCDH10")

plot_list <- list()
for (g in sel8) {
  df <- data.frame(Methylation = as.numeric(meth_sub[g, ]),
                   Expression  = as.numeric(expr_sub[g, ]))
  df <- df[complete.cases(df), ]
  rho <- cor_results$rho[cor_results$gene == g][1]
  fdr <- cor_results$fdr[cor_results$gene == g][1]
  lab <- if (is.na(fdr)) "FDR = NA" else if (fdr < 0.001) "FDR < 0.001" else sprintf("FDR = %.3f", fdr)
  p <- ggplot(df, aes(Methylation, Expression)) +
    geom_point(alpha = 0.4, size = 0.8, color = "gray30") +
    geom_smooth(method = "lm", se = FALSE, color = "#D55E00", linewidth = 1) +
    labs(title = g,
         x = "Promoter Methylation (beta)",
         y = "Expression (log2 TPM)") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
             label = paste0("rho = ", sprintf("%.2f", rho), "\n", lab),
             size = 3.5) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  plot_list[[g]] <- p
}
fig2 <- wrap_plots(plot_list, ncol = 4) +
  plot_annotation(title = "Promoter Methylation vs Gene Expression (TCGA NSCLC)",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16)))
ggsave("Figure2_Correlation_Scatter.pdf", fig2, width = 16, height = 8)

## ---- Figure S1: all genes ----------------------------------------------------
pdf("FigureS1_Correlation_All_Genes.pdf", width = 14, height = 12)
par(mfrow = c(5, 5), mar = c(4, 4, 3, 1))
for (g in rownames(meth_sub)) {
  m <- as.numeric(meth_sub[g, ]); e <- as.numeric(expr_sub[g, ])
  ok <- is.finite(m) & is.finite(e)
  if (sum(ok) < 10) next
  plot(m[ok], e[ok],
       xlab = "Promoter Methylation (beta)",
       ylab = "Expression (log2 TPM)",
       main = g,
       pch = 19, col = rgb(0,0,0,0.3), cex = 0.5)
  abline(lm(e[ok] ~ m[ok]), col = "red", lwd = 1.5)
  rho <- cor_results$rho[cor_results$gene == g][1]
  fdr <- cor_results$fdr[cor_results$gene == g][1]
  rho_lab <- if (is.na(rho)) "rho = NA" else sprintf("rho = %+.2f", rho)
  fdr_lab <- if (is.na(fdr)) "FDR = NA" else if (fdr < 0.001) "FDR < 0.001" else sprintf("FDR = %.3f", fdr)
  legend("topright", legend = paste0(rho_lab, "\n", fdr_lab), bty = "n", cex = 0.7)
}
dev.off()

## ---- Figure 3: stratified correlation heatmap (6 strata) --------------------
cov_meta <- covariates_clean[, c("patient_id","stage_group","histology")]
surv_meta <- cov_meta[cov_meta$patient_id %in% common_patients, ]

strata <- list(
  "Overall"           = surv_meta$patient_id,
  "LUAD"              = surv_meta$patient_id[surv_meta$histology == "LUAD"],
  "LUSC"              = surv_meta$patient_id[surv_meta$histology == "LUSC"],
  "Early (I/II)"      = surv_meta$patient_id[surv_meta$stage_group == "Early"],
  "Advanced (III/IV)" = surv_meta$patient_id[surv_meta$stage_group == "Advanced"],
  "Unknown (missing)" = surv_meta$patient_id[surv_meta$stage_group == "Unknown"]
)

cat("\n=== Stratum sizes ===\n")
for (s in names(strata)) cat(sprintf("%-20s : %d\n", s, length(strata[[s]])))

cat("\nSum of stage strata (Early + Advanced + Unknown):",
    length(strata[["Early (I/II)"]]) +
      length(strata[["Advanced (III/IV)"]]) +
      length(strata[["Unknown (missing)"]]), "\n")

rho_mat <- matrix(NA_real_, nrow = length(common_genes), ncol = length(strata),
                  dimnames = list(common_genes, names(strata)))
fdr_mat <- rho_mat
for (s in names(strata)) {
  pats <- strata[[s]]
  if (length(pats) < 10) next
  pvals <- numeric(length(common_genes)); pvals[] <- NA
  for (i in seq_along(common_genes)) {
    g <- common_genes[i]
    m <- as.numeric(meth_sub[g, pats]); e <- as.numeric(expr_sub[g, pats])
    ok <- is.finite(m) & is.finite(e)
    if (sum(ok) < 10) next
    if (sd(m[ok]) == 0 || sd(e[ok]) == 0) next
    tt <- tryCatch(cor.test(m[ok], e[ok], method = "spearman", exact = FALSE),
                   error = function(e) NULL)
    if (is.null(tt)) next
    rho_mat[i, s] <- tt$estimate
    pvals[i]      <- tt$p.value
  }
  ok <- !is.na(pvals)
  if (any(ok)) fdr_mat[ok, s] <- p.adjust(pvals[ok], method = "BH")
}

## Save Table S2 with both rho and FDR
write.csv(data.frame(gene = rownames(rho_mat), rho_mat, check.names = FALSE),
          "TableS2_Stratified_Correlations.csv", row.names = FALSE)
write.csv(data.frame(gene = rownames(fdr_mat), fdr_mat, check.names = FALSE),
          "TableS2b_Stratified_FDR.csv", row.names = FALSE)

## Print key genes
cat("\n=== Key genes across strata (rho values) ===\n")
key_genes <- c("CHRNB4","DAPK1","KDM1A","SPARC","WIF1","SHOX2","KDM2A",
               "PCDH10","RASSF1","BRMS1","EGLN2","TGFBI","LRRC3B")
key_genes <- intersect(key_genes, rownames(rho_mat))
print(round(rho_mat[key_genes, , drop = FALSE], 3))

cat("\n=== Key genes across strata (FDR values) ===\n")
print(signif(fdr_mat[key_genes, , drop = FALSE], 3))

## ---- Heatmap ---------------------------------------------------------------
plot_mat <- rho_mat; plot_mat[is.na(plot_mat)] <- 0
lim <- max(abs(range(plot_mat)))
palette_colors <- colorRampPalette(c("#2166AC","#67A9CF","#F7F7F7","#EF8A62","#B2182B"))(100)
breaks <- seq(-lim, lim, length.out = 101)

pdf("Figure3_Stratified_Correlation_Heatmap.pdf", width = 9, height = 10)
pheatmap(plot_mat, cluster_rows = TRUE, cluster_cols = FALSE,
         color = palette_colors, breaks = breaks,
         display_numbers = TRUE, number_format = "%.2f", number_color = "black",
         fontsize = 11, fontsize_number = 8,
         main = "Spearman Correlation: Promoter Methylation vs Expression",
         angle_col = 45, border_color = "gray80", legend = TRUE)
dev.off()

cat("\nScript 01 complete.\n")