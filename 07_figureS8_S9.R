# ==============================================================================
# Script 07: Figure S8 (rescaled HRs) and Figure S9 (20-gene forest)
#
# Author:   Getachew Derbew Belay
# Date:     October 2026
# Platform: R 4.5.2
# Contact:  derbew1953@gmail.com / d610113006@tmu.edu.tw
#
# Purpose:  Generate the combined forest plot (Figure S8) comparing the adjusted
#           hazard ratios of expression vs methylation signatures, and the full
#           multivariable Cox forest plot (Figure S9) showing all 20 LASSO genes.
#
# Produces: FigureS8_Forest_Expression_vs_Methylation.pdf
#           FigureS9_Multivariable_Forest.pdf
#
# Save as:  07_figureS8_S9.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/
#
# Depends:  nsclc_base.rds, nsclc_script03_signature.rds,
#           nsclc_script04b_hr.rds
#
# Seeds:    456
# ==============================================================================

setwd("~/TCGA_Local_Storage/TCGA NSCLC Project")

## ---- Load inputs -------------------------------------------------------------
base <- readRDS("nsclc_base.rds");               list2env(base, envir = environment())
sig  <- readRDS("nsclc_script03_signature.rds"); list2env(sig,  envir = environment())
hr   <- readRDS("nsclc_script04b_hr.rds")

library(survival); library(glmnet); library(ggplot2); library(dplyr)

## ---- Shared utility (self-contained) -----------------------------------------
calc_risk <- function(data, genes, coefs) {
  genes <- intersect(genes, colnames(data))
  if (length(genes) == 0) return(rep(0, nrow(data)))
  cf <- as.numeric(coefs[genes])
  as.numeric(as.matrix(data[, genes, drop = FALSE]) %*% cf)
}

## =============================================================================
## FIGURE S8 — Corrected HRs: Expression vs Methylation
## =============================================================================
cat("\n=== FIGURE S8: Corrected HR Forest Plot ===\n")

# ---- 1. Extract the 4 adjusted rows per signature ---------------------------
# 04b produces these Model labels:
#   "Univariate (whole cohort)"
#   "Adjusted (whole cohort, in-sample)"
#   "Adjusted (test set only, out-of-sample)"
# We keep only the two "Adjusted" models, and only per IQR / per SD.
extract_risk_adj <- function(tbl, signature_label) {
  d <- tbl[tbl$Scale %in% c("per IQR", "per SD") &
             grepl("^Adjusted", tbl$Model), ]
  if (nrow(d) == 0) {
    stop("No adjusted rows matched for ", signature_label)
  }
  data.frame(
    Signature = signature_label,
    Model     = d$Model,
    Scale     = d$Scale,
    HR        = d$HR,
    Lower     = d$CI_lower,
    Upper     = d$CI_upper,
    P         = d$P_value,
    stringsAsFactors = FALSE
  )
}

s8_data <- rbind(
  extract_risk_adj(hr$expr$table, "Expression"),
  extract_risk_adj(hr$meth$table, "Methylation")
)

# ---- 2. Annotate cohort, scale, significance --------------------------------
s8_data$Cohort_Label <- factor(
  ifelse(grepl("in-sample", s8_data$Model),
         "A. Whole cohort (in-sample)",
         "B. Test set (out-of-sample)"),
  levels = c("A. Whole cohort (in-sample)",
             "B. Test set (out-of-sample)")
)

s8_data$Sig <- ifelse(s8_data$P < 0.001, "***",
                      ifelse(s8_data$P < 0.01,  "**",
                             ifelse(s8_data$P < 0.05,  "*", "")))

s8_data$Scale    <- factor(s8_data$Scale,    levels = c("per IQR", "per SD"))
s8_data$Signature <- factor(s8_data$Signature, levels = c("Expression", "Methylation"))

s8_data$Row_Label <- paste0(s8_data$Signature, " (", s8_data$Scale, ") ", s8_data$Sig)
s8_data$HR_text   <- sprintf("%.2f (%.2f-%.2f)",
                             s8_data$HR, s8_data$Lower, s8_data$Upper)

# ---- 3. Plot ----------------------------------------------------------------
text_x_pos <- max(s8_data$Upper, na.rm = TRUE) * 2.5

p_s8 <- ggplot(s8_data,
               aes(x = HR, y = reorder(Row_Label, HR),
                   color = Signature, shape = Scale)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_point(size = 5, position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper),
                 height = 0.2, linewidth = 1.2,
                 position = position_dodge(width = 0.5)) +
  geom_text(aes(x = text_x_pos, label = HR_text),
            position = position_dodge(width = 0.5),
            hjust = 0, size = 4, color = "black", show.legend = FALSE) +
  facet_grid(Cohort_Label ~ ., scales = "free_y", space = "free_y") +
  scale_x_log10(limits = c(0.5, text_x_pos * 3)) +
  scale_color_manual(values = c("Expression"  = "#E7B800",
                                "Methylation" = "#2E86AB")) +
  labs(
    x        = "Adjusted Hazard Ratio (95% CI, log scale)",
    y        = "",
    title    = "Figure S8. Adjusted HRs: Expression vs. Methylation Risk Scores",
    subtitle = "Adjusted for age, gender, pathologic stage, and histology; HR per IQR (primary) and per SD"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "gray95"),
        strip.text       = element_text(face = "bold", size = 12),
        plot.title       = element_text(face = "bold", size = 14),
        plot.subtitle    = element_text(size = 11))

ggsave("FigureS8_Forest_Expression_vs_Methylation.pdf",
       p_s8, width = 12, height = 8)
cat("Saved: FigureS8_Forest_Expression_vs_Methylation.pdf\n")

## =============================================================================
## FIGURE S9 — Multivariable Cox including the 20 expression LASSO genes
## (REVISED: Standardized genes for interpretable HRs)
## =============================================================================
cat("\n=== FIGURE S9: Multivariable Forest Plot ===\n")

# ---- 1. Re-fit expression LASSO to obtain exact coefficients ----
available_genes <- intersect(candidate_genes, colnames(train_data))
x_tr <- as.matrix(train_data[, available_genes])
y_tr <- Surv(train_data$os_time, train_data$os_event)

set.seed(456)
cv_expr <- cv.glmnet(x_tr, y_tr, family = "cox", alpha = 1,
                     nfolds = 10, type.measure = "C")
cf_min  <- as.matrix(coef(cv_expr, s = "lambda.min"))
sel_min <- rownames(cf_min)[cf_min[, 1] != 0]

# ---- 2. Prepare cohort and STANDARDIZE the selected genes ----
cohort$stage_group <- factor(cohort$stage_group, levels = c("Early", "Advanced", "Unknown"))
cohort$histology   <- factor(cohort$histology, levels = c("LUAD", "LUSC"))
cohort$gender      <- factor(cohort$gender)

# Standardize genes to Z-scores (mean=0, SD=1) for interpretable HRs (per 1 SD)
cohort_scaled <- cohort
cohort_scaled[, sel_min] <- scale(cohort_scaled[, sel_min, drop = FALSE])

# ---- 3. Fit multivariable Cox (clinical + standardized genes) ----
if (length(sel_min) > 0) {
  # CRITICAL FIX: Added backticks around sel_min to handle names like MIR124-2
  f_multi <- as.formula(paste(
    "Surv(os_time, os_event) ~",
    paste(c("age", "gender", "stage_group", "histology", 
            paste0("`", sel_min, "`")), collapse = " + ")
  ))
} else {
  f_multi <- Surv(os_time, os_event) ~ age + gender + stage_group + histology
}

multi_expr <- coxph(f_multi, data = cohort_scaled)
cat("\nMultivariable Cox summary:\n")
print(summary(multi_expr))

# ---- 4. Extract coefficients for forest plot ----
coef_df <- data.frame(
  Variable = rownames(summary(multi_expr)$coefficients),
  HR       = exp(coef(multi_expr)),
  Lower    = exp(confint(multi_expr))[, 1],
  Upper    = exp(confint(multi_expr))[, 2],
  P        = summary(multi_expr)$coefficients[, "Pr(>|z|)"],
  stringsAsFactors = FALSE
)

# Clean up names (remove backticks so they look nice in the plot)
coef_df$Variable <- gsub("`", "", coef_df$Variable)
coef_df$Variable <- gsub("stage_group", "Stage: ", coef_df$Variable)
coef_df$Variable <- gsub("histology", "Histology: ", coef_df$Variable)
coef_df$Variable <- gsub("gendermale", "Male Gender", coef_df$Variable)

coef_df$Sig <- ifelse(coef_df$P < 0.001, "***",
                      ifelse(coef_df$P < 0.01,  "**",
                             ifelse(coef_df$P < 0.05,  "*", "")))
coef_df$Label   <- paste0(coef_df$Variable, " ", coef_df$Sig)

# Use 3 decimal places so they don't all look like 1.00
coef_df$HR_text <- sprintf("%.3f (%.3f-%.3f)", coef_df$HR, coef_df$Lower, coef_df$Upper)

# ---- 5. Plot ----
text_x_s9 <- max(coef_df$Upper, na.rm = TRUE) * 2.0

p_s9 <- ggplot(coef_df, aes(x = HR, y = reorder(Label, HR))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_point(size = 4, color = "#2E86AB") +
  geom_errorbar(aes(xmin = Lower, xmax = Upper), 
                width = 0.25, linewidth = 1.0, color = "#2E86AB", orientation = "y") +
  geom_text(aes(x = text_x_s9, label = HR_text),
            hjust = 0, size = 3.5, color = "black") +
  scale_x_log10(limits = c(0.1, text_x_s9 * 2.5)) +
  labs(
    x        = "Hazard Ratio (95% CI, log scale)",
    y        = "",
    title    = "Figure S9. Multivariable Cox Regression (Clinical + LASSO Genes)",
    subtitle = sprintf("n = %d; events = %d. Clinical covariates + %d genes (standardized).", 
                       nrow(cohort), sum(cohort$os_event), length(sel_min))
  ) +
  theme_bw(base_size = 13) +
  theme(panel.grid.major.y = element_blank(),
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10))

ggsave("FigureS9_Multivariable_Forest.pdf", p_s9, width = 11, height = 9)
cat("Saved: FigureS9_Multivariable_Forest.pdf\n")

cat("\nScript 07 complete: Figure S8 and S9 saved.\n")