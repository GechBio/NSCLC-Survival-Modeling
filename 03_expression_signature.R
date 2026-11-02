# ==============================================================================
# Script 03: Expression Signature, Clinical Comparison, and Negative Findings
#
# Author:   Getachew Derbew Belay
# Date:     October 2026
# Platform: R 4.5.2
#
# Purpose:  Develop the expression LASSO signature (lambda.min/lambda.1se);
#           demonstrate that the gene panel has no prognostic signal by
#           comparing Clinical-Only vs Integrated (Clinical + Genes) models;
#           generate replacement figures (Clinical KM, ROC, LASSO paths).
#
# Save as:  03_expression_signature.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/
#
# Output:   TableS3_Univariate_Cox.csv
#           Table_TimeROC_AUC.csv
#           Table_Multivariable_Cox_Expression.csv
#           Figure5_KM_Clinical.pdf             
#           FigureS2_TimeROC_Clinical_vs_Int.pdf
#           FigureS3_LASSO_Paths.pdf            (Proves 0 genes at lambda.min)
#           FigureS6_KM_Clinical_Stratified.pdf 
#           Figure_Forest_Plot_Expression.pdf
#           nsclc_script03_signature.rds
#
# Depends:  nsclc_base.rds  (from Script 00)
# ==============================================================================

setwd("~/TCGA_Local_Storage/TCGA NSCLC Project")
base <- readRDS("nsclc_base.rds")
list2env(base, envir = environment())
library(survival); library(glmnet); library(timeROC); library(ggplot2); library(dplyr)

calc_risk <- function(data, genes, coefs) {
  genes <- intersect(genes, colnames(data))
  if (length(genes) == 0) return(rep(0, nrow(data))) # Return 0 if no genes
  cf <- as.numeric(coefs[genes])
  as.numeric(as.matrix(data[, genes, drop = FALSE]) %*% cf)
}
safe_cindex <- function(lp, time, event) {
  if (all(is.na(lp)) || length(unique(lp)) == 1) return(0.5)
  unname(summary(survival::coxph(survival::Surv(time, event) ~ lp))$concordance[1])
}

available_genes <- intersect(candidate_genes, colnames(train_data))
x_tr <- as.matrix(train_data[, available_genes])
y_tr <- Surv(train_data$os_time, train_data$os_event)
x_te <- as.matrix(test_data[,  available_genes])
y_te <- Surv(test_data$os_time,  test_data$os_event)

## =============================================================================
## BLOCK 1 — TABLE S3: UNIVARIATE COX SCREENING FOR ALL 25 GENES
## =============================================================================
cat("\n=== UNIVARIATE COX SCREENING (Table S3) ===\n")
univ_res <- data.frame(Gene = character(), HR = numeric(), Lower = numeric(),
                       Upper = numeric(), P_value = numeric(), FDR = numeric(),
                       stringsAsFactors = FALSE)

for (g in available_genes) {
  f <- as.formula(sprintf("Surv(os_time, os_event) ~ `%s`", g))
  fit <- tryCatch(coxph(f, data = train_data), error = function(e) NULL)
  if (!is.null(fit)) {
    s <- summary(fit)
    univ_res <- rbind(univ_res, data.frame(
      Gene = g, HR = s$conf.int[1, 1], Lower = s$conf.int[1, 3],
      Upper = s$conf.int[1, 4], P_value = s$coefficients[1, 5]
    ))
  }
}
univ_res$FDR <- p.adjust(univ_res$P_value, method = "BH")
univ_res <- univ_res[order(univ_res$P_value), ]
write.csv(univ_res, "TableS3_Univariate_Cox.csv", row.names = FALSE)
print(univ_res[1:5, ]) # Print top 5

## =============================================================================
## BLOCK 2 — LASSO (EXPRESSION) & FIGURE S3
## =============================================================================
cat("\n=== LASSO EXPRESSION ===\n")
set.seed(456)
cv_expr <- cv.glmnet(x_tr, y_tr, family = "cox", alpha = 1,
                     nfolds = 10, type.measure = "C")

# Check number of genes at lambda.min
cf_min <- as.matrix(coef(cv_expr, s = "lambda.min"))
sel_min <- rownames(cf_min)[cf_min[,1] != 0]
cat("Genes selected at lambda.min:", length(sel_min), "\n")
if (length(sel_min) > 0) cat(paste(sel_min, collapse = ", "), "\n")

# Figure S3: LASSO Paths and CV Curve
pdf("FigureS3_LASSO_Paths.pdf", width = 12, height = 6)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))
plot(cv_expr$glmnet.fit, xvar = "lambda", label = FALSE)
abline(v = log(cv_expr$lambda.min), lty = 2, col = "red", lwd = 2)
title(main = "A: LASSO Coefficient Paths", cex.main = 1.4)
plot(cv_expr)
title(main = "B: CV Partial Likelihood Deviance", cex.main = 1.4)
dev.off()
cat("Saved: FigureS3_LASSO_Paths.pdf\n")

## =============================================================================
## BLOCK 3 — CLINICAL VS INTEGRATED MODELS
## =============================================================================
cat("\n=== CLINICAL VS INTEGRATED (CLINICAL + GENES) ===\n")

# 1. Clinical-only model
cox_clin <- coxph(Surv(os_time, os_event) ~ age + gender + stage_group + histology,
                  data = train_data)
lp_clin_te <- predict(cox_clin, newdata = test_data, type = "lp")
c_clin_te  <- safe_cindex(lp_clin_te, test_data$os_time, test_data$os_event)

# 2. Integrated Penalized Cox (clinical forced in, genes penalized)
X_clin_tr <- model.matrix(~ age + gender + stage_group + histology - 1, data = train_data)
X_clin_te <- model.matrix(~ age + gender + stage_group + histology - 1, data = test_data)
X_full_tr <- cbind(X_clin_tr, x_tr)
X_full_te <- cbind(X_clin_te, x_te)

pf <- c(rep(0, ncol(X_clin_tr)), rep(1, ncol(x_tr)))
set.seed(456)
cv_int <- cv.glmnet(X_full_tr, y_tr, family = "cox", alpha = 1,
                    nfolds = 10, penalty.factor = pf, type.measure = "C")
lp_int_te <- as.numeric(predict(cv_int, newx = X_full_te, s = "lambda.min"))
c_int_te  <- safe_cindex(lp_int_te, test_data$os_time, test_data$os_event)

# 3. Pure Gene LASSO (if any genes were selected)
if (length(sel_min) > 0) {
  lp_expr_te <- as.numeric(predict(cv_expr, newx = x_te, s = "lambda.min"))
  c_expr_te  <- safe_cindex(lp_expr_te, test_data$os_time, test_data$os_event)
} else {
  lp_expr_te <- rep(0, nrow(test_data))
  c_expr_te  <- 0.500
}

cat("Clinical-only C-index        :", round(c_clin_te, 3), "\n")
cat("Integrated (Clin + Genes)    :", round(c_int_te, 3), "\n")
cat("Expression-only C-index      :", round(c_expr_te, 3), "\n")
cat("Delta C (Integrated - Clin)  :", round(c_int_te - c_clin_te, 4), "\n")

## =============================================================================
## BLOCK 4 — FIGURE 5: KM CLINICAL MODEL (REPLACES GENE KM)
## =============================================================================
cat("\n=== FIGURE 5 (CLINICAL KM) ===\n")

# Use clinical risk score for grouping
train_data$risk_clin <- predict(cox_clin, newdata = train_data, type = "lp")
test_data$risk_clin  <- predict(cox_clin, newdata = test_data,  type = "lp")

risk_cutoff <- median(train_data$risk_clin, na.rm = TRUE)
train_data$risk_group <- factor(ifelse(train_data$risk_clin > risk_cutoff, "High","Low"), levels = c("Low","High"))
test_data$risk_group  <- factor(ifelse(test_data$risk_clin > risk_cutoff, "High","Low"), levels = c("Low","High"))

pdf("Figure5_KM_Clinical.pdf", width = 14, height = 7)
par(mfrow = c(1,2), mar = c(5.5,5.5,7,3), cex.lab = 1.6, cex.axis = 1.3)
for (d in list(list(data=train_data, title="A: Training Set"), list(data=test_data, title="B: Validation Set"))) {
  fit <- survfit(Surv(os_time, os_event) ~ risk_group, data = d$data)
  lr  <- survdiff(Surv(os_time, os_event) ~ risk_group, data = d$data)$pvalue
  c_idx <- safe_cindex(d$data$risk_clin, d$data$os_time, d$data$os_event)
  plot(fit, col = c("#2E9FDF","#E7B800"), lwd = 4, ylim = c(0,1.05),
       xlab = "Time (days)", ylab = "Survival Probability", main = "")
  grid(col = "gray90", lty = 2)
  mtext(d$title, side = 3, line = 5, cex = 1.8, font = 2, adj = 0)
  mtext("Clinical Model (Age + Gender + Stage + Histology)", side = 3, line = 3.5, cex = 1.1, font = 3, adj = 0)
  mtext(sprintf("N = %d | C-index = %.3f | Log-rank p = %.4f", nrow(d$data), c_idx, lr),
        side = 3, line = 2, cex = 1.2, adj = 0)
  legend("topright", c("Low Risk","High Risk"), col = c("#2E9FDF","#E7B800"), lwd = 4, bty = "n", cex = 1.4)
}
dev.off()
cat("Saved: Figure5_KM_Clinical.pdf\n")

## =============================================================================
## BLOCK 5 — FIGURE S2: TIME-ROC CLINICAL VS INTEGRATED
## =============================================================================
cat("\n=== FIGURE S2 (TIME-ROC) ===\n")

# We need valid positive markers. For timeROC, higher marker = higher risk.
# Both lp_clin_te and lp_int_te are hazard scores.
roc_clin <- timeROC(T = test_data$os_time, delta = test_data$os_event,
                    marker = lp_clin_te, cause = 1, times = c(365, 1095, 1825), iid = FALSE)
roc_int  <- timeROC(T = test_data$os_time, delta = test_data$os_event,
                    marker = lp_int_te, cause = 1, times = c(365, 1095, 1825), iid = FALSE)

pdf("FigureS2_TimeROC_Clinical_vs_Int.pdf", width = 10, height = 8)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))
# Clinical ROC
plot(roc_clin, time = 1095, col = "#2E86AB", lwd = 2, title = FALSE)
abline(0, 1, lty = 3, col = "gray50")
title(main = "A: Clinical Model (3-year ROC)")
legend("bottomright", legend = paste0("AUC = ", round(roc_clin$AUC[2], 3)),
       col = "#2E86AB", lwd = 2, bty = "n", cex = 1.2)

# Integrated ROC
plot(roc_int, time = 1095, col = "#D55E00", lwd = 2, title = FALSE)
abline(0, 1, lty = 3, col = "gray50")
title(main = "B: Integrated (Clinical + Genes, 3-year ROC)")
legend("bottomright", legend = paste0("AUC = ", round(roc_int$AUC[2], 3)),
       col = "#D55E00", lwd = 2, bty = "n", cex = 1.2)
dev.off()

write.csv(data.frame(
  Model = c(rep("Clinical",3), rep("Integrated",3)),
  Time = rep(c("1y","3y","5y"), 2),
  AUC = c(roc_clin$AUC, roc_int$AUC)
), "Table_TimeROC_AUC.csv", row.names = FALSE)
cat("Saved: FigureS2_TimeROC_Clinical_vs_Int.pdf\n")

## =============================================================================
## BLOCK 6 — FIGURE S6: STRATIFIED KM FOR CLINICAL MODEL (FIXED)
## =============================================================================
cat("\n=== FIGURE S6 (STRATIFIED CLINICAL KM) ===\n")

cohort$stage_group <- factor(cohort$stage_group, levels = c("Early","Advanced","Unknown"))
cohort$histology   <- factor(cohort$histology, levels = c("LUAD","LUSC"))
cohort$risk_clin   <- predict(cox_clin, newdata = cohort, type = "lp")

# FIX: Recalculate cutoff within each subgroup to avoid "only 1 group" error
plot_km_strat <- function(d, title, panel) {
  if (nrow(d) < 10) return() # Skip if group is too small
  
  # Use within-subgroup median to guarantee two groups for plotting
  sub_cutoff <- median(d$risk_clin, na.rm = TRUE)
  d$risk_group_sub <- factor(ifelse(d$risk_clin > sub_cutoff, "High","Low"), levels = c("Low","High"))
  
  fit <- survfit(Surv(os_time, os_event) ~ risk_group_sub, data = d)
  
  # Safe log-rank test
  if (length(unique(d$risk_group_sub)) > 1) {
    lr  <- survdiff(Surv(os_time, os_event) ~ risk_group_sub, data = d)$pvalue
  } else {
    lr <- NA
  }
  
  cidx <- safe_cindex(d$risk_clin, d$os_time, d$os_event)
  
  plot(fit, col = c("#2E86AB","#A23B72"), lwd = 3.5, ylim = c(0,1.05),
       xlab = "Time (days)", ylab = "Survival Probability", main = "")
  grid(col = "gray90", lty = 2)
  mtext(paste0(panel, ": ", title), side = 3, line = 5, cex = 1.5, font = 2, adj = 0)
  
  pval_txt <- ifelse(is.na(lr), "NA", sprintf("%.4f", lr))
  mtext(sprintf("N = %d | Events = %d | C-index = %.3f | Log-rank p = %s",
                nrow(d), sum(d$os_event), cidx, pval_txt), side = 3, line = 3, cex = 1.1, adj = 0)
  legend("topright", c("Low","High"), col = c("#2E86AB","#A23B72"), lwd = 3.5, bty = "n", cex = 1.2)
}

pdf("FigureS6_KM_Clinical_Stratified.pdf", width = 20, height = 14)
par(mfrow = c(2, 3), mar = c(5.5, 5.5, 7, 3), cex.lab = 1.6, cex.axis = 1.4)

plot_km_strat(cohort, "All Patients", "A")
plot_km_strat(cohort[cohort$histology == "LUAD", ], "LUAD", "B")
plot_km_strat(cohort[cohort$histology == "LUSC", ], "LUSC", "C")
plot_km_strat(cohort[cohort$stage_group == "Early", ], "Early Stage", "D")
plot_km_strat(cohort[cohort$stage_group == "Advanced", ], "Advanced Stage", "E")
plot_km_strat(cohort[cohort$stage_group == "Unknown", ], "Unknown Stage", "F")

dev.off()
cat("Saved: FigureS6_KM_Clinical_Stratified.pdf\n")

## =============================================================================
## BLOCK 7 — MULTIVARIABLE COX & FOREST PLOT
## =============================================================================
cat("\n=== MULTIVARIABLE COX ===\n")
# Fit clinical + whatever genes survived (if any)
if (length(sel_min) > 0) {
  f_multi <- as.formula(paste("Surv(os_time, os_event) ~",
                              paste(c("age","gender","stage_group","histology",
                                      paste0("`", sel_min, "`")), collapse=" + ")))
} else {
  f_multi <- Surv(os_time, os_event) ~ age + gender + stage_group + histology
}

multi_expr <- coxph(f_multi, data = cohort)
print(summary(multi_expr))

# Extract coefficients for forest plot
coef_df <- data.frame(
  Variable = gsub("`", "", rownames(summary(multi_expr)$coefficients)),
  HR = exp(coef(multi_expr)),
  Lower = exp(confint(multi_expr))[,1],
  Upper = exp(confint(multi_expr))[,2],
  P = summary(multi_expr)$coefficients[, "Pr(>|z|)"]
)
coef_df$Sig <- ifelse(coef_df$P < 0.001, "***", ifelse(coef_df$P < 0.01, "**", ifelse(coef_df$P < 0.05, "*", "")))
coef_df$Label <- paste0(coef_df$Variable, " ", coef_df$Sig)

p_forest <- ggplot(coef_df, aes(x = HR, y = reorder(Label, HR))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_point(size = 5, color = "#2E86AB") +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.25, linewidth = 1.2, color = "#2E86AB") +
  scale_x_log10() +
  labs(x = "Hazard Ratio (95% CI, log scale)", y = "",
       title = "Multivariable Cox Regression",
       subtitle = "Clinical variables + LASSO selected genes (if any)") +
  theme_bw(base_size = 14) +
  geom_text(aes(label = sprintf("%.2f", HR)), vjust = -1, size = 4, fontface = "bold")
ggsave("Figure_Forest_Plot_Expression.pdf", p_forest, width = 10, height = 6)

write.csv(coef_df, "Table_Multivariable_Cox_Expression.csv", row.names = FALSE)

# Save RDS for Script 05
saveRDS(list(sel_min = sel_min, c_expr_te = c_expr_te, c_int_te = c_int_te,
             c_clin_te = c_clin_te, cv_integrated = cv_int, cox_expr_tr = cox_clin),
        "nsclc_script03_signature.rds")

cat("\nScript 03 complete.\n")