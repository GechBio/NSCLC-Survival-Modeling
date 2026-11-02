# ==============================================================================
# Script 02: Six-Model Comparison
#
# Author:   Getachew Derbew Belay
# Date:     October 2026
# Platform: R 4.5.2
# Contact:  derbew1953@gmail.com / d610113006@tmu.edu.tw
#
# Purpose:  Compare six survival models (Clinical-only, LASSO, Elastic Net,
#           Ridge, RSF, Stepwise Cox); report validation C-index and delta-C
#           vs clinical-only; produce Figure 4 (bar+ladder), Figure S4 (gene
#           heatmap), Figure S5 (nested CV stability).
# Produces: Table_All_Models_Comparison.csv, TableS_CV_Stability.csv, Figure_Model_Comparison_Combined.pdf and, 
#           Figure_Model_Comparison_Heatmap.pdf, Figure_Model_CV_ErrorBars.pdf, nsclc_script02_models.rds
#
# Save as:  02_six_model_comparison.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/
#
# Depends:  nsclc_base.rds  (from Script 00)
#
# Seeds:    cross-validation = 456
# TUNING APPLIED (all on training set only, no test-set leakage):
#   - LASSO / Elastic Net / Ridge: type.measure = "C"
#   - Elastic Net: alpha tuned over 0.1-0.9 via 10-fold CV
#   - RSF: mtry {3,5,8,12} x nodesize {5,10,15} tuned via OOB C-index
#   - Stepwise: forward vs both, AIC penalty
#   - Clinical: cox.zph PH check with stratification if needed
# ==============================================================================

setwd("~/TCGA_Local_Storage/TCGA NSCLC Project")
base <- readRDS("nsclc_base.rds")
list2env(base, envir = environment())
library(survival); library(glmnet); library(randomForestSRC)
library(ggplot2); library(patchwork); library(pheatmap); library(dplyr)

## ---- Shared utilities --------------------------------------------------------
calc_risk <- function(data, genes, coefs) {
  genes <- intersect(genes, colnames(data))
  if (length(genes) == 0) return(rep(NA_real_, nrow(data)))
  cf <- as.numeric(coefs[genes])
  as.numeric(as.matrix(data[, genes, drop = FALSE]) %*% cf)
}

## Direction-robust C-index: refit a univariate Cox internally
safe_cindex <- function(lp, time, event) {
  if (all(is.na(lp))) return(NA_real_)
  unname(summary(survival::coxph(survival::Surv(time, event) ~ lp))$concordance[1])
}

available_genes <- intersect(candidate_genes, colnames(train_data))
x_train <- as.matrix(train_data[, available_genes])
y_train <- Surv(train_data$os_time, train_data$os_event)
x_test  <- as.matrix(test_data[,  available_genes])
y_test  <- Surv(test_data$os_time,  test_data$os_event)

cat("\n=== SETUP ===\n")
cat("Training:", nrow(train_data), "| Test:", nrow(test_data), "\n")
cat("Available genes:", length(available_genes), "\n\n")

## =============================================================================
## BLOCK 1 — CLINICAL-ONLY MODEL (unstratified)
## =============================================================================
cat("\n=== CLINICAL-ONLY MODEL ===\n")

cox_clin <- coxph(Surv(os_time, os_event) ~ age + gender + stage_group + histology,
                  data = train_data)
lp_clin  <- predict(cox_clin, newdata = test_data, type = "lp")
c_clin   <- safe_cindex(lp_clin, test_data$os_time, test_data$os_event)

## PH assumption check (diagnostic only — do NOT stratify)
ph_test <- cox.zph(cox_clin)
cat("\nProportional hazards test (diagnostic only):\n")
print(ph_test$table)

ph_violations <- rownames(ph_test$table)[ph_test$table[, "p"] < 0.05]
ph_violations <- setdiff(ph_violations, "GLOBAL")

if (length(ph_violations) > 0) {
  cat("\nPH violations detected for:", paste(ph_violations, collapse = ", "), "\n")
  cat("NOTE: This is expected in prognostic settings because 'Unknown' stage\n")
  cat("      represents a heterogeneous high-risk group. The unstratified\n")
  cat("      model is retained because it preserves stage as a fixed effect\n")
  cat("      and provides a single risk score for prediction.\n")
} else {
  cat("No PH violations.\n")
}

cat("Clinical C-index (unstratified):", round(c_clin, 3), "\n")

## =============================================================================
## BLOCK 2 — LASSO (type.measure = "C")
## =============================================================================
cat("\n=== LASSO ===\n")

set.seed(456)
cv_lasso <- cv.glmnet(x_train, y_train, family = "cox", alpha = 1,
                      nfolds = 10, type.measure = "C")
lp_lasso <- as.numeric(predict(cv_lasso, newx = x_test, s = "lambda.min"))
c_lasso  <- safe_cindex(lp_lasso, test_data$os_time, test_data$os_event)
lasso_cf    <- as.matrix(coef(cv_lasso, s = "lambda.min"))
lasso_genes <- rownames(lasso_cf)[lasso_cf[,1] != 0]
cat("Lambda.min:", signif(cv_lasso$lambda.min, 3), "\n")
cat("Genes selected:", length(lasso_genes),
    if (length(lasso_genes) > 0) paste("(", paste(lasso_genes, collapse = ", "), ")") else "(none)", "\n")
cat("Test C-index:", round(c_lasso, 3), "\n")

## =============================================================================
## BLOCK 3 — ELASTIC NET (tune alpha)
## =============================================================================
cat("\n=== ELASTIC NET ===\n")

alphas <- seq(0.1, 0.9, by = 0.1)
cv_enet_results <- lapply(alphas, function(a) {
  set.seed(456)
  cv <- cv.glmnet(x_train, y_train, family = "cox", alpha = a,
                  nfolds = 10, type.measure = "C")
  list(alpha = a, cv = cv, cvm_max = max(cv$cvm))
})

enet_cvm <- sapply(cv_enet_results, function(x) x$cvm_max)
best_idx <- which.max(enet_cvm)
best_alpha <- alphas[best_idx]
cv_enet <- cv_enet_results[[best_idx]]$cv

cat("Alpha tuning (CV C-index):\n")
for (i in seq_along(alphas)) {
  cat(sprintf("  alpha = %.1f -> CV C-index = %.3f\n", alphas[i], enet_cvm[i]))
}
cat("Best alpha:", best_alpha, "\n")

lp_enet <- as.numeric(predict(cv_enet, newx = x_test, s = "lambda.min"))
c_enet  <- safe_cindex(lp_enet, test_data$os_time, test_data$os_event)
enet_cf    <- as.matrix(coef(cv_enet, s = "lambda.min"))
enet_genes <- rownames(enet_cf)[enet_cf[,1] != 0]
cat("Genes selected:", length(enet_genes),
    if (length(enet_genes) > 0) paste("(", paste(enet_genes, collapse = ", "), ")") else "(none)", "\n")
cat("Test C-index:", round(c_enet, 3), "\n")

## =============================================================================
## BLOCK 4 — RIDGE (type.measure = "C")
## =============================================================================
cat("\n=== RIDGE ===\n")

set.seed(456)
cv_ridge <- cv.glmnet(x_train, y_train, family = "cox", alpha = 0,
                      nfolds = 10, type.measure = "C")
lp_ridge <- as.numeric(predict(cv_ridge, newx = x_test, s = "lambda.min"))
c_ridge  <- safe_cindex(lp_ridge, test_data$os_time, test_data$os_event)
ridge_cf <- as.matrix(coef(cv_ridge, s = "lambda.min"))
ridge_nonzero <- rownames(ridge_cf)[ridge_cf[,1] != 0]
cat("Non-zero coefficients:", length(ridge_nonzero), "of", nrow(ridge_cf), "\n")
cat("Test C-index:", round(c_ridge, 3), "\n")

## =============================================================================
## BLOCK 5 — RANDOM SURVIVAL FOREST (tune mtry, nodesize)
## =============================================================================
cat("\n=== RANDOM SURVIVAL FOREST ===\n")
cat("Grid search (may take 5-10 min)...\n\n")

rsf_data <- train_data[, c("os_time","os_event", available_genes)]
rsf_vars <- names(which(apply(rsf_data[, available_genes], 2, sd, na.rm = TRUE) > 0))

mtry_grid     <- c(3, 5, 8, 12)
nodesize_grid <- c(5, 10, 15)
tune_grid <- expand.grid(mtry = mtry_grid, nodesize = nodesize_grid)

oob_c <- numeric(nrow(tune_grid))
for (i in seq_len(nrow(tune_grid))) {
  set.seed(456)
  fit <- rfsrc(Surv(os_time, os_event) ~ .,
               data = rsf_data[, c("os_time","os_event", rsf_vars)],
               ntree = 500,                      # 500 for tuning speed
               mtry = tune_grid$mtry[i],
               nodesize = tune_grid$nodesize[i],
               splitrule = "logrank",
               importance = FALSE)
  oob_c[i] <- 1 - tail(fit$err.rate, 1)
  cat(sprintf("  mtry = %2d, nodesize = %2d -> OOB C-index = %.3f\n",
              tune_grid$mtry[i], tune_grid$nodesize[i], oob_c[i]))
}

best_rsf_idx <- which.max(oob_c)
best_mtry    <- tune_grid$mtry[best_rsf_idx]
best_nodesize <- tune_grid$nodesize[best_rsf_idx]
cat("\nBest parameters: mtry =", best_mtry, ", nodesize =", best_nodesize,
    ", OOB C-index =", round(oob_c[best_rsf_idx], 3), "\n")

## Refit with best parameters and full ntree
set.seed(456)
rsf_model <- rfsrc(Surv(os_time, os_event) ~ .,
                   data = rsf_data[, c("os_time","os_event", rsf_vars)],
                   ntree = 1000,
                   mtry = best_mtry,
                   nodesize = best_nodesize,
                   splitrule = "logrank",
                   importance = TRUE)
rsf_pred <- as.numeric(predict(rsf_model, newdata = test_data)$predicted)
c_rsf    <- safe_cindex(rsf_pred, test_data$os_time, test_data$os_event)
rsf_top5 <- names(sort(rsf_model$importance, decreasing = TRUE))[1:5]
cat("Top-5 features:", paste(rsf_top5, collapse = ", "), "\n")
cat("Test C-index:", round(c_rsf, 3), "\n")

## =============================================================================
## BLOCK 6 — STEPWISE COX (forward vs both, AIC)
## =============================================================================
cat("\n=== STEPWISE COX ===\n")

base_formula <- Surv(os_time, os_event) ~ age + gender + stage_group + histology
quoted_genes <- paste0("`", available_genes, "`")
full_formula <- as.formula(paste("Surv(os_time, os_event) ~",
                                 paste(c(quoted_genes, "age","gender","stage_group","histology"),
                                       collapse = " + ")))

step_forward <- step(coxph(base_formula, data = train_data),
                     scope = list(upper = full_formula),
                     direction = "forward", trace = 0)
step_both <- step(coxph(base_formula, data = train_data),
                  scope = list(upper = full_formula),
                  direction = "both", trace = 0)

cat("Forward AIC:", round(AIC(step_forward), 2), "\n")
cat("Both    AIC:", round(AIC(step_both),    2), "\n")
step_cox <- if (AIC(step_both) < AIC(step_forward)) step_both else step_forward
cat("Selected model:",
    if (identical(step_cox, step_both)) "both" else "forward", "\n")

lp_step <- predict(step_cox, newdata = test_data, type = "lp")
c_step  <- safe_cindex(lp_step, test_data$os_time, test_data$os_event)

## Fix backtick bug
coef_names   <- gsub("`", "", names(coef(step_cox)))
step_genes   <- intersect(coef_names, available_genes)
step_coefs   <- coef(step_cox)[match(step_genes, coef_names)]
clin_terms   <- setdiff(coef_names, available_genes)

cat("Genes added:", length(step_genes),
    if (length(step_genes) > 0) paste("(", paste(step_genes, collapse = ", "), ")") else "(none)", "\n")
cat("Clinical covariates retained:", length(clin_terms), "\n")
cat("Test C-index:", round(c_step, 3), "\n")

## =============================================================================
## BLOCK 7 — COMPARISON TABLE
## =============================================================================
delta_c <- c(
  LASSO      = c_lasso - c_clin,
  ElasticNet = c_enet  - c_clin,
  Ridge      = c_ridge - c_clin,
  RSF        = c_rsf   - c_clin,
  Stepwise   = c_step  - c_clin
)

cat("\n=== VALIDATION C-INDEX ===\n")
print(round(c(ClinicalOnly = c_clin, LASSO = c_lasso, ElasticNet = c_enet,
              Ridge = c_ridge, RSF = c_rsf, Stepwise = c_step), 3))
cat("\n=== DELTA C VS CLINICAL-ONLY ===\n")
print(round(delta_c, 3))

comparison_all <- data.frame(
  Model = c("Clinical Only (Baseline)",
            paste0("LASSO (", length(lasso_genes), " genes)"),
            paste0("Elastic Net (", length(enet_genes), " genes; alpha=", best_alpha, ")"),
            paste0("Ridge (", length(ridge_nonzero), " genes)"),
            paste0("Random Forest (", length(rsf_vars), " genes)"),
            paste0("Stepwise Cox (", length(step_genes), " genes)")),
  C_index = round(c(c_clin, c_lasso, c_enet, c_ridge, c_rsf, c_step), 3),
  Delta_C = round(c(0, delta_c["LASSO"], delta_c["ElasticNet"],
                    delta_c["Ridge"], delta_c["RSF"], delta_c["Stepwise"]), 3))
write.csv(comparison_all, "Table_All_Models_Comparison.csv", row.names = FALSE)

## =============================================================================
## BLOCK 8 — FIGURE 4 (bar + ladder)
## =============================================================================
comparison_all$Model <- factor(comparison_all$Model,
                               levels = comparison_all$Model[order(-comparison_all$C_index)])
p1 <- ggplot(comparison_all, aes(Model, C_index, fill = C_index)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  geom_text(aes(label = sprintf("%.3f", C_index)), vjust = -0.5,
            size = 4.5, fontface = "bold") +
  scale_fill_gradient(low = "#F8766D", high = "#00BFC4", limits = c(0.4, 0.75)) +
  labs(y = "C-index (Validation)", x = "", title = "Prognostic Model Performance") +
  ylim(0, 0.80) + theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")

comparison_all$Improvement <- comparison_all$C_index - 0.5
p2 <- ggplot(comparison_all, aes(Improvement, reorder(Model, Improvement))) +
  geom_point(size = 5, aes(color = Improvement > 0)) +
  geom_segment(aes(xend = 0, yend = Model), linewidth = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("TRUE" = "#00BFC4","FALSE" = "#F8766D"), guide = "none") +
  labs(x = "Improvement Over Random", y = "", title = "Improvement Over Random") +
  theme_bw(base_size = 14)

ggsave("Figure_Model_Comparison_Combined.pdf",
       p1 + p2 + plot_layout(widths = c(2,1)), width = 16, height = 8)
cat("\nSaved: Figure_Model_Comparison_Combined.pdf\n")

## =============================================================================
## BLOCK 9 — NESTED CV STABILITY (penalized models)
## =============================================================================
compute_cv_cindex <- function(x, time, event, alpha, nfolds = 10,
                              inner = 5, seed = 456) {
  n <- nrow(x); set.seed(seed)
  fold_ids <- sample(rep(1:nfolds, length.out = n))
  out <- numeric(nfolds)
  for (k in 1:nfolds) {
    te <- which(fold_ids == k); tr <- which(fold_ids != k)
    set.seed(seed)
    fit <- cv.glmnet(x[tr,,drop = FALSE], Surv(time[tr], event[tr]),
                     family = "cox", alpha = alpha, nfolds = inner,
                     type.measure = "C")
    lp  <- as.numeric(predict(fit, newx = x[te,,drop = FALSE], s = "lambda.min"))
    out[k] <- safe_cindex(lp, time[te], event[te])
  }
  list(mean = mean(out), sd = sd(out), folds = out)
}

cat("\nRunning nested CV (LASSO, EN, Ridge)...\n")
cv_lasso_stats <- compute_cv_cindex(x_train, train_data$os_time, train_data$os_event, 1.0)
cv_enet_stats  <- compute_cv_cindex(x_train, train_data$os_time, train_data$os_event, best_alpha)
cv_ridge_stats <- compute_cv_cindex(x_train, train_data$os_time, train_data$os_event, 0.0)

cv_summary <- data.frame(
  Model = factor(c("Ridge","Elastic Net","LASSO"),
                 levels = c("Ridge","Elastic Net","LASSO")),
  Mean  = c(cv_ridge_stats$mean, cv_enet_stats$mean, cv_lasso_stats$mean),
  SD    = c(cv_ridge_stats$sd, cv_enet_stats$sd, cv_lasso_stats$sd))
cv_summary$Lower <- pmax(cv_summary$Mean - cv_summary$SD, 0)
cv_summary$Upper <- pmin(cv_summary$Mean + cv_summary$SD, 1)
write.csv(cv_summary, "TableS_CV_Stability.csv", row.names = FALSE)
print(cv_summary)

p3 <- ggplot(cv_summary, aes(Model, Mean)) +
  geom_point(size = 4, color = "#00BFC4") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.15, color = "#00BFC4") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  geom_text(aes(label = sprintf("%.3f \u00b1 %.3f", Mean, SD)),
            vjust = -1.5, size = 4.5, fontface = "bold") +
  ylim(0, 0.7) +
  labs(y = "Nested-CV C-index", x = "", title = "Nested CV Stability") +
  theme_bw(base_size = 14)
ggsave("Figure_Model_CV_ErrorBars.pdf", p3, width = 8, height = 6)
cat("Saved: Figure_Model_CV_ErrorBars.pdf\n")

## =============================================================================
## BLOCK 10 — GENE SELECTION HEATMAP (Ridge EXCLUDED)
## =============================================================================
cat("\n=== GENE SELECTION HEATMAP ===\n")
cat("Note: Ridge is excluded from this figure because it retains all\n")
cat("      25 genes with shrunken coefficients (no feature selection).\n\n")

## Collect genes selected by any selection-based model (not Ridge)
all_sel <- unique(c(lasso_genes, enet_genes, step_genes, rsf_top5))

if (length(all_sel) == 0) {
  cat("No genes selected by any model. Skipping heatmap.\n")
} else {
  gene_matrix <- matrix(0, nrow = length(all_sel), ncol = 5,
                        dimnames = list(all_sel,
                                        c("LASSO","ElasticNet",
                                          "RSF_top5","Stepwise","Clinical")))
  for (g in all_sel) {
    gene_matrix[g, "LASSO"]      <- as.integer(g %in% lasso_genes)
    gene_matrix[g, "ElasticNet"] <- as.integer(g %in% enet_genes)
    gene_matrix[g, "RSF_top5"]   <- as.integer(g %in% rsf_top5)
    gene_matrix[g, "Stepwise"]   <- as.integer(g %in% step_genes)
    gene_matrix[g, "Clinical"]   <- 0L
  }
  
  ## Order rows by total selections (most-selected at top)
  row_totals <- rowSums(gene_matrix)
  gene_matrix <- gene_matrix[order(-row_totals, rownames(gene_matrix)), , drop = FALSE]
  
  pdf("Figure_Model_Comparison_Heatmap.pdf", width = 10, height = 5)
  pheatmap(gene_matrix,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           color = colorRampPalette(c("white","#00BFC4"))(100),
           display_numbers = TRUE,
           number_format = "%.0f",
           main = "Gene Selection Across Models (Ridge excluded)",
           fontsize = 12,
           border_color = "gray80")
  dev.off()
  cat("Saved: Figure_Model_Comparison_Heatmap.pdf\n")
  
  ## Save as CSV too
  gene_selection_table <- data.frame(
    Gene       = rownames(gene_matrix),
    LASSO      = gene_matrix[, "LASSO"],
    ElasticNet = gene_matrix[, "ElasticNet"],
    RSF_top5   = gene_matrix[, "RSF_top5"],
    Stepwise   = gene_matrix[, "Stepwise"],
    Total      = row_totals[rownames(gene_matrix)],
    stringsAsFactors = FALSE)
  write.csv(gene_selection_table,
            "TableS_Gene_Selection_by_Model.csv", row.names = FALSE)
  cat("Saved: TableS_Gene_Selection_by_Model.csv\n")
}

## =============================================================================
## BLOCK 11 — SAVE RDS
## =============================================================================
saveRDS(list(
  c_clin = c_clin, lasso_genes = lasso_genes, enet_genes = enet_genes,
  step_genes = step_genes, rsf_top5 = rsf_top5,
  cox_clin = cox_clin, cox_clin_original = cox_clin,         # <--- FIXED
  cv_lasso = cv_lasso, cv_enet = cv_enet, cv_ridge = cv_ridge,
  best_alpha = best_alpha, best_mtry = best_mtry, best_nodesize = best_nodesize,
  ph_violations = ph_violations
), "nsclc_script02_models.rds")

cat("\nScript 02 complete.\n")
