# ==============================================================================
# Script 04b: Methylation and Expression HR Reporting per SD, and per IQR
#
# Author:   Getachew Derbew Belay
# Date:     October 2026
# Platform: R 4.5.2
#
# Purpose:   Compute risk scores of the cohorts, fit univariate and adjusted Cox models, and report
#           HR per raw unit, per SD, and per IQR. Produce  Table S6 and data for revised Figure S8.
# ==============================================================================
# Outputs:
#   TableS6_Corrected_HR.csv               (whole cohort + test set)
#   TableS6_Corrected_HR_formatted.csv     (display version)
#   nsclc_script04b_hr.rds
#
# Depends:  nsclc_base.rds, nsclc_script03_signature.rds,
#           nsclc_script04_methylation.rds
#
# Seeds:    456 (for consistency with Script 03 and 04)
# ==============================================================================

setwd("~/TCGA_Local_Storage/TCGA NSCLC Project")
base <- readRDS("nsclc_base.rds");               list2env(base, envir = environment())
sig  <- readRDS("nsclc_script03_signature.rds"); list2env(sig,  envir = environment())
library(survival); library(glmnet); library(dplyr)

## ---- Shared utilities --------------------------------------------------------
calc_risk <- function(data, genes, coefs) {
  genes <- intersect(genes, colnames(data))
  if (length(genes) == 0) return(rep(NA_real_, nrow(data)))
  cf <- as.numeric(coefs[genes])
  as.numeric(as.matrix(data[, genes, drop = FALSE]) %*% cf)
}

rescale_hr <- function(beta, se, scale_factor, scale_label = "per unit") {
  hr    <- exp(beta * scale_factor)
  ci_lo <- exp((beta - 1.96 * se) * scale_factor)
  ci_hi <- exp((beta + 1.96 * se) * scale_factor)
  p     <- 2 * pnorm(-abs(beta / se))
  data.frame(Scale = scale_label, Factor = scale_factor,
             HR = hr, CI_lower = ci_lo, CI_upper = ci_hi, P_value = p,
             stringsAsFactors = FALSE)
}

## =============================================================================
## PART 0 — Recompute expression LASSO (fix for coef_v not found)
## =============================================================================
available_genes <- intersect(candidate_genes, colnames(train_data))
x_tr <- as.matrix(train_data[, available_genes])
y_tr <- Surv(train_data$os_time, train_data$os_event)

set.seed(456)
cv_expr <- cv.glmnet(x_tr, y_tr, family = "cox", alpha = 1,
                     nfolds = 10, type.measure = "C")
cf_min  <- as.matrix(coef(cv_expr, s = "lambda.min"))
sel_min <- rownames(cf_min)[cf_min[,1] != 0]
coef_v  <- as.numeric(cf_min[sel_min, 1]); names(coef_v) <- sel_min

cat("Expression genes retained at lambda.min:", length(sel_min), "\n")
if (length(sel_min) > 0) cat(paste(sel_min, collapse = ", "), "\n\n")

## =============================================================================
## PART A — METHYLATION RISK SCORE (correct cohort: meth_cohort, n = 701)
## =============================================================================
cat("\n=== METHYLATION RISK SCORE ===\n")

meth_cohort$risk_meth <- calc_risk(meth_cohort, sel_meth, coef_meth_v)
meth_cohort$stage_group <- factor(meth_cohort$stage_group,
                                  levels = c("Early","Advanced","Unknown"))
meth_cohort$histology   <- factor(meth_cohort$histology, levels = c("LUAD","LUSC"))
meth_cohort$gender      <- factor(meth_cohort$gender)

# Also compute on test_meth (out-of-sample)
test_meth$risk_meth <- calc_risk(test_meth, sel_meth, coef_meth_v)
test_meth$stage_group <- factor(test_meth$stage_group,
                                levels = c("Early","Advanced","Unknown"))
test_meth$histology   <- factor(test_meth$histology, levels = c("LUAD","LUSC"))
test_meth$gender      <- factor(test_meth$gender)

r <- meth_cohort$risk_meth[is.finite(meth_cohort$risk_meth)]
sd_meth <- sd(r); iqr_meth <- IQR(r)
cat(sprintf("Score: mean=%.4f | SD=%.4f | IQR=%.4f | range=[%.4f, %.4f]\n",
            mean(r), sd_meth, iqr_meth, min(r), max(r)))
cat(sprintf("Whole cohort: n=%d | events=%d\n",
            nrow(meth_cohort), sum(meth_cohort$os_event)))
cat(sprintf("Test set   : n=%d | events=%d\n",
            nrow(test_meth), sum(test_meth$os_event)))

# ---- A1. Univariate (whole cohort) ----
cox_meth_uni <- coxph(Surv(os_time, os_event) ~ risk_meth, data = meth_cohort)
b_u <- coef(cox_meth_uni)["risk_meth"]; se_u <- sqrt(vcov(cox_meth_uni)["risk_meth","risk_meth"])
hr_meth_uni <- rbind(
  rescale_hr(b_u, se_u, 1,        "per raw unit"),
  rescale_hr(b_u, se_u, sd_meth,  "per SD"),
  rescale_hr(b_u, se_u, iqr_meth, "per IQR")
)
hr_meth_uni$Model <- "Univariate (whole cohort)"; hr_meth_uni$Signature <- "Methylation"

# ---- A2. Adjusted (whole cohort, IN-SAMPLE) ----
cox_meth_adj <- coxph(
  Surv(os_time, os_event) ~ risk_meth + age + gender + stage_group + histology,
  data = meth_cohort)
b_a <- coef(cox_meth_adj)["risk_meth"]; se_a <- sqrt(vcov(cox_meth_adj)["risk_meth","risk_meth"])
hr_meth_adj <- rbind(
  rescale_hr(b_a, se_a, 1,        "per raw unit"),
  rescale_hr(b_a, se_a, sd_meth,  "per SD"),
  rescale_hr(b_a, se_a, iqr_meth, "per IQR")
)
hr_meth_adj$Model <- "Adjusted (whole cohort, in-sample)"; hr_meth_adj$Signature <- "Methylation"

# ---- A3. Adjusted (test set only, OUT-OF-SAMPLE) ----
cox_meth_adj_test <- coxph(
  Surv(os_time, os_event) ~ risk_meth + age + gender + stage_group + histology,
  data = test_meth)
b_t <- coef(cox_meth_adj_test)["risk_meth"]; se_t <- sqrt(vcov(cox_meth_adj_test)["risk_meth","risk_meth"])
hr_meth_test <- rbind(
  rescale_hr(b_t, se_t, 1,        "per raw unit"),
  rescale_hr(b_t, se_t, sd_meth,  "per SD"),
  rescale_hr(b_t, se_t, iqr_meth, "per IQR")
)
hr_meth_test$Model <- "Adjusted (test set only, out-of-sample)"; hr_meth_test$Signature <- "Methylation"

print(rbind(hr_meth_uni, hr_meth_adj, hr_meth_test))

## =============================================================================
## PART B — EXPRESSION RISK SCORE (cohort, n = 867)
## =============================================================================
cat("\n=== EXPRESSION RISK SCORE ===\n")

cohort$risk_expr <- calc_risk(cohort, sel_min, coef_v)
cohort$stage_group <- factor(cohort$stage_group,
                             levels = c("Early","Advanced","Unknown"))
cohort$histology   <- factor(cohort$histology, levels = c("LUAD","LUSC"))
cohort$gender      <- factor(cohort$gender)

test_data$risk_expr <- calc_risk(test_data, sel_min, coef_v)
test_data$stage_group <- factor(test_data$stage_group,
                                levels = c("Early","Advanced","Unknown"))
test_data$histology   <- factor(test_data$histology, levels = c("LUAD","LUSC"))
test_data$gender      <- factor(test_data$gender)

r_e <- cohort$risk_expr[is.finite(cohort$risk_expr)]
sd_expr <- sd(r_e); iqr_expr <- IQR(r_e)
cat(sprintf("Score: mean=%.4f | SD=%.4f | IQR=%.4f | range=[%.4f, %.4f]\n",
            mean(r_e), sd_expr, iqr_expr, min(r_e), max(r_e)))

cox_expr_uni <- coxph(Surv(os_time, os_event) ~ risk_expr, data = cohort)
b_e_u <- coef(cox_expr_uni)["risk_expr"]; se_e_u <- sqrt(vcov(cox_expr_uni)["risk_expr","risk_expr"])
hr_expr_uni <- rbind(
  rescale_hr(b_e_u, se_e_u, 1,        "per raw unit"),
  rescale_hr(b_e_u, se_e_u, sd_expr,  "per SD"),
  rescale_hr(b_e_u, se_e_u, iqr_expr, "per IQR")
)
hr_expr_uni$Model <- "Univariate (whole cohort)"; hr_expr_uni$Signature <- "Expression"

cox_expr_adj <- coxph(
  Surv(os_time, os_event) ~ risk_expr + age + gender + stage_group + histology,
  data = cohort)
b_e_a <- coef(cox_expr_adj)["risk_expr"]; se_e_a <- sqrt(vcov(cox_expr_adj)["risk_expr","risk_expr"])
hr_expr_adj <- rbind(
  rescale_hr(b_e_a, se_e_a, 1,        "per raw unit"),
  rescale_hr(b_e_a, se_e_a, sd_expr,  "per SD"),
  rescale_hr(b_e_a, se_e_a, iqr_expr, "per IQR")
)
hr_expr_adj$Model <- "Adjusted (whole cohort, in-sample)"; hr_expr_adj$Signature <- "Expression"

cox_expr_adj_test <- coxph(
  Surv(os_time, os_event) ~ risk_expr + age + gender + stage_group + histology,
  data = test_data)
b_e_t <- coef(cox_expr_adj_test)["risk_expr"]; se_e_t <- sqrt(vcov(cox_expr_adj_test)["risk_expr","risk_expr"])
hr_expr_test <- rbind(
  rescale_hr(b_e_t, se_e_t, 1,        "per raw unit"),
  rescale_hr(b_e_t, se_e_t, sd_expr,  "per SD"),
  rescale_hr(b_e_t, se_e_t, iqr_expr, "per IQR")
)
hr_expr_test$Model <- "Adjusted (test set only, out-of-sample)"; hr_expr_test$Signature <- "Expression"

print(rbind(hr_expr_uni, hr_expr_adj, hr_expr_test))

## =============================================================================
## PART C — Assemble and save Table S6
## =============================================================================
tableS6 <- rbind(hr_meth_uni, hr_meth_adj, hr_meth_test,
                 hr_expr_uni, hr_expr_adj, hr_expr_test)
tableS6 <- tableS6[, c("Signature","Model","Scale","Factor",
                       "HR","CI_lower","CI_upper","P_value")]

tableS6_display <- tableS6
tableS6_display$HR_95CI <- sprintf("%.3f (%.3f-%.3f)",
                                   tableS6$HR, tableS6$CI_lower, tableS6$CI_upper)
tableS6_display$P_value <- signif(tableS6$P_value, 3)
tableS6_display$Factor  <- signif(tableS6$Factor, 3)

write.csv(tableS6,          "TableS6_Corrected_HR.csv",            row.names = FALSE)
write.csv(tableS6_display,  "TableS6_Corrected_HR_formatted.csv",  row.names = FALSE)

cat("\n=== CORRECTED TABLE S6 ===\n")
print(tableS6_display[, c("Signature","Model","Scale","HR_95CI","P_value")],
      row.names = FALSE)

saveRDS(list(meth = list(sd = sd_meth, iqr = iqr_meth,
                         cox_uni = cox_meth_uni, cox_adj = cox_meth_adj,
                         cox_adj_test = cox_meth_adj_test,
                         table = rbind(hr_meth_uni, hr_meth_adj, hr_meth_test)),
             expr = list(sd = sd_expr, iqr = iqr_expr,
                         cox_uni = cox_expr_uni, cox_adj = cox_expr_adj,
                         cox_adj_test = cox_expr_adj_test,
                         table = rbind(hr_expr_uni, hr_expr_adj, hr_expr_test)),
             tableS6 = tableS6),
        "nsclc_script04b_hr.rds")

cat("\nScript 04b complete: TableS6_Corrected_HR.csv and .rds saved.\n")