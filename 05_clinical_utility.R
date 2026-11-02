# ==============================================================================
# Script 05: Clinical Utility (Nomogram, Calibration, DCA, Figure 6)
#
# Author:   Getachew Derbew Belay
# Date:     October 2026
# Platform: R 4.5.2
#
# Purpose:  Fit multivariable Cox models; build nomogram using basehaz();
#           bootstrap calibration; decision curve analysis using dcurves.
#           Assemble composite Figure 6 (A: nomogram, B: calibration, C: DCA).
#
# Save as:  05_clinical_utility.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/
#
# Output:   Table_LRT.csv
#           Figure6_PanelA_Nomogram.pdf
#           Figure6_PanelB_Calibration_rms.pdf
#           Figure6_PanelC_DCA.pdf
#           Figure6_Clinical_Utility.pdf
#           nsclc_script05_clinical.rds
#
# Depends:  nsclc_base.rds, nsclc_script03_signature.rds
# ==============================================================================

setwd("~/TCGA_Local_Storage/TCGA NSCLC Project")
base <- readRDS("nsclc_base.rds"); list2env(base, envir = environment())
sig  <- readRDS("nsclc_script03_signature.rds"); list2env(sig, envir = environment())
library(survival); library(rms); library(dcurves); library(ggplot2)
library(glmnet); library(dplyr); library(magick)

# 1. Re-fit Expression LASSO to get the exact coefficient vector (coef_v)
# (Ensures this script runs independently even if coef_v wasn't saved)
available_genes <- intersect(candidate_genes, colnames(train_data))
x_tr <- as.matrix(train_data[, available_genes])
y_tr <- Surv(train_data$os_time, train_data$os_event)
set.seed(456)
cv_expr <- cv.glmnet(x_tr, y_tr, family = "cox", alpha = 1, nfolds = 10, type.measure = "C")
cf_min <- as.matrix(coef(cv_expr, s = "lambda.min"))
sel_min <- rownames(cf_min)[cf_min[,1] != 0]
coef_v <- as.numeric(cf_min[sel_min, 1]); names(coef_v) <- sel_min

calc_risk <- function(data, genes, coefs) {
  genes <- intersect(genes, colnames(data))
  if (length(genes) == 0) return(rep(0, nrow(data)))
  cf <- as.numeric(coefs[genes])
  as.numeric(as.matrix(data[, genes, drop = FALSE]) %*% cf)
}

# 2. Setup Cohort and Risk Scores
cohort$stage_group <- factor(cohort$stage_group, levels = c("Early","Advanced","Unknown"))
cohort$histology   <- factor(cohort$histology,   levels = c("LUAD","LUSC"))
cohort$gender      <- factor(cohort$gender)
cohort$risk_score  <- calc_risk(cohort, sel_min, coef_v)

## =============================================================================
## BLOCK 1 — MULTIVARIABLE COX + LRT
## =============================================================================
cat("\n=== LIKELIHOOD-RATIO TEST ===\n")
cox_clin <- coxph(Surv(os_time, os_event) ~ age + gender + stage_group + histology,
                  data = cohort)
cox_full <- coxph(Surv(os_time, os_event) ~ risk_score + age + gender + stage_group + histology,
                  data = cohort)
lrt <- anova(cox_clin, cox_full, test = "LRT")
print(lrt)
write.csv(data.frame(model = c("Clinical only","Clinical + risk"),
                     logLik = c(as.numeric(logLik(cox_clin)), as.numeric(logLik(cox_full))),
                     LRT_df = c(NA, lrt[2, "Df"]),
                     LRT_p  = c(NA, lrt[2, "Pr(>|Chi|)"])),
          "Table_LRT.csv", row.names = FALSE)

## =============================================================================
## BLOCK 2 — NOMOGRAM (CORRECT BASELINE)
## =============================================================================
cat("\n=== NOMOGRAM ===\n")
dd <- datadist(cohort); options(datadist = "dd")
cph_model <- cph(Surv(os_time, os_event) ~ risk_score + age + gender + stage_group + histology,
                 data = cohort, surv = TRUE, x = TRUE, y = TRUE)

# True X=0 baseline survival (centered = FALSE)
bh <- survival::basehaz(cph_model, centered = FALSE)
S0_fun <- function(t) exp(-approx(bh$time, bh$hazard, xout = t, rule = 2)$y)
S0_1y <- S0_fun(1 * 365.25)
S0_3y <- S0_fun(3 * 365.25)
S0_5y <- S0_fun(5 * 365.25)
surv_fun <- function(lp, t, S0) S0 ^ exp(lp)

nom <- nomogram(cph_model,
                fun = list(function(x) surv_fun(x, 1, S0_1y),
                           function(x) surv_fun(x, 3, S0_3y),
                           function(x) surv_fun(x, 5, S0_5y)),
                funlabel = c("1-Year Survival","3-Year Survival","5-Year Survival"),
                lp = FALSE)

pdf("Figure6_PanelA_Nomogram.pdf", width = 10, height = 8)
plot(nom, xfrac = 0.3, cex.axis = 0.85, cex.var = 0.9)
title("Nomogram: 1-, 3-, 5-Year Overall Survival")
dev.off()
cat("Saved: Figure6_PanelA_Nomogram.pdf\n")

## =============================================================================
## BLOCK 3 — CALIBRATION
## =============================================================================
cat("\n=== CALIBRATION ===\n")
test_sub <- cohort[cohort$patient_id %in% test_data$patient_id, ]
cal_3y <- calibrate(cph_model, method = "boot", B = 500, u = 3 * 365.25, cmethod = "hare", data = test_sub)

pdf("Figure6_PanelB_Calibration_rms.pdf", width = 8, height = 8)
plot(cal_3y, xlab = "Predicted 3-year survival", ylab = "Observed 3-year survival",
     main = "Calibration: 3-Year Survival", subtitles = FALSE)
dev.off()
cat("Saved: Figure6_PanelB_Calibration_rms.pdf\n")

## =============================================================================
## BLOCK 4 — DECISION CURVE ANALYSIS (dcurves)
## =============================================================================
cat("\n=== DCA ===\n")
test_sub$prob_full <- 1 - surv_fun(predict(cph_model, newdata = test_sub, type = "lp"), 3, S0_3y)
test_sub$prob_clin <- 1 - surv_fun(predict(cox_clin, newdata = test_sub, type = "lp"), 3, S0_3y)

dca_obj <- dca(Surv(os_time, os_event) ~ prob_full + prob_clin,
               data = test_sub, time = 3 * 365.25, thresholds = seq(0, 0.5, by = 0.01))

pdf("Figure6_PanelC_DCA.pdf", width = 8, height = 6)
print(plot(dca_obj, smooth = TRUE, main = "Decision Curve Analysis (3-year OS)"))
dev.off()
cat("Saved: Figure6_PanelC_DCA.pdf\n")

## =============================================================================
## BLOCK 5 — COMPOSITE FIGURE 6
## =============================================================================
cat("\n=== MERGING FIGURE 6 ===\n")
img_nom <- magick::image_read_pdf("Figure6_PanelA_Nomogram.pdf", pages = 1, density = 200)[1]
img_cal <- magick::image_read_pdf("Figure6_PanelB_Calibration_rms.pdf", pages = 1, density = 200)[1]
img_dca <- magick::image_read_pdf("Figure6_PanelC_DCA.pdf", pages = 1, density = 200)[1]

target_w <- 1600
img_nom <- magick::image_resize(img_nom, paste0(target_w, "x"))
img_cal <- magick::image_resize(img_cal, paste0(target_w, "x"))
img_dca <- magick::image_resize(img_dca, paste0(target_w, "x"))

img_nom <- magick::image_annotate(img_nom, "A", size = 40, weight = 700, location = "+30+30")
img_cal <- magick::image_annotate(img_cal, "B", size = 40, weight = 700, location = "+30+30")
img_dca <- magick::image_annotate(img_dca, "C", size = 40, weight = 700, location = "+30+30")

bottom  <- magick::image_append(c(img_cal, img_dca), stack = FALSE)
bottom_w <- magick::image_info(bottom)$width
img_nom  <- magick::image_resize(img_nom, paste0(bottom_w, "x"))
combined <- magick::image_flatten(magick::image_append(c(img_nom, bottom), stack = TRUE))

magick::image_write(combined, "Figure6_Clinical_Utility.pdf", format = "pdf", density = 200)
cat("Saved: Figure6_Clinical_Utility.pdf\n")

options(datadist = NULL)

saveRDS(list(cox_full = cox_full, cox_clin = cox_clin, cph_model = cph_model, lrt = lrt), 
        "nsclc_script05_clinical.rds")

cat("\nScript 05 complete.\n")