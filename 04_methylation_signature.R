# ==============================================================================
# Script 04: Methylation Signature (Gene-level and CpG-level)
#
# Author:   Getachew Derbew Belay
# Date:     October 2026
# Platform: R 4.5.2
#
# Purpose:  Build gene-level promoter methylation LASSO signature; build a
#           CpG-level signature via top-1 CpG per gene; evaluate validation
#           C-index against clinical baseline; produce stratified KM panels.
#
# Save as:  04_methylation_signature.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/
#
# Output:   TableS4_Methylation_Signature.csv
#           TableS5_CpG_Signature.csv
#           FigureS7_KM_Methylation_Stratified.pdf
#           nsclc_script04_methylation.rds
#
# Depends:  nsclc_base.rds  (from Script 00)
# ==============================================================================

setwd("~/TCGA_Local_Storage/TCGA NSCLC Project")
base <- readRDS("nsclc_base.rds")
list2env(base, envir = environment())
library(survival); library(glmnet); library(ggplot2); library(dplyr); library(caret)

## Shared utilities
calc_risk <- function(data, genes, coefs) {
  genes <- intersect(genes, colnames(data))
  if (length(genes) == 0) return(rep(0, nrow(data)))
  cf <- as.numeric(coefs[genes])
  as.numeric(as.matrix(data[, genes, drop = FALSE]) %*% cf)
}
safe_cindex <- function(lp, time, event) {
  if (all(is.na(lp)) || length(unique(lp)) == 1) return(0.5)
  unname(summary(survival::coxph(survival::Surv(time, event) ~ lp))$concordance[1])
}
# NEW: C-index with 95% CI from the concordance object
safe_cindex_ci <- function(lp, time, event) {
  if (all(is.na(lp)) || length(unique(lp)) == 1)
    return(c(cindex = 0.5, lower = NA, upper = NA))
  cc <- summary(survival::coxph(survival::Surv(time, event) ~ lp))$concordance
  c(cindex = unname(cc[1]),
    lower  = unname(cc[1] - 1.96 * cc[2]),
    upper  = unname(cc[1] + 1.96 * cc[2]))
}

meth_genes <- rownames(gene_meth)
x_tr <- as.matrix(train_meth[, meth_genes])
y_tr <- Surv(train_meth$os_time, train_meth$os_event)
x_te <- as.matrix(test_meth[,  meth_genes])
y_te <- Surv(test_meth$os_time,  test_meth$os_event)

## =============================================================================
## BLOCK 0 (NEW) — Verify gene-level split sizes
## =============================================================================
cat("\n=== GENE-LEVEL SPLIT VERIFICATION ===\n")
cat(sprintf("meth_cohort (total) : n = %d | events = %d\n",
            nrow(meth_cohort), sum(meth_cohort$os_event)))
cat(sprintf("train_meth          : n = %d | events = %d\n",
            nrow(train_meth), sum(train_meth$os_event)))
cat(sprintf("test_meth           : n = %d | events = %d\n",
            nrow(test_meth),  sum(test_meth$os_event)))
cat(sprintf("Split fraction      : %.4f / %.4f\n",
            nrow(train_meth) / nrow(meth_cohort),
            nrow(test_meth)  / nrow(meth_cohort)))
cat("Expected (manuscript): train = 494, test = 207\n")
if (nrow(train_meth) == 494 && nrow(test_meth) == 207) {
  cat(">>> Match: split is 494/207, as expected from current pipeline.\n")
} else {
  cat(">>> WARNING: split differs from expected 494/207. Update manuscript.\n")
}

## =============================================================================
## BLOCK 1 — GENE-LEVEL METHYLATION LASSO
## =============================================================================
cat("\n=== GENE-LEVEL METHYLATION LASSO ===\n")
set.seed(456)
cv_meth <- cv.glmnet(x_tr, y_tr, family = "cox", alpha = 1,
                     nfolds = 10, type.measure = "C")
cf_meth <- as.matrix(coef(cv_meth, s = "lambda.min"))
sel_meth <- rownames(cf_meth)[cf_meth[,1] != 0]
coef_meth_v <- as.numeric(cf_meth[sel_meth, 1]); names(coef_meth_v) <- sel_meth

cat("Genes selected at lambda.min:", length(sel_meth), "\n")
if (length(sel_meth) > 0) cat(paste(sel_meth, collapse = ", "), "\n")

train_meth$risk_meth   <- calc_risk(train_meth, sel_meth, coef_meth_v)
test_meth$risk_meth    <- calc_risk(test_meth,  sel_meth, coef_meth_v)
meth_cohort$risk_meth  <- calc_risk(meth_cohort, sel_meth, coef_meth_v)

# Clinical-only baseline on the SAME 701-patient cohort
cox_clin_meth <- coxph(Surv(os_time, os_event) ~ age + gender + stage_group + histology,
                       data = train_meth)
lp_clin_meth  <- predict(cox_clin_meth, newdata = test_meth, type = "lp")
c_clin_meth   <- safe_cindex(lp_clin_meth, test_meth$os_time, test_meth$os_event)
c_clin_meth_ci <- safe_cindex_ci(lp_clin_meth, test_meth$os_time, test_meth$os_event)

cox_meth_tr <- coxph(Surv(os_time, os_event) ~ risk_meth, data = train_meth)
lp_meth_te  <- predict(cox_meth_tr, newdata = test_meth, type = "lp")
c_meth_te   <- safe_cindex(lp_meth_te, test_meth$os_time, test_meth$os_event)
c_meth_te_ci <- safe_cindex_ci(lp_meth_te, test_meth$os_time, test_meth$os_event)

cat(sprintf("\nGene-level Methylation C-index (Validation): %.4f (95%% CI %.3f-%.3f)\n",
            c_meth_te, c_meth_te_ci["lower"], c_meth_te_ci["upper"]))
cat(sprintf("Clinical-Only C-index (same cohort)        : %.4f (95%% CI %.3f-%.3f)\n",
            c_clin_meth, c_clin_meth_ci["lower"], c_clin_meth_ci["upper"]))

# NEW: comparison to manuscript value
cat("\n--- Comparison to manuscript ---\n")
cat(sprintf("Manuscript (Phase 2 Results): gene-level = 0.545 ; CpG-level = 0.557\n"))
cat(sprintf("Current run               : gene-level = %.3f\n", c_meth_te))
delta_gene <- abs(c_meth_te - 0.545)
cat(sprintf("Absolute difference       : %.4f\n", delta_gene))
if (delta_gene < 0.01) {
  cat(">>> Match: gene-level C-index within 0.01 of manuscript. No update needed.\n")
} else {
  cat(">>> WARNING: gene-level C-index differs by >= 0.01. Update manuscript.\n")
}

## =============================================================================
## BLOCK 2 — CpG-LEVEL METHYLATION LASSO
## =============================================================================
cat("\n=== CpG-LEVEL METHYLATION LASSO ===\n")
promoter_probes <- intersect(anno_prom$Name, rownames(meth_patient))
meth_prom    <- meth_patient[promoter_probes, , drop = FALSE]
common_pats  <- intersect(colnames(meth_prom), meth_cohort$patient_id)
meth_prom_df <- as.data.frame(t(meth_prom[, common_pats]))
meth_prom_df$patient_id <- rownames(meth_prom_df)

meth_cpg <- merge(meth_prom_df,
                  meth_cohort[, c("patient_id","os_time","os_event",
                                  "age","gender","stage_group","histology")],
                  by = "patient_id")

# Impute missing CpG values safely
for (j in promoter_probes) {
  if (!j %in% colnames(meth_cpg)) next
  miss <- is.na(meth_cpg[[j]])
  if (any(miss)) {
    med <- median(meth_cpg[[j]], na.rm = TRUE)
    if (!is.na(med)) {
      meth_cpg[[j]][miss] <- med
    } else {
      meth_cpg[[j]][miss] <- 0
    }
  }
}

meth_cpg <- meth_cpg[complete.cases(meth_cpg[, c("patient_id", "os_time", "os_event",
                                                 "age", "gender", "stage_group", "histology")]), ]

cat("CpG cohort patients:", nrow(meth_cpg), "\n")

set.seed(123)
strata_meth <- interaction(meth_cpg$os_event, meth_cpg$stage_group, meth_cpg$histology)
tr_i <- caret::createDataPartition(strata_meth, p = 0.7, list = FALSE)
train_cpg <- meth_cpg[ tr_i, ]
test_cpg  <- meth_cpg[-tr_i, ]

# NEW: report CpG split sizes
cat(sprintf("CpG split: train = %d | test = %d\n",
            nrow(train_cpg), nrow(test_cpg)))

# Safe univariate screening
probe_to_gene <- sapply(probe_to_genes[match(promoter_probes, anno_prom$Name)],
                        function(x) x[1])
probe_to_gene[is.na(probe_to_gene)] <- "unknown"

univ_p <- sapply(promoter_probes, function(pp) {
  if (!pp %in% colnames(train_cpg)) return(NA)
  vals <- train_cpg[[pp]]
  if (length(vals) == 0 || all(is.na(vals))) return(NA)
  sdev <- sd(vals, na.rm = TRUE)
  if (is.na(sdev) || sdev == 0) return(NA)
  f <- as.formula(sprintf("Surv(os_time, os_event) ~ `%s`", pp))
  tryCatch(summary(coxph(f, data = train_cpg))$coefficients[1, "Pr(>|z|)"],
           error = function(e) NA)
})

top_df <- data.frame(probe = promoter_probes, gene = probe_to_gene, p = univ_p)
top_df <- top_df[!is.na(top_df$p) & top_df$gene != "unknown", ]
top_probes <- top_df %>% group_by(gene) %>% slice_min(p, n = 1, with_ties = FALSE) %>% pull(probe)
top_probes <- intersect(top_probes, colnames(train_cpg))

x_tr_cpg <- as.matrix(train_cpg[, top_probes])
y_tr_cpg <- Surv(train_cpg$os_time, train_cpg$os_event)
x_te_cpg <- as.matrix(test_cpg[,  top_probes])
y_te_cpg <- Surv(test_cpg$os_time,  test_cpg$os_event)

set.seed(456)
cv_cpg <- cv.glmnet(x_tr_cpg, y_tr_cpg, family = "cox", alpha = 1,
                    nfolds = 10, type.measure = "C")
cf_cpg <- as.matrix(coef(cv_cpg, s = "lambda.min"))
sel_cpg <- rownames(cf_cpg)[cf_cpg[,1] != 0]
coef_cpg <- as.numeric(cf_cpg[sel_cpg, 1]); names(coef_cpg) <- sel_cpg

cat("CpGs entering LASSO:", length(top_probes), "\n")
cat("CpGs retained at lambda.min:", length(sel_cpg), "\n")

train_cpg$risk_cpg <- calc_risk(train_cpg, sel_cpg, coef_cpg)
test_cpg$risk_cpg  <- calc_risk(test_cpg,  sel_cpg, coef_cpg)

cox_cpg_tr <- coxph(Surv(os_time, os_event) ~ risk_cpg, data = train_cpg)
lp_cpg_te  <- predict(cox_cpg_tr, newdata = test_cpg, type = "lp")
c_cpg_te   <- safe_cindex(lp_cpg_te, test_cpg$os_time, test_cpg$os_event)
c_cpg_te_ci <- safe_cindex_ci(lp_cpg_te, test_cpg$os_time, test_cpg$os_event)

cat(sprintf("\nCpG-level C-index (Validation): %.4f (95%% CI %.3f-%.3f)\n",
            c_cpg_te, c_cpg_te_ci["lower"], c_cpg_te_ci["upper"]))

# NEW: comparison to manuscript value
cat("\n--- Comparison to manuscript ---\n")
delta_cpg <- abs(c_cpg_te - 0.557)
cat(sprintf("Current run               : CpG-level = %.3f\n", c_cpg_te))
cat(sprintf("Absolute difference       : %.4f\n", delta_cpg))
if (delta_cpg < 0.01) {
  cat(">>> Match: CpG-level C-index within 0.01 of manuscript. No update needed.\n")
} else {
  cat(">>> WARNING: CpG-level C-index differs by >= 0.01. Update manuscript.\n")
}

## =============================================================================
## BLOCK 3 — FIGURE S7: STRATIFIED KM FOR METHYLATION
## =============================================================================
cat("\n=== FIGURE S7 (STRATIFIED METHYLATION KM) ===\n")

meth_cohort$stage_group <- factor(meth_cohort$stage_group, levels = c("Early","Advanced","Unknown"))
meth_cohort$histology   <- factor(meth_cohort$histology, levels = c("LUAD","LUSC"))

plot_km_meth <- function(d, title, panel) {
  if (nrow(d) < 10) return()
  sub_cutoff <- median(d$risk_meth, na.rm = TRUE)
  d$grp <- factor(ifelse(d$risk_meth > sub_cutoff, "High","Low"), levels = c("Low","High"))
  
  fit <- survfit(Surv(os_time, os_event) ~ grp, data = d)
  lr  <- if(length(unique(d$grp)) > 1) survdiff(Surv(os_time, os_event) ~ grp, data = d)$pvalue else NA
  cidx <- safe_cindex(d$risk_meth, d$os_time, d$os_event)
  
  plot(fit, col = c("#2E86AB", "#A23B72"), lwd = 3.5, ylim = c(0, 1.05),
       xlab = "Time (days)", ylab = "Survival Probability", main = "")
  grid(col = "gray90", lty = 2)
  mtext(paste0(panel, ": ", title), side = 3, line = 5, cex = 1.5, font = 2, adj = 0)
  pval_txt <- ifelse(is.na(lr), "NA", sprintf("%.4f", lr))
  mtext(sprintf("N = %d | Events = %d | C-index = %.3f | Log-rank p = %s",
                nrow(d), sum(d$os_event), cidx, pval_txt), side = 3, line = 3, cex = 1.1, adj = 0)
  legend("topright", c("Low","High"), col = c("#2E86AB", "#A23B72"), lwd = 3.5, bty = "n", cex = 1.2)
}

pdf("FigureS7_KM_Methylation_Stratified.pdf", width = 20, height = 14)
par(mfrow = c(2, 3), mar = c(5.5, 5.5, 7, 3), cex.lab = 1.6, cex.axis = 1.4)

plot_km_meth(meth_cohort, "All Patients", "A")
plot_km_meth(meth_cohort[meth_cohort$histology == "LUAD", ], "LUAD", "B")
plot_km_meth(meth_cohort[meth_cohort$histology == "LUSC", ], "LUSC", "C")
plot_km_meth(meth_cohort[meth_cohort$stage_group == "Early", ], "Early Stage", "D")
plot_km_meth(meth_cohort[meth_cohort$stage_group == "Advanced", ], "Advanced Stage", "E")
plot_km_meth(meth_cohort[meth_cohort$stage_group == "Unknown", ], "Unknown Stage", "F")

dev.off()
cat("Saved: FigureS7_KM_Methylation_Stratified.pdf\n")

## =============================================================================
## BLOCK 4 — SAVE TABLES
## =============================================================================
write.csv(data.frame(Gene = sel_meth, Coef = as.numeric(coef_meth_v)),
          "TableS4_Methylation_Signature.csv", row.names = FALSE)
write.csv(data.frame(Probe = sel_cpg, Coef = as.numeric(coef_cpg)),
          "TableS5_CpG_Signature.csv", row.names = FALSE)

# NEW: verification CSV
verif <- data.frame(
  Analysis = c("Gene-level methylation", "Clinical-only (meth cohort)",
               "CpG-level methylation"),
  Train_n  = c(nrow(train_meth), nrow(train_meth), nrow(train_cpg)),
  Test_n   = c(nrow(test_meth),  nrow(test_meth),  nrow(test_cpg)),
  Cindex   = c(round(c_meth_te, 4), round(c_clin_meth, 4), round(c_cpg_te, 4)),
  CI_lower = c(round(c_meth_te_ci["lower"], 4),
               round(c_clin_meth_ci["lower"], 4),
               round(c_cpg_te_ci["lower"], 4)),
  CI_upper = c(round(c_meth_te_ci["upper"], 4),
               round(c_clin_meth_ci["upper"], 4),
               round(c_cpg_te_ci["upper"], 4)),
  Manuscript_value = c(0.545, 0.640, 0.557),
  Delta_from_manuscript = c(round(c_meth_te - 0.545, 4),
                            round(c_clin_meth - 0.640, 4),
                            round(c_cpg_te - 0.557, 4)),
  stringsAsFactors = FALSE
)
write.csv(verif, "Table_Cindex_Verification.csv", row.names = FALSE)
cat("\nSaved: Table_Cindex_Verification.csv\n")
print(verif, row.names = FALSE)

saveRDS(list(sel_meth = sel_meth, coef_meth_v = coef_meth_v,
             c_meth_te = c_meth_te, c_clin_meth = c_clin_meth,
             c_meth_te_ci = c_meth_te_ci, c_clin_meth_ci = c_clin_meth_ci,
             sel_cpg = sel_cpg, coef_cpg = coef_cpg,
             c_cpg_te = c_cpg_te, c_cpg_te_ci = c_cpg_te_ci),
        "nsclc_script04_methylation.rds")

cat("\nScript 04 complete.\n")