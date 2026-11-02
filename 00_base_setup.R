# ==============================================================================
# Script 00: Shared Setup, Utilities, and Cohort Assembly
#
# Author:   Getachew Derbew Belay
# Date:     October 2026
# Platform: R 4.5.2
# Contact:  derbew1953@gmail.com / d610113006@tmu.edu.tw
#
# Purpose:  Build the analytic NSCLC cohort from TCGA RNA-seq and clinical
#           data; construct promoter-region gene-averaged methylation matrices;
#           define shared utilities; perform stratified 70/30 split; save all
#           base objects to disk.
# Produces: nsclc_base.rds
#
# Save as:  00_base_setup.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/
#
# Output:   nsclc_base.rds
#
# Depends:  master_nsclc_matrix.rds
#           clinical_final.rds
#           clinical_covariates.rds
#           NSCLC_GDC_Sample_Sheet.tsv
#           Methylation beta files under ~/TCGA_Local_Storage/nsclc_raw_data
#
# Seeds:    split = 123
#
# CHANGES FROM ORIGINAL (Issue 1 fix):
#   - Added rescale_hr() to the shared utilities block. Required by
#     Script 04b (04b_methylation_hr_rescaling.R) for reporting HRs per
#     IQR and per SD without refitting the model.
# ==============================================================================

## ---- 0. Paths ----------------------------------------------------------------
PROJ_DIR     <- "~/TCGA_Local_Storage/TCGA NSCLC Project"
METH_DIR     <- "~/TCGA_Local_Storage/nsclc_raw_data"
SAMPLE_SHEET <- file.path(PROJ_DIR, "NSCLC_GDC_Sample_Sheet.tsv")
EXPR_RDS     <- file.path(PROJ_DIR, "master_nsclc_matrix.rds")
SURV_RDS     <- file.path(PROJ_DIR, "clinical_final.rds")
COV_RDS      <- file.path(PROJ_DIR, "clinical_covariates.rds")

if (!dir.exists(PROJ_DIR)) stop("Project directory not found: ", PROJ_DIR)
setwd(PROJ_DIR)

## ---- 1. Packages -------------------------------------------------------------
pkgs <- c("glmnet","survival","survminer","timeROC","rms","dcurves","caret",
          "dplyr","tidyr","ggplot2","patchwork","pheatmap","data.table","tibble",
          "randomForestSRC","purrr","magick",
          "IlluminaHumanMethylation450kanno.ilmn12.hg19")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (p == "IlluminaHumanMethylation450kanno.ilmn12.hg19") {
      BiocManager::install(p, update = FALSE, ask = FALSE)
    } else install.packages(p)
  }
}
invisible(lapply(pkgs, library, character.only = TRUE))

## ---- 2. Shared utilities -----------------------------------------------------
calc_risk <- function(data, genes, coefs) {
  genes <- intersect(genes, colnames(data))
  if (length(genes) == 0) return(rep(NA_real_, nrow(data)))
  cf <- as.numeric(coefs[genes])
  as.numeric(as.matrix(data[, genes, drop = FALSE]) %*% cf)
}

safe_cindex <- function(train_fit, newdata, time, event) {
  lp <- tryCatch(predict(train_fit, newdata = newdata, type = "lp"),
                 error = function(e) rep(NA_real_, nrow(newdata)))
  if (all(is.na(lp))) return(NA_real_)
  cc <- tryCatch(
    unname(summary(survival::coxph(survival::Surv(time, event) ~ lp))$concordance[1]),
    error = function(e) NA_real_)
  as.numeric(cc)
}

# Rescale a fitted Cox coefficient to interpretable units without refitting.
# Returns a one-row data.frame with HR, 95% CI, and p-value at the chosen scale.
# Usage: rescale_hr(beta, se, scale_factor, "per IQR")
rescale_hr <- function(beta, se, scale_factor, scale_label = "per unit") {
  if (!is.numeric(beta) || length(beta) != 1 || is.na(beta)) {
    stop("rescale_hr: beta must be a single non-NA numeric value")
  }
  if (!is.numeric(se) || length(se) != 1 || is.na(se)) {
    stop("rescale_hr: se must be a single non-NA numeric value")
  }
  if (!is.numeric(scale_factor) || length(scale_factor) != 1 ||
      is.na(scale_factor) || scale_factor == 0) {
    stop("rescale_hr: scale_factor must be a single non-zero numeric value")
  }
  hr    <- exp(beta * scale_factor)
  ci_lo <- exp((beta - 1.96 * se) * scale_factor)
  ci_hi <- exp((beta + 1.96 * se) * scale_factor)
  p     <- 2 * pnorm(-abs(beta / se))
  data.frame(
    Scale    = scale_label,
    Factor   = scale_factor,
    HR       = hr,
    CI_lower = ci_lo,
    CI_upper = ci_hi,
    P_value  = p,
    stringsAsFactors = FALSE
  )
}

parse_stage <- function(x) {
  xu <- toupper(trimws(as.character(x)))
  out <- dplyr::case_when(
    grepl("IV",  xu) ~ "Advanced",
    grepl("III", xu) ~ "Advanced",
    grepl("II",  xu) ~ "Early",
    grepl("I",   xu) ~ "Early",
    TRUE ~ "Unknown"
  )
  factor(out, levels = c("Early","Advanced","Unknown"))
}

## ---- 3. Load raw objects -----------------------------------------------------
master_nsclc_matrix <- readRDS(EXPR_RDS)
clinical_final      <- readRDS(SURV_RDS)
clinical_covariates <- readRDS(COV_RDS)

candidate_genes <- c(
  "APC","BRMS1","BTG2","CDKN2A","CHRNB4","DAPK1","EGLN2",
  "HOXA9","HDGFL3","KDM1A","KDM2A","LRRC3B","MGMT",
  "PCDH10","PITX2","RASSF1","RUNX3","SHOX2","SPARC",
  "TGFBI","WIF1","MIR34B","MIR124-1","MIR124-2","MIR124-3"
)

cat("Expression matrix:", nrow(master_nsclc_matrix), "x",
    ncol(master_nsclc_matrix), "\n")

## ---- 4. Expression cohort ---------------------------------------------------
surv_endpoints <- data.frame(
  patient_id             = clinical_final$cases.submitter_id,
  vital_status           = clinical_final$demographic.vital_status,
  days_to_death          = suppressWarnings(as.numeric(clinical_final$demographic.days_to_death)),
  days_to_last_follow_up = suppressWarnings(as.numeric(clinical_final$diagnoses.days_to_last_follow_up)),
  stringsAsFactors = FALSE)
surv_endpoints$os_time  <- ifelse(surv_endpoints$vital_status == "Dead",
                                  surv_endpoints$days_to_death,
                                  surv_endpoints$days_to_last_follow_up)
surv_endpoints$os_event <- as.integer(surv_endpoints$vital_status == "Dead")
surv_endpoints <- surv_endpoints[!is.na(surv_endpoints$os_time) &
                                   surv_endpoints$os_time >= 30, ]

covariates_clean <- data.frame(
  patient_id = clinical_covariates$cases.submitter_id,
  age        = suppressWarnings(as.numeric(clinical_covariates$demographic.age_at_index)),
  gender     = tolower(clinical_covariates$demographic.gender),
  stage      = clinical_covariates$diagnoses.ajcc_pathologic_stage,
  project    = clinical_covariates$project.project_id,
  stringsAsFactors = FALSE)
covariates_clean <- covariates_clean[grepl("TCGA-LUAD|TCGA-LUSC", covariates_clean$project), ]
covariates_clean <- covariates_clean[!duplicated(covariates_clean$patient_id), ]
covariates_clean$stage_group <- parse_stage(covariates_clean$stage)
covariates_clean$histology   <- ifelse(grepl("LUAD", covariates_clean$project), "LUAD",
                                       ifelse(grepl("LUSC", covariates_clean$project), "LUSC", NA))
covariates_clean$histology   <- factor(covariates_clean$histology, levels = c("LUAD","LUSC"))
covariates_clean$gender      <- factor(covariates_clean$gender, levels = c("female","male"))

present_genes <- intersect(candidate_genes, rownames(master_nsclc_matrix))
sample_ids    <- colnames(master_nsclc_matrix)
expr_t        <- as.data.frame(t(master_nsclc_matrix[present_genes, , drop = FALSE]))
expr_t$patient_id <- sample_ids

cohort <- merge(expr_t, surv_endpoints, by = "patient_id")
cohort <- merge(cohort, covariates_clean, by = "patient_id", all.x = TRUE)
cohort <- cohort[!is.na(cohort$stage_group) &
                   !is.na(cohort$histology) &
                   !is.na(cohort$age) &
                   !is.na(cohort$gender), ]
cohort <- cohort[order(cohort$patient_id, -cohort$os_time), ]
cohort <- cohort[!duplicated(cohort$patient_id), ]

cat("Analytic cohort (expression):", nrow(cohort), "\n")
cat("Stage distribution:\n"); print(table(cohort$stage_group))
cat("Histology distribution:\n"); print(table(cohort$histology))

## ---- 5. Methylation matrices (promoter-region only) --------------------------
meth_files <- list.files(METH_DIR,
                         pattern = "\\.methylation_array\\.sesame\\.level3betas\\.txt$",
                         recursive = TRUE, full.names = TRUE)
if (length(meth_files) == 0) stop("No methylation files in ", METH_DIR)
cat("Methylation files found:", length(meth_files), "\n")

sheet <- read.table(SAMPLE_SHEET, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sheet_meth <- sheet[sheet$Data.Type == "Methylation Beta Value", ]

## ---- CRITICAL FIX: use FOLDER UUID (basename of dirname) for mapping --------
folder_uuids <- basename(dirname(meth_files))

## Diagnostic — check sheet ID format
cat("First 3 folder UUIDs:", head(folder_uuids, 3), "\n")
cat("First 3 sheet File.ID:", head(sheet_meth$File.ID, 3), "\n")
cat("Sheet File.ID %in% folder_uuids:", sum(sheet_meth$File.ID %in% folder_uuids), "\n")

## Build folder-UUID -> Case.ID map
uuid_to_patient <- setNames(sheet_meth$Case.ID, sheet_meth$File.ID)

first_meth <- data.table::fread(meth_files[1], data.table = FALSE)
probe_ids  <- first_meth[[1]]

meth_beta <- matrix(NA_real_, nrow = length(probe_ids), ncol = length(meth_files),
                    dimnames = list(probe_ids, folder_uuids))
for (i in seq_along(meth_files)) {
  dat <- data.table::fread(meth_files[i], data.table = FALSE)
  meth_beta[, i] <- dat[[2]]
}

patient_ids <- uuid_to_patient[folder_uuids]
cat("Patient IDs non-NA:", sum(!is.na(patient_ids)), "of", length(patient_ids), "\n")

colnames(meth_beta) <- patient_ids
meth_beta <- meth_beta[, !is.na(colnames(meth_beta)), drop = FALSE]

unique_patients <- unique(colnames(meth_beta))
meth_patient <- matrix(NA_real_, nrow = nrow(meth_beta), ncol = length(unique_patients),
                       dimnames = list(rownames(meth_beta), unique_patients))
for (p in unique_patients) {
  cols <- which(colnames(meth_beta) == p)
  meth_patient[, p] <- if (length(cols) > 1)
    rowMeans(meth_beta[, cols, drop = FALSE], na.rm = TRUE)
  else meth_beta[, cols]
}
cat("meth_patient:", nrow(meth_patient), "probes x", ncol(meth_patient), "patients\n")

if (ncol(meth_patient) == 0) {
  stop("methylation matrix has 0 patients. The sample sheet's File.ID column ",
       "does not match the folder UUIDs. Verify by comparing ",
       "head(sheet_meth$File.ID) to head(folder_uuids).")
}

## ---- Promoter probe filter ----
anno    <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno_df <- as.data.frame(anno)
promoter_groups <- c("TSS1500","TSS200","5'UTR","1stExon")
keep_region <- vapply(anno_df$UCSC_RefGene_Group, function(g) {
  if (is.na(g)) return(FALSE)
  any(strsplit(g, ";")[[1]] %in% promoter_groups)
}, logical(1))
anno_prom <- anno_df[keep_region, c("Name","UCSC_RefGene_Name","UCSC_RefGene_Group")]
anno_prom <- anno_prom[anno_prom$Name %in% rownames(meth_patient), ]
cat("Promoter probes present in meth_patient:", nrow(anno_prom), "\n")

## ---- Exact symbol matching (primary) ----
probe_to_genes <- lapply(strsplit(anno_prom$UCSC_RefGene_Name, ";"),
                         function(x) {
                           x <- x[!is.na(x)]
                           intersect(x, candidate_genes)
                         })
keep_match <- lengths(probe_to_genes) > 0
anno_prom  <- anno_prom[keep_match, ]
probe_to_genes <- probe_to_genes[keep_match]
cat("Exact-match probes mapped to candidate genes:", length(probe_to_genes), "\n")

## ---- Substring fallback (only if exact match returned nothing) ----
if (length(probe_to_genes) == 0) {
  cat("Fallback: substring matching across candidate genes...\n")
  anno_prom <- anno_df[keep_region, c("Name","UCSC_RefGene_Name","UCSC_RefGene_Group")]
  anno_prom <- anno_prom[anno_prom$Name %in% rownames(meth_patient), ]
  pattern    <- paste(candidate_genes, collapse = "|")
  keep_match <- grepl(pattern, anno_prom$UCSC_RefGene_Name)
  anno_prom  <- anno_prom[keep_match, ]
  probe_to_genes <- lapply(anno_prom$UCSC_RefGene_Name, function(name) {
    if (is.na(name)) return(character(0))
    hits <- vapply(candidate_genes,
                   function(g) grepl(g, name, fixed = TRUE),
                   logical(1))
    candidate_genes[hits]
  })
  cat("Fallback probes mapped:", length(probe_to_genes), "\n")
}

## ---- Gene-averaged promoter methylation ----
gene_meth <- matrix(NA_real_, nrow = length(candidate_genes), ncol = ncol(meth_patient),
                    dimnames = list(candidate_genes, colnames(meth_patient)))
for (g in candidate_genes) {
  idx <- which(vapply(probe_to_genes, function(x) g %in% x, logical(1)))
  probes <- anno_prom$Name[idx]
  probes <- intersect(probes, rownames(meth_patient))
  if (length(probes) >= 1) {
    gene_meth[g, ] <- colMeans(meth_patient[probes, , drop = FALSE], na.rm = TRUE)
  }
}
gene_meth <- gene_meth[rowSums(!is.na(gene_meth)) > 0, , drop = FALSE]
cat("Genes with promoter methylation:", nrow(gene_meth), "\n")
if (nrow(gene_meth) > 0) {
  cat("Genes with methylation data:", paste(rownames(gene_meth), collapse = ", "), "\n")
}

if (nrow(gene_meth) == 0) {
  stop("No genes with promoter methylation could be constructed.")
}

## ---- 6. Methylation survival cohort ------------------------------------------
cov_meta <- covariates_clean[, c("patient_id","age","gender","stage_group","histology")]
meth_t   <- as.data.frame(t(gene_meth))
meth_t$patient_id <- colnames(gene_meth)
meth_cohort <- merge(meth_t, surv_endpoints, by = "patient_id")
meth_cohort <- merge(meth_cohort, cov_meta, by = "patient_id", all.x = TRUE)
meth_cohort <- meth_cohort[
  complete.cases(meth_cohort[, c(rownames(gene_meth), "os_time", "os_event",
                                 "age","gender","stage_group","histology")]), ]
cat("Analytic cohort (methylation):", nrow(meth_cohort), "\n")

## ---- 7. Stratified split (event x stage x histology) -------------------------
set.seed(123)
strata_var   <- interaction(cohort$os_event, cohort$stage_group, cohort$histology, drop = TRUE)
train_index  <- caret::createDataPartition(strata_var, p = 0.7, list = FALSE)
train_data   <- cohort[ train_index, ]
test_data    <- cohort[-train_index, ]

set.seed(123)
strata_meth  <- interaction(meth_cohort$os_event, meth_cohort$stage_group,
                            meth_cohort$histology, drop = TRUE)
train_meth_i <- caret::createDataPartition(strata_meth, p = 0.7, list = FALSE)
train_meth   <- meth_cohort[ train_meth_i, ]
test_meth    <- meth_cohort[-train_meth_i, ]

cat("\nExpression split: train =", nrow(train_data), " test =", nrow(test_data), "\n")
cat("Methylation split: train =", nrow(train_meth), " test =", nrow(test_meth), "\n")

## ---- 8. Save ----------------------------------------------------------------
saveRDS(list(
  cohort = cohort, meth_cohort = meth_cohort,
  train_data = train_data, test_data = test_data,
  train_meth = train_meth, test_meth = test_meth,
  gene_meth = gene_meth, meth_patient = meth_patient,
  anno_prom = anno_prom, probe_to_genes = probe_to_genes,
  candidate_genes = candidate_genes, present_genes = present_genes,
  covariates_clean = covariates_clean, surv_endpoints = surv_endpoints
), "nsclc_base.rds")

cat("\nSaved nsclc_base.rds\n")