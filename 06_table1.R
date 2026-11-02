# ==============================================================================
# Script 06: Table 1 - Cohort Characteristics
#
# Author:   Getachew Derbew Belay
# Date:     October 2026
# Platform: R 4.5.2
# Contact:  derbew1953@gmail.com / d610113006@tmu.edu.tw
#
# Purpose:  Produce Table 1 for the overall, training, and validation sets.
# Produces: Table1_cohort_characteristics.csv
# Save as:  06_table1.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/
#
# Depends:  nsclc_base.rds  (from Script 00)
#
# Seeds:    none
# ==============================================================================

setwd("~/TCGA_Local_Storage/TCGA NSCLC Project")
base <- readRDS("nsclc_base.rds"); list2env(base, envir = environment())
library(dplyr)

summ_cont <- function(x) sprintf("%.1f (%.1f) | %.0f [%.0f]",
                                 mean(x, na.rm = TRUE), sd(x, na.rm = TRUE),
                                 median(x, na.rm = TRUE), IQR(x, na.rm = TRUE))
summ_cat  <- function(x, lvls = NULL) {
  if (!is.null(lvls)) x <- factor(x, levels = lvls)
  t <- table(x); p <- round(prop.table(t) * 100, 1)
  setNames(paste0(t, " (", p, "%)"), names(t))
}

rows <- list()
add_row <- function(var, ov, tr, te) {
  rows[[length(rows) + 1]] <<- data.frame(Variable = var, Overall = ov,
                                          Training = tr, Testing = te,
                                          stringsAsFactors = FALSE)
}

add_row("Age - Mean (SD) | Median [IQR]",
        summ_cont(cohort$age), summ_cont(train_data$age), summ_cont(test_data$age))
add_row("OS time (days) - Mean (SD) | Median [IQR]",
        summ_cont(cohort$os_time), summ_cont(train_data$os_time),
        summ_cont(test_data$os_time))

h_o <- summ_cat(cohort$histology, c("LUAD","LUSC"))
h_t <- summ_cat(train_data$histology, c("LUAD","LUSC"))
h_v <- summ_cat(test_data$histology, c("LUAD","LUSC"))
add_row("LUAD; n (%)", h_o["LUAD"], h_t["LUAD"], h_v["LUAD"])
add_row("LUSC; n (%)", h_o["LUSC"], h_t["LUSC"], h_v["LUSC"])

s_o <- summ_cat(cohort$stage_group, c("Early","Advanced","Unknown"))
s_t <- summ_cat(train_data$stage_group, c("Early","Advanced","Unknown"))
s_v <- summ_cat(test_data$stage_group, c("Early","Advanced","Unknown"))
add_row("Stage Early (I/II); n (%)",      s_o["Early"],    s_t["Early"],    s_v["Early"])
add_row("Stage Advanced (III/IV); n (%)", s_o["Advanced"], s_t["Advanced"], s_v["Advanced"])
add_row("Stage Unknown; n (%)",           s_o["Unknown"],  s_t["Unknown"],  s_v["Unknown"])

g_o <- summ_cat(cohort$gender, c("male","female"))
g_t <- summ_cat(train_data$gender, c("male","female"))
g_v <- summ_cat(test_data$gender, c("male","female"))
add_row("Male; n (%)",   g_o["male"],   g_t["male"],   g_v["male"])
add_row("Female; n (%)", g_o["female"], g_t["female"], g_v["female"])

add_row("Deaths; n (%)",
        sprintf("%d (%.1f%%)", sum(cohort$os_event), 100*mean(cohort$os_event)),
        sprintf("%d (%.1f%%)", sum(train_data$os_event), 100*mean(train_data$os_event)),
        sprintf("%d (%.1f%%)", sum(test_data$os_event), 100*mean(test_data$os_event)))

table1 <- do.call(rbind, rows)
write.csv(table1, "Table1_cohort_characteristics.csv", row.names = FALSE)
print(table1, row.names = FALSE)
cat("Script 06 complete.\n")