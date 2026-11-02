# ==============================================================================
# Script 08: Provenance for Figure 1
#
# Author:   Getachew Derbew Belay
# Date:     October 2026
# Platform: R 4.5.2
# Contact:  derbew1953@gmail.com / d610113006@tmu.edu.tw
#
# Purpose:  Copy the hand-drawn Figure 1 (CONSORT flowchart) from the 
#           provenance folder into the main project directory so that it is 
#           included when the reproducibility package is zipped.
# Produces: Figure1_final.pdf (copied to project root)
# Save as:  08_provenance_figure1.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/
#
# Depends:  provenance/Figure1_consort/Figure1_final.pdf
#
# Seeds:    none
# ==============================================================================

setwd("~/TCGA_Local_Storage/TCGA NSCLC Project")

# If the image is in the provenance folder, copy it to the main folder
if (file.exists("provenance/Figure1_consort/Figure1_final.pdf")) {
  file.copy("provenance/Figure1_consort/Figure1_final.pdf", "Figure1_final.pdf", overwrite = TRUE)
  cat("Figure 1 successfully copied to main directory.\n")
} else {
  cat("WARNING: Figure1_final.pdf not found in provenance folder. Please ensure it is provided.\n")
}