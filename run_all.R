# =============================================================================
# run_all.R : regenerate every table and figure of the SAGA manuscript
#
# Usage (VS Code / R terminal, working directory = this folder):
#   source("run_all.R")
#
# Input : ../../raw_data/packages/<SITE>.zip  (CohortMethod result packages)
# Output: ../../manuscript/tables, ../../manuscript/figures
# Adding a new site (e.g. KYUH.zip) only requires dropping the zip into
# raw_data/packages and re-running this file.
# =============================================================================

t0 <- Sys.time()
source("R/02_meta_analysis.R")   # pooled HRs, site estimates, leave-one-out
source("R/03_table1.R")          # baseline characteristics
source("R/04_figures.R")         # figures
source("R/05_summary_numbers.R") # numbers quoted in the manuscript text
message(sprintf("All done in %.1f min.", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
