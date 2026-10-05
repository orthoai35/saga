# =============================================================================
# SAGA — Spinal versus General Anesthesia for lower-extremity surgery
# 00_config.R : paths, site list, study IDs, labels, analysis options
#
# Every other script starts with source("R/00_config.R").
# Working directory must be the code folder (in_progress/code).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
})

# ---- Paths ------------------------------------------------------------------
# The code folder lives at <project>/in_progress/code. We walk upward from the
# working directory until we find the folder that contains raw_data/.
find_project_root <- function(start = getwd()) {
  d <- normalizePath(start, winslash = "/")
  repeat {
    if (dir.exists(file.path(d, "raw_data"))) return(d)
    parent <- dirname(d)
    if (parent == d) stop("Could not find a parent folder containing raw_data/. ",
                          "Open in_progress/code as the working directory.")
    d <- parent
  }
}
PROJECT_ROOT <- Sys.getenv("SAGA_ROOT", find_project_root())
PACKAGE_DIR  <- file.path(PROJECT_ROOT, "raw_data", "packages")   # <SITE>.zip
TABLE_DIR    <- file.path(PROJECT_ROOT, "manuscript", "tables")
FIGURE_DIR   <- file.path(PROJECT_ROOT, "manuscript", "figures")
dir.create(TABLE_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)

# Path of the CohortMethod export inside each site zip
EXPORT_PATH <- "script/extras/result/export"

# ---- Sites ------------------------------------------------------------------
# Every <SITE>.zip in raw_data/packages is analysed. Adding KYUH.zip later is
# enough: re-run run_all.R and all tables/figures are regenerated.
ALL_SITES <- sort(sub("\\.zip$", "", list.files(PACKAGE_DIR, pattern = "\\.zip$")))

# KCCH contributed very few spinal-anesthesia patients (n = 57 after trimming).
# Main analyses are run with and without it (see 02_meta_analysis.R).
EXCLUDE_SITES_MAIN <- c("KCCH")

SITE_NAMES <- c(
  AUMC  = "Ajou University Medical Center",
  DCMC  = "Daegu Catholic University Medical Center",
  DSMC  = "Keimyung University Dongsan Medical Center",
  EUMC  = "Ewha Womans University Medical Center",
  GNUH  = "Gyeongsang National University Hospital",
  HUMC  = "Hallym University Medical Center",          # to be confirmed
  ISH   = "International St. Mary's Hospital",
  KCCH  = "Korea Cancer Center Hospital",               # to be confirmed
  KDH   = "Kyung Hee University Hospital at Gangdong",
  KHMC  = "Kyung Hee University Medical Center",
  KWMC  = "Kangwon National University Hospital",
  KYUH  = "Konyang University Hospital",
  PNUH  = "Pusan National University Hospital",
  SCHBC = "Soonchunhyang University Bucheon Hospital",
  SCHCA = "Soonchunhyang University Cheonan Hospital",
  SCHGM = "Soonchunhyang University Gumi Hospital",
  SCHSU = "Soonchunhyang University Seoul Hospital",
  WKUH  = "Wonkwang University Hospital"
)

# ---- Cohort IDs (ATLAS) -------------------------------------------------------
COMPARISONS <- tibble::tribble(
  ~comparison,        ~target_id, ~comparator_id, ~comparison_label,
  "main",             8045,       8047,           "Lower-extremity surgery",
  "tjr",              8000,       7998,           "Total joint replacement",
  "hip_fracture",     7829,       7828,           "Hip fracture surgery"
)

OUTCOMES <- tibble::tribble(
  ~outcome_id, ~outcome,      ~outcome_label,                 ~role,
  8077,        "mortality",   "All-cause in-hospital death",  "primary",
  7639,        "delirium",    "Postoperative delirium",       "secondary",
  7996,        "mace",        "MACE",                         "secondary",
  7761,        "pneumonia",   "Pneumonia or respiratory failure", "secondary",
  7739,        "vte",         "DVT/PTE",                      "secondary"
)

NEGATIVE_CONTROLS <- tibble::tribble(
  ~outcome_id, ~outcome_label,
  380038,      "Viral conjunctivitis",
  4291005,     "Viral hepatitis",
  440193,      "Wristdrop",
  4115367,     "Pain of joint of wrist"
)

# CohortMethod analysis IDs (from cohort_method_analysis.csv)
ANALYSES <- tibble::tribble(
  ~analysis_id, ~method, ~tar_days, ~analysis_label,
  1,            "IPTW",  60,        "IPTW, 60-day",
  2,            "IPTW",  7,         "IPTW, 7-day",
  3,            "IPTW",  30,        "IPTW, 30-day",
  4,            "PSM",   60,        "1:1 PSM, 60-day",
  5,            "PSM",   7,         "1:1 PSM, 7-day",
  6,            "PSM",   30,        "1:1 PSM, 30-day"
)

# Pre-specified time-at-risk per outcome for the main analysis:
# delirium at 7 days, all other outcomes at 60 days.
MAIN_TAR <- c(mortality = 60, delirium = 7, mace = 60, pneumonia = 60, vte = 60)

# Which analysis_id answers (method, time-at-risk)
analysis_id_for <- function(method, tar) {
  ANALYSES$analysis_id[ANALYSES$method == method & ANALYSES$tar_days == tar]
}

# ---- Meta-analysis options ------------------------------------------------------
META_METHOD_TAU <- "REML"   # between-site variance estimator
META_MODEL      <- "random" # "random" (primary) ; fixed effect reported as sensitivity

# ---- Plot theme -------------------------------------------------------------------
COL_SA <- "#1F6FB2"  # spinal
COL_GA <- "#C8553D"  # general
theme_saga <- function(base_size = 9) {
  theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")
}

save_figure <- function(plot, name, width, height) {
  ggsave(file.path(FIGURE_DIR, paste0(name, ".pdf")), plot,
         width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(file.path(FIGURE_DIR, paste0(name, ".tiff")), plot,
         width = width, height = height, units = "in", dpi = 600,
         compression = "lzw")
  message("  saved figures/", name, ".pdf/.tiff")
}

message("SAGA config: ", length(ALL_SITES), " site packages found: ",
        paste(ALL_SITES, collapse = ", "))
