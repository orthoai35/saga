# =============================================================================
# 01_read_results.R : read the CohortMethod export tables of every site
#
# Files are read straight out of raw_data/packages/<SITE>.zip, without
# unzipping to disk, so no temporary folders or caches are left behind.
# Only aggregate export tables are read (no patient-level files).
# =============================================================================

read_site_table <- function(site, table) {
  zip_file <- file.path(PACKAGE_DIR, paste0(site, ".zip"))
  inner    <- file.path(EXPORT_PATH, paste0(table, ".csv"))
  con <- unz(zip_file, inner)
  df  <- tryCatch(utils::read.csv(con, stringsAsFactors = FALSE),
                  error = function(e) NULL)
  if (is.null(df)) {
    warning(site, ": ", table, ".csv not found")
    return(NULL)
  }
  # database_id in the packages is the CDM snapshot name; replace with site code
  df$database_id <- NULL
  df$site <- site
  df
}

read_all_sites <- function(table, sites = ALL_SITES) {
  bind_rows(lapply(sites, read_site_table, table = table))
}

# Counts below 5 are exported as -5 (minimum cell count). Keep them as NA and
# flag them so that totals can be footnoted.
unmask <- function(x) ifelse(!is.na(x) & x < 0, NA_real_, as.numeric(x))

load_results <- function(sites = ALL_SITES) {
  message("Reading cohort_method_result ...")
  res <- read_all_sites("cohort_method_result", sites) |>
    mutate(across(c(target_subjects, comparator_subjects,
                    target_outcomes, comparator_outcomes), unmask,
                  .names = "{.col}"),
           masked_outcomes = is.na(target_outcomes) | is.na(comparator_outcomes))
  message("Reading attrition ...")
  att <- read_all_sites("attrition", sites) |> mutate(subjects = unmask(subjects))
  list(result = res, attrition = att)
}

load_balance <- function(sites = ALL_SITES, target_id, comparator_id,
                         outcome_id, analysis_id) {
  message("Reading covariate_balance (", length(sites), " sites) ...")
  bind_rows(lapply(sites, function(s) {
    b <- read_site_table(s, "covariate_balance")
    if (is.null(b)) return(NULL)
    b[b$target_id == target_id & b$comparator_id == comparator_id &
        b$outcome_id == outcome_id & b$analysis_id == analysis_id, ]
  }))
}

load_covariate_names <- function(sites = ALL_SITES) {
  read_all_sites("covariate", sites) |>
    distinct(covariate_id, covariate_name, covariate_analysis_id)
}
