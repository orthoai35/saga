# =============================================================================
# 02_meta_analysis.R : pool site-specific hazard ratios
#
# Site estimates (log HR, SE) come from the CohortMethod outcome models fitted
# locally at each site. They are pooled with a random-effects model (REML).
# Fixed-effect estimates are kept as a sensitivity analysis.
#
# Outputs (manuscript/tables):
#   meta_results.csv       pooled HR for every comparison x outcome x analysis,
#                          for two site sets: all sites / excluding KCCH
#   site_estimates.csv     site-level HRs used in the pooling
#   leave_one_out.csv      leave-one-site-out for the main IPTW analysis
# =============================================================================

source("R/00_config.R")
source("R/01_read_results.R")
library(meta)

dat <- load_results()

site_est <- dat$result |>
  inner_join(COMPARISONS, by = c("target_id", "comparator_id")) |>
  left_join(OUTCOMES |> select(outcome_id, outcome, outcome_label), by = "outcome_id") |>
  left_join(NEGATIVE_CONTROLS |> rename(nc_label = outcome_label), by = "outcome_id") |>
  mutate(outcome = coalesce(outcome, paste0("nc_", outcome_id)),
         outcome_label = coalesce(outcome_label, nc_label),
         negative_control = !is.na(nc_label)) |>
  select(-nc_label) |>
  left_join(ANALYSES, by = "analysis_id")

# ---- Pool one set of site estimates ------------------------------------------
pool <- function(d) {
  d <- d |> filter(is.finite(log_rr), is.finite(se_log_rr), se_log_rr > 0)
  out <- tibble(k = nrow(d),
                sites = paste(d$site, collapse = ";"),
                target_n = sum(d$target_subjects, na.rm = TRUE),
                comparator_n = sum(d$comparator_subjects, na.rm = TRUE),
                target_events = sum(d$target_outcomes, na.rm = TRUE),
                comparator_events = sum(d$comparator_outcomes, na.rm = TRUE),
                any_masked = any(d$masked_outcomes))
  if (nrow(d) == 0) return(out)
  if (nrow(d) == 1) {
    return(mutate(out, hr = exp(d$log_rr),
                  lower = exp(d$log_rr - 1.96 * d$se_log_rr),
                  upper = exp(d$log_rr + 1.96 * d$se_log_rr)))
  }
  m <- tryCatch(
    metagen(TE = d$log_rr, seTE = d$se_log_rr, studlab = d$site, sm = "HR",
            common = TRUE, random = TRUE, method.tau = META_METHOD_TAU,
            method.random.ci = "classic", prediction = TRUE,
            control = list(stepadj = 0.5, maxiter = 1000)),
    error = function(e) NULL)
  if (is.null(m)) return(out)
  mutate(out,
         hr = exp(m$TE.random), lower = exp(m$lower.random),
         upper = exp(m$upper.random), p = m$pval.random,
         hr_fixed = exp(m$TE.common), lower_fixed = exp(m$lower.common),
         upper_fixed = exp(m$upper.common), p_fixed = m$pval.common,
         i2 = m$I2, tau2 = m$tau2, q_p = m$pval.Q,
         pi_lower = exp(m$lower.predict), pi_upper = exp(m$upper.predict))
}

site_sets <- list(
  all       = ALL_SITES,
  excl_kcch = setdiff(ALL_SITES, EXCLUDE_SITES_MAIN)
)

meta_results <- imap_dfr(site_sets, function(sites, set_name) {
  site_est |>
    filter(site %in% sites) |>
    group_by(comparison, comparison_label, outcome_id, outcome, outcome_label,
             negative_control, analysis_id, method, tar_days, analysis_label) |>
    group_modify(~ pool(.x)) |>
    ungroup() |>
    mutate(site_set = set_name, .before = 1)
})

# ---- Leave-one-site-out (main comparison, IPTW, pre-specified TAR) ------------
loo <- map_dfr(names(MAIN_TAR), function(o) {
  a <- analysis_id_for("IPTW", MAIN_TAR[[o]])
  d <- site_est |> filter(comparison == "main", outcome == o, analysis_id == a,
                          site %in% site_sets$excl_kcch)
  map_dfr(unique(d$site), function(s) {
    pool(filter(d, site != s)) |> mutate(outcome = o, omitted_site = s)
  })
})

write.csv(site_est, file.path(TABLE_DIR, "site_estimates.csv"), row.names = FALSE)
write.csv(meta_results, file.path(TABLE_DIR, "meta_results.csv"), row.names = FALSE)
write.csv(loo, file.path(TABLE_DIR, "leave_one_out.csv"), row.names = FALSE)

# ---- Console summary: main analysis, with vs without KCCH -----------------------
main_summary <- meta_results |>
  filter(comparison == "main", !negative_control, method == "IPTW") |>
  filter(tar_days == MAIN_TAR[outcome]) |>
  transmute(site_set, outcome, k,
            `HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", hr, lower, upper),
            p = signif(p, 2), I2 = round(i2 * 100))
print(as.data.frame(main_summary))
message("02_meta_analysis.R done.")
