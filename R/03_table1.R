# =============================================================================
# 03_table1.R : pooled baseline characteristics (Table 1, Supplementary tables)
#
# Each site exports covariate means (proportions) for the target and comparator
# groups before and after propensity-score adjustment. Pooled proportions are
# sum(site proportion x site N) / sum(site N). A covariate missing from a site's
# balance file had zero prevalence at that site, so that site still counts in
# the denominator (full-cohort denominators).
#
# Outputs (manuscript/tables):
#   table1_<comparison>.csv / .docx    pooled characteristics
#   supp_medications_main.csv          all pre-specified medications
# =============================================================================

source("R/00_config.R")
source("R/01_read_results.R")
library(flextable)
library(officer)

SITES_TABLE1 <- setdiff(ALL_SITES, EXCLUDE_SITES_MAIN)
# Covariates are identical across outcomes (no prior-outcome exclusion), so the
# primary-outcome population (60-day IPTW) defines the cohort description.
T1_OUTCOME  <- 8077
T1_ANALYSIS <- analysis_id_for("IPTW", 60)

# Pre-specified comorbidities in the propensity model (condition group eras)
CONDITIONS <- tibble::tribble(
  ~concept_id, ~label,
  320128,   "Hypertension",
  201826,   "Type 2 diabetes mellitus",
  317576,   "Coronary artery disease",
  316139,   "Heart failure",
  313217,   "Atrial fibrillation",
  314054,   "Aortic valve disorder",
  381591,   "Cerebrovascular disease",
  321052,   "Peripheral vascular disease",
  255573,   "Chronic obstructive pulmonary disease",
  46271022, "Chronic kidney disease",
  4064161,  "Liver cirrhosis",
  439777,   "Anemia",
  432585,   "Coagulation disorder",
  443392,   "Malignant neoplasm",
  4182210,  "Dementia",
  373995,   "Delirium (prior)",
  440383,   "Depressive disorder",
  381270,   "Parkinson disease"
) |> mutate(covariate_id = concept_id * 1000 + 209)

AGE_GROUPS <- tibble(covariate_id = (10:19) * 1000 + 3) |>
  mutate(lower = (covariate_id %/% 1000) * 5,
         label = case_when(lower < 60 ~ "50-59", lower < 70 ~ "60-69",
                           lower < 80 ~ "70-79", lower < 90 ~ "80-89",
                           TRUE ~ "≥90"))

# Index (calendar) year, FeatureExtraction covariate id = year * 1000 + 6
YEAR_GROUPS <- tibble(year = 1990:2026, covariate_id = (1990:2026) * 1000 + 6) |>
  mutate(label = case_when(year <= 2009 ~ "Before 2010", year <= 2014 ~ "2010-2014",
                           year <= 2019 ~ "2015-2019", TRUE ~ "2020-2025"))

# Note: sex is not available. The pre-specified covariate list
# (includedCovariateConceptIds) did not contain the gender concepts, so no sex
# covariate was constructed or exported, and sex was not in the PS model.

# Medications shown in Table 1 (all pre-specified agents go to the supplement)
TABLE1_MEDS <- c("aspirin", "clopidogrel", "atorvastatin", "rosuvastatin",
                 "metformin", "insulin glargine", "furosemide", "losartan",
                 "carvedilol", "prednisolone", "tramadol", "haloperidol",
                 "quetiapine", "alprazolam", "lorazepam", "escitalopram")

pooled_prop <- function(bal, ids, n) {
  # bal: balance rows for one comparison; n: per-site N (before/after, T/C)
  # Sites export a cell with fewer than 5 subjects as a negative mean
  # (-minCellCount / N). Such cells are counted as 0 and flagged.
  b <- bal |> filter(covariate_id %in% ids) |>
    mutate(masked = target_mean_before < 0 | comparator_mean_before < 0 |
             target_mean_after < 0 | comparator_mean_after < 0) |>
    group_by(site) |>
    summarise(across(c(target_mean_before, comparator_mean_before,
                       target_mean_after, comparator_mean_after), ~ sum(pmax(.x, 0))),
              masked = any(masked), .groups = "drop")
  n |> left_join(b, by = "site") |>
    mutate(across(c(target_mean_before, comparator_mean_before,
                    target_mean_after, comparator_mean_after), ~ coalesce(.x, 0)),
           masked = coalesce(masked, FALSE))
}

smd_prop <- function(p1, p0) {
  s <- sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
  ifelse(s > 0, (p1 - p0) / s, 0)
}

make_table1 <- function(comp) {
  cmp <- COMPARISONS |> filter(comparison == comp)
  dat <- load_results(SITES_TABLE1)

  # Site N before adjustment = population after the study-population step
  # ("Have at least 1 days at risk"); after = trimmed population
  att <- dat$attrition |>
    filter(target_id == cmp$target_id, comparator_id == cmp$comparator_id,
           outcome_id == T1_OUTCOME, analysis_id == T1_ANALYSIS,
           grepl("days at risk", description)) |>
    mutate(arm = ifelse(exposure_id == cmp$target_id, "t_before", "c_before")) |>
    select(site, arm, subjects) |>
    pivot_wider(names_from = arm, values_from = subjects)
  aft <- dat$result |>
    filter(target_id == cmp$target_id, comparator_id == cmp$comparator_id,
           outcome_id == T1_OUTCOME, analysis_id == T1_ANALYSIS) |>
    transmute(site, t_after = target_subjects, c_after = comparator_subjects)
  n <- inner_join(att, aft, by = "site")

  bal <- load_balance(SITES_TABLE1, cmp$target_id, cmp$comparator_id,
                      T1_OUTCOME, T1_ANALYSIS) |>
    filter(site %in% n$site)

  row_for <- function(ids, label, group) {
    p <- pooled_prop(bal, ids, n)
    tb <- sum(p$target_mean_before * p$t_before) / sum(p$t_before)
    cb <- sum(p$comparator_mean_before * p$c_before) / sum(p$c_before)
    ta <- sum(p$target_mean_after * p$t_after) / sum(p$t_after)
    ca <- sum(p$comparator_mean_after * p$c_after) / sum(p$c_after)
    tibble(group = group, label = label,
           t_before = tb, c_before = cb, smd_before = smd_prop(tb, cb),
           t_after = ta, c_after = ca, smd_after = smd_prop(ta, ca),
           n_t_before = sum(p$t_before), n_c_before = sum(p$c_before),
           n_t_after = sum(p$t_after), n_c_after = sum(p$c_after), masked = any(p$masked))
  }

  rows <- bind_rows(
    map_dfr(unique(AGE_GROUPS$label), function(g)
      row_for(AGE_GROUPS$covariate_id[AGE_GROUPS$label == g], g, "Age group, y")),
    map_dfr(unique(YEAR_GROUPS$label), function(g)
      row_for(YEAR_GROUPS$covariate_id[YEAR_GROUPS$label == g], g, "Index year")),
    map2_dfr(CONDITIONS$covariate_id, CONDITIONS$label,
             function(id, lab) row_for(id, lab, "Comorbidities"))
  )

  # Charlson index: continuous; pooled mean and N-weighted site SMD
  ch <- bal |> filter(covariate_id == 1901) |> inner_join(n, by = "site")
  charlson <- tibble(
    group = "Comorbidity burden", label = "Charlson comorbidity index, mean",
    t_before = sum(ch$target_mean_before * ch$t_before) / sum(ch$t_before),
    c_before = sum(ch$comparator_mean_before * ch$c_before) / sum(ch$c_before),
    smd_before = weighted.mean(ch$std_diff_before, ch$t_before + ch$c_before),
    t_after = sum(ch$target_mean_after * ch$t_after) / sum(ch$t_after),
    c_after = sum(ch$comparator_mean_after * ch$c_after) / sum(ch$c_after),
    smd_after = weighted.mean(ch$std_diff_after, ch$t_after + ch$c_after))

  # Medications: every pre-specified ingredient (drug group era, analysis 409)
  meds_ids <- bal |> filter(covariate_id %% 1000 == 409) |> distinct(covariate_id)
  covnames <- load_covariate_names(SITES_TABLE1) |>
    filter(covariate_analysis_id == 409) |> distinct(covariate_id, .keep_all = TRUE)
  meds <- map_dfr(meds_ids$covariate_id, function(id) {
    nm <- covnames$covariate_name[covnames$covariate_id == id][1]
    row_for(id, sub(".*relative to index: ", "", nm), "Prior medication")
  }) |> arrange(desc(t_before + c_before))

  list(n = n, rows = rows, charlson = charlson, meds = meds,
       n_tot = c(tb = sum(n$t_before), cb = sum(n$c_before),
                 ta = sum(n$t_after), ca = sum(n$c_after)))
}

# ---- Formatting --------------------------------------------------------------
fmt_n_pct <- function(p, n) {
  cnt <- round(p * n)
  pct <- 100 * p
  ifelse(pct > 0 & pct < 0.05, sprintf("%s (<0.1)", format(cnt, big.mark = ",")),
         sprintf("%s (%.1f)", format(cnt, big.mark = ","), pct))
}

build_docx <- function(t1, comp_label, file, top_meds = 12) {
  r <- t1$rows
  fmt_rows <- function(r) tibble(
    Characteristic = ifelse(r$masked, paste0(r$label, "†"), r$label),
    `SA before` = fmt_n_pct(r$t_before, t1$n_tot["tb"]),
    `GA before` = fmt_n_pct(r$c_before, t1$n_tot["cb"]),
    `SMD before` = sprintf("%.3f", abs(r$smd_before)),
    `SA after` = sprintf("%.1f", 100 * r$t_after),
    `GA after` = sprintf("%.1f", 100 * r$c_after),
    `SMD after` = sprintf("%.3f", abs(r$smd_after)))
  head_row <- function(x) tibble(Characteristic = x, `SA before` = "", `GA before` = "",
                                 `SMD before` = "", `SA after` = "", `GA after` = "",
                                 `SMD after` = "")
  ch <- t1$charlson
  body <- bind_rows(
    tibble(Characteristic = "No. of patients",
           `SA before` = format(t1$n_tot["tb"], big.mark = ","),
           `GA before` = format(t1$n_tot["cb"], big.mark = ","), `SMD before` = "",
           `SA after` = format(t1$n_tot["ta"], big.mark = ","),
           `GA after` = format(t1$n_tot["ca"], big.mark = ","), `SMD after` = ""),
    head_row("Age group, y, No. (%)"), fmt_rows(filter(r, group == "Age group, y")),
    head_row("Index year, No. (%)"), fmt_rows(filter(r, group == "Index year")),
    tibble(Characteristic = "Charlson comorbidity index, mean",
           `SA before` = sprintf("%.2f", ch$t_before), `GA before` = sprintf("%.2f", ch$c_before),
           `SMD before` = sprintf("%.3f", abs(ch$smd_before)),
           `SA after` = sprintf("%.2f", ch$t_after), `GA after` = sprintf("%.2f", ch$c_after),
           `SMD after` = sprintf("%.3f", abs(ch$smd_after))),
    head_row("Comorbidities, No. (%)"), fmt_rows(filter(r, group == "Comorbidities")),
    head_row("Prior medications (selected), No. (%)"),
    fmt_rows(t1$meds |> filter(label %in% TABLE1_MEDS) |>
               mutate(label = str_to_sentence(label)))
  )
  indent_rows <- which(!body$Characteristic %in% c("No. of patients",
                        "Charlson comorbidity index, mean") & body$`SA before` != "")
  ft <- flextable(body) |>
    set_header_labels(`SA before` = "Spinal", `GA before` = "General", `SMD before` = "SMD",
                      `SA after` = "Spinal, %", `GA after` = "General, %", `SMD after` = "SMD") |>
    add_header_row(values = c("", "Before adjustment", "After equipoise trimming and IPTW"),
                   colwidths = c(1, 3, 3)) |>
    padding(i = indent_rows, j = 1, padding.left = 12) |>
    bold(i = which(body$`SA before` == ""), j = 1) |>
    bold(i = which(suppressWarnings(as.numeric(body$`SMD after`)) >= 0.1), j = 7) |>
    font(fontname = "Arial", part = "all") |> fontsize(size = 8, part = "all") |>
    align(j = 2:7, align = "right", part = "body") |>
    align(align = "center", part = "header") |>
    theme_booktabs() |> autofit()
  doc <- read_docx() |>
    body_add_par(sprintf("Table 1. Baseline characteristics before and after propensity-score adjustment — %s",
                         comp_label), style = "heading 2") |>
    body_add_flextable(ft) |>
    body_add_par(paste(
      "Data are No. (%) unless otherwise indicated. Before adjustment: all eligible patients;",
      "after adjustment: patients within the preference-score equipoise region (0.25-0.75), with",
      "covariate means after propensity-score stratification and inverse probability of treatment",
      "weighting. Proportions are pooled across sites with the full site population as denominator;",
      "a covariate absent from a site's export was counted as 0 at that site. Comorbidities and",
      "medications were assessed any time before the index date. SMD, absolute standardized mean",
      "difference (values >= 0.1 in bold). †Sites suppress counts below 5; suppressed cells were counted as 0. Sex was not among the exported covariates.",
      "The complete list of pre-specified medications is given in the Supplement. Abbreviations: IPTW, inverse probability of treatment weighting;",
      "SMD, standardized mean difference."), style = "Normal")
  print(doc, target = file)
}

for (comp in COMPARISONS$comparison) {
  label <- COMPARISONS$comparison_label[COMPARISONS$comparison == comp]
  message("Table 1: ", label)
  t1 <- make_table1(comp)
  out <- bind_rows(t1$rows, t1$charlson, t1$meds)
  write.csv(out, file.path(TABLE_DIR, paste0("table1_", comp, ".csv")), row.names = FALSE)
  write.csv(t1$n, file.path(TABLE_DIR, paste0("table1_site_n_", comp, ".csv")), row.names = FALSE)
  build_docx(t1, label, file.path(TABLE_DIR, paste0("table1_", comp, ".docx")))
}
message("03_table1.R done.")
