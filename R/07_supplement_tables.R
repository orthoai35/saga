# =============================================================================
# 07_supplement_tables.R : supplementary tables (Word)
#
#   eTable 1  participating hospitals: data period and cohort sizes
#   eTable 2  target trial specification and its emulation
#   eTable 3  pooled results for every cohort, outcome and analysis
#   eTable 4  pre-specified medications before and after adjustment (main cohort)
# Output: manuscript/tables/supplement_etables.docx
# =============================================================================

source("R/00_config.R")
library(flextable)
library(officer)

meta_res <- read.csv(file.path(TABLE_DIR, "meta_results.csv"))
sites    <- read.csv(file.path(TABLE_DIR, "etable_sites.csv"))
t1       <- read.csv(file.path(TABLE_DIR, "table1_main.csv"))

style_ft <- function(ft) {
  ft |> font(fontname = "Arial", part = "all") |> fontsize(size = 8, part = "all") |>
    theme_booktabs() |> align(align = "center", part = "header") |> autofit()
}
n_fmt <- function(x) ifelse(is.na(x), "<5", formatC(x, big.mark = ",", format = "d"))

# ---- eTable 1 -----------------------------------------------------------------------
e1 <- sites |>
  transmute(Code = site, Hospital = name,
            `Data period` = paste(format(as.Date(first_index), "%Y-%m"),
                                  format(as.Date(last_index), "%Y-%m"), sep = " to "),
            `Eligible, spinal` = n_fmt(spinal_eligible),
            `Eligible, general` = n_fmt(general_eligible),
            `Analyzed, spinal` = n_fmt(spinal_analyzed),
            `Analyzed, general` = n_fmt(general_analyzed),
            `Main analysis` = ifelse(in_main_analysis, "Yes", "No (sparse spinal cohort)"))

# ---- eTable 2 -----------------------------------------------------------------------
e2 <- tibble::tribble(
  ~`Protocol component`, ~`Target trial`, ~`Emulation in FEEDERNET`,
  "Eligibility", "Adults >=50 y scheduled for lower-extremity surgery; no operation in the prior 60 days",
    "Age >=50 y; lower-extremity surgery procedure on the day of anesthesia; no spinal or general anesthesia in the prior 60 days; >=1 day of prior observation",
  "Treatment strategies", "(1) Spinal anesthesia; (2) general anesthesia for the index operation",
    "Spinal anesthesia procedure concept vs general anesthesia concepts (spinal and epidural concepts excluded from the general cohort)",
  "Assignment", "Random assignment",
    "Patients in the preference-score equipoise region (0.25-0.75); IPTW within PS-quintile strata; covariates prespecified from the literature",
  "Time zero", "Randomization (day of surgery)", "Day of surgery; outcomes counted from day 1",
  "Follow-up", "60 days (7 days for delirium)",
    "Until outcome, end of the outcome window, or end of observation",
  "Outcomes", "Death; delirium; MACE; pneumonia; DVT/PTE",
    "All-cause in-hospital death (60 d); delirium (7 d); MACE, pneumonia/respiratory failure, DVT/PTE (60 d)",
  "Causal contrast", "Intention-to-treat effect", "Observational analogue of the intention-to-treat effect",
  "Unit of assignment", "One assignment per patient",
    "First qualifying operation; patients in both cohorts kept in the earlier; same-day dual entry excluded",
  "Analysis", "Cox model",
    "Site-level IPTW Cox models stratified by PS quintile; random-effects (REML) pooling across hospitals")

# ---- eTable 3 -----------------------------------------------------------------------
e3 <- meta_res |>
  filter(!negative_control) |>
  mutate(Cohort = COMPARISONS$comparison_label[match(comparison, COMPARISONS$comparison)],
         `Site set` = ifelse(site_set == "all", "All hospitals", "Main"),
         Outcome = outcome_label, Analysis = analysis_label,
         `Spinal, events/N` = paste0(n_fmt(target_events), "/", n_fmt(target_n)),
         `General, events/N` = paste0(n_fmt(comparator_events), "/", n_fmt(comparator_n)),
         `HR (95% CI), random` = ifelse(is.na(hr), "NA", sprintf("%.2f (%.2f-%.2f)", hr, lower, upper)),
         `HR (95% CI), fixed` = ifelse(is.na(hr_fixed), "NA",
                                       sprintf("%.2f (%.2f-%.2f)", hr_fixed, lower_fixed, upper_fixed)),
         `I2, %` = ifelse(is.na(i2), "", sprintf("%.0f", 100 * i2)),
         `95% PI` = ifelse(is.na(pi_lower), "", sprintf("%.2f-%.2f", pi_lower, pi_upper)),
         k = k) |>
  arrange(match(comparison, COMPARISONS$comparison), desc(site_set),
          match(outcome, OUTCOMES$outcome), analysis_id) |>
  select(Cohort, `Site set`, Outcome, Analysis, k, `Spinal, events/N`, `General, events/N`,
         `HR (95% CI), random`, `HR (95% CI), fixed`, `I2, %`, `95% PI`)

# ---- eTable 4 -----------------------------------------------------------------------
e4 <- t1 |> filter(group == "Prior medication") |>
  transmute(Medication = str_to_sentence(label),
            `Spinal before, %` = sprintf("%.1f", 100 * t_before),
            `General before, %` = sprintf("%.1f", 100 * c_before),
            `SMD before` = sprintf("%.3f", abs(smd_before)),
            `Spinal after, %` = sprintf("%.1f", 100 * t_after),
            `General after, %` = sprintf("%.1f", 100 * c_after),
            `SMD after` = sprintf("%.3f", abs(smd_after)))

doc <- read_docx() |>
  body_add_par("eTable 1. Participating Hospitals", style = "heading 2") |>
  body_add_flextable(style_ft(flextable(e1))) |>
  body_add_par("Eligible: cohort entry; analyzed: preference-score equipoise region in the 60-day IPTW analysis of the primary outcome. Hospital names should be confirmed by each site.", style = "Normal") |>
  body_add_break() |>
  body_add_par("eTable 2. Specification of the Target Trial and Its Emulation", style = "heading 2") |>
  body_add_flextable(style_ft(flextable(e2)) |> width(j = 2:3, width = 3)) |>
  body_add_break() |>
  body_add_par("eTable 3. Pooled Hazard Ratios for All Cohorts, Outcomes, and Analyses", style = "heading 2") |>
  body_add_flextable(style_ft(flextable(e3))) |>
  body_add_par("Event counts below 5 at a hospital were suppressed and are excluded from totals. k, number of hospitals contributing an estimate; PI, prediction interval.", style = "Normal") |>
  body_add_break() |>
  body_add_par("eTable 4. Prespecified Prior Medications Before and After Adjustment, Lower-Extremity Surgery Cohort", style = "heading 2") |>
  body_add_flextable(style_ft(flextable(e4)))
print(doc, target = file.path(TABLE_DIR, "supplement_etables.docx"))
message("07_supplement_tables.R done.")
