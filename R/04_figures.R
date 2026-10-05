# =============================================================================
# 04_figures.R : main and supplementary figures
#
# Needs manuscript/tables/meta_results.csv and site_estimates.csv
# (run 02_meta_analysis.R first).
#
#   figure1_flow              study flow (attrition summed over sites)
#   figure2_pooled            pooled HRs, main IPTW analysis, 3 surgical cohorts
#   efigure_site_forest       site-specific HRs for each outcome (main cohort)
#   efigure_sensitivity       pooled HRs across sensitivity analyses
#   efigure_preference        preference-score distributions by site
#   efigure_balance           covariate balance before vs after, by site
# =============================================================================

source("R/00_config.R")
source("R/01_read_results.R")
library(patchwork)

SITES_MAIN <- setdiff(ALL_SITES, EXCLUDE_SITES_MAIN)
meta_res <- read.csv(file.path(TABLE_DIR, "meta_results.csv"))
site_est <- read.csv(file.path(TABLE_DIR, "site_estimates.csv"))

outcome_levels <- OUTCOMES$outcome_label
main_rows <- function(d) {
  d |> filter(!negative_control, method == "IPTW",
              tar_days == MAIN_TAR[outcome])
}
fmt_hr <- function(hr, lo, hi) ifelse(is.na(hr), "NA",
                                      sprintf("%.2f (%.2f-%.2f)", hr, lo, hi))
fmt_ev <- function(e, n) ifelse(is.na(n), NA_character_, sprintf("%s/%s", ifelse(is.na(e), "NA", formatC(e, big.mark = ",", format = "d")), formatC(n, big.mark = ",", format = "d")))

# ---- Forest-plot helper: text columns on the left, HR axis on the right ------------
forest_panel <- function(d, title = NULL, xlim = c(0.25, 2.5)) {
  d <- d |> mutate(row = rev(seq_len(n())))
  txt <- ggplot(d, aes(y = row)) +
    geom_text(aes(x = 0, label = label), hjust = 0, size = 2.7,
              fontface = ifelse(d$is_header, "bold", "plain")) +
    geom_text(aes(x = 2.35, label = ev_t), hjust = 1, size = 2.6) +
    geom_text(aes(x = 3.25, label = ev_c), hjust = 1, size = 2.6) +
    geom_text(aes(x = 4.4, label = hr_txt), hjust = 1, size = 2.6) +
    annotate("text", x = c(0, 2.35, 3.25, 4.4), y = max(d$row) + 1,
             label = c("Outcome", "Spinal, events/N", "General, events/N", "HR (95% CI)"),
             hjust = c(0, 1, 1, 1), size = 2.6, fontface = "bold") +
    scale_x_continuous(limits = c(0, 4.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, max(d$row) + 1.5)) +
    theme_void()
  fp <- ggplot(d, aes(y = row, x = hr, xmin = lower, xmax = upper)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    geom_errorbar(orientation = "y", width = 0.25, colour = "grey20", na.rm = TRUE) +
    geom_point(aes(size = is_pooled), shape = 15, colour = COL_SA, na.rm = TRUE) +
    scale_size_manual(values = c(`TRUE` = 2.4, `FALSE` = 1.6), guide = "none") +
    scale_x_log10(limits = xlim, breaks = c(0.25, 0.5, 1, 2),
                  oob = scales::squish) +
    scale_y_continuous(limits = c(0.5, max(d$row) + 1.5)) +
    annotate("text", x = c(0.6, 1.6), y = max(d$row) + 1, size = 2.5,
             label = c("Favors spinal", "Favors general")) +
    labs(x = "Hazard ratio (log scale)", y = NULL) +
    theme_saga() +
    theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank())
  p <- txt + fp + plot_layout(widths = c(2.3, 1.2))
  if (!is.null(title)) p <- p + plot_annotation(title = title)
  p
}

# ---- Figure 2: pooled HRs, main analysis (sites excluding KCCH) --------------------
fig2_dat <- meta_res |>
  filter(site_set == "excl_kcch") |> main_rows() |>
  mutate(comparison_label = factor(comparison_label, COMPARISONS$comparison_label),
         outcome_label = factor(outcome_label, outcome_levels)) |>
  arrange(comparison_label, outcome_label)

fig2_rows <- fig2_dat |>
  group_by(comparison_label) |>
  group_modify(~ bind_rows(
    tibble(label = as.character(.y$comparison_label), is_header = TRUE),
    tibble(label = paste0("   ", .x$outcome_label,
                          ifelse(.x$tar_days == 7, " (7 d)", " (60 d)")),
           ev_t = fmt_ev(.x$target_events, .x$target_n),
           ev_c = fmt_ev(.x$comparator_events, .x$comparator_n),
           hr = .x$hr, lower = .x$lower, upper = .x$upper,
           hr_txt = fmt_hr(.x$hr, .x$lower, .x$upper), is_header = FALSE))) |>
  ungroup() |>
  mutate(across(c(ev_t, ev_c, hr_txt), ~ coalesce(.x, "")),
         is_pooled = TRUE)
save_figure(forest_panel(fig2_rows), "figure2_pooled", 8.5, 5.2)

# ---- eFigure: site-specific forest plots (main cohort, main analysis) ---------------
site_forest <- map(OUTCOMES$outcome, function(o) {
  a <- analysis_id_for("IPTW", MAIN_TAR[[o]])
  s <- site_est |> filter(comparison == "main", outcome == o, analysis_id == a,
                          site %in% SITES_MAIN) |>
    arrange(site)
  m <- meta_res |> filter(site_set == "excl_kcch", comparison == "main",
                          outcome == o, analysis_id == a)
  rows <- bind_rows(
    tibble(label = s$site,
           ev_t = fmt_ev(s$target_outcomes, s$target_subjects),
           ev_c = fmt_ev(s$comparator_outcomes, s$comparator_subjects),
           hr = s$rr, lower = s$ci_95_lb, upper = s$ci_95_ub, is_pooled = FALSE),
    tibble(label = sprintf("Random effects (I² = %.0f%%)", 100 * m$i2),
           ev_t = fmt_ev(m$target_events, m$target_n),
           ev_c = fmt_ev(m$comparator_events, m$comparator_n),
           hr = m$hr, lower = m$lower, upper = m$upper, is_pooled = TRUE),
    tibble(label = sprintf("95%% prediction interval: %.2f-%.2f", m$pi_lower, m$pi_upper),
           is_pooled = TRUE)) |>
    mutate(ev_t = gsub("NA", "<5", coalesce(ev_t, "")),
           ev_c = gsub("NA", "<5", coalesce(ev_c, "")),
           hr_txt = ifelse(is.na(hr), "", fmt_hr(hr, lower, upper)),
           is_header = FALSE)
  forest_panel(rows, title = sprintf("%s, %d-day",
                                     OUTCOMES$outcome_label[OUTCOMES$outcome == o],
                                     MAIN_TAR[[o]]), xlim = c(0.05, 10))
})
pdf(file.path(FIGURE_DIR, "efigure_site_forest.pdf"), width = 8.5, height = 5.5)
invisible(lapply(site_forest, print)); dev.off()
message("  saved figures/efigure_site_forest.pdf (one page per outcome)")

# ---- eFigure: sensitivity analyses ---------------------------------------------------
sens <- meta_res |>
  filter(comparison == "main", !negative_control) |>
  mutate(spec = case_when(
    site_set == "excl_kcch" & method == "IPTW" & tar_days == MAIN_TAR[outcome] ~ "Main: IPTW, random effects",
    site_set == "excl_kcch" & method == "IPTW" & tar_days == 30 ~ "IPTW, common 30-day window",
    site_set == "excl_kcch" & method == "PSM" & tar_days == MAIN_TAR[outcome] ~ "1:1 PS matching",
    site_set == "all" & method == "IPTW" & tar_days == MAIN_TAR[outcome] ~ "IPTW, including KCCH",
    TRUE ~ NA_character_)) |>
  filter(!is.na(spec))
sens <- bind_rows(sens,
  sens |> filter(spec == "Main: IPTW, random effects") |>
    mutate(spec = "IPTW, fixed effect", hr = hr_fixed, lower = lower_fixed, upper = upper_fixed))
sens_levels <- c("Main: IPTW, random effects", "IPTW, fixed effect", "IPTW, including KCCH",
                 "IPTW, common 30-day window", "1:1 PS matching")
sens_rows <- sens |>
  mutate(outcome_label = factor(outcome_label, outcome_levels),
         spec = factor(spec, sens_levels)) |>
  arrange(outcome_label, spec) |>
  group_by(outcome_label) |>
  group_modify(~ bind_rows(
    tibble(label = as.character(.y$outcome_label), is_header = TRUE),
    tibble(label = paste0("   ", .x$spec),
           ev_t = fmt_ev(.x$target_events, .x$target_n),
           ev_c = fmt_ev(.x$comparator_events, .x$comparator_n),
           hr = .x$hr, lower = .x$lower, upper = .x$upper,
           hr_txt = fmt_hr(.x$hr, .x$lower, .x$upper), is_header = FALSE))) |>
  ungroup() |>
  mutate(across(c(ev_t, ev_c, hr_txt), ~ coalesce(.x, "")), is_pooled = TRUE)
save_figure(forest_panel(sens_rows), "efigure_sensitivity", 8.5, 8.5)

# ---- eFigure: preference-score distributions ------------------------------------------
ps <- read_all_sites("preference_score_dist", SITES_MAIN) |>
  filter(target_id == 8045, comparator_id == 8047,
         analysis_id == analysis_id_for("IPTW", 60)) |>
  pivot_longer(c(target_density, comparator_density), names_to = "group",
               values_to = "density") |>
  mutate(group = ifelse(group == "target_density", "Spinal", "General"))
p_ps <- ggplot(ps, aes(preference_score, density, fill = group, colour = group)) +
  annotate("rect", xmin = 0.25, xmax = 0.75, ymin = -Inf, ymax = Inf,
           fill = "grey92") +
  geom_area(alpha = 0.35, position = "identity", linewidth = 0.3) +
  facet_wrap(~ site, ncol = 4, scales = "free_y") +
  scale_fill_manual(values = c(Spinal = COL_SA, General = COL_GA)) +
  scale_colour_manual(values = c(Spinal = COL_SA, General = COL_GA)) +
  labs(x = "Preference score", y = "Density", fill = NULL, colour = NULL) +
  theme_saga() + theme(axis.text.y = element_blank())
save_figure(p_ps, "efigure_preference", 8, 8)

# ---- eFigure: covariate balance, by site ----------------------------------------------
bal <- load_balance(SITES_MAIN, 8045, 8047, 8077, analysis_id_for("IPTW", 60))
bal_sum <- bal |> group_by(site) |>
  summarise(n_cov = n(), n_imb = sum(abs(std_diff_after) >= 0.1, na.rm = TRUE),
            .groups = "drop") |>
  mutate(strip = sprintf("%s\n%d covariates; %d with |SMD| ≥ 0.1", site, n_cov, n_imb))
p_bal <- bal |> left_join(bal_sum, by = "site") |>
  ggplot(aes(abs(std_diff_before), abs(std_diff_after))) +
  geom_point(alpha = 0.3, size = 0.6, colour = COL_SA) +
  geom_hline(yintercept = 0.1, linetype = 3) + geom_vline(xintercept = 0.1, linetype = 3) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60") +
  facet_wrap(~ strip, ncol = 4) +
  coord_cartesian(xlim = c(0, 0.8), ylim = c(0, 0.8)) +
  labs(x = "|SMD| before adjustment", y = "|SMD| after trimming and IPTW") +
  theme_saga(8)
save_figure(p_bal, "efigure_balance", 8, 9)
write.csv(bal_sum |> select(-strip), file.path(TABLE_DIR, "balance_summary_by_site.csv"),
          row.names = FALSE)

# ---- Figure 1: study flow -------------------------------------------------------------
att <- read_all_sites("attrition", SITES_MAIN) |>
  filter(target_id == 8045, comparator_id == 8047, outcome_id == 8077,
         analysis_id == analysis_id_for("IPTW", 60)) |>
  mutate(arm = ifelse(exposure_id == 8045, "Spinal", "General"),
         subjects = unmask(subjects)) |>
  group_by(arm, sequence_number, description) |>
  summarise(n = sum(subjects, na.rm = TRUE), .groups = "drop")
write.csv(att, file.path(TABLE_DIR, "flow_counts.csv"), row.names = FALSE)

flow_box <- function(arm, x) {
  a <- att |> filter(arm == !!arm) |> arrange(sequence_number)
  steps <- c("Eligible patients (cohort entry)",
             "First qualifying surgery; not in both cohorts;\n≥1 day of prior observation",
             "≥1 day at risk",
             "Within preference-score equipoise (0.25-0.75)\nanalyzed with IPTW")
  tibble(x = x, y = rev(seq_along(steps)) * 1.6,
         label = sprintf("%s anesthesia\n%s\nn = %s", arm, steps,
                         format(a$n, big.mark = ",")))
}
boxes <- bind_rows(flow_box("Spinal", 0), flow_box("General", 3))
p_flow <- ggplot(boxes, aes(x, y)) +
  geom_segment(data = boxes |> group_by(x) |> arrange(desc(y)) |>
                 mutate(yend = lead(y)) |> filter(!is.na(yend)),
               aes(xend = x, yend = yend + 0.55, y = y - 0.55),
               arrow = arrow(length = unit(2, "mm")), colour = "grey40") +
  geom_label(aes(label = label), size = 2.6, label.padding = unit(2.5, "mm"),
             lineheight = 0.95) +
  annotate("text", x = 1.5, y = max(boxes$y) + 1.1, size = 3, fontface = "bold",
           label = sprintf("%d FEEDERNET hospitals, adults ≥50 y undergoing lower-extremity surgery",
                           length(SITES_MAIN))) +
  scale_x_continuous(limits = c(-1.3, 4.3)) +
  scale_y_continuous(limits = c(0.6, max(boxes$y) + 1.5)) +
  theme_void()
save_figure(p_flow, "figure1_flow", 7, 6)

message("04_figures.R done.")
