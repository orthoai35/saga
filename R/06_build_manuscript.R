# =============================================================================
# 06_build_manuscript.R : fill the manuscript template with current results
#
# The manuscript text lives outside this repository
# (../../manuscript/template/*.md). Numbers are written as {{key}} placeholders
# and filled here from the csv outputs of scripts 02-05, so the Word files are
# regenerated whenever a site is added or the analysis is re-run.
#
# Markdown subset understood by the converter:
#   # / ## / ### headings, blank-line separated paragraphs,
#   **bold**, *italic*, ^superscript^ (citations), "- " bullets,
#   a line "\pagebreak" for a page break.
# Outputs: manuscript/<name>.docx for every template file.
# =============================================================================

source("R/00_config.R")
library(officer)

TEMPLATE_DIR <- file.path(PROJECT_ROOT, "manuscript", "template")
OUT_DIR      <- file.path(PROJECT_ROOT, "manuscript")

meta_res <- read.csv(file.path(TABLE_DIR, "meta_results.csv"))
loo      <- read.csv(file.path(TABLE_DIR, "leave_one_out.csv"))
flow     <- read.csv(file.path(TABLE_DIR, "flow_counts.csv"))
t1       <- read.csv(file.path(TABLE_DIR, "table1_main.csv"))
bal_sum  <- read.csv(file.path(TABLE_DIR, "balance_summary_by_site.csv"))
sites    <- read.csv(file.path(TABLE_DIR, "etable_sites.csv"))

# ---- Formatting helpers (JAMA style) -----------------------------------------------
n_fmt  <- function(x) formatC(round(x), big.mark = ",", format = "d")
p1     <- function(x) sprintf("%.1f", x)
hr_ci  <- function(d) sprintf("HR, %.2f; 95%% CI, %.2f-%.2f", d$hr, d$lower, d$upper)
hr_ci_short <- function(d) sprintf("%.2f (95%% CI, %.2f-%.2f)", d$hr, d$lower, d$upper)
p_fmt  <- function(p) ifelse(p < 0.001, "P < .001", paste0("P = ", sub("^0\\.", ".", sprintf("%.3f", p))))

main_spec <- function(d) d |> filter(!negative_control, tar_days == MAIN_TAR[outcome])
pick <- function(set, cmp, method) {
  meta_res |> filter(site_set == set, comparison == cmp, method == !!method) |> main_spec()
}

V <- list()
V$n_sites_total <- length(ALL_SITES)
V$n_sites_main  <- length(setdiff(ALL_SITES, EXCLUDE_SITES_MAIN))
V$excluded_site <- paste(EXCLUDE_SITES_MAIN, collapse = ", ")
V$excl_site_sa  <- n_fmt(sum(sites$spinal_analyzed[sites$site %in% EXCLUDE_SITES_MAIN]))
main_sites <- sites |> filter(in_main_analysis)
V$year_start <- format(min(as.Date(main_sites$first_index)), "%Y")
V$year_end   <- format(max(as.Date(main_sites$last_index)), "%Y")

fl <- function(arm, step) flow$n[flow$arm == arm & flow$sequence_number == step]
V$sa_elig <- n_fmt(fl("Spinal", 1));  V$ga_elig <- n_fmt(fl("General", 1))
V$n_elig  <- n_fmt(fl("Spinal", 1) + fl("General", 1))
V$sa_first <- n_fmt(fl("Spinal", 3)); V$ga_first <- n_fmt(fl("General", 3))
V$n_first <- n_fmt(fl("Spinal", 3) + fl("General", 3))
V$sa_an <- n_fmt(fl("Spinal", 4));    V$ga_an <- n_fmt(fl("General", 4))
V$n_an  <- n_fmt(fl("Spinal", 4) + fl("General", 4))
V$pct_retained <- p1(100 * (fl("Spinal", 4) + fl("General", 4)) /
                       (fl("Spinal", 3) + fl("General", 3)))

row1 <- function(lab) t1[t1$label == lab, ][1, ]
a80 <- t1 |> filter(label %in% c("80-89", "≥90"))
V$age80_sa <- p1(100 * sum(a80$t_before)); V$age80_ga <- p1(100 * sum(a80$c_before))
V$age80_sa_after <- p1(100 * sum(a80$t_after)); V$age80_ga_after <- p1(100 * sum(a80$c_after))
ch <- row1("Charlson comorbidity index, mean")
V$cci_sa <- sprintf("%.2f", ch$t_before); V$cci_ga <- sprintf("%.2f", ch$c_before)
V$cci_smd_before <- sprintf("%.2f", abs(ch$smd_before))
V$cci_sa_after <- sprintf("%.2f", ch$t_after); V$cci_ga_after <- sprintf("%.2f", ch$c_after)
V$cci_smd_after <- sprintf("%.2f", abs(ch$smd_after))
for (lab in c("Hypertension", "Dementia", "Chronic kidney disease", "Heart failure")) {
  r <- row1(lab); k <- tolower(gsub("[^A-Za-z]", "", lab))
  V[[paste0(k, "_sa")]] <- p1(100 * r$t_before); V[[paste0(k, "_ga")]] <- p1(100 * r$c_before)
}
V$cov_min <- n_fmt(min(bal_sum$n_cov)); V$cov_max <- n_fmt(max(bal_sum$n_cov))
V$sites_imb <- sum(bal_sum$n_imb > 0); V$max_imb <- max(bal_sum$n_imb)
V$pct_cov_balanced <- sprintf("%.1f", 100 * (1 - sum(bal_sum$n_imb) / sum(bal_sum$n_cov)))

fill_outcomes <- function(prefix, d) {
  for (i in seq_len(nrow(d))) {
    r <- d[i, ]; k <- paste0(prefix, "_", r$outcome)
    V[[paste0(k, "_hr")]]    <<- hr_ci(r)
    V[[paste0(k, "_hrs")]]   <<- hr_ci_short(r)
    V[[paste0(k, "_hrnum")]] <<- sprintf("%.2f", r$hr)
    V[[paste0(k, "_p")]]     <<- p_fmt(r$p)
    V[[paste0(k, "_i2")]]    <<- sprintf("%.0f%%", 100 * r$i2)
    V[[paste0(k, "_pi")]]    <<- sprintf("%.2f-%.2f", r$pi_lower, r$pi_upper)
    V[[paste0(k, "_k")]]     <<- r$k
    V[[paste0(k, "_ev_sa")]] <<- n_fmt(r$target_events)
    V[[paste0(k, "_ev_ga")]] <<- n_fmt(r$comparator_events)
    V[[paste0(k, "_n_sa")]]  <<- n_fmt(r$target_n)
    V[[paste0(k, "_n_ga")]]  <<- n_fmt(r$comparator_n)
    V[[paste0(k, "_pct_sa")]] <<- p1(100 * r$target_events / r$target_n)
    V[[paste0(k, "_pct_ga")]] <<- p1(100 * r$comparator_events / r$comparator_n)
    V[[paste0(k, "_fixed")]] <<- sprintf("%.2f (95%% CI, %.2f-%.2f)", r$hr_fixed,
                                         r$lower_fixed, r$upper_fixed)
    # E-value (VanderWeele & Ding) for the point estimate and the CI limit closest
    # to 1; outcomes are rare (<15%), so the HR approximates the risk ratio.
    ev <- function(rr) { rr <- ifelse(rr < 1, 1 / rr, rr); rr + sqrt(rr * (rr - 1)) }
    ci_lim <- if (r$upper < 1) r$upper else if (r$lower > 1) r$lower else 1
    V[[paste0(k, "_evalue")]] <<- sprintf("%.2f", ev(r$hr))
    V[[paste0(k, "_evalue_ci")]] <<- sprintf("%.2f", ev(ci_lim))
  }
}
fill_outcomes("main", pick("excl_kcch", "main", "IPTW"))
fill_outcomes("tjr",  pick("excl_kcch", "tjr", "IPTW"))
fill_outcomes("hip",  pick("excl_kcch", "hip_fracture", "IPTW"))
fill_outcomes("psm",  pick("excl_kcch", "main", "PSM"))
fill_outcomes("all",  pick("all", "main", "IPTW"))
fill_outcomes("d30",  meta_res |> filter(site_set == "excl_kcch", comparison == "main",
                                         method == "IPTW", !negative_control, tar_days == 30))
nc <- meta_res |> filter(site_set == "excl_kcch", comparison == "main", negative_control,
                         analysis_id == analysis_id_for("IPTW", 60), k > 1)
V$nc_hep_hr <- hr_ci(nc[nc$outcome_id == 4291005, ])
V$nc_hep_k  <- nc$k[nc$outcome_id == 4291005]
for (o in unique(loo$outcome)) {
  l <- loo[loo$outcome == o, ]
  V[[paste0("loo_", o)]] <- sprintf("%.2f to %.2f", min(l$hr), max(l$hr))
}

# ---- Template rendering ------------------------------------------------------------
render <- function(txt) {
  keys <- unique(regmatches(txt, gregexpr("\\{\\{[a-z0-9_]+\\}\\}", txt))[[1]])
  for (k in keys) {
    nm <- gsub("[{}]", "", k)
    val <- V[[nm]]
    if (is.null(val)) { warning("No value for placeholder ", k); val <- paste0("[", nm, "?]") }
    txt <- gsub(k, as.character(val), txt, fixed = TRUE)
  }
  txt
}

# Inline markup -> officer fpar
inline_fpar <- function(s, style_base = fp_text(font.family = "Arial", font.size = 11)) {
  tokens <- regmatches(s, gregexpr("\\*\\*[^*]+\\*\\*|\\*[^*]+\\*|\\^[^^]+\\^|[^*^]+", s))[[1]]
  runs <- lapply(tokens, function(tk) {
    if (grepl("^\\*\\*", tk)) ftext(gsub("\\*\\*", "", tk), update(style_base, bold = TRUE))
    else if (grepl("^\\*", tk)) ftext(gsub("\\*", "", tk), update(style_base, italic = TRUE))
    else if (grepl("^\\^", tk)) ftext(gsub("\\^", "", tk), update(style_base, vertical.align = "superscript"))
    else ftext(tk, style_base)
  })
  do.call(fpar, c(runs, list(fp_p = fp_par(line_spacing = 2, padding.bottom = 6))))
}

# Citations: [@key] or [@key1;@key2] -> superscript numbers in order of first
# appearance; the line {{references}} becomes the numbered reference list.
REFS_FILE <- file.path(TEMPLATE_DIR, "references.csv")
REFS <- if (file.exists(REFS_FILE)) {
  read.csv(REFS_FILE, stringsAsFactors = FALSE, encoding = "UTF-8")
} else data.frame(key = character(), citation = character())
collapse_nums <- function(n) {
  n <- sort(unique(n)); out <- c(); i <- 1
  while (i <= length(n)) {
    j <- i
    while (j < length(n) && n[j + 1] == n[j] + 1) j <- j + 1
    out <- c(out, if (j - i >= 2) paste0(n[i], "-", n[j]) else if (j > i) paste(n[i], n[j], sep = ",") else n[i])
    i <- j + 1
  }
  paste(out, collapse = ",")
}
cite <- function(txt) {
  groups <- regmatches(txt, gregexpr("\\[@[^]]+\\]", txt))[[1]]
  order <- c()
  for (g in groups) for (k in strsplit(gsub("[][@ ]", "", g), ";")[[1]]) order <- union(order, k)
  missing <- setdiff(order, REFS$key)
  if (length(missing)) warning("Unknown reference keys: ", paste(missing, collapse = ", "))
  for (g in unique(groups)) {
    ks <- strsplit(gsub("[][@ ]", "", g), ";")[[1]]
    txt <- gsub(g, paste0("^", collapse_nums(match(ks, order)), "^"), txt, fixed = TRUE)
  }
  ref_lines <- sprintf("%d. %s", seq_along(order), REFS$citation[match(order, REFS$key)])
  sub("{{references}}", paste(ref_lines, collapse = "\n\n"), txt, fixed = TRUE)
}

md_to_docx <- function(md_file, out_file) {
  txt <- paste(readLines(md_file, encoding = "UTF-8"), collapse = "\n")
  txt <- render(cite(txt))
  blocks <- strsplit(txt, "\n[ \t]*\n")[[1]]
  doc <- read_docx()
  for (b in blocks) {
    b <- trimws(b)
    if (b == "") next
    if (b == "\\pagebreak") { doc <- body_add_break(doc); next }
    if (grepl("^#{1,3} ", b)) {
      lvl <- nchar(sub("^(#+).*", "\\1", b))
      doc <- body_add_par(doc, sub("^#+ ", "", b), style = paste("heading", lvl))
      next
    }
    lines <- strsplit(b, "\n")[[1]]
    if (all(grepl("^- ", lines))) {
      for (l in lines) doc <- body_add_fpar(doc, inline_fpar(paste0("• ", sub("^- ", "", l))))
      next
    }
    doc <- body_add_fpar(doc, inline_fpar(paste(lines, collapse = " ")))
  }
  print(doc, target = out_file)
  message("  wrote ", basename(out_file))
}

templates <- list.files(TEMPLATE_DIR, pattern = "\\.md$", full.names = TRUE)
for (f in templates) {
  md_to_docx(f, file.path(OUT_DIR, sub("\\.md$", ".docx", basename(f))))
}
writeLines(sprintf("%s = %s", names(V), unlist(V)),
           file.path(TABLE_DIR, "manuscript_values.txt"), useBytes = TRUE)
message("06_build_manuscript.R done.")
