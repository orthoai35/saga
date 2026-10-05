# =============================================================================
# fetch_references.R : format references (AMA style) from DOIs via Crossref
#
# Reads  ../../manuscript/template/references_keys.csv  (key, doi, manual)
# Writes ../../manuscript/template/references.csv       (key, doi, citation)
# Run once, or whenever a reference is added. Needs internet access.
# Every formatted entry should still be checked against the PDF in literature/.
# =============================================================================

source("R/00_config.R")
library(httr2)

keys_file <- file.path(PROJECT_ROOT, "manuscript", "template", "references_keys.csv")
out_file  <- file.path(PROJECT_ROOT, "manuscript", "template", "references.csv")
refs <- read.csv(keys_file, stringsAsFactors = FALSE, encoding = "UTF-8")

initials <- function(given) {
  if (is.null(given) || is.na(given)) return("")
  parts <- strsplit(gsub("\\.", " ", given), "[ -]+")[[1]]
  paste0(substr(parts[nchar(parts) > 0], 1, 1), collapse = "")
}

# NLM abbreviations for journals whose Crossref short title is missing or long
JOURNAL_ABBREV <- c(
  "British Journal of Anaesthesia" = "Br J Anaesth",
  "The Journal of Bone and Joint Surgery" = "J Bone Joint Surg Am",
  "The Journal of Arthroplasty" = "J Arthroplasty",
  "JCM" = "J Clin Med",
  "Journal of the American Medical Informatics Association" = "J Am Med Inform Assoc",
  "International Journal of Epidemiology" = "Int J Epidemiol",
  "Journal of Biomedical Informatics" = "J Biomed Inform",
  "Statistics in Medicine" = "Stat Med",
  "Pharmacoepidemiology and Drug" = "Pharmacoepidemiol Drug Saf",
  "Regional Anesthesia and Pain Medicine" = "Reg Anesth Pain Med",
  "J American Geriatrics Society" = "J Am Geriatr Soc",
  "The Lancet" = "Lancet",
  "Proc Natl Acad Sci USA" = "Proc Natl Acad Sci U S A",
  "Medicina" = "Medicina (Kaunas)")

ama <- function(m) {
  au <- vapply(m$author, function(a) {
    if (!is.null(a$family)) paste(a$family, initials(a$given)) else a$name
  }, character(1))
  au <- if (length(au) > 6) paste0(paste(au[1:3], collapse = ", "), ", et al") else paste(au, collapse = ", ")
  title <- gsub("<[^>]+>", "", unlist(m$title)[1])
  title <- trimws(gsub("[[:space:]]+", " ", sub(": Table 1$", "", title)))
  jr <- c(unlist(m$`short-container-title`), unlist(m$`container-title`))[1]
  jr <- gsub("\\.", "", jr)
  if (jr %in% names(JOURNAL_ABBREV)) jr <- JOURNAL_ABBREV[[jr]]
  date_src <- if (!is.null(m$`published-print`)) m$`published-print` else m$issued
  dp <- unlist(date_src$`date-parts`[[1]])
  vol <- if (!is.null(m$volume)) m$volume else NULL
  iss <- if (!is.null(m$issue)) paste0("(", m$issue, ")") else ""
  page <- unlist(m$page)[1]
  pg  <- if (!is.null(page)) paste0(":", page) else if (!is.null(m$`article-number`)) paste0(":", m$`article-number`) else ""
  loc <- if (is.null(vol)) "" else paste0(";", vol, iss, pg)
  sprintf("%s. %s. *%s*. %s%s. doi:%s", au, sub("\\.$", "", title), jr, dp[1], loc, m$DOI)
}

refs$citation <- NA_character_
for (i in seq_len(nrow(refs))) {
  if (!is.na(refs$manual[i]) && nzchar(refs$manual[i])) { refs$citation[i] <- refs$manual[i]; next }
  m <- tryCatch(request(paste0("https://api.crossref.org/works/", refs$doi[i])) |>
                  req_user_agent("saga-refs (mailto:jongmin.lee.35@gmail.com)") |>
                  req_perform() |> resp_body_json(),
                error = function(e) NULL)
  refs$citation[i] <- if (is.null(m)) paste("[Crossref lookup failed]", refs$doi[i]) else ama(m$message)
}
write.csv(refs[, c("key", "doi", "citation")], out_file, row.names = FALSE, fileEncoding = "UTF-8")
cat(paste0(refs$key, ": ", refs$citation), sep = "\n")
