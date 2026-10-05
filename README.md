# SAGA: Spinal versus General Anesthesia for lower-extremity surgery

Analysis code for a multicenter target trial emulation of spinal versus general
anesthesia in adults aged 50 years or older undergoing lower-extremity surgery,
run across hospitals of FEEDERNET, the Korean distributed research network
whose data are standardized to the OMOP Common Data Model (CDM) v5.3.

## Design in brief

| Element | Specification |
|---|---|
| Data | 17 FEEDERNET hospital EHR databases (OMOP CDM v5.3), each analysed locally |
| Population | Age ≥50 y; lower-extremity surgery on the day of spinal or general anesthesia; no operation in the prior 60 days |
| Strategies | Spinal anesthesia vs general anesthesia (general cohort excludes spinal/epidural concepts); patients entering both cohorts are kept only in the first; same-day dual entry removed |
| Time zero | Day of surgery; time at risk starts on day 1 |
| Primary outcome | All-cause death recorded during the hospital episode, 60 days |
| Secondary outcomes | Delirium (7 days); MACE, pneumonia/respiratory failure, DVT/PTE (60 days) |
| Confounding control | L1-regularized propensity model over pre-specified covariates (age group, index year, Charlson index, pre-specified conditions and drugs); preference-score equipoise trimming (0.25–0.75); IPTW within a Cox model stratified on PS quintiles |
| Sensitivity | 1:1 PS matching (caliper 0.2 SD of logit); common 30-day window; fixed-effect pooling; including the smallest site; leave-one-site-out |
| Subgroups | Total joint replacement; hip fracture surgery |
| Pooling | Random-effects meta-analysis (REML) of site log hazard ratios |

Site analyses used OHDSI [CohortMethod](https://github.com/OHDSI/CohortMethod).
The cohort definitions (ATLAS JSON) and the full CohortMethod analysis
specification are in `study_specification/`.

## Repository layout

```
run_all.R                    regenerate every table and figure
R/00_config.R                paths, site list, cohort/outcome IDs, options
R/01_read_results.R          read site result packages (zip, no unzip to disk)
R/02_meta_analysis.R         random-effects pooling, sensitivity, leave-one-out
R/03_table1.R                pooled baseline characteristics (Table 1)
R/04_figures.R               flow diagram, forest plots, PS and balance figures
R/05_summary_numbers.R       numbers quoted in the manuscript
R/06_build_manuscript.R      fill manuscript templates with current results
R/07_supplement_tables.R     supplementary eTables
R/fetch_references.R         AMA references from DOIs (Crossref)
study_specification/         ATLAS cohort JSON and CohortMethod settings
```

## Running

Requirements: R ≥ 4.3 with `dplyr tidyr purrr stringr ggplot2 meta officer flextable patchwork scales httr2`.

```r
# working directory = this folder
source("run_all.R")
```

The scripts expect the site result packages at `../../raw_data/packages/<SITE>.zip`
(each a CohortMethod `export/` folder). Outputs go to `../../manuscript/tables`
and `../../manuscript/figures`. To add a site, drop its zip into the packages
folder and re-run.

## Data availability

Patient-level data cannot leave the participating hospitals and are **not** in
this repository. Only aggregate statistics were shared by each site.
