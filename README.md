# Replication Package

**Herpes Zoster Vaccination and Cognitive Decline: A Quasi-Experimental Study**

Xi Lin¹˒²†, Lilac Zihui Zhao¹†, Jonathan Rodean¹, Matthew Leipzig³, Reid Dale³,
Nicole M. Gatto⁴, Anna S. Tiefenbacher¹˒², and Pascal Geldsetzer¹˒⁵˒⁶˒⁷*

¹ Division of Primary Care and Population Health, Department of Medicine, Stanford University, Stanford, CA, USA
² Heidelberg Institute of Global Health, Heidelberg University Hospital, Heidelberg, Germany
³ Department of Cardiothoracic Surgery, Stanford University, Stanford, CA, USA
⁴ Kaiser Permanente Washington Health Research Institute, Seattle, WA, USA
⁵ Stanford Center for Digital Health, Stanford University, Stanford, CA, USA
⁶ The Phil and Penny Knight Initiative for Brain Resilience at the Wu Tsai Neurosciences Institute, Stanford University, Stanford, CA, USA
⁷ Biohub, San Francisco, CA, USA

† These authors contributed equally.
\* Corresponding author: Pascal Geldsetzer ([pgeldsetzer@stanford.edu](mailto:pgeldsetzer@stanford.edu))

This repository reproduces every figure and table in the manuscript and its
Supplement, starting only from raw ACT Study data in
`data/raw/`. `R/00_data_prep.R` builds the analysis panel; every other
script consumes that panel (or, for two Stata scripts, its `.dta` export)
and needs nothing else.

## Table of Contents

- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo](#demo)
- [Instructions for Use](#instructions-for-use)
- [Directory Structure](#directory-structure)
- [Raw Data](#raw-data)
- [Figure and Table Correspondence](#figure-and-table-correspondence)
- [Data Availability](#data-availability)
- [License](#license)
- [Citation](#citation)

## System Requirements

### Hardware

No non-standard hardware is required. The pipeline was developed and
tested on a standard laptop (macOS); it does not require a GPU, a compute
cluster, or more than a few GB of RAM.

### Software

Tested on macOS (Darwin), but all three languages below are
cross-platform and the pipeline is expected to run unmodified on Linux or
Windows.

| Software | Version | Notes |
|---|---|---|
| R | >= 4.3 | packages below installed automatically by `run_all.R` if missing |
| Python | >= 3.9 | only for the LASSO risk-score step |
| Stata | 16 or later | only for the Sun & Abraham estimator and the HTE panels |

**R packages:** `dplyr`, `tidyr`, `readr`, `lubridate`, `stringr`, `purrr`,
`fixest`, `ggplot2`, `broom`, `MASS`, `patchwork`, `glmnet`, `haven`,
`HonestDiD`, `splines`.

**Python packages:**
```bash
python3 -m pip install scikit-learn pandas numpy matplotlib
```

**Stata packages:**
```stata
ssc install reghdfe
ssc install coefplot
ssc install parmest
```

## Installation Guide

1. Install R (>= 4.3), Python (>= 3.9), and Stata (16+) if not already present.
2. Clone or download this repository.
3. From the package root, install the Python and Stata dependencies listed
   above. R dependencies install automatically on first run.

No compilation step is required; typical install time for the R and
Python dependencies is under 5 minutes on a standard internet connection.

## Demo

From the package root:

```bash
Rscript run_all.R
```

This runs the full R pipeline (21 scripts) and then, if Stata is found on
the machine, automatically runs `stata/SunAbraham.do` and
`stata/HTE_all.do` and finalizes Table S22. **Expected run time: about 2
minutes** end to end (R pipeline + both Stata scripts) on a standard
laptop. **Expected output:** every figure (PNG) in `output/figures/` and
every table (CSV) in `output/tables/`, matching the numbering used in the
manuscript and Supplement (see [Figure and Table
Correspondence](#figure-and-table-correspondence)). The console log
reports each script's sample sizes and key estimates as it runs, so a
successful run can be checked against the numbers reported in the paper
without opening the output files.

If Stata isn't auto-detected, the script says so and exits after the R
portion; run the two Stata scripts and the finalization step yourself:

```bash
stata -b do stata/SunAbraham.do
stata -b do stata/HTE_all.do
Rscript R/27_finalize_table_s22.R
```

## Instructions for Use

- **Reproduce the full paper:** run `run_all.R` as above (see
  [Demo](#demo)). This is the intended way to use this package and
  requires no arguments or configuration.
- **Reproduce a single figure or table:** each R script runs standalone
  and only needs `data/analysis_panel.rds` to already exist (built by
  `R/00_data_prep.R`, which runs once as part of `run_all.R`). Use the
  [correspondence table](#figure-and-table-correspondence) below to find
  the script for any given figure or table, then run e.g.
  `Rscript R/08_sunab_estimator.R` directly.
- **Run on different data:** the pipeline expects the raw ACT data
  described under [Raw Data](#raw-data) in `data/raw/`, in the exact
  column format the Data Coordinating Center supplies. Substituting
  differently-formatted data will require adapting `R/00_data_prep.R`.
- **Additional robustness checks**, not part of the main pipeline, are
  described under [Standalone Scripts](#standalone-scripts) below.

## Directory Structure

```
data/
  raw/                            raw ACT data (see "Raw data" below)
  analysis_panel.rds              built by R -- analysis panel, R format
  analysis_panel.dta               built by R -- same panel, Stata format;
                                   the input to both Stata scripts
  risk_scores.csv                 LASSO-predicted dementia-risk score
R/
  00_data_prep.R                  raw ACT data -> analysis panel
  00b_validate_data_prep.R        sanity checks vs. the manuscript's reported sample statistics
  00c_lasso_risk_score.R          calls python/lasso_risk_score.py; merges the result into the panel
  01_event_study_utils.R          shared event-time / regression / plot helpers
  02-08, 17-20, 22-24, 26, 27     one script per figure/table group (see below)
  09-14                           standalone RZV-vs-ZVL comparison analyses; run independently
python/
  lasso_risk_score.py             LASSO dementia-risk score (scikit-learn ElasticNetCV)
stata/
  SunAbraham.do                   Sun & Abraham (2021) estimator (Figure S20)
  HTE_all.do                      every HTE result, via one reusable `hte_panel` program
output/
  figures/                        every figure, PNG
  tables/                         every regression/descriptive result, CSV
run_all.R                         master script
```

## Raw Data

**These data are not public and are not included in this repository.**
`data/raw/` is where the approved raw ACT data must be placed by a user
who has independently obtained access (see [Data
Availability](#data-availability)); the pipeline will not run without it.

The data used here were released under **ACT Proposal #296**
("Association of herpes zoster episodes and cognitive decline in the
Adult Changes in Thought (ACT) Study"), approved by the ACT Proposals and
Publications Committee, covering cognitive-assessment and electronic
health record data. Separate data released under **ACT Proposal #398**
("Association of herpes zoster vaccination and neuroimaging markers in
the Adult Changes in Thought (ACT) Study") cover the neuroimaging data
used in the MRI dose-response analysis (Table S26, Figure S29). Neither
dataset, nor any file derived from it, is distributed with this
repository; `R/00_data_prep.R` documents the exact input files and
columns it expects for a user with approved access to
reconstruct.

## Figure and Table Correspondence

| Script | Produces |
|---|---|
| `R/00c_lasso_risk_score.R` + `python/lasso_risk_score.py` | Figure S22; risk score used in S23 and the HTE Table S22/S21 |
| `R/02_table1_baseline_characteristics.R` | Table 1 |
| `R/03_domains_and_dementia_outcomes.R` | Figure 3, Tables S8-S16; Figure S26, Tables S23-S25; Figures S27-S28 |
| `R/04_robustness_checks.R` | Figures S4, S7, S9, S9f, S9g, S10-S12 (extended window, negative controls, falsification, attrition) |
| `R/05_descriptive_figures.R` | Figure 2 |
| `R/06_figure1_plot.R` | Figure 1 |
| `R/07_conceptual_illustration.R` | Figure S41 |
| `R/08_sunab_estimator.R` | Figure S20 (Sun & Abraham 2021 estimator) |
| `R/17_figS38_age_at_vaccination.R` | Figure S38 |
| `R/18_figS39_S40_age_groups.R` | Figures S39-S40 |
| `R/19_rambachan_roth_sensitivity.R` | Figure S19 (Rambachan & Roth 2023 sensitivity) |
| `R/20_prevalence_ranking.R` | Figure S8 |
| `R/21_cem_matching.R` | Figures S13-S15 (Coarsened Exact Matching) |
| `R/22_chronic_hospitalization_adjusted.R` | Figure S16 |
| `R/23_subdomain_decomposition.R` | Figure S21 |
| `R/24_descriptive_figures_S1_S42.R` | Figures S1, S42 |
| `R/26_additional_figures_tables.R` | Figures S2, S3, S5, S6, S17, S18, S23, S24; Tables S4-S7, S22 (partial) |
| `R/28_table_s3_generalizability.R` | Table S3 |
| `R/29_table_s26_mri_doseresponse.R` | Table S26, Figure S29 |
| `stata/HTE_all.do` + `R/27_finalize_table_s22.R` | Figure 4 (panels A-E), Figures S30-S37, Tables S17-S21, Table S22 (final) |

Negative-control and robustness figures in `R/04` write a companion CSV (`TableS_negctrl_*`, `TableS_falsification_*`, `TableS_attrition_*`) alongside each figure.

### Standalone Scripts

- **`R/09`-`R/14`** (6 scripts): RZV-vs-ZVL comparison designs (direct interaction test, COVID-era confound check, and covariate/propensity-score matching). Each script's header explains what it covers; run standalone, in order, once `data/analysis_panel.rds` exists.
- **`cluster_se_robustness/run_cluster_se_check.R`**: re-estimates Figure 1, the 9 CASI domain subscores, and the 3 dementia-outcome regressions with standard errors clustered at the subject level instead of the heteroskedasticity-robust ("hetero") SEs used throughout the main pipeline. Run from the package root: `Rscript cluster_se_robustness/run_cluster_se_check.R`.

## Data Availability

**The data analyzed in this study are not publicly available and are not
included in this repository.** The data were obtained from the Adult
Changes in Thought (ACT) Study, which is managed by the Kaiser Permanente
Washington Health Research Institute (KPWHRI); this includes electronic
health record, utilization, and cognitive-assessment data, as well as
ACT neuroimaging data. **The authors do not have permission to
redistribute these data, under any license, to any third party.**

Access to ACT Study data requires separate, independent approval and is
granted only through the ACT Study itself, not through the authors of
this paper:

- Requests must be submitted through the **ACT Proposals and Publications
  Committee**.
- General information about the ACT Study, and the current data-access
  request process, is available at
  [https://www.actagingresearch.org/](https://www.actagingresearch.org/).

This repository is provided so that a researcher who has independently
obtained approved access to the ACT data described above can
regenerate every number, figure, and table in the paper from raw data.
Nothing in `data/raw/` is shipped with this repository (see [Raw
Data](#raw-data)); only code, and the fully derived, de-identified
summary statistics it produces (every CSV under `output/tables/`), are
distributed here.

Analytic code is additionally deposited at
[https://github.com/jrodean/ACT-HZ-Vaccination](https://github.com/jrodean/ACT-HZ-Vaccination),
per the manuscript's Data and Code Availability statement.

## License

No license file is currently set for this repository or the deposited
GitHub copy above. [Authors to specify a license before public release.]

## Citation

Lin X†, Zhao LZ†, Rodean J, Leipzig M, Dale R, Gatto NM, Tiefenbacher AS,
Geldsetzer P. **Herpes Zoster Vaccination and Cognitive Decline: A
Quasi-Experimental Study.** [Journal, volume/issue, DOI -- to be filled
in on acceptance.]

For questions about this repository, contact the lead contact: Pascal
Geldsetzer ([pgeldsetzer@stanford.edu](mailto:pgeldsetzer@stanford.edu)).
