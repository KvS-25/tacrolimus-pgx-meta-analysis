# Deviations from the original plan

This review's protocol was registered retrospectively on the Open Science Framework ([10.17605/OSF.IO/XAPBR](https://doi.org/10.17605/OSF.IO/XAPBR)), after searching, screening, data extraction and the main analyses were complete. It is not a preregistration. It lists every known difference between the original project plan and the methods actually carried out, so that readers can judge the risk of bias from analytic flexibility. This page repeats that list as it applies to the code in this repository.

| # | Original plan | What was done | Note |
|---|---|---|---|
| 1 | Pooled combined multi-gene model (CYP3A4 + CYP3A5 + ABCB1) as the novelty claim | Combined multi-gene evidence is reported narratively only. An earlier pooled combined model was run, found invalid and discarded | No mean/SD for joint genotype cells |
| 2 | DerSimonian-Laird random effects, weighted mean difference of C0/D | Multilevel REML model on log C0/D with fold-changes and cluster-robust CIs | Fold-changes are within-study ratios, so comparable across differing dose and assay scales |
| 3 | Weighting by sample size | Weights come from reported SD only and ignore sample size | Per-genotype n not available |
| 4 | Funnel plot and Egger test | Funnel plots with total sample size as precision proxy; exploratory regression and rank correlation; no formal Egger or trim-and-fill | No per-genotype SE or n |
| 5 | Ethnicity: Asian vs Caucasian (and others) | Ancestry classified by keyword from the free-text ethnicity field; analysis added after the main analyses | Post hoc; no African cohorts available |
| 6 | Time window: early under 1 month vs stable over 3 months | Early (under 3 months) vs stable (6 months or more), classified by keyword from the reported time-point labels; time points that fit neither are left unspecified | Reason not recorded at the time |
| 7 | Median/IQR converted to mean/SD (Wan 2014, Luo 2018) | 25 median/IQR rows entered as mean/SD without conversion | Flagged as a limitation |
| 8 | Donor-genotype confound in liver transplants flagged separately | Donor-genotype rows excluded from pooled models; no donor vs recipient subgroup | Liver estimate less precise |
| 9 | Consortium and duplicate cohorts excluded | Master publication used; overlap-risk studies kept in the primary analysis, with a sensitivity run excluding them | Results unchanged materially |
| 10 | Acute rejection and nephrotoxicity as outcomes | Extracted where reported, not meta-analysed | Sparse and inconsistent reporting |
| 11 | Full in-silico toolset | MutationTaster, PROVEAN and I-Mutant not run | Listed as a limitation |

## Where each deviation lives in the code

| # | Script behaviour |
|---|---|
| 1 | The combined model is off by default. Set `run_combined_model: true` in `config/config.yaml` to fit it as a diagnostic. Its output is labelled "diagnostic only", is kept out of the pooled tables and figures, and should not be cited |
| 2, 3 | `rma.mv(yi, vi, mods = ~ grp, random = ~ 1 | study_id, method = "REML")` with `yi = log(C0/D mean)` and `vi = (SD/mean)^2`; cluster-robust CIs |
| 4 | Figure F4 and `small_study_tests.csv` |
| 5 | Section 20 of `scripts/01_meta_analysis.R` and the `ethnicity_*.csv` files |
| 6 | The `time_window` keyword rule in the cleaning step |
| 7 | Median-type rows are counted in the log and pooled as mean/SD |
| 8 | Donor and graft genotype labels are dropped before pooling |
| 9 | `OVERLAP-RISK` Study IDs are kept in the `main` run and dropped in the `sensitivity` run |
| 10, 11 | Not part of this code |
