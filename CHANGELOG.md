# Changelog

## 1.1.0 (2026-10-22)
- Changed: the combined multi-gene model no longer runs by default, matching the registered protocol (multi-gene evidence is narrative only). It can be switched on with `run_combined_model: true` in the config as a diagnostic; its output is kept out of the pooled tables and figures.
- Added: `docs/deviations.md`, the list of deviations from the original plan and where each appears in the code.

## 1.0.1 (2026-10-22)
- Added the OSF protocol registration (DOI 10.17605/OSF.IO/XAPBR) to the README (badge and a Protocol registration section) and to `CITATION.cff` as a related reference.

## 1.0.0 (2026-10-17)
- First public version. Analysis script is the v6 + plots + ancestry version used for the MSc project (Umeå University, 5MO134), plus the F7 and F8 figures ported to R.
- Added: config file, `run_all.R` entry point, package installer, Excel table builder, synthetic example workbook and generator.
- Changed: hard-coded Windows paths and the manual `drop_overlap_risk` switch replaced by `config/config.yaml` and a `main` / `sensitivity` command-line argument; outputs go to `results/main` and `results/sensitivity`.
- Model, filters and exclusion rules unchanged from the v6 script.
