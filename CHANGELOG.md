# Changelog

## 1.0.0 (2026-10-17)
- First public version. Analysis script is the v6 + plots + ancestry version used for the MSc project (Umeå University, 5MO134), plus the F7 and F8 figures ported to R.
- Added: config file, `run_all.R` entry point, package installer, Excel table builder, synthetic example workbook and generator.
- Changed: hard-coded Windows paths and the manual `drop_overlap_risk` switch replaced by `config/config.yaml` and a `main` / `sensitivity` command-line argument; outputs go to `results/main` and `results/sensitivity`.
- Model, filters and exclusion rules unchanged from the v6 script.
