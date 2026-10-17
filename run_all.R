# One-and-done: installs packages, runs main + sensitivity analyses, builds the Excel tables.
# Usage (from the repo root):  Rscript run_all.R [config/config.yaml]
args <- commandArgs(trailingOnly = TRUE)
cfg <- if (length(args) >= 1) args[1] else "config/config.yaml"
if (!file.exists("scripts/01_meta_analysis.R")) stop("Run this from the repo root.")
if (!file.exists(cfg)) { file.copy("config/config.yaml.template", cfg); message("Created ", cfg, " from the template.") }
rscript <- file.path(R.home("bin"), "Rscript")
step <- function(label, ...) {
  message("\n=== ", label, " ===")
  status <- system2(rscript, c(...))
  if (status != 0) stop("Step failed: ", label)
}
step("0/4 packages",    "scripts/00_install_packages.R")
step("1/4 main run",    "scripts/01_meta_analysis.R", "main", cfg)
step("2/4 sensitivity", "scripts/01_meta_analysis.R", "sensitivity", cfg)
step("3/4 main again (adds the sensitivity overlay to F1)", "scripts/01_meta_analysis.R", "main", cfg)
step("4/4 tables",      "scripts/02_build_tables.R", cfg)
message("\nDone. Results in: ", yaml::read_yaml(cfg)$output_dir)
