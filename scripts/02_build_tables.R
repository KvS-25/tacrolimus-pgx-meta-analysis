# Collects the CSV outputs of the main and sensitivity runs into one Excel file.
# Usage (repo root): Rscript scripts/02_build_tables.R [config/config.yaml]
args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1) args[1] else "config/config.yaml"
cfg <- yaml::read_yaml(cfg_path)
out_root <- cfg$output_dir
suppressWarnings(library(writexl))

rd <- function(run, f) {
  p <- file.path(out_root, run, f)
  if (file.exists(p)) read.csv(p, check.names = FALSE) else NULL
}
fc_main <- rd("main", "all_foldchanges.csv")
fc_sens <- rd("sensitivity", "all_foldchanges.csv")
if (is.null(fc_main)) stop("No main-run results in ", out_root, "/main. Run 01_meta_analysis.R first.")

is_sub <- function(d) if (is.null(d)) NULL else d[grepl("\\|", d$model), ]
not_sub <- function(d) if (is.null(d)) NULL else d[!grepl("\\||Combined", d$model), ]

sheets <- list(
  Pooled_main        = not_sub(fc_main),
  Pooled_sensitivity = not_sub(fc_sens),
  Combined_diagnostic = if (!is.null(fc_main)) fc_main[grepl("^Combined", fc_main$model), ] else NULL,
  CYP3A5_subgroups   = is_sub(fc_main)[grepl("^CYP3A5", is_sub(fc_main)$model) & !grepl("ancestry", is_sub(fc_main)$model), ],
  Ancestry_CYP3A5    = rd("main", "ethnicity_cyp3a5_stratified.csv"),
  Ancestry_study_lvl = rd("main", "ethnicity_study_level_summary.csv"),
  Ancestry_tests     = rd("main", "ethnicity_study_level_tests.csv"),
  Ancestry_ABCB1_3435 = rd("main", "ethnicity_abcb1_3435_stratified.csv"),
  Leave_one_out      = rd("main", "loo_cyp3a5.csv"),
  Small_study_tests  = rd("main", "small_study_tests.csv"),
  Study_level        = rd("main", "study_level_contrasts.csv"),
  Studies_excluded   = rd("main", "studies_excluded.csv")
)
sheets <- sheets[!vapply(sheets, is.null, logical(1))]
sheets <- sheets[vapply(sheets, nrow, integer(1)) > 0]
out <- file.path(out_root, "Meta_Results_Tables.xlsx")
write_xlsx(sheets, out)
cat("Wrote", out, "with sheets:", paste(names(sheets), collapse = ", "), "\n")
