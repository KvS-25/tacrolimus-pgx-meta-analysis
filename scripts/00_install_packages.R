# Installs the CRAN packages the pipeline needs (skips ones already present).
pkgs <- c("readxl", "dplyr", "tidyr", "stringr", "metafor", "purrr", "tibble",
          "yaml", "writexl")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  message("Installing: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
} else message("All packages already installed.")
cat("R", as.character(getRversion()), "| metafor", as.character(packageVersion("metafor")), "\n")
