# =============================================================================
# TACROLIMUS PHARMACOGENETICS META-ANALYSIS - main analysis script
# CYP3A4 (rs2242480, rs35599367), CYP3A5 (rs776746), ABCB1 (rs1045642, rs2032582, rs1128503)
# Models: rma.mv, yi = log(C0/D mean), vi = (SD/mean)^2, random = ~1 | study_id,
#         cluster-robust CIs. Based on v6 + plots + ethnicity (2026-10-15).
# Writes: <output_dir>/<main|sensitivity>/ (log, CSVs, plots/)
# =============================================================================

# ---- 0. SETUP ---------------------------------------------------------------
# Usage (from the repo root):
#   Rscript scripts/01_meta_analysis.R main        [config/config.yaml]
#   Rscript scripts/01_meta_analysis.R sensitivity [config/config.yaml]
args <- commandArgs(trailingOnly = TRUE)
run_mode <- if (length(args) >= 1) args[1] else "main"
cfg_path <- if (length(args) >= 2) args[2] else "config/config.yaml"
if (!run_mode %in% c("main", "sensitivity")) stop("First argument must be 'main' or 'sensitivity'")
if (!file.exists(cfg_path)) stop("Config not found: ", cfg_path, " (run from the repo root)")
cfg <- yaml::read_yaml(cfg_path)

wb_path  <- normalizePath(cfg$workbook, mustWork = FALSE)
wb_sheet <- if (!is.null(cfg$sheet)) cfg$sheet else "Meta_Data_Extraction"
if (!file.exists(wb_path)) stop("Workbook not found: ", wb_path, "\nEdit 'workbook' in ", cfg_path)

drop_overlap_risk <- (run_mode == "sensitivity")  # TRUE = drop OVERLAP-RISK studies
manual_exclude <- as.character(unlist(cfg$manual_exclude))   # Study IDs to drop
run_tag  <- run_mode
out_root <- normalizePath(cfg$output_dir, mustWork = FALSE)
out_dir  <- file.path(out_root, run_mode)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
out_root <- normalizePath(out_root)
sink(file.path(out_dir, paste0("meta_run_", run_tag, ".txt")), split = TRUE)

required_packages <- c("readxl", "dplyr", "tidyr", "stringr",
                       "metafor", "purrr", "tibble")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Package '", pkg, "' missing. Run: Rscript scripts/00_install_packages.R")
  library(pkg, character.only = TRUE)
}
set.seed(if (!is.null(cfg$seed)) cfg$seed else 20260913)
options(stringsAsFactors = FALSE, width = 140)
cat("Run:", run_tag, "| outputs go to:", out_dir, "\n")
cat("Workbook:", wb_path, "\n")
cat("metafor version:", as.character(packageVersion("metafor")), "\n")

# ---- 1. COLUMN NAME CLEANER -------------------------------------------------
safe_clean <- function(df) {
  nm <- names(df)
  nm <- str_replace_all(nm, "[\u00A0\u2007\u202F\uFEFF]", " ")
  nm <- str_squish(nm)
  nm <- str_replace_all(nm, "[^A-Za-z0-9]+", "_")
  nm <- str_to_lower(nm)
  nm <- str_replace_all(nm, "^_+|_+$", "")
  nm <- ifelse(nm == "", paste0("col_", seq_along(nm)), nm)
  nm <- make.unique(nm, sep = "_")
  names(df) <- nm
  df
}

# ---- 2. IMPORT --------------------------------------------------------------
extraction <- read_excel(wb_path, sheet = wb_sheet) %>% safe_clean()
cat("Rows imported:", nrow(extraction), "\n")

# ---- 3. DEFENSIVE RENAME ----------------------------------------------------
extraction <- extraction %>%
  rename_with(~ "snp_rsid", .cols = any_of(c(
    "snp_rsid","snp_rs_id","snp_id","rsid","snp"))) %>%
  rename_with(~ "gene", .cols = any_of(c("gene","gene_name","target_gene"))) %>%
  rename_with(~ "genotype_groups", .cols = any_of(c(
    "genotype_groups","genotype_group","genotypes","genotype"))) %>%
  rename_with(~ "sample_size", .cols = any_of(c(
    "sample_size","n","n_total"))) %>%
  rename_with(~ "population_adult_pediatric", .cols = any_of(c(
    "population_adult_pediatric","population_adult_paediatric",
    "population","age_group"))) %>%
  rename_with(~ "organ_type", .cols = any_of(c(
    "organ_type","organ","transplant_type"))) %>%
  rename_with(~ "tacrolimus_dose_mg_kg_day", .cols = any_of(c(
    "tacrolimus_dose_mg_kg_day","tacrolimus_dose","dose_mg_kg_day","dose"))) %>%
  rename_with(~ "c0_mean",   .cols = any_of(c("c0_mean","c0","trough_mean"))) %>%
  rename_with(~ "c0_sd",     .cols = any_of(c("c0_sd","sd_c0","trough_sd"))) %>%
  rename_with(~ "c0_d_mean", .cols = any_of(c("c0_d_mean","c0d_mean","c0_d"))) %>%
  rename_with(~ "c0_d_sd",   .cols = any_of(c("c0_d_sd","c0d_sd","sd_c0_d"))) %>%
  rename_with(~ "nos_quality_score", .cols = any_of(c(
    "nos_quality_score","nos_score","nos"))) %>%
  rename_with(~ "study_id", .cols = any_of(c("study_id","studyid","id"))) %>%
  rename_with(~ "tacrolimus_timepoint", .cols = any_of(c(
    "tacrolimus_timepoint","timepoint","time_point")))

# ---- 4. VALIDATE ------------------------------------------------------------
required_cols <- c(
  "study_id","first_author","year","gene","snp_rsid","genotype_groups",
  "sample_size","organ_type","population_adult_pediatric",
  "c0_mean","c0_sd","c0_d_mean","c0_d_sd","nos_quality_score","notes",
  "data_type_mean_sd_or_median_iqr"
)
missing_cols <- setdiff(required_cols, names(extraction))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "),
       "\nActual: ", paste(names(extraction), collapse = ", "))
}

# ---- 5. CLEAN + FLAGS -------------------------------------------------------
num_cols <- c("tacrolimus_dose_mg_kg_day","c0_mean","c0_sd",
              "c0_d_mean","c0_d_sd","sample_size","nos_quality_score")

extraction <- extraction %>%
  mutate(
    across(any_of(num_cols), ~ suppressWarnings(as.numeric(.x))),
    year  = suppressWarnings(as.integer(year)),
    notes = coalesce(as.character(notes), ""),   # blank notes must not drop rows
    genotype_groups = str_squish(as.character(genotype_groups)),
    is_combined  = str_detect(str_to_lower(notes), "combined|multi-gene|cluster|novelty"),
    is_duplicate = str_detect(str_to_lower(notes), "duplicate|overlap|supersed"),
    is_pediatric = str_detect(str_to_lower(coalesce(population_adult_pediatric, "")),
                              "pediatric|paediatric"),
    is_liver     = str_detect(str_to_lower(coalesce(organ_type, "")), "liver"),
    is_kidney    = str_detect(str_to_lower(coalesce(organ_type, "")), "kidney"),
    is_heart     = str_detect(str_to_lower(coalesce(organ_type, "")), "heart"),
    is_lung      = str_detect(str_to_lower(coalesce(organ_type, "")), "lung"),
    organ_clean  = case_when(
      is_kidney ~ "Kidney",
      is_liver  ~ "Liver",
      is_heart  ~ "Heart",
      is_lung   ~ "Lung",
      TRUE      ~ "Mixed"
    ),
    time_window = case_when(
      str_detect(str_to_lower(coalesce(tacrolimus_timepoint, "")),
                 "day|1\\s*mo|1\\s*month|early|discharge|post-?op") ~ "Early (<3 mo)",
      str_detect(str_to_lower(coalesce(tacrolimus_timepoint, "")),
                 "6\\s*mo|6\\s*month|12\\s*mo|year|stable|maintenance") ~ "Stable (>=6 mo)",
      TRUE ~ "Unspecified"
    ),
    is_median = str_detect(
      str_to_lower(coalesce(as.character(data_type_mean_sd_or_median_iqr), "")),
      "median")
  )

# ---- 6. GENE / SNP LABELS ---------------------------------------------------
extraction <- extraction %>%
  mutate(
    gene_clean = case_when(
      str_detect(gene, regex("CYP3A5", ignore_case = TRUE)) ~ "CYP3A5",
      str_detect(gene, regex("CYP3A4", ignore_case = TRUE)) ~ "CYP3A4",
      str_detect(gene, regex("ABCB1|MDR", ignore_case = TRUE)) ~ "ABCB1",
      str_detect(gene, regex("combined", ignore_case = TRUE)) ~ "Combined",
      TRUE ~ gene
    )
  )

# FIX 1: rows that name more than one target SNP are set aside as "MULTI"
targets <- c("rs776746","rs2242480","rs35599367","rs1045642","rs2032582","rs1128503")
hits  <- sapply(targets, function(t) str_detect(coalesce(as.character(extraction$snp_rsid), ""), t))
n_hit <- rowSums(hits)
extraction$snp_clean <- case_when(
  n_hit == 1 ~ targets[max.col(hits * 1, ties.method = "first")],
  n_hit > 1  ~ "MULTI",
  TRUE       ~ as.character(extraction$snp_rsid)
)
cat("Rows naming >1 target SNP (set aside as MULTI):", sum(n_hit > 1), "\n")
print(table(extraction$snp_clean)[order(-table(extraction$snp_clean))][1:10])

# ---- 7. POOLING DATASET -----------------------------------------------------
# v5: exclusions are decided at STUDY level from the Study ID (not by words in
# the Notes column, which also matched things like "LD overlap" or "not an overlap").
extraction <- extraction %>%
  mutate(excl_reason = case_when(
    str_detect(study_id, "DUPLICATE-RISK") ~ "Study ID flag: DUPLICATE-RISK",
    drop_overlap_risk & str_detect(study_id, "OVERLAP-RISK") ~
      "Study ID flag: OVERLAP-RISK (sensitivity run)",
    study_id %in% manual_exclude ~ "manual exclusion",
    TRUE ~ NA_character_))

excl_studies <- extraction %>% filter(!is.na(excl_reason)) %>%
  count(study_id, excl_reason, name = "rows")
write.csv(excl_studies, file.path(out_dir, "studies_excluded.csv"), row.names = FALSE)
cat("Rows excluded:", sum(!is.na(extraction$excl_reason)),
    "| studies:", n_distinct(excl_studies$study_id),
    " (list in the outputs folder: studies_excluded.csv)\n")
cat("drop_overlap_risk =", drop_overlap_risk, "\n")

pooling_data <- extraction %>%
  filter(is.na(excl_reason)) %>%
  mutate(
    effect_metric = if_else(!is.na(c0_d_mean), "C0/D", "C0"),
    effect_mean   = coalesce(c0_d_mean, c0_mean),
    effect_sd     = coalesce(c0_d_sd,   c0_sd),
    n             = coalesce(sample_size, 1)
  ) %>%
  filter(!is.na(effect_mean), !is.na(effect_sd),
         effect_mean > 0, effect_sd > 0, n > 0) %>%
  mutate(
    yi  = log(effect_mean),
    sei = effect_sd / effect_mean,
    vi  = sei^2
  ) %>%
  filter(is.finite(yi), is.finite(vi), vi > 0)

cat("Rows available for pooling:", nrow(pooling_data),
    "| studies:", n_distinct(pooling_data$study_id), "\n")

# ---- 8. GENOTYPE CLASSIFICATION ---------------------------------------------
pooling_data <- pooling_data %>%
  mutate(
    gl = str_to_lower(coalesce(genotype_groups, "")),
    genotype_class = case_when(
      # CYP3A5: group labels (expresser/non-expresser) are tested BEFORE the
      # single genotypes, so "*1/*1+*1/*3 (expressers)" is no longer read as *1/*1
      snp_clean == "rs776746" & str_detect(gl, "non-?\\s*express") ~ "CYP3A5 non-expresser",
      snp_clean == "rs776746" & (str_detect(gl, "express") |
          str_detect(gl, "\\*1\\s*carrier") |
          str_detect(gl, "\\*1/\\*1\\s*(\\+|or|and|,|/)\\s*\\*1/\\*3") |
          str_detect(gl, "\\*1/\\*3\\s*(\\+|or|and|,|/)\\s*\\*1/\\*1")) ~ "CYP3A5 expresser",
      snp_clean == "rs776746" & str_detect(gl, "\\*1/\\*1") ~ "CYP3A5 *1/*1",
      snp_clean == "rs776746" & str_detect(gl, "\\*1/\\*3") ~ "CYP3A5 *1/*3",
      snp_clean == "rs776746" & str_detect(gl, "\\*3/\\*3") ~ "CYP3A5 *3/*3",
      snp_clean == "rs35599367" & str_detect(genotype_groups, "\\*1/\\*1|GG|CC") ~ "CYP3A4*22 WT",
      snp_clean == "rs35599367" & str_detect(genotype_groups, "\\*1/\\*22|\\*22|GA|CT") ~ "CYP3A4*22 carrier",
      snp_clean == "rs2242480" & str_detect(genotype_groups, "\\*1/\\*1|CC|GG") ~ "CYP3A4*1G WT",
      snp_clean == "rs2242480" & str_detect(genotype_groups, "\\*1/\\*1G|\\*1G|CT|TT|GA|AA") ~ "CYP3A4*1G carrier",
      snp_clean == "rs1045642" & str_detect(genotype_groups, "CC") ~ "ABCB1 3435 CC",
      snp_clean == "rs1045642" & str_detect(genotype_groups, "CT") ~ "ABCB1 3435 CT",
      snp_clean == "rs1045642" & str_detect(genotype_groups, "TT") ~ "ABCB1 3435 TT",
      snp_clean == "rs2032582" & str_detect(genotype_groups, "GG") ~ "ABCB1 2677 GG",
      snp_clean == "rs2032582" & str_detect(genotype_groups, "GT|GA") ~ "ABCB1 2677 GT/GA",
      snp_clean == "rs2032582" & str_detect(genotype_groups, "TT|TA|AA") ~ "ABCB1 2677 TT/TA/AA",
      snp_clean == "rs1128503" & str_detect(genotype_groups, "CC") ~ "ABCB1 1236 CC",
      snp_clean == "rs1128503" & str_detect(genotype_groups, "CT") ~ "ABCB1 1236 CT",
      snp_clean == "rs1128503" & str_detect(genotype_groups, "TT") ~ "ABCB1 1236 TT",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(genotype_class)) %>%
  mutate(row_id = row_number())
cat("Rows after genotype classification:", nrow(pooling_data), "\n")

# FIX 2a: drop rows whose genotype label involves a donor / graft
n_before <- nrow(pooling_data)
pooling_data <- pooling_data %>%
  filter(!str_detect(gl, "donor|graft"))
cat("Rows dropped (donor/graft genotype labels):", n_before - nrow(pooling_data), "\n")

# FIX 2b: ABCB1 labels that are grouped/ambiguous ("CT+TT", "carrier", "vs") are
# set aside and listed for your review, not guessed
amb <- pooling_data %>%
  filter(snp_clean %in% c("rs1045642","rs2032582","rs1128503"),
         str_detect(gl, "\\+|carrier|\\bor\\b|\\bvs\\b|&"))
write.csv(amb %>% select(study_id, snp_clean, genotype_groups, genotype_class),
          file.path(out_dir, "abcb1_ambiguous_labels_set_aside.csv"), row.names = FALSE)
cat("ABCB1 rows with grouped/ambiguous labels set aside:", nrow(amb),
    " (list in the outputs folder: abcb1_ambiguous_labels_set_aside.csv)\n")
pooling_data <- pooling_data %>% filter(!(row_id %in% amb$row_id))

cat("Median-type rows still in pooling:", sum(pooling_data$is_median),
    "of", nrow(pooling_data), "\n")
cat("Final pooling rows:", nrow(pooling_data),
    "| studies:", n_distinct(pooling_data$study_id), "\n")
print(table(pooling_data$snp_clean))

# ---- 9. HELPERS -------------------------------------------------------------
fc_store <- list()

fc_table <- function(fit) {
  b  <- as.numeric(coef(fit))
  lb <- as.numeric(fit$ci.lb)
  ub <- as.numeric(fit$ci.ub)
  pv <- as.numeric(fit$pval)
  nm <- sub("^grp", "", names(coef(fit)))
  data.frame(
    contrast    = nm[-1],
    fold_change = round(exp(b[-1]), 3),
    ci_lower    = round(exp(lb[-1]), 3),
    ci_upper    = round(exp(ub[-1]), 3),
    p_value     = signif(pv[-1], 3),
    inverse_fc  = round(exp(-b[-1]), 3),   # reference / group
    row.names   = NULL
  )
}

build_grp <- function(dat, snp, mapping, levels_vec) {
  dat %>%
    filter(snp_clean == snp) %>%
    mutate(grp = factor(unname(mapping[genotype_class]), levels = levels_vec)) %>%
    filter(!is.na(grp))
}

fit_model <- function(dat, label, test = "z", verbose = TRUE) {
  dat <- droplevels(dat)
  if (verbose) {
    cat("\n-- ", label, " --\n", sep = "")
    cat("Rows:", nrow(dat), "| studies:", n_distinct(dat$study_id),
        "| levels:", nlevels(dat$grp), "\n")
    print(table(dat$grp))
  }
  if (!("Ref" %in% levels(dat$grp)) || nlevels(dat$grp) < 2 ||
      nrow(dat) < 3 || n_distinct(dat$study_id) < 2) {
    if (verbose) cat("   -> NOT FITTED (needs a Ref level, >=1 other level, >=3 rows, >=2 studies)\n")
    return(NULL)
  }
  fit <- tryCatch(
    rma.mv(yi, vi, mods = ~ grp, random = ~ 1 | study_id,
           data = dat, method = "REML", test = test),
    error = function(e) { cat("   -> error:", conditionMessage(e), "\n"); NULL })
  if (is.null(fit)) return(NULL)
  tab <- fc_table(fit)
  # v5: cluster-robust CIs (rows from the same study are not independent)
  fit_r <- tryCatch(robust(fit, cluster = dat$study_id), error = function(e) NULL)
  if (!is.null(fit_r)) {
    tab_r <- fc_table(fit_r)
    tab$ci_lower_robust <- tab_r$ci_lower
    tab$ci_upper_robust <- tab_r$ci_upper
    tab$p_value_robust  <- tab_r$p_value
  }
  if (verbose) {
    print(summary(fit))
    cat("\nFold-changes vs reference:\n"); print(tab)
    cat("Reference mean log-C0/D:", round(coef(fit)[1], 4),
        "-> C0/D =", round(exp(coef(fit)[1]), 2), "\n")
  }
  fc_store[[label]] <<- cbind(model = label, tab,
                              k_rows = nrow(dat), k_studies = n_distinct(dat$study_id))
  fit
}

cat("\n=============================================\n")
cat("  MAIN MODELS (contrast-coded, Ref = first level)\n")
cat("=============================================\n")

# ---- 10. CYP3A5 rs776746 (Ref = *3/*3 or non-expresser) ----------------------
cyp3a5_map <- c("CYP3A5 *3/*3" = "Ref", "CYP3A5 non-expresser" = "Ref",
                "CYP3A5 *1/*1" = "CYP3A5 *1/*1", "CYP3A5 *1/*3" = "CYP3A5 *1/*3",
                "CYP3A5 expresser" = "CYP3A5 expresser")
cyp3a5_final <- build_grp(pooling_data, "rs776746", cyp3a5_map,
                          c("Ref","CYP3A5 *1/*1","CYP3A5 *1/*3","CYP3A5 expresser"))
ma_cyp3a5_final <- fit_model(cyp3a5_final, "CYP3A5 rs776746", test = "t")

# ---- 11. CYP3A4 *22 and *1G ---------------------------------------------------
cyp3a4_22_final <- build_grp(pooling_data, "rs35599367",
  c("CYP3A4*22 WT" = "Ref", "CYP3A4*22 carrier" = "CYP3A4*22 carrier"),
  c("Ref","CYP3A4*22 carrier"))
ma_cyp3a4_22_final <- fit_model(cyp3a4_22_final, "CYP3A4*22 rs35599367")

cyp3a4_1g_final <- build_grp(pooling_data, "rs2242480",
  c("CYP3A4*1G WT" = "Ref", "CYP3A4*1G carrier" = "CYP3A4*1G carrier"),
  c("Ref","CYP3A4*1G carrier"))
ma_cyp3a4_1g_final <- fit_model(cyp3a4_1g_final, "CYP3A4*1G rs2242480")

# ---- 12. ABCB1 ----------------------------------------------------------------
abcb1_3435 <- build_grp(pooling_data, "rs1045642",
  c("ABCB1 3435 CC" = "Ref", "ABCB1 3435 CT" = "ABCB1 3435 CT", "ABCB1 3435 TT" = "ABCB1 3435 TT"),
  c("Ref","ABCB1 3435 CT","ABCB1 3435 TT"))
ma_abcb1_3435 <- fit_model(abcb1_3435, "ABCB1 rs1045642 (3435C>T)")

abcb1_2677 <- build_grp(pooling_data, "rs2032582",
  c("ABCB1 2677 GG" = "Ref", "ABCB1 2677 GT/GA" = "ABCB1 2677 GT/GA",
    "ABCB1 2677 TT/TA/AA" = "ABCB1 2677 TT/TA/AA"),
  c("Ref","ABCB1 2677 GT/GA","ABCB1 2677 TT/TA/AA"))
ma_abcb1_2677 <- fit_model(abcb1_2677, "ABCB1 rs2032582 (2677G>T/A)")

abcb1_1236 <- build_grp(pooling_data, "rs1128503",
  c("ABCB1 1236 CC" = "Ref", "ABCB1 1236 CT" = "ABCB1 1236 CT", "ABCB1 1236 TT" = "ABCB1 1236 TT"),
  c("Ref","ABCB1 1236 CT","ABCB1 1236 TT"))
ma_abcb1_1236 <- fit_model(abcb1_1236, "ABCB1 rs1128503 (1236C>T)")

# ---- 13. COMBINED MULTI-GENE (PROVISIONAL - to be rebuilt together) -----------
cat("\n=============================================\n")
cat("  COMBINED MODEL - PROVISIONAL, DO NOT CITE YET\n")
cat("=============================================\n")
combined_map <- c("CYP3A5 *3/*3" = "Ref", "CYP3A5 *1/*1" = "CYP3A5 *1/*1",
                  "CYP3A5 *1/*3" = "CYP3A5 *1/*3", "CYP3A5 expresser" = "CYP3A5 expresser",
                  "CYP3A5 non-expresser" = "CYP3A5 non-expresser", "CYP3A4*1G WT" = "CYP3A4*1G WT")
combined_final <- pooling_data %>%
  filter(gene_clean == "Combined" | is_combined) %>%
  mutate(grp = factor(unname(combined_map[genotype_class]),
                      levels = c("Ref","CYP3A5 *1/*1","CYP3A5 *1/*3","CYP3A5 expresser",
                                 "CYP3A5 non-expresser","CYP3A4*1G WT"))) %>%
  filter(!is.na(grp))
ma_combined_final <- fit_model(combined_final, "Combined multi-gene (provisional)")

# (v5: the single-gene vs multi-gene moderator test was removed. It compared
#  absolute C0/D levels between different rows, which is not a valid comparison.)

# ---- 14. STRATIFIED CYP3A5 ESTIMATES (organ, age, time window) --------------
cat("\n=============================================\n")
cat("  CYP3A5 STRATIFIED BY SUBGROUP\n")
cat("=============================================\n")
cyp3a5_final <- cyp3a5_final %>%
  mutate(age_group = if_else(is_pediatric, "Pediatric", "Adult"))

strat_fit <- function(dat, var) {
  for (lv in sort(unique(as.character(dat[[var]])))) {
    sub <- dat[dat[[var]] == lv & !is.na(dat[[var]]), ]
    fit <- fit_model(sub, paste0("CYP3A5 | ", var, " = ", lv), test = "t", verbose = FALSE)
    if (is.null(fit)) {
      cat(sprintf("%-40s : not fitted (rows %d, studies %d)\n",
                  paste0(var, " = ", lv), nrow(sub), n_distinct(sub$study_id)))
    } else {
      cat("\n", var, "=", lv, "(rows", nrow(sub), ", studies", n_distinct(sub$study_id), ")\n")
      print(fc_table(fit))
    }
  }
}
strat_fit(cyp3a5_final, "organ_clean")
strat_fit(cyp3a5_final, "age_group")
strat_fit(cyp3a5_final, "time_window")

# ---- 15. NOS META-REGRESSION (CYP3A5) -----------------------------------------
nos_dat <- cyp3a5_final %>% filter(!is.na(nos_quality_score)) %>% droplevels()
cat("\nNOS meta-regression rows:", nrow(nos_dat), "| studies:", n_distinct(nos_dat$study_id), "\n")
if (nrow(nos_dat) >= 10 && nlevels(nos_dat$grp) >= 2) {
  ma_nos <- tryCatch(
    rma.mv(yi, vi, mods = ~ grp + nos_quality_score, random = ~ 1 | study_id,
           data = nos_dat, method = "REML"),
    error = function(e) { cat("NOS error:", conditionMessage(e), "\n"); NULL })
  if (!is.null(ma_nos)) print(summary(ma_nos))
}

# ---- 16. LEAVE-ONE-OUT (CYP3A5) - this one is slow, let it run -----------------
if (!is.null(ma_cyp3a5_final)) {
  cat("\nLeave-one-out over", n_distinct(cyp3a5_final$study_id), "studies...\n")
  loo_results <- purrr::map_dfr(unique(cyp3a5_final$study_id), function(sid) {
    d <- droplevels(dplyr::filter(cyp3a5_final, study_id != sid))
    if (nlevels(d$grp) < 2) return(NULL)
    fit <- tryCatch(
      rma.mv(yi, vi, mods = ~ grp, random = ~ 1 | study_id, data = d, method = "REML"),
      error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    tb <- fc_table(fit); tb$omitted_study <- sid; tb
  })
  write.csv(loo_results, file.path(out_dir, "loo_cyp3a5.csv"), row.names = FALSE)
  loo_summary <- loo_results %>% group_by(contrast) %>%
    summarise(fc_min = min(fold_change), fc_max = max(fold_change),
              n_p_ge_0.05 = sum(p_value >= 0.05), n_fits = n(), .groups = "drop")
  print(loo_summary)
}

# ---- 17. HETEROGENEITY + PUBLICATION BIAS (CYP3A5) -----------------------------
cat("\n=============================================\n")
cat("  HETEROGENEITY & PUBLICATION BIAS (CYP3A5)\n")
cat("=============================================\n")
if (!is.null(ma_cyp3a5_final)) {
  fit <- ma_cyp3a5_final
  k    <- fit$k
  tau2 <- sum(as.numeric(fit$sigma2), na.rm = TRUE)
  vbar <- k / sum(1 / fit$vi, na.rm = TRUE)
  I2   <- tau2 / (tau2 + vbar) * 100
  cat(sprintf("k = %d | tau2 = %.4f | tau = %.4f | typical v = %.4f | I2 = %.1f %%\n",
              k, tau2, sqrt(tau2), vbar, I2))
  cat(sprintf("Q = %.2f (df = %d, p = %s)\n", fit$QE, k - fit$p, format(fit$QEp, digits = 3)))

  # v5: Egger test and funnel plot removed for this model. With between-study
  # variance (tau2) far larger than the within-row variance, every point sits at
  # the same height and the test does not measure publication bias.
}
# (Forest plot intentionally left out of this run: with hundreds of rows it
#  needs a per-study layout, we will design it once the numbers are settled.)

# ---- 18. RESULTS SUMMARY ------------------------------------------------------
cat("\n=============================================\n")
cat("  RESULTS SUMMARY (all fitted models)\n")
cat("=============================================\n")
if (length(fc_store) > 0) {
  all_fc <- dplyr::bind_rows(fc_store)
  print(as.data.frame(all_fc), row.names = FALSE)
  write.csv(all_fc, file.path(out_dir, "all_foldchanges.csv"), row.names = FALSE)
}
cat("\nNOT FITTED models:",
    paste(c(if (is.null(ma_cyp3a5_final)) "CYP3A5",
            if (is.null(ma_cyp3a4_22_final)) "CYP3A4*22",
            if (is.null(ma_cyp3a4_1g_final)) "CYP3A4*1G",
            if (is.null(ma_abcb1_3435)) "rs1045642",
            if (is.null(ma_abcb1_2677)) "rs2032582",
            if (is.null(ma_abcb1_1236)) "rs1128503"), collapse = ", "), "\n")
# ---- 19. PLOTS + STUDY-LEVEL TABLES (added after v6; base R only, no new packages) ----
# Writes PNGs to <out_dir>/plots. Study-level estimates are the within-study log
# fold-change of each non-reference genotype group vs the reference group (matched on
# timepoint + metric where possible). Per-study intervals are NOT drawn: per-genotype n
# was not extracted, so vi = (SD/mean)^2 is not a standard error.
tryCatch({
plot_dir <- file.path(out_dir, "plots"); dir.create(plot_dir, showWarnings = FALSE)
col_main <- "#2A78D6"; col_sens <- "#EB6834"; col_grey <- "#8A8983"
png_open <- function(name, w, h) png(file.path(plot_dir, name), width = w, height = h, units = "in", res = 300)
logx_axis <- function(ticks, lo, hi) { t <- ticks[ticks >= lo & ticks <= hi]; axis(1, at = log(t), labels = t, las = 1) }

# sibling run (the other drop_overlap_risk setting), if its csv exists, for the overlay
other_tag <- if (drop_overlap_risk) "main" else "sensitivity"
other_csv <- file.path(out_root, other_tag, "all_foldchanges.csv")
other_fc  <- if (file.exists(other_csv)) read.csv(other_csv, check.names = FALSE) else NULL
this_lab  <- if (drop_overlap_risk) "Sensitivity (OVERLAP-RISK removed)" else "Primary (all studies)"
other_lab <- if (drop_overlap_risk) "Primary (all studies)" else "Sensitivity (OVERLAP-RISK removed)"
all_fc <- as.data.frame(all_fc)

# ---- F1: summary forest ----
single <- all_fc[!grepl("\\||Combined", all_fc$model), ]
single$label <- paste0(sub(" rs.*$", "", single$model), ": ", single$contrast)
n <- nrow(single); ypos <- rev(seq_len(n))
png_open("F1_summary_forest.png", 9, 6)
par(mar = c(4.5, 15, 2.5, 8), las = 1)
plot(NA, xlim = log(c(0.2, 5)), ylim = c(0.5, n + 1.6), axes = FALSE, xlab = "", ylab = "")
abline(v = 0, col = "grey50"); abline(v = log(c(0.3, 0.5, 2, 3)), col = "grey90")
for (i in seq_len(n)) {
  r <- single[i, ]
  segments(log(r$ci_lower_robust), ypos[i] + 0.15, log(r$ci_upper_robust), ypos[i] + 0.15, col = col_main, lwd = 2)
  points(log(r$fold_change), ypos[i] + 0.15, pch = 16, col = col_main, cex = 1.1)
  if (!is.null(other_fc)) {
    o <- other_fc[other_fc$model == r$model & other_fc$contrast == r$contrast, ]
    if (nrow(o)) { segments(log(o$ci_lower_robust), ypos[i] - 0.15, log(o$ci_upper_robust), ypos[i] - 0.15, col = col_sens, lwd = 2)
                   points(log(o$fold_change), ypos[i] - 0.15, pch = 15, col = col_sens, cex = 1.0) }
  }
  mtext(sprintf("%.2f (%.2f-%.2f)", r$fold_change, r$ci_lower_robust, r$ci_upper_robust), side = 4, at = ypos[i], line = 0.5, cex = 0.7, adj = 0)
}
axis(2, at = ypos, labels = single$label, tick = FALSE, cex.axis = 0.75); logx_axis(c(0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5), 0.2, 5); box(bty = "l")
title(xlab = "Fold-change in C0/D vs reference genotype (log scale), cluster-robust 95% CI")
title(main = "Pooled single-gene effects on tacrolimus C0/D", adj = 0, cex.main = 1)
legend("topright", legend = c(this_lab, if (!is.null(other_fc)) other_lab), col = c(col_main, if (!is.null(other_fc)) col_sens),
       pch = c(16, if (!is.null(other_fc)) 15), lwd = 2, bty = "n", cex = 0.75)
invisible(dev.off())

# ---- F2: CYP3A5 subgroup forest ----
sub <- all_fc[grepl("^CYP3A5 \\|", all_fc$model), ]
sub$dim <- sub("^CYP3A5 \\| (\\w+) = .*$", "\\1", sub$model); sub$stratum <- sub("^.* = ", "", sub$model)
base3 <- all_fc[all_fc$model == "CYP3A5 rs776746", ]; base3$dim <- "all"; base3$stratum <- "All CYP3A5 studies"
sub <- rbind(base3, sub)
ord <- c("All CYP3A5 studies","Kidney","Liver","Heart","Lung","Adult","Pediatric","Early (<3 mo)","Stable (>=6 mo)","Unspecified")
strata <- unique(sub[, c("dim", "stratum")]); strata <- strata[match(ord, strata$stratum), ]; strata <- strata[!is.na(strata$stratum), ]; ys <- rev(seq_len(nrow(strata)))
png_open("F2_CYP3A5_subgroups.png", 11, 5.2)
par(mfrow = c(1, 3), mar = c(4.5, 1, 3, 1), oma = c(2, 11, 0, 3), las = 1)
for (ci_ in 1:3) {
  con <- c("CYP3A5 expresser", "CYP3A5 *1/*3", "CYP3A5 *1/*1")[ci_]
  plot(NA, xlim = log(c(0.1, 10)), ylim = c(0.5, nrow(strata) + 0.5), axes = FALSE, xlab = "Fold-change in C0/D (log)", ylab = "")
  abline(v = 0, col = "grey50"); abline(v = log(c(0.3, 3)), col = "grey90")
  for (i in seq_len(nrow(strata))) {
    r <- sub[sub$dim == strata$dim[i] & sub$stratum == strata$stratum[i] & sub$contrast == con, ]
    if (!nrow(r)) { text(log(0.2), ys[i], "not fitted", cex = 0.7, col = "grey40"); next }
    cl <- if (r$k_studies < 6) col_grey else col_main
    segments(log(max(r$ci_lower_robust, 0.1)), ys[i], log(min(r$ci_upper_robust, 10)), ys[i], col = cl, lwd = 2)
    points(log(r$fold_change), ys[i], pch = if (strata$dim[i] == "all") 18 else 16, col = cl, cex = 1.2)
    if (ci_ == 3) mtext(paste0("k=", r$k_studies), side = 4, at = ys[i], line = 0.3, cex = 0.7, las = 1)
  }
  logx_axis(c(0.1, 0.3, 1, 3, 10), 0.1, 10); box(bty = "l")
  title(main = c("Expresser vs non-expresser", "*1/*3 vs *3/*3", "*1/*1 vs *3/*3")[ci_], adj = 0, cex.main = 0.95)
  if (ci_ == 1) mtext(strata$stratum, side = 2, at = ys, line = 0.5, cex = 0.7, las = 1, adj = 1, outer = FALSE)
}
mtext("Grey = fewer than 6 studies (indicative only). k = studies in the expresser model.", side = 1, outer = TRUE, line = 0.5, cex = 0.65, adj = 0.01, col = "grey30")
invisible(dev.off())

# ---- study-level contrasts (CYP3A5, ABCB1, CYP3A4) ----
mods <- list("CYP3A5 rs776746" = cyp3a5_final, "CYP3A4*22" = cyp3a4_22_final, "CYP3A4*1G" = cyp3a4_1g_final,
             "ABCB1 3435" = abcb1_3435, "ABCB1 2677" = abcb1_2677, "ABCB1 1236" = abcb1_1236)
sl <- do.call(rbind, lapply(names(mods), function(nm) {
  d <- mods[[nm]]
  do.call(rbind, lapply(unique(d$study_id), function(sid) {
    x <- d[d$study_id == sid, ]; r <- x[x$grp == "Ref", ]
    if (!nrow(r)) return(NULL)
    do.call(rbind, lapply(setdiff(unique(as.character(x$grp)), "Ref"), function(g) {
      o <- x[x$grp == g, ]
      key <- merge(data.frame(tp = o$tacrolimus_timepoint, em = o$effect_metric, yo = o$yi),
                   data.frame(tp = r$tacrolimus_timepoint, em = r$effect_metric, yr = r$yi), by = c("tp", "em"))
      if (nrow(key)) { est <- mean(key$yo - key$yr); how <- "timepoint-matched"; np <- nrow(key) }
      else {
        mo <- tapply(o$yi, o$effect_metric, mean); mr <- tapply(r$yi, r$effect_metric, mean)
        cm <- intersect(names(mo), names(mr)); if (!length(cm)) return(NULL)
        est <- mean(mo[cm] - mr[cm]); how <- "study-mean"; np <- length(cm)
      }
      ns <- suppressWarnings(max(x$sample_size, na.rm = TRUE)); if (!is.finite(ns)) ns <- NA
      data.frame(model = nm, contrast = g, study_id = sid, log_fc = est, fc = exp(est), pairing = how,
                 n_pairs = np, n_study = ns, organ = x$organ_clean[1], stringsAsFactors = FALSE)
    }))
  }))
}))
write.csv(sl, file.path(out_dir, "study_level_contrasts.csv"), row.names = FALSE)

# ---- F3: per-study CYP3A5 forest (point estimates, size = sample size) ----
c5 <- sl[sl$model == "CYP3A5 rs776746", ]
car <- do.call(rbind, lapply(split(c5, c5$study_id), function(z)
  data.frame(study_id = z$study_id[1], log_fc = mean(z$log_fc), n_study = suppressWarnings(max(z$n_study, na.rm = TRUE)), organ = z$organ[1])))
car$n_study[!is.finite(car$n_study)] <- NA; car <- car[order(car$log_fc), ]; k <- nrow(car)
pooled <- all_fc[all_fc$model == "CYP3A5 rs776746" & all_fc$contrast == "CYP3A5 expresser", ]
png_open("F3_CYP3A5_study_level.png", 7.4, 0.17 * k + 1.8)
par(mar = c(4.5, 13, 2.5, 1), las = 1)
plot(NA, xlim = log(c(0.03, 12)), ylim = c(0, k + 1), axes = FALSE, xlab = "Study-level fold-change in C0/D, CYP3A5 expresser vs *3/*3 (log)", ylab = "")
rect(log(pooled$ci_lower_robust), 0, log(pooled$ci_upper_robust), k + 1, col = adjustcolor(col_main, 0.12), border = NA)
abline(v = log(pooled$fold_change), col = col_main, lwd = 1.5); abline(v = 0, col = "grey50")
cx <- ifelse(is.na(car$n_study), 0.7, 0.6 + sqrt(pmin(car$n_study, 300)) / 8)
points(car$log_fc, seq_len(k), pch = ifelse(is.na(car$n_study), 1, 16), cex = cx, col = ifelse(car$organ == "Liver", col_sens, col_main))
axis(2, at = seq_len(k), labels = gsub("_(CYP3A5|CYP3A4|CYP3A|Kidney|Liver|Heart|Lung)", "", car$study_id), tick = FALSE, cex.axis = 0.38)
logx_axis(c(0.03, 0.1, 0.2, 0.3, 0.5, 1, 2, 5, 10), 0.03, 12); box(bty = "l")
title(main = sprintf("CYP3A5 effect in each study (k = %d)", k), adj = 0, cex.main = 0.95)
legend("bottomright", c(sprintf("Pooled %.2f (%.2f-%.2f)", pooled$fold_change, pooled$ci_lower_robust, pooled$ci_upper_robust), "Liver (donor graft)", "n not extracted (open)"),
       pch = c(15, 16, 1), col = c(adjustcolor(col_main, 0.4), col_sens, "grey30"), bty = "n", cex = 0.6)
invisible(dev.off())

# ---- F4: funnel plots (precision proxy = study sample size) + small-study tests ----
tt <- function(z) { z <- z[!is.na(z$n_study) & z$n_study > 0, ]
  if (nrow(z) < 8) return(NULL)
  f <- summary(lm(log_fc ~ I(1 / sqrt(n_study)), data = z))$coefficients; ct <- suppressWarnings(cor.test(z$log_fc, z$n_study, method = "spearman", exact = FALSE))
  data.frame(k_with_n = nrow(z), slope = f[2, 1], slope_p = f[2, 4], rho = unname(ct$estimate), rho_p = ct$p.value) }
t3435 <- sl[sl$model == "ABCB1 3435" & sl$contrast == "ABCB1 3435 TT", ]
fun <- list("CYP3A5 expresser vs non-expresser" = car[, c("log_fc", "n_study")], "ABCB1 3435 TT vs CC" = t3435[, c("log_fc", "n_study")])
pf <- c(log(pooled$fold_change), log(all_fc$fold_change[all_fc$contrast == "ABCB1 3435 TT" & all_fc$model == "ABCB1 rs1045642 (3435C>T)"]))
png_open("F4_funnel_plots.png", 10, 4.6); par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1), las = 1); tests <- list()
for (i in 1:2) {
  z <- fun[[i]]; z <- z[!is.na(z$n_study) & z$n_study > 0, ]
  plot(z$log_fc, log(z$n_study), ylim = c(min(log(z$n_study)) - 1.3, max(log(z$n_study)) + 0.1), pch = 16, col = col_main, axes = FALSE, xlab = "Study-level fold-change in C0/D (log)", ylab = "Study sample size (log)")
  abline(v = 0, col = "grey50"); abline(v = pf[i], lty = 2)
  logx_axis(c(0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5), exp(min(z$log_fc)) * 0.8, exp(max(z$log_fc)) * 1.25)
  axis(2, at = log(c(2, 5, 10, 20, 50, 100, 200)), labels = c(2, 5, 10, 20, 50, 100, 200)); box(bty = "l")
  title(main = names(fun)[i], adj = 0, cex.main = 0.9)
  r_ <- tt(fun[[i]]); if (!is.null(r_)) { tests[[i]] <- cbind(analysis = names(fun)[i], r_)
    legend("bottomleft", c(sprintf("%d of %d studies report n", r_$k_with_n, nrow(fun[[i]])), sprintf("slope on 1/sqrt(n): p = %.2f", r_$slope_p), sprintf("Spearman rho = %.2f, p = %.2f", r_$rho, r_$rho_p)), bty = "n", cex = 0.75) }
}
invisible(dev.off())
write.csv(do.call(rbind, tests), file.path(out_dir, "small_study_tests.csv"), row.names = FALSE)

# ---- F5: leave-one-out, most influential studies ----
if (exists("loo_results") && nrow(loo_results)) {
  png_open("F5_leave_one_out.png", 13.5, 4.8); par(mfrow = c(1, 3), mar = c(4.5, 10, 4, 2), las = 1)
  for (con in c("CYP3A5 expresser", "CYP3A5 *1/*3", "CYP3A5 *1/*1")) {
    l <- loo_results[loo_results$contrast == con, ]; full <- all_fc$fold_change[all_fc$model == "CYP3A5 rs776746" & all_fc$contrast == con]
    top <- head(l[order(-abs(log(l$fold_change / full))), ], 8); top <- top[nrow(top):1, ]
    rg <- range(l$fold_change); pad <- diff(rg) * 0.3
    plot(NA, xlim = c(rg[1] - pad, rg[2] + pad), ylim = c(0.5, nrow(top) + 0.5), axes = FALSE, xlab = "Pooled fold-change with this study omitted", ylab = "")
    rect(rg[1], 0, rg[2], nrow(top) + 1, col = adjustcolor(col_main, 0.10), border = NA); abline(v = full, lty = 2)
    segments(full, seq_len(nrow(top)), top$fold_change, seq_len(nrow(top)), col = col_main, lwd = 2); points(top$fold_change, seq_len(nrow(top)), pch = 16, col = col_main)
    axis(2, at = seq_len(nrow(top)), labels = gsub("_(CYP3A5|CYP3A4|CYP3A|Kidney|Liver|Heart|Lung)", "", top$omitted_study), tick = FALSE, cex.axis = 0.7); axis(1); box(bty = "l")
    title(main = sprintf("%s\nfull %.3f; range %.3f-%.3f", con, full, rg[1], rg[2]), adj = 0, cex.main = 0.85)
  }
  invisible(dev.off())
}

# ---- F6: ABCB1 genotype gradients ----
png_open("F6_ABCB1_gradients.png", 10, 3.8); par(mfrow = c(1, 3), mar = c(4, 4.5, 3, 1), las = 1)
for (m in list(list("ABCB1 rs1045642 (3435C>T)", c("CC", "CT", "TT"), c("ABCB1 3435 CT", "ABCB1 3435 TT")),
               list("ABCB1 rs2032582 (2677G>T/A)", c("GG", "GT/GA", "TT/TA/AA"), c("ABCB1 2677 GT/GA", "ABCB1 2677 TT/TA/AA")),
               list("ABCB1 rs1128503 (1236C>T)", c("CC", "CT", "TT"), c("ABCB1 1236 CT", "ABCB1 1236 TT")))) {
  r <- do.call(rbind, lapply(m[[3]], function(cn) all_fc[all_fc$model == m[[1]] & all_fc$contrast == cn, ]))
  plot(NA, xlim = c(0.6, 3.4), ylim = log(c(0.6, 2)), axes = FALSE, xlab = "", ylab = "Fold-change in C0/D (log)")
  abline(h = 0, col = "grey50"); lines(1:3, c(0, log(r$fold_change)), col = adjustcolor(col_main, 0.5))
  segments(2:3, log(r$ci_lower_robust), 2:3, log(r$ci_upper_robust), col = col_main, lwd = 2); points(1:3, c(0, log(r$fold_change)), pch = 16, col = col_main, cex = 1.2)
  axis(1, at = 1:3, labels = m[[2]]); axis(2, at = log(c(0.7, 1, 1.5, 2)), labels = c(0.7, 1, 1.5, 2)); box(bty = "l")
  title(main = sprintf("%s (k = %d)", sub(" rs\\d+ ", " ", m[[1]]), r$k_studies[1]), adj = 0, cex.main = 0.85)
}
invisible(dev.off())
# ---- F7: CYP3A5 gradient across studies (one dot per study, bar = pooled estimate) ----
cons <- c("CYP3A5 *1/*3", "CYP3A5 *1/*1", "CYP3A5 expresser")
cons_lab <- c("*1/*3 vs *3/*3", "*1/*1 vs *3/*3", "Expresser (grouped)\nvs non-expresser")
png_open("F7_CYP3A5_gradient_by_study.png", 6.8, 4.4)
par(mar = c(5.5, 4.8, 3, 3), las = 1)
plot(NA, xlim = c(0.5, 3.7), ylim = log(c(0.1, 3)), axes = FALSE, xlab = "", ylab = "")
abline(h = log(c(0.1, 0.2, 0.5, 2)), col = "grey90"); abline(h = 0, col = "grey50")
set.seed(7)
for (i in seq_along(cons)) {
  z <- sl[sl$model == "CYP3A5 rs776746" & sl$contrast == cons[i], ]
  if (!nrow(z)) next
  points(i + runif(nrow(z), -0.18, 0.18), log(z$fc), pch = 16, col = adjustcolor(col_main, 0.65), cex = 0.9)
  r <- all_fc[all_fc$model == "CYP3A5 rs776746" & all_fc$contrast == cons[i], ]
  if (nrow(r)) { segments(i - 0.3, log(r$fold_change), i + 0.3, log(r$fold_change), lwd = 3)
                 text(i + 0.34, log(r$fold_change), sprintf("%.2f", r$fold_change), adj = 0, cex = 0.8) }
  text(i, log(0.115), paste0("k = ", nrow(z)), cex = 0.75, col = "grey40")
}
axis(2, at = log(c(0.1, 0.2, 0.5, 1, 2)), labels = c(0.1, 0.2, 0.5, 1, 2))
axis(1, at = seq_along(cons), labels = cons_lab, tick = TRUE, padj = 0.6, cex.axis = 0.8)
title(ylab = "Study-level fold-change in C0/D (log)")
title(main = "CYP3A5 gradient across studies", adj = 0, cex.main = 1)
mtext("Dots = individual studies; black bar = pooled estimate.", side = 1, line = 4.3, adj = 0, cex = 0.65, col = "grey30")
box(bty = "l")
invisible(dev.off())

# ---- F8: evidence base per target SNP ----
ev <- all_fc[!grepl("\\||Combined", all_fc$model), c("model", "k_studies", "k_rows")]
ev <- ev[!duplicated(ev$model), ]
ev$lab <- c("CYP3A5 rs776746 (*3)", "CYP3A4 rs35599367 (*22)", "CYP3A4 rs2242480 (*1G)",
            "ABCB1 rs1045642 (3435)", "ABCB1 rs2032582 (2677)", "ABCB1 rs1128503 (1236)")[
  match(sub(" \\(.*$", "", ev$model),
        c("CYP3A5 rs776746", "CYP3A4*22 rs35599367", "CYP3A4*1G rs2242480",
          "ABCB1 rs1045642", "ABCB1 rs2032582", "ABCB1 rs1128503"))]
ev <- ev[!is.na(ev$lab), ]; ev <- ev[nrow(ev):1, ]
xmax <- max(150, max(ev$k_studies) * 1.4)
png_open("F8_evidence_base.png", 6.6, 3.2)
par(mar = c(4, 11, 3, 1), las = 1)
bp <- barplot(ev$k_studies, horiz = TRUE, names.arg = ev$lab, col = col_main, border = NA,
              xlim = c(0, xmax), xlab = "Studies contributing to the model", cex.names = 0.75, width = 0.55, space = 0.8)
text(ev$k_studies + xmax * 0.01, bp, sprintf("%d studies, %d rows", ev$k_studies, ev$k_rows), adj = 0, cex = 0.7)
title(main = "Evidence base per target SNP", adj = 0, cex.main = 1)
invisible(dev.off())
cat("\nPlots and study-level tables written to:", plot_dir, "\n")

}, error = function(e) cat("\nPLOT SECTION FAILED:", conditionMessage(e), "\n(models and CSVs above are unaffected; run sink() before retrying)\n"))


# ---- 20. ETHNICITY (added 2026-10-15) ------------------------------------------
# Harmonises the free-text Ethnicity field into ancestry groups with a keyword
# rule (no guessing: "Not stated/reported/specified" stays Unspecified), then
# (a) fits the CYP3A5 model within each group, (b) tests genotype x ancestry
# interaction, (c) repeats the ABCB1 3435 model by ancestry (exploratory).
tryCatch({
cat("\n=============================================\n")
cat("  ETHNICITY ANALYSIS\n")
cat("=============================================\n")
eth_group <- function(x) {
  s <- str_to_lower(coalesce(as.character(x), ""))
  dplyr::case_when(
    s == "" | str_detect(s, "^(not stated|not reported|not specified)")           ~ "Unspecified",
    str_detect(s, "mixed|multiethnic|black and white|black african|native american|/black") ~ "Mixed",
    str_detect(s, "^(caucasian|white|european|predominantly|mostly|\\d+(\\.\\d+)?% european)") ~ "European",
    str_detect(s, "malay|thai|vietnam|myanmar|burmese")                           ~ "Southeast Asian",
    str_detect(s, "indian|south asian")                                           ~ "South Asian",
    str_detect(s, "chinese|\\bhan\\b|japan|korea|taiwan|(?<!cauc)asian")           ~ "East Asian",
    str_detect(s, "egypt|iran|tunisia|turk|jordan|middle east|north africa|arab")  ~ "MENA",
    str_detect(s, "mexic|hispanic|brazil|argentin")                               ~ "Latin American",
    str_detect(s, "african")                                                      ~ "African",
    str_detect(s, "caucasian|european|white|polish|spanish|italian|hungarian|serbian|german|dutch|scandinavian|sicilian|scottish|french") ~ "European",
    TRUE ~ "Unspecified")
}
cy <- cyp3a5_final %>% mutate(eth = eth_group(ethnicity))
eth_map <- cy %>% distinct(ethnicity, eth) %>% arrange(eth, ethnicity)
write.csv(eth_map, file.path(out_dir, "ethnicity_classification_CYP3A5.csv"), row.names = FALSE)
cat("\nCYP3A5 rows / studies by ancestry group:\n")
print(cy %>% group_by(eth) %>% summarise(rows = n(), studies = n_distinct(study_id)) %>% as.data.frame())

eth_out <- list()
for (g in setdiff(sort(unique(cy$eth)), c("Unspecified", "Mixed"))) {
  d <- droplevels(filter(cy, eth == g))
  if (n_distinct(d$study_id) < 2) { cat("\n", g, ": <2 studies, not pooled\n"); next }
  f <- fit_model(d, paste("CYP3A5 | ancestry =", g), test = "t", verbose = FALSE)
  if (is.null(f)) { cat("\n", g, ": not fitted\n"); next }
  tb <- fc_store[[paste("CYP3A5 | ancestry =", g)]]
  tb$tau2 <- round(f$sigma2, 3); eth_out[[g]] <- tb
}
eth_tab <- dplyr::bind_rows(eth_out)
cat("\nCYP3A5 fold-changes by ancestry (robust 95% CI):\n"); print(eth_tab, row.names = FALSE)
write.csv(eth_tab, file.path(out_dir, "ethnicity_cyp3a5_stratified.csv"), row.names = FALSE)

# (b) does the CYP3A5 effect differ by ancestry?  The row-level model has no residual
# between-study variation in the genotype effect, so a Wald or likelihood-ratio test
# of genotype x ancestry terms comes out far too small (checked: p about 1e-19).
# The test below therefore works on one number per study: the mean within-study log
# fold-change for *1 carriers vs *3/*3 (same study-level effects as the funnel plots).
cat("\nStudy-level comparison of the CYP3A5 effect across ancestry groups:\n")
slc <- sl[sl$model == "CYP3A5 rs776746", ]
st <- slc %>% group_by(study_id) %>% summarise(log_fc = mean(log_fc), .groups = "drop") %>%
  left_join(cy %>% group_by(study_id) %>% summarise(eth = first(eth), .groups = "drop"), by = "study_id")
sg <- st %>% filter(!eth %in% c("Unspecified", "Mixed")) %>% group_by(eth) %>% filter(n() >= 3) %>% ungroup()
sg_tab <- sg %>% group_by(eth) %>% summarise(studies = n(), median_fc = round(exp(median(log_fc)), 3),
                                             geo_mean_fc = round(exp(mean(log_fc)), 3),
                                             q25_fc = round(exp(quantile(log_fc, .25)), 3),
                                             q75_fc = round(exp(quantile(log_fc, .75)), 3), .groups = "drop")
print(as.data.frame(sg_tab), row.names = FALSE)
write.csv(sg_tab, file.path(out_dir, "ethnicity_study_level_summary.csv"), row.names = FALSE)
an <- anova(lm(log_fc ~ eth, data = sg)); kw <- kruskal.test(log_fc ~ factor(eth), data = sg)
cat(sprintf("One-way ANOVA on study-level log fold-change: F(%d, %d) = %.2f, p = %.3f\n", an$Df[1], an$Df[2], an$`F value`[1], an$`Pr(>F)`[1]))
cat(sprintf("Kruskal-Wallis: chi2(%d) = %.2f, p = %.3f\n", kw$parameter, kw$statistic, kw$p.value))
write.csv(data.frame(studies = nrow(sg), groups = paste(sort(unique(sg$eth)), collapse = "; "),
                     anova_F = round(an$`F value`[1], 3), df1 = an$Df[1], df2 = an$Df[2], anova_p = signif(an$`Pr(>F)`[1], 3),
                     kruskal_chi2 = round(unname(kw$statistic), 3), kruskal_p = signif(kw$p.value, 3)),
          file.path(out_dir, "ethnicity_study_level_tests.csv"), row.names = FALSE)

# (c) ABCB1 3435 by ancestry (exploratory; not adjusted for CYP3A5)
ab <- abcb1_3435 %>% mutate(eth = eth_group(ethnicity))
cat("\nABCB1 3435 rows / studies by ancestry group:\n")
print(ab %>% group_by(eth) %>% summarise(rows = n(), studies = n_distinct(study_id)) %>% as.data.frame())
ab_out <- list()
for (g in setdiff(sort(unique(ab$eth)), c("Unspecified", "Mixed"))) {
  d <- droplevels(filter(ab, eth == g)); if (n_distinct(d$study_id) < 3) next
  f <- fit_model(d, paste("ABCB1 3435 | ancestry =", g), test = "t", verbose = FALSE)
  if (!is.null(f)) ab_out[[g]] <- fc_store[[paste("ABCB1 3435 | ancestry =", g)]]
}
ab_tab <- dplyr::bind_rows(ab_out); print(ab_tab, row.names = FALSE)
write.csv(ab_tab, file.path(out_dir, "ethnicity_abcb1_3435_stratified.csv"), row.names = FALSE)
}, error = function(e) cat("\nETHNICITY SECTION FAILED:", conditionMessage(e), "\n"))

cat("\nDONE.\n")
sink()
