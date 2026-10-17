# Generates test_data/example_workbook.xlsx: a SYNTHETIC dataset in the layout the pipeline expects. No real study is represented; author names are invented.
# Effect sizes are simulated around plausible values so the pipeline has something to find.
# Usage (repo root): Rscript test_data/make_example_data.R
suppressWarnings(library(writexl))
set.seed(42)

cols <- c("Study ID","Zotero Key","First Author","Year","Country","Ethnicity","Study Design",
          "Organ Type","Population (Adult/Pediatric)","Sample Size","Gene","SNP rsID",
          "Allele Nomenclature","Genotype Groups","Tacrolimus Timepoint",
          "Tacrolimus Dose (mg/kg/day)","C0 Mean","C0 SD","C0/D Mean","C0/D SD",
          "Data Type (Mean/SD or Median/IQR)","Acute Rejection","Nephrotoxicity",
          "Genotyping Method","Hardy-Weinberg Equilibrium","NOS Quality Score","Notes")

authors <- c("Alder","Berg","Castell","Dube","Eklund","Fonseca","Gupta","Haddad","Ito","Jansen",
             "Kowal","Lindqvist","Moreau","Nakamura","Okafor","Park","Quintana","Rossi","Sato","Tan",
             "Ueda","Varga","Wang","Xu","Yilmaz","Zhou","Abbas","Bianchi","Chen","Dimitrov",
             "Erdem","Fischer","Gomez","Hassan","Ivanov","Jung","Khan","Li","Martin","Novak",
             "Ortiz","Patel","Qian","Rahman","Silva","Tanaka","Usman","Vega","Wu","Yamada")
eth_pool <- list(
  "East Asian"  = c("Chinese (Han)", "Korean", "Japanese", "Taiwanese"),
  "European"    = c("Caucasian", "Polish", "Spanish", "Italian"),
  "South Asian" = c("Indian", "South Indian"),
  "MENA"        = c("Egyptian", "Iranian", "Tunisian"),
  "SE Asian"    = c("Malay", "Thai")
)
eth_w <- c(0.35, 0.30, 0.12, 0.13, 0.10)
countries <- c("China","Poland","India","Egypt","Thailand")
organ_w <- c(Kidney = 0.70, Liver = 0.12, Heart = 0.10, Lung = 0.08)

n_studies <- 50
studies <- data.frame(
  author = authors[1:n_studies],
  year   = sample(2008:2025, n_studies, replace = TRUE),
  grp    = sample(names(eth_pool), n_studies, replace = TRUE, prob = eth_w),
  organ  = sample(names(organ_w), n_studies, replace = TRUE, prob = organ_w),
  ped    = runif(n_studies) < 0.15,
  stringsAsFactors = FALSE)
studies$eth <- mapply(function(g) sample(eth_pool[[g]], 1), studies$grp)
studies$country <- countries[match(studies$grp, names(eth_pool))]
studies$base <- exp(rnorm(n_studies, log(95), 0.35))      # study baseline C0/D in non-expressers
studies$id <- sprintf("%s_%d_CYP3A5_%s", studies$author, studies$year, studies$organ)
studies$id[c(7, 19)]  <- paste0(studies$id[c(7, 19)],  "_OVERLAP-RISK")   # exercises the sensitivity run

fc <- c(`*1/*1` = 0.34, `*1/*3` = 0.45, `*3/*3` = 1)
rows <- list()
add <- function(s, gene, rs, allele, geno, tp, mu, n, note = "", dtype = "Mean/SD", use_c0 = FALSE) {
  sd_ <- mu * runif(1, 0.25, 0.45)
  r <- setNames(rep(NA, length(cols)), cols)
  r[c("Study ID","Zotero Key","First Author","Year","Country","Ethnicity","Study Design","Organ Type",
      "Population (Adult/Pediatric)","Sample Size","Gene","SNP rsID","Allele Nomenclature",
      "Genotype Groups","Tacrolimus Timepoint","Data Type (Mean/SD or Median/IQR)",
      "Genotyping Method","Hardy-Weinberg Equilibrium","NOS Quality Score","Notes")] <- c(
    s$id, paste0("SYN", sprintf("%05d", sample(1e5, 1))), s$author, s$year, s$country, s$eth,
    "Retrospective cohort", s$organ, if (s$ped) "Pediatric" else "Adult", n, gene, rs, allele, geno, tp,
    dtype, "TaqMan", "Yes", sample(c(5:8, NA), 1), note)
  if (use_c0) { r[["C0 Mean"]] <- round(mu * 0.08, 2); r[["C0 SD"]] <- round(sd_ * 0.08, 2) }
  else { r[["C0/D Mean"]] <- round(mu, 2); r[["C0/D SD"]] <- round(sd_, 2) }
  rows[[length(rows) + 1]] <<- r
}

tps <- c("Day 7 post-tx", "1 month post-tx", "6 mo post-tx", "12 mo post-tx")
for (i in seq_len(n_studies)) {
  s <- studies[i, ]
  tp_set <- if (runif(1) < 0.35) sample(tps, 2) else sample(tps, 1)
  # CYP3A5 (all studies): either 3 genotype rows or expresser/non-expresser rows
  use_c0 <- runif(1) < 0.1
  for (tp in tp_set) {
    tp_eff <- if (grepl("Day|1 month", tp)) 1.1 else 0.95
    liver_blur <- if (s$organ == "Liver") 0.6 else 1       # weaker apparent effect (donor graft)
    f <- function(g) { x <- fc[[g]]; exp(log(x) * liver_blur) * exp(rnorm(1, 0, 0.12)) }
    if (runif(1) < 0.2) {
      add(s, "CYP3A5", "rs776746", "*3", "Expressers (*1/*1 + *1/*3)", tp, s$base * tp_eff * 0.42 * exp(rnorm(1,0,.1)), sample(15:60, 1), use_c0 = use_c0)
      add(s, "CYP3A5", "rs776746", "*3", "Non-expressers (*3/*3)", tp, s$base * tp_eff, sample(30:150, 1), use_c0 = use_c0)
    } else {
      add(s, "CYP3A5", "rs776746", "*3", "*3/*3", tp, s$base * tp_eff, sample(30:150, 1), use_c0 = use_c0)
      add(s, "CYP3A5", "rs776746", "*3", "*1/*3", tp, s$base * tp_eff * f("*1/*3"), sample(15:80, 1), use_c0 = use_c0)
      if (runif(1) < 0.5) add(s, "CYP3A5", "rs776746", "*3", "*1/*1", tp, s$base * tp_eff * f("*1/*1"), sample(3:20, 1), use_c0 = use_c0)
    }
  }
  tp <- tp_set[1]
  # ABCB1 3435: mild TT effect in ~45% of studies
  if (runif(1) < 0.45) for (g in c("CC","CT","TT")) {
    m <- c(CC = 1, CT = 1.1, TT = 1.25)[[g]]
    add(s, "ABCB1", "rs1045642", "3435C>T", g, tp, s$base * m * exp(rnorm(1, 0, 0.1)), sample(15:90, 1))
  }
  # ABCB1 2677 and 1236: null
  if (runif(1) < 0.25) for (g in c("GG","GT","TT")) add(s, "ABCB1", "rs2032582", "2677G>T/A", g, tp, s$base * exp(rnorm(1, 0, 0.08)), sample(15:90, 1))
  if (runif(1) < 0.15) for (g in c("CC","CT","TT")) add(s, "ABCB1", "rs1128503", "1236C>T", g, tp, s$base * exp(rnorm(1, 0, 0.08)), sample(15:90, 1))
  # CYP3A4*22 (carrier effect ~1.5x) and *1G: only a handful of studies
  if (i %in% c(3, 11, 23, 31, 40, 47)) {
    add(s, "CYP3A4", "rs35599367", "*22", "*1/*1", tp, s$base * 0.9, sample(60:200, 1))
    add(s, "CYP3A4", "rs35599367", "*22", "*1/*22", tp, s$base * 0.9 * 1.5 * exp(rnorm(1, 0, .15)), sample(5:25, 1))
  }
  if (i %in% c(5, 14, 28, 36, 44)) {
    add(s, "CYP3A4", "rs2242480", "*1G", "CC", tp, s$base * 0.9, sample(40:120, 1))
    add(s, "CYP3A4", "rs2242480", "*1G", "CT", tp, s$base * 0.9 * 0.8 * exp(rnorm(1, 0, .15)), sample(10:40, 1))
  }
}
dat <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE, check.names = FALSE)
for (nm in c("Year","Sample Size","C0 Mean","C0 SD","C0/D Mean","C0/D SD","NOS Quality Score"))
  dat[[nm]] <- suppressWarnings(as.numeric(dat[[nm]]))

# One 'median/IQR' row, one row naming a donor genotype, one grouped ABCB1 label: the cleaning rules should handle these
dat[["Data Type (Mean/SD or Median/IQR)"]][5] <- "Median/IQR"
about <- data.frame(
  README = c("SYNTHETIC EXAMPLE DATA. Author names, Study IDs and values are invented and do not come from any publication.",
             "Layout expected by the pipeline: one row per genotype group x timepoint, sheet 'Meta_Data_Extraction'.",
             "Study IDs ending in _OVERLAP-RISK are dropped in the sensitivity run; IDs containing DUPLICATE-RISK are always dropped.",
             "Regenerate with: Rscript test_data/make_example_data.R"))
write_xlsx(list(Meta_Data_Extraction = dat, README = about), "test_data/example_workbook.xlsx")
cat("Wrote test_data/example_workbook.xlsx:", nrow(dat), "rows,", length(unique(dat[["Study ID"]])), "studies\n")
