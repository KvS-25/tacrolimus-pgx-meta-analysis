# Input workbook format

One sheet (default name `Meta_Data_Extraction`), **one row per genotype group per timepoint per study**. Column names are cleaned automatically (lower-cased, non-alphanumerics become `_`), so small spelling changes in headers are tolerated, but these must exist:

| Column | Used for |
|---|---|
| Study ID | Random-effects cluster. Convention `FirstAuthor_Year_PrimaryGene_Organ` (add `a`/`b` for same author+year). |
| First Author, Year | Labels |
| Ethnicity | Free text; keyword-classified into ancestry groups (East Asian, European, South Asian, MENA, Southeast Asian, Latin American, African, Mixed, Unspecified) |
| Organ Type | Kidney / Liver / Heart / Lung / other (other becomes "Mixed") |
| Population (Adult/Pediatric) | Contains "Pediatric" for the age subgroup |
| Sample Size | Study size (funnel plots); may be text, non-numeric becomes missing |
| Gene, SNP rsID | One of the six target rsIDs. A row naming more than one target SNP is set aside |
| Genotype Groups | e.g. `*1/*3`, `*3/*3`, `Expressers (*1/*1 + *1/*3)`, `CC`, `CT`, `TT`, `GG`, `*1/*22` |
| Tacrolimus Timepoint | Free text; matched for the early (<3 mo) vs stable (>=6 mo) split |
| C0 Mean, C0 SD, C0/D Mean, C0/D SD | Numeric. C0/D is used when present, otherwise C0 |
| Data Type (Mean/SD or Median/IQR) | Rows containing "median" are flagged (counted in the log, still pooled as mean/SD) |
| NOS Quality Score | Optional, used for the quality meta-regression |
| Notes | Free text |

Other columns (Country, Study Design, Zotero Key, Acute Rejection, ...) are ignored by the analysis.

## Special Study ID flags
- `DUPLICATE-RISK` in the Study ID: always excluded (same cohort as another included study).
- `OVERLAP-RISK` in the Study ID: kept in the `main` run, dropped in the `sensitivity` run.
- Anything in `manual_exclude` (config) is dropped in both.

## Cleaning rules applied before pooling
- Rows without positive mean and SD are dropped (`vi = (SD/mean)^2` needs both).
- Genotype labels mentioning a donor or graft are dropped (liver graft genotype is a different exposure).
- Grouped or ambiguous ABCB1 labels (`CT+TT`, `carrier`, `vs`) are written to `abcb1_ambiguous_labels_set_aside.csv` and not pooled.
- CYP3A5 group labels (expresser / non-expresser) are recognised before single genotypes, so `*1/*1 + *1/*3` is read as "expresser".
