DILI Discovery Proteomics with Hotgenes
================

## Building a Hotgenes Object from DILI Discovery Proteomics Data

This vignette demonstrates how to build a Hotgenes object from the
Drug-Induced Liver Injury (DILI) discovery proteomics dataset published
in:

> Ravindra, K.C., Vaidya, V.S., Wang, Z. et al. *Tandem mass tag-based
> quantitative proteomic profiling identifies candidate serum biomarkers
> of drug-induced liver injury in humans* **Nat Commun**, 14, 1215
> (2023). <https://doi.org/10.1038/s41467-023-36858-6>

The raw data are publicly deposited in MassIVE under accession
**MSV000089782**. A pre-parsed snapshot is shipped with this package as
`inst/extdata/dili_raw.RDS` to avoid a runtime download.

### A note on private metadata

The original analysis used two variables from the private metadata file
`1219_p1-p5_protein_KEY.xlsx` that are **not** available in the public
deposit:

1.  **subject** – patient identifier, used with `duplicateCorrelation()`
    to control for repeated measures (paired DO/DF samples from the same
    patient). Without subject IDs this blocking step is omitted here.

2.  **channel** – exact TMT channel label (e.g. `126C`, `127N`) within
    each pool, used with `voomaByGroup()` to compensate for
    channel-specific variance. `Pool` (P1–P5) is used as a proxy because
    it captures the dominant source of TMT batch variance.

All other analytical steps — VSN normalization, `voomaByGroup`, robust
`lmFit`, and a BH-adjusted p-value threshold of 0.1 — reproduce those
described in Federspiel et al.

## Step 1 — Load data

`dili_raw.RDS` is a pre-parsed snapshot of the public MassIVE deposit
(MSV000089782), saved to `inst/extdata` to avoid a runtime download.

``` r
library(Hotgenes)

raw_df <- readRDS(
  system.file("extdata", "dili_raw.RDS",
              package = "Hotgenes",
              mustWork = TRUE)
)
```

To re-download and refresh `dili_raw.RDS` directly from MassIVE, run the
following block once (set `eval = FALSE` to skip during package build):

``` r
url <- paste0(
  "https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?",
  "file=f.MSV000089782%2Fupdates%2F2023-02-28_jfederspiel_1bc96582",
  "%2Fother%2FDILI_discovery_data.xlsx&forceDownload=true"
)

tmp <- tempfile(fileext = ".xlsx")
download.file(url, tmp, mode = "wb")

# Requires readxl and janitor packages (install if needed):
# install.packages(c("readxl", "janitor"))
raw_df <- readxl::read_excel(tmp) |>
  janitor::clean_names()

saveRDS(
  raw_df,
  file = file.path(getwd(), "inst", "extdata", "dili_raw.RDS")
)
```

## Step 2 — Identify sample columns

Sample columns follow the pattern `<condition>_p<pool>_<number>` where
condition is one of `do`, `df`, `hv`, `nafld`, `ndo`, or `ndf`.

``` r
sample_cols <- colnames(raw_df)[
  grepl("^(do|df|hv|nafld|ndo|ndf)_", colnames(raw_df))
]

length(sample_cols)
## [1] 50
head(sample_cols)
## [1] "do_p1_1"  "do_p1_2"  "do_p2_19" "do_p2_20" "do_p3_24" "do_p3_27"
```

## Step 3 — Filter proteins

Retain only:

- Reviewed UniProt entries (accession starts with `sp|`)
- Proteins that are not contaminants
- Proteins quantified with at least one peptide (`min_peps > 0`)

Unique gene symbols are generated as the `Feature` column, which is
required downstream for GSEA.

``` r
filtered_exps <- raw_df |>
  dplyr::filter(!grepl("contaminant", .data$protein)) |>
  dplyr::filter(grepl("^sp[|]", .data$protein)) |>
  dplyr::filter(.data$min_peps > 0) |>
  dplyr::mutate(Feature = make.names(.data$gene_symbol, unique = TRUE))

nrow(filtered_exps)
## [1] 1463
```

## Step 4 — Build the expression matrix

``` r
expr_matrix <- filtered_exps |>
  dplyr::select("Feature", dplyr::any_of(sample_cols)) |>
  tibble::column_to_rownames("Feature") |>
  as.matrix()

dim(expr_matrix)
## [1] 1463   50
```

## Step 5 — Build the protein ID mapper

The mapper links `Feature` names back to the original gene symbol, full
protein accession, and description. Hotgenes functions use the mapper to
display human-readable labels in plots and tables.

``` r
mapper_df <- filtered_exps |>
  dplyr::select(
    "Feature",
    "Gene"        = "gene_symbol",
    "Protein"     = "protein",
    "Description" = "description"
  )

head(mapper_df)
## # A tibble: 6 × 4
##   Feature  Gene     Protein                   Description                                                        
##   <chr>    <chr>    <chr>                     <chr>                                                              
## 1 IGLV4.69 IGLV4-69 sp|A0A075B6H9|LV469_HUMAN Immunoglobulin lambda variable 4-69 OS=Homo sapiens GN=IGLV4-69 PE…
## 2 IGLV4.60 IGLV4-60 sp|A0A075B6I1|LV460_HUMAN Immunoglobulin lambda variable 4-60 OS=Homo sapiens GN=IGLV4-60 PE…
## 3 IGLV3.16 IGLV3-16 sp|A0A075B6K0|LV316_HUMAN Immunoglobulin lambda variable 3-16 OS=Homo sapiens GN=IGLV3-16 PE…
## 4 IGLV3.10 IGLV3-10 sp|A0A075B6K4|LV310_HUMAN Immunoglobulin lambda variable 3-10 OS=Homo sapiens GN=IGLV3-10 PE…
## 5 IGLV3.9  IGLV3-9  sp|A0A075B6K5|LV39_HUMAN  Immunoglobulin lambda variable 3-9 OS=Homo sapiens GN=IGLV3-9 PE=3…
## # ℹ 1 more row
```

## Step 6 — Build sample metadata (coldata)

Condition and pool are parsed from the sample column names.  
`Condition` levels are ordered so that healthy volunteers (`HV`) appear
first.

``` r
coldata <- data.frame(
  Sample    = sample_cols,
  Condition = toupper(sub("_p[0-9]+_[0-9]+$", "", sample_cols)),
  Pool      = toupper(sub(".*_(p[0-9]+)_.*", "\\1", sample_cols)),
  row.names = sample_cols,
  stringsAsFactors = TRUE
)

coldata[["Condition"]] <- factor(
  coldata[["Condition"]],
  levels = c("HV", "DO", "DF", "NDO", "NDF", "NAFLD")
)

table(coldata$Condition)
## 
##    HV    DO    DF   NDO   NDF NAFLD 
##    10    10    10     5     5    10
```

## Step 7 — VSN normalization

Variance-stabilising normalisation (`limma::normalizeVSN`) is applied to
correct for intensity-dependent variance across the dynamic range of the
TMT measurements. Missing values are replaced with 0 before
normalisation.

``` r
expr_matrix[is.na(expr_matrix)] <- 0
expr_matrix_vsn <- limma::normalizeVSN(expr_matrix)
## vsn2: 1463 x 50 matrix (1 stratum).
## Please use 'meanSdPlot' to verify the fit.
```

## Step 8 — Design matrix

A no-intercept design matrix is constructed so that each coefficient
directly represents the mean VSN-normalised expression for one
condition, making contrasts straightforward to specify.

``` r
design <- model.matrix(~ 0 + Condition, data = coldata)
colnames(design) <- gsub("Condition", "", colnames(design))

head(design)
##          HV DO DF NDO NDF NAFLD
## do_p1_1   0  1  0   0   0     0
## do_p1_2   0  1  0   0   0     0
## do_p2_19  0  1  0   0   0     0
## do_p2_20  0  1  0   0   0     0
## do_p3_24  0  1  0   0   0     0
## do_p3_27  0  1  0   0   0     0
```

## Step 9 — voomaByGroup transformation

`limma::voomaByGroup` estimates mean–variance trends separately for each
TMT pool and converts the VSN-normalised intensities into precision
weights for downstream linear modelling. `Pool` is used as the grouping
variable because it captures the dominant source of batch variance in
the public deposit.

``` r
vm_exp <- limma::voomaByGroup(
  y      = expr_matrix_vsn,
  group  = coldata[["Pool"]],
  design = design,
  plot   = FALSE
)
```

## Step 10 — Robust linear model fit

`limma::lmFit` with `method = "robust"` downweights outlier samples when
estimating per-protein coefficients, improving robustness against
occasional poorly quantified channels.

Note: `duplicateCorrelation()` is omitted here because the patient
subject IDs required to model repeated measures (paired DO/DF samples)
are not available in the public deposit.

``` r
fit <- limma::lmFit(
  vm_exp,
  design,
  method = "robust"
)
## Warning in rlm.default(x = X, y = y, weights = w, ...): 'rlm' failed to converge in 20 steps
## Warning in rlm.default(x = X, y = y, weights = w, ...): 'rlm' failed to converge in 20 steps
## Warning in rlm.default(x = X, y = y, weights = w, ...): 'rlm' failed to converge in 20 steps
## Warning in rlm.default(x = X, y = y, weights = w, ...): 'rlm' failed to converge in 20 steps
## Warning in rlm.default(x = X, y = y, weights = w, ...): 'rlm' failed to converge in 20 steps
## Warning in rlm.default(x = X, y = y, weights = w, ...): 'rlm' failed to converge in 20 steps
## Warning in rlm.default(x = X, y = y, weights = w, ...): 'rlm' failed to converge in 20 steps
## Warning in rlm.default(x = X, y = y, weights = w, ...): 'rlm' failed to converge in 20 steps
## Warning in rlm.default(x = X, y = y, weights = w, ...): 'rlm' failed to converge in 20 steps
## Warning in rlm.default(x = X, y = y, weights = w, ...): 'rlm' failed to converge in 20 steps
## Warning in rlm.default(x = X, y = y, weights = w, ...): 'rlm' failed to converge in 20 steps
## Warning in rlm.default(x = X, y = y, weights = w, ...): 'rlm' failed to converge in 20 steps
## Warning in rlm.default(x = X, y = y, weights = w, ...): 'rlm' failed to converge in 20 steps
## Warning in rlm.default(x = X, y = y, weights = w, ...): 'rlm' failed to converge in 20 steps
```

## Step 11 — Define contrasts and run eBayes

Seven clinically motivated contrasts are tested, all relative to healthy
volunteers (`HV`), plus one intra-DILI comparison (`DF_vs_DO`) and one
cross-group comparison (`NDO_vs_DO`).

``` r
contrasts_mat <- limma::makeContrasts(
  DO_vs_HV    = DO    - HV,
  DF_vs_HV    = DF    - HV,
  NDO_vs_HV   = NDO   - HV,
  NDF_vs_HV   = NDF   - HV,
  NAFLD_vs_HV = NAFLD - HV,
  DF_vs_DO    = DF    - DO,
  NDO_vs_DO   = NDO   - DO,
  levels = design
)

fit2 <- limma::contrasts.fit(fit, contrasts_mat)
fit2 <- limma::eBayes(fit2)
```

## Step 12 — Create the Hotgenes object

All upstream objects are bundled into a single `Hotgenes` object using
`Hotgeneslimma()`. The raw log2-transformed matrix is stored as an
auxiliary assay so it is available for alternative visualizations.

``` r
dili_hotgenes <- Hotgeneslimma(
  limmafit        = fit2,
  coldata         = coldata,
  Expression      = vm_exp,
  Expression_name = "VSN",
  Exps_list       = list(log2 = log2(expr_matrix + 1)),
  Mapper          = mapper_df
)

dili_hotgenes
## class: Hotgenes 
## Original class/package:  EList/limma
## 
## Differential expression (default thresholds): 
## |contrast    | total|
## |:-----------|-----:|
## |DF_vs_DO    |    61|
## |DF_vs_HV    |    73|
## |DO_vs_HV    |   212|
## |NAFLD_vs_HV |   428|
## |NDF_vs_HV   |    78|
## |NDO_vs_DO   |    73|
## |NDO_vs_HV   |   303|
## 
## Available feature mapping:  Feature, Gene, Protein, Description 
## ExpressionSlots:  VSN, log2 
## Total auxiliary assays:  0 
## Total samples:  50
```

## What can you do with the DILI Hotgenes object?
