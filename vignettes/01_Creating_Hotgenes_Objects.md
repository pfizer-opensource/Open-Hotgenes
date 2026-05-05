Creating Hotgenes Objects
================

## Introduction

A **Hotgenes object** is a container that bundles together all the key
outputs of a differential expression (DE) analysis into a single,
consistently structured S4 object. Once you have a Hotgenes object, you
can use the same downstream functions—plots, reports, and the Shiny
app—regardless of which upstream analysis tool produced your DE results.

The core slots of a Hotgenes object are:

| Slot | Description | Accessor |
|----|----|----|
| `Output_DE` | Named list of DE result tables | `Output_DE_()` |
| `Normalized_Expression` | Named list of normalized expression matrices | `Normalized_Data_()` |
| `coldata` | Sample metadata data.frame | `coldata_()` |
| `designMatrix` | Model matrix used for DE | `designMatrix_()` |
| `contrastMatrix` | Contrast matrix (where applicable) | `contrastMatrix_()` |
| `Mapper` | Feature-alias mapping table | `Mapper_()` |
| `auxiliary_assays` | Optional sample-level auxiliary data | `auxiliary_assays_()` |
| `Original_Object` | The original DE object | `O_()` |

Three dedicated constructor functions handle the most common platforms:

- `HotgenesDRomics()` — DRomics benchmark dose modeling
- `HotgenesDEseq2()` — DESeq2 negative-binomial DE
- `Hotgeneslimma()` — limma/voom linear modeling

For any other platform, `HotgenesUniversal()` accepts pre-formatted
inputs directly.

------------------------------------------------------------------------

## 1. Creating a Hotgenes Object from DRomics

DRomics models dose-response relationships using benchmark dose (BMD)
methodology. `HotgenesDRomics()` converts the output of
`DRomics::bmdcalc()` into a Hotgenes object.

### Step 1 — Run the DRomics pipeline on the built-in RNAseq sample

``` r
library(Hotgenes)
library(DRomics)

# Load the example RNAseq count data bundled with DRomics
datafilename <- system.file("extdata", "RNAseq_sample.txt",
                            package = "DRomics")

# Pre-process: VST normalisation + filter
o <- DRomics::RNAseqdata(datafilename,
                         check = TRUE,
                         transfo.method = "vst")

# Select items with a significant dose-response (quadratic test)
s_quad <- DRomics::itemselect(o,
                               select.method = "quadratic",
                               FDR = 0.05)

# Fit dose-response models
f <- DRomics::drcfit(itemselect = s_quad, parallel = "no")
## The fitting may be long if the number of selected items is high.
##   |                                                                                       |                                                                               |   0%  |                                                                                       |==                                                                             |   2%  |                                                                                       |====                                                                           |   5%  |                                                                                       |======                                                                         |   7%  |                                                                                       |=======                                                                        |   9%  |                                                                                       |=========                                                                      |  12%  |                                                                                       |===========                                                                    |  14%  |                                                                                       |=============                                                                  |  16%  |                                                                                       |===============                                                                |  19%  |                                                                                       |=================                                                              |  21%  |                                                                                       |==================                                                             |  23%  |                                                                                       |====================                                                           |  26%  |                                                                                       |======================                                                         |  28%  |                                                                                       |========================                                                       |  30%  |                                                                                       |==========================                                                     |  33%  |                                                                                       |============================                                                   |  35%  |                                                                                       |=============================                                                  |  37%  |                                                                                       |===============================                                                |  40%  |                                                                                       |=================================                                              |  42%  |                                                                                       |===================================                                            |  44%  |                                                                                       |=====================================                                          |  47%  |                                                                                       |=======================================                                        |  49%  |                                                                                       |========================================                                       |  51%  |                                                                                       |==========================================                                     |  53%  |                                                                                       |============================================                                   |  56%  |                                                                                       |==============================================                                 |  58%  |                                                                                       |================================================                               |  60%  |                                                                                       |==================================================                             |  63%  |                                                                                       |===================================================                            |  65%  |                                                                                       |=====================================================                          |  67%  |                                                                                       |=======================================================                        |  70%  |                                                                                       |=========================================================                      |  72%  |                                                                                       |===========================================================                    |  74%  |                                                                                       |=============================================================                  |  77%  |                                                                                       |==============================================================                 |  79%  |                                                                                       |================================================================               |  81%  |                                                                                       |==================================================================             |  84%  |                                                                                       |====================================================================           |  86%  |                                                                                       |======================================================================         |  88%  |                                                                                       |========================================================================       |  91%  |                                                                                       |=========================================================================      |  93%  |                                                                                       |===========================================================================    |  95%  |                                                                                       |=============================================================================  |  98%  |                                                                                       |===============================================================================| 100%

# Compute benchmark doses
bmdcalc_out <- DRomics::bmdcalc(f)
```

### Step 2 — Convert to a Hotgenes object

``` r
# Convert bmdcalc result to a Hotgenes object
# The Mapper argument is optional; omit it if you have no alias table
hotDR <- HotgenesDRomics(bmdcalc = bmdcalc_out)

# Print summary
hotDR
## class: Hotgenes 
## Original class/package:  bmdcalc/
## 
## Differential expression (default thresholds): 
## |contrast | total|
## |:--------|-----:|
## |DRomics  |    38|
## 
## Available feature mapping:  Feature 
## ExpressionSlots:  Normalized_data 
## Total auxiliary assays:  0 
## Total samples:  14
```

The single contrast is named `"DRomics"` and contains one row per
dose-responsive feature. The `log2FoldChange` is computed as the log2
ratio between the expression at the maximum (or extreme) dose and the
baseline, using a small pseudocount for count data.

### Optional: add a feature alias mapper

If you have a gene symbol ↔ transcript ID table, pass it as `Mapper`.
The mapper must contain a `"Feature"` column whose values match the
rownames of the expression matrix.

``` r
# Example: build mapper from org.Mm.eg.db (mouse gene IDs)
if (requireNamespace("org.Mm.eg.db", quietly = TRUE) &&
    requireNamespace("DBI", quietly = TRUE)) {

  dbCon <- org.Mm.eg.db::org.Mm.eg_dbconn()
  sqlQuery <- "SELECT * FROM ensembl, refseq, gene_info
               WHERE ensembl._id == gene_info._id
               AND   refseq._id  == gene_info._id;"

  mapper <- DBI::dbGetQuery(dbCon, sqlQuery) |>
    dplyr::select(
      "Feature"     = "accession",
      "ensembl_id",
      "gene_name",
      "symbol"
    ) |>
    tibble::as_tibble()

  hotDR_mapped <- HotgenesDRomics(
    bmdcalc = bmdcalc_out,
    Mapper  = mapper
  )
}
```

------------------------------------------------------------------------

## 2. Creating a Hotgenes Object from DESeq2

`HotgenesDEseq2()` wraps `DESeq2::lfcShrink()` and extracts normalised
expression data in one step.

### Step 1 — Run the DESeq2 pipeline

The example uses `dds_con.Rdata`, an RNA-seq count dataset shipped with
the Hotgenes package.

``` r
library(DESeq2)

# Load example count data and coldata
dds_con_dir <- system.file("extdata", "dds_con.Rdata",
                            package = "Hotgenes",
                            mustWork = TRUE)
load(dds_con_dir)

# Extract raw counts and sample metadata
cts    <- DESeq2::counts(dds_con) |> as.data.frame()
Design <- SummarizedExperiment::colData(dds_con) |>
  base::as.data.frame() |>
  dplyr::select_if(is.factor) |>
  dplyr::mutate(Time = as.numeric(levels(.data$Hrs))[.data$Hrs])

# Build a DESeqDataSet
model_DESeq <- eval(~ sh * Hrs)

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = cts,
  colData   = Design,
  design    = model_DESeq
)

# Fit the model
dds <- DESeq2::DESeq(dds)
```

### Step 2 — Optional: build an ENSEMBL → gene symbol mapper

``` r
# Requires org.Hs.eg.db (in Suggests)
if (requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
    requireNamespace("DBI",          quietly = TRUE)) {

  dbCon    <- org.Hs.eg.db::org.Hs.eg_dbconn()
  sqlQuery <- "SELECT * FROM ENSEMBL, gene_info
               WHERE ENSEMBL._id == gene_info._id;"

  ensembl_Symbol <- DBI::dbGetQuery(dbCon, sqlQuery) |>
    dplyr::select("Feature" = "symbol", "ensembl_id")
}
```

### Step 3 — Convert to a Hotgenes object

``` r
# Convert to Hotgenes – request only one contrast for speed
dds_Hotgenes <- HotgenesDEseq2(
  DEseq2_object  = dds,
  lfcShrink_type = "apeglm",
  contrasts      = "sh_EWS_vs_Ctrl",
  ExpressionData = "vsd"
)

dds_Hotgenes
## class: Hotgenes 
## Original class/package:  DESeqDataSet/DESeq2
## 
## Differential expression (default thresholds): 
## |contrast       | total|
## |:--------------|-----:|
## |sh_EWS_vs_Ctrl |    52|
## 
## Available feature mapping:  Feature 
## ExpressionSlots:  vsd 
## Total auxiliary assays:  0 
## Total samples:  12
```

Key parameters:

| Parameter | Purpose |
|----|----|
| `lfcShrink_type` | Shrinkage estimator: `"normal"`, `"apeglm"`, or `"ashr"` |
| `contrasts` | Character vector of contrast names (from `DESeq2::resultsNames(dds)`) |
| `contrast_vectors` | Alternative: named list of `c("factor","level1","level2")` vectors |
| `ExpressionData` | Normalization: `"vsd"`, `"rld"`, or `"log2"` |
| `Mapper` | Optional feature alias table |

------------------------------------------------------------------------

## 3. Creating a Hotgenes Object from limma

`Hotgeneslimma()` accepts a limma `MArrayLM` fit object together with
expression data and sample metadata. This example reproduces the
Drug-Induced Liver Injury (DILI) proteomics analysis that is fully
documented in the [DILI Discovery Proteomics
vignette](DILI_Hotgenes.md).

``` r
# Load the pre-parsed DILI proteomics snapshot
raw_df <- readRDS(
  system.file("extdata", "dili_raw.RDS",
              package = "Hotgenes",
              mustWork = TRUE)
)
```

``` r
# Identify sample columns
sample_cols <- colnames(raw_df)[
  grepl("^(do|df|hv|nafld|ndo|ndf)_", colnames(raw_df))
]

# Filter proteins and build expression matrix
filtered_exps <- raw_df |>
  dplyr::filter(!grepl("contaminant", .data$protein)) |>
  dplyr::filter(grepl("^sp[|]", .data$protein)) |>
  dplyr::filter(.data$min_peps > 0) |>
  dplyr::mutate(Feature = make.names(.data$gene_symbol, unique = TRUE))

expr_matrix <- filtered_exps |>
  dplyr::select("Feature", dplyr::any_of(sample_cols)) |>
  tibble::column_to_rownames("Feature") |>
  as.matrix()

# Feature alias mapper
mapper_df <- filtered_exps |>
  dplyr::select(
    "Feature",
    "Gene"        = "gene_symbol",
    "Protein"     = "protein",
    "Description" = "description"
  )

# Sample metadata
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
```

``` r
# VSN normalisation
expr_matrix[is.na(expr_matrix)] <- 0
expr_matrix_vsn <- limma::normalizeVSN(expr_matrix)

# Design and voomaByGroup
design <- model.matrix(~ 0 + Condition, data = coldata)
colnames(design) <- gsub("Condition", "", colnames(design))

vm_exp <- limma::voomaByGroup(
  y      = expr_matrix_vsn,
  group  = coldata[["Pool"]],
  design = design,
  plot   = FALSE
)

# Robust linear model
fit <- limma::lmFit(vm_exp, design, method = "robust")

# Define contrasts and run eBayes
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

``` r
# Bundle everything into a Hotgenes object
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

Key parameters for `Hotgeneslimma()`:

| Parameter | Purpose |
|----|----|
| `limmafit` | `MArrayLM` object from `limma::lmFit()` + `eBayes()` |
| `coldata` | Sample metadata (rownames must match `colnames(Expression)`) |
| `Expression` | `EList` from `limma::voom()` or `voomaByGroup()` |
| `Expression_name` | Name for the primary expression slot (e.g., `"VSN"`) |
| `Exps_list` | Optional named list of additional expression matrices |
| `Mapper` | Optional feature alias table |
| `contrasts` | Subset of contrasts to include (default: all) |

For the complete DILI walkthrough with detailed commentary, see the
companion vignette **DILI Discovery Proteomics with Hotgenes**.

------------------------------------------------------------------------

## 4. Summary: Common Features of All Hotgenes Objects

Regardless of origin, every Hotgenes object exposes the same accessors:

``` r
# Using the DESeq2 example from Section 2
htgs <- readRDS(
  system.file("extdata", "dds_Hotgenes.RDS",
              package = "Hotgenes",
              mustWork = TRUE)
) |> update_object()

# Available contrasts
contrasts_(htgs)
## [1] "sh_EWS_vs_Ctrl" "Hrs_2_vs_0"     "Hrs_6_vs_0"     "shEWS.Hrs2"     "shEWS.Hrs6"

# Sample IDs
SampleIDs_(htgs) |> head()
## [1] "shCON_0hrs_1" "shCON_0hrs_2" "shCON_2hrs_1" "shCON_2hrs_2" "shCON_6hrs_1"
## [6] "shCON_6hrs_2"

# Available expression slots
ExpressionSlots_(htgs)
## [1] "rld"  "vsd"  "log2"

# Available feature aliases
Mapper_(htgs) |> names()
## [1] "Feature"    "ensembl_id"
```

The object summary (printed automatically on call) shows the number of
DE features per contrast at the default thresholds:

``` r
htgs
## class: Hotgenes 
## Original class/package:  DESeqDataSet/DESeq2
## 
## Differential expression (default thresholds): 
## |contrast       | total|
## |:--------------|-----:|
## |sh_EWS_vs_Ctrl |    52|
## |Hrs_2_vs_0     |    40|
## |Hrs_6_vs_0     |    63|
## |shEWS.Hrs2     |     5|
## |shEWS.Hrs6     |     1|
## 
## Available feature mapping:  Feature, ensembl_id 
## ExpressionSlots:  rld, vsd, log2 
## Total auxiliary assays:  0 
## Total samples:  12
```

For details on how to use these objects for analysis and visualization,
see the companion vignettes:

- **02 API and Methods** — accessor functions and DE reports
- **03 Visualization and Exploration** — plots and heatmaps
- **04 Interactive Exploration with Shiny** — the Shiny app
