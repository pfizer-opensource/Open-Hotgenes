Hotgenes API and Methods
================

## Overview

This vignette covers the core API for interacting with a Hotgenes object
after it has been created. All examples use the pre-built `dds_Hotgenes`
object that ships with the package (created from a DESeq2 analysis of
Ewing sarcoma RNAseq data).

For instructions on creating a Hotgenes object from scratch, see **01
Creating Hotgenes Objects**.

``` r
library(Hotgenes)

# Load the pre-built DESeq2-based example
htgs <- readRDS(
  system.file("extdata", "dds_Hotgenes.RDS",
              package = "Hotgenes",
              mustWork = TRUE)
) |> update_object()
```

------------------------------------------------------------------------

## 1. Object Summary

Calling a Hotgenes object prints a summary of its contents:

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

------------------------------------------------------------------------

## 2. Accessing Slots

### Sample metadata

``` r
coldata_(htgs)
##                sh Bio_Rep Hrs Time
## shCON_0hrs_1 Ctrl       1   0    0
## shCON_0hrs_2 Ctrl       2   0    0
## shCON_2hrs_1 Ctrl       1   2    2
## shCON_2hrs_2 Ctrl       2   2    2
## shCON_6hrs_1 Ctrl       1   6    6
## shCON_6hrs_2 Ctrl       2   6    6
## shEWS_0hrs_1  EWS       1   0    0
## shEWS_0hrs_2  EWS       2   0    0
## shEWS_2hrs_1  EWS       1   2    2
## shEWS_2hrs_2  EWS       2   2    2
## shEWS_6hrs_1  EWS       1   6    6
## shEWS_6hrs_2  EWS       2   6    6
```

### Normalized expression data

`Normalized_Data_()` returns the named list of expression matrices. Use
`ExpressionSlots_()` to see available slot names.

``` r
ExpressionSlots_(htgs)
## [1] "rld"  "vsd"  "log2"
```

``` r
# Access a specific slot by name
Normalized_Data_(htgs, slot = "vsd")[1:4, 1:4]
##        shCON_0hrs_1 shCON_0hrs_2 shCON_2hrs_1 shCON_2hrs_2
## AGER      10.062846    10.125494    10.092091    10.040882
## ALOX12     9.927860     9.929179     9.931316     9.928084
## ALOX15    10.062846    10.072337    10.047724    10.022307
## ALOX5      9.969462     9.971327     9.900875    10.008864
```

### Sample IDs

``` r
SampleIDs_(htgs)
##  [1] "shCON_0hrs_1" "shCON_0hrs_2" "shCON_2hrs_1" "shCON_2hrs_2" "shCON_6hrs_1"
##  [6] "shCON_6hrs_2" "shEWS_0hrs_1" "shEWS_0hrs_2" "shEWS_2hrs_1" "shEWS_2hrs_2"
## [11] "shEWS_6hrs_1" "shEWS_6hrs_2"
```

### Feature alias mapper

``` r
Mapper_(htgs) |> head()
## # A tibble: 6 × 2
##   Feature ensembl_id     
##   <chr>   <chr>          
## 1 AGER    ENSG00000204305
## 2 AGER    ENSG00000231268
## 3 AGER    ENSG00000229058
## 4 AGER    ENSG00000237405
## 5 AGER    ENSG00000234729
## 6 AGER    ENSG00000206320
```

### Design matrix

``` r
designMatrix_(htgs) |> head()
##              Intercept sh_EWS_vs_Ctrl Hrs_2_vs_0 Hrs_6_vs_0 shEWS.Hrs2 shEWS.Hrs6
## shCON_0hrs_1         1              0          0          0          0          0
## shCON_0hrs_2         1              0          0          0          0          0
## shCON_2hrs_1         1              0          1          0          0          0
## shCON_2hrs_2         1              0          1          0          0          0
## shCON_6hrs_1         1              0          0          1          0          0
## shCON_6hrs_2         1              0          0          1          0          0
```

### Original DE object

`O_()` returns the original analysis object stored at creation time (a
`DESeqDataSet` for DESeq2-derived objects, an `EList` with an embedded
limma fit for limma-derived objects).

``` r
O_(htgs) |> class()
## [1] "DESeqDataSet"
## attr(,"package")
## [1] "DESeq2"
```

### Auxiliary assays

Auxiliary assays store per-sample data that do not fit into the
expression matrix—clinical scores, QC metrics, low-throughput assay
values, etc. The slot starts empty (only a `SampleIDs` column).

``` r
auxiliary_assays_(htgs)
## data frame with 0 columns and 12 rows
```

Add data with the assignment form of the accessor:

``` r
set.seed(42)
n_samples <- length(SampleIDs_(htgs))

new_aux <- auxiliary_assays_default(htgs) |>
  dplyr::mutate(
    score_A = rnorm(n_samples, mean = 5, sd = 1),
    score_B = runif(n_samples, min = 0, max = 10)
  )

auxiliary_assays_(htgs) <- new_aux
auxiliary_assays_(htgs)
##               score_A    score_B
## shCON_0hrs_1 6.370958 0.82437558
## shCON_0hrs_2 4.435302 5.14211784
## shCON_2hrs_1 5.363128 3.90203467
## shCON_2hrs_2 5.632863 9.05738131
## shCON_6hrs_1 5.404268 4.46969628
## shCON_6hrs_2 4.893875 8.36004260
## shEWS_0hrs_1 6.511522 7.37595618
## shEWS_0hrs_2 4.905341 8.11055141
## shEWS_2hrs_1 7.018424 3.88108283
## shEWS_2hrs_2 4.937286 6.85169729
## shEWS_6hrs_1 6.304870 0.03948339
## shEWS_6hrs_2 7.286645 8.32916080
```

------------------------------------------------------------------------

## 3. The `DE()` Function

`DE()` is the main entry point for querying differential expression
results. It supports filtering by contrast, adjusted p-value, log2
fold-change, and gene name.

### Default: all significant features across all contrasts

``` r
DE(htgs, 
   # limit to top 10
   Topn = 10) 
## $sh_EWS_vs_Ctrl
## # A tibble: 10 × 11
##    Feature contrast_dir       baseMean log2FoldChange    FC   stat   pvalue     padj
##    <chr>   <chr>                 <dbl>          <dbl> <dbl>  <dbl>    <dbl>    <dbl>
##  1 C1R     sh_EWS_vs_Ctrl_do…    5967.         -0.799 0.575 -15.1  1.16e-51 1.99e-49
##  2 C1R     sh_EWS_vs_Ctrl_do…    5967.         -0.799 0.575 -15.1  1.16e-51 1.99e-49
##  3 CCL2    sh_EWS_vs_Ctrl_do…    3119.         -0.945 0.519 -14.8  1.61e-49 1.37e-47
##  4 C1S     sh_EWS_vs_Ctrl_do…    3514.         -0.594 0.662 -11.7  1.36e-31 7.76e-30
##  5 MMP3    sh_EWS_vs_Ctrl_up     3986.          0.578 1.49    8.80 1.31e-18 5.60e-17
##  6 HMGB2   sh_EWS_vs_Ctrl_up     1118.          0.621 1.54    7.92 2.29e-15 7.83e-14
##  7 STAT1   sh_EWS_vs_Ctrl_do…    4448.         -0.355 0.782  -7.24 4.52e-13 1.29e-11
##  8 CXCL6   sh_EWS_vs_Ctrl_do…     130.         -1.65  0.318  -7.02 2.27e-12 5.54e-11
##  9 IL1R1   sh_EWS_vs_Ctrl_do…    1519.         -0.414 0.751  -6.65 2.88e-11 6.16e-10
## 10 RAC1    sh_EWS_vs_Ctrl_up     9136.          0.265 1.20    6.38 1.80e-10 3.43e- 9
## # ℹ 3 more variables: lfcSE <dbl>, ensembl_id <chr>, significant <chr>
## 
## $Hrs_2_vs_0
## # A tibble: 10 × 11
##    Feature contrast_dir  baseMean log2FoldChange    FC  stat    pvalue      padj
##    <chr>   <chr>            <dbl>          <dbl> <dbl> <dbl>     <dbl>     <dbl>
##  1 IL6     Hrs_2_vs_0_up     914.          2.30   4.94 26.7  7.56e-157 1.93e-154
##  2 TNFAIP3 Hrs_2_vs_0_up     780.          1.92   3.79 21.3  5.75e-101 7.34e- 99
##  3 CXCL8   Hrs_2_vs_0_up    1409.          2.88   7.35 19.4  8.35e- 84 7.10e- 82
##  4 IL11    Hrs_2_vs_0_up     659.          1.80   3.49 18.3  1.93e- 74 1.23e- 72
##  5 PTGS2   Hrs_2_vs_0_up    1907.          1.12   2.17 15.5  3.41e- 54 1.74e- 52
##  6 SMAD7   Hrs_2_vs_0_up    3674.          0.901  1.87 14.9  1.98e- 50 8.40e- 49
##  7 CXCL1   Hrs_2_vs_0_up    1533.          1.87   3.67 14.7  1.08e- 48 3.92e- 47
##  8 CCL2    Hrs_2_vs_0_up    3119.          0.808  1.75 13.5  1.63e- 41 5.19e- 40
##  9 JUN     Hrs_2_vs_0_up    1619.          0.742  1.67  9.44 3.84e- 21 1.09e- 19
## 10 CXCL3   Hrs_2_vs_0_up     132.          1.30   2.47  8.27 1.37e- 16 3.50e- 15
## # ℹ 3 more variables: lfcSE <dbl>, ensembl_id <chr>, significant <chr>
## 
## $Hrs_6_vs_0
## # A tibble: 10 × 11
##    Feature contrast_dir baseMean log2FoldChange    FC  stat   pvalue     padj  lfcSE
##    <chr>   <chr>           <dbl>          <dbl> <dbl> <dbl>    <dbl>    <dbl>  <dbl>
##  1 CCL2    Hrs_6_vs_0_…    3119.          1.25   2.38  20.8 3.57e-96 7.14e-94 0.0507
##  2 PTGS2   Hrs_6_vs_0_…    1907.          1.42   2.68  20.2 1.35e-90 1.35e-88 0.0576
##  3 IL6     Hrs_6_vs_0_…     914.          1.73   3.32  20.0 8.93e-89 5.95e-87 0.0724
##  4 CXCL8   Hrs_6_vs_0_…    1409.          3.04   8.24  19.8 4.19e-87 2.09e-85 0.122 
##  5 TNFAIP3 Hrs_6_vs_0_…     780.          1.67   3.18  17.9 6.41e-72 2.56e-70 0.0741
##  6 CXCL1   Hrs_6_vs_0_…    1533.          2.23   4.69  16.4 1.46e-60 4.86e-59 0.101 
##  7 RELB    Hrs_6_vs_0_…     436.          1.83   3.55  16.0 9.66e-58 2.76e-56 0.0874
##  8 IL11    Hrs_6_vs_0_…     659.          1.26   2.40  13.0 1.21e-38 3.03e-37 0.0823
##  9 BIRC2   Hrs_6_vs_0_…    5074.          0.702  1.63  12.1 1.42e-33 3.15e-32 0.0509
## 10 NFKB1   Hrs_6_vs_0_…     705.          0.943  1.92  10.6 3.94e-26 7.88e-25 0.0724
## # ℹ 2 more variables: ensembl_id <chr>, significant <chr>
## 
## $shEWS.Hrs2
## # A tibble: 5 × 11
##   Feature contrast_dir    baseMean log2FoldChange    FC  stat   pvalue   padj  lfcSE
##   <chr>   <chr>              <dbl>          <dbl> <dbl> <dbl>    <dbl>  <dbl>  <dbl>
## 1 TNFAIP3 shEWS.Hrs2_down     780.        -0.272  0.828 -4.07  4.62e-5 0.0118 0.0680
## 2 IL6     shEWS.Hrs2_down     914.        -0.204  0.868 -3.57  3.61e-4 0.0330 0.0672
## 3 TGFBR1  shEWS.Hrs2_down    3125.        -0.194  0.874 -3.55  3.88e-4 0.0330 0.0570
## 4 CXCL3   shEWS.Hrs2_down     132.        -0.0801 0.946 -3.18  1.46e-3 0.0932 0.0625
## 5 CXCL8   shEWS.Hrs2_down    1409.        -0.0896 0.940 -3.10  1.92e-3 0.0980 0.0675
## # ℹ 2 more variables: ensembl_id <chr>, significant <chr>
## 
## $shEWS.Hrs6
## # A tibble: 1 × 11
##   Feature contrast_dir  baseMean log2FoldChange    FC  stat   pvalue   padj  lfcSE
##   <chr>   <chr>            <dbl>          <dbl> <dbl> <dbl>    <dbl>  <dbl>  <dbl>
## 1 FOS     shEWS.Hrs6_up     263.          0.155  1.11  3.57 0.000359 0.0916 0.0545
## # ℹ 2 more variables: ensembl_id <chr>, significant <chr>
```

### Select specific contrasts

``` r
# List available contrasts
contrasts_(htgs)
## [1] "sh_EWS_vs_Ctrl" "Hrs_2_vs_0"     "Hrs_6_vs_0"     "shEWS.Hrs2"    
## [5] "shEWS.Hrs6"
```

``` r
DE(htgs,
   contrasts      = "sh_EWS_vs_Ctrl",
   padj_cut       = 0.05,
   .log2FoldChange = 1) 
## $sh_EWS_vs_Ctrl
## # A tibble: 2 × 11
##   Feature contrast_dir   baseMean log2FoldChange    FC  stat   pvalue     padj lfcSE
##   <chr>   <chr>             <dbl>          <dbl> <dbl> <dbl>    <dbl>    <dbl> <dbl>
## 1 CXCL6   sh_EWS_vs_Ctr…    130.           -1.65 0.318 -7.02 2.27e-12 5.54e-11 0.107
## 2 C3      sh_EWS_vs_Ctr…     49.3          -1.23 0.428 -4.78 1.76e- 6 1.67e- 5 0.138
## # ℹ 2 more variables: ensembl_id <chr>, significant <chr>
```

### Look up specific genes with `hotList`

When `hotList` is provided, all matches are returned regardless of
significance thresholds. A `significant` column is added to flag
features that pass the thresholds.

``` r
DE(htgs,
   hotList        = c("CSF2", "IL6", "CCL2"),
   padj_cut       = 0.05,
   .log2FoldChange = 0.5)
## $sh_EWS_vs_Ctrl
## # A tibble: 3 × 11
##   Feature contrast_dir baseMean log2FoldChange    FC   stat   pvalue     padj  lfcSE
##   <chr>   <chr>           <dbl>          <dbl> <dbl>  <dbl>    <dbl>    <dbl>  <dbl>
## 1 CCL2    sh_EWS_vs_C…   3119.          -0.945 0.519 -14.8  1.61e-49 1.37e-47 0.0474
## 2 IL6     sh_EWS_vs_C…    914.          -0.462 0.726  -2.41 1.59e- 2 5.83e- 2 0.0596
## 3 CSF2    sh_EWS_vs_C…     38.4          0.784 1.72    1.65 9.95e- 2 2.37e- 1 0.147 
## # ℹ 2 more variables: ensembl_id <chr>, significant <chr>
## 
## $Hrs_2_vs_0
## # A tibble: 3 × 11
##   Feature contrast_dir  baseMean log2FoldChange    FC  stat    pvalue      padj
##   <chr>   <chr>            <dbl>          <dbl> <dbl> <dbl>     <dbl>     <dbl>
## 1 IL6     Hrs_2_vs_0_up    914.           2.30   4.94 26.7  7.56e-157 1.93e-154
## 2 CCL2    Hrs_2_vs_0_up   3119.           0.808  1.75 13.5  1.63e- 41 5.19e- 40
## 3 CSF2    Hrs_2_vs_0_up     38.4          0.966  1.95  2.94 3.31e-  3 2.72e-  2
## # ℹ 3 more variables: lfcSE <dbl>, ensembl_id <chr>, significant <chr>
## 
## $Hrs_6_vs_0
## # A tibble: 3 × 11
##   Feature contrast_dir  baseMean log2FoldChange    FC  stat   pvalue     padj  lfcSE
##   <chr>   <chr>            <dbl>          <dbl> <dbl> <dbl>    <dbl>    <dbl>  <dbl>
## 1 CCL2    Hrs_6_vs_0_up   3119.            1.25  2.38 20.8  3.57e-96 7.14e-94 0.0507
## 2 IL6     Hrs_6_vs_0_up    914.            1.73  3.32 20.0  8.93e-89 5.95e-87 0.0724
## 3 CSF2    Hrs_6_vs_0_up     38.4           2.15  4.43  4.99 6.13e- 7 4.09e- 6 0.232 
## # ℹ 2 more variables: ensembl_id <chr>, significant <chr>
## 
## $shEWS.Hrs2
## # A tibble: 3 × 11
##   Feature contrast_dir    baseMean log2FoldChange    FC   stat  pvalue   padj  lfcSE
##   <chr>   <chr>              <dbl>          <dbl> <dbl>  <dbl>   <dbl>  <dbl>  <dbl>
## 1 IL6     shEWS.Hrs2_down    914.          -0.204 0.868 -3.57  3.61e-4 0.0330 0.0672
## 2 CCL2    shEWS.Hrs2_down   3119.          -0.107 0.929 -1.19  2.33e-1 0.983  0.0616
## 3 CSF2    shEWS.Hrs2_up       38.4          0.121 1.09   0.437 6.62e-1 0.983  0.0421
## # ℹ 2 more variables: ensembl_id <chr>, significant <chr>
## 
## $shEWS.Hrs6
## # A tibble: 3 × 11
##   Feature contrast_dir    baseMean log2FoldChange    FC   stat pvalue  padj  lfcSE
##   <chr>   <chr>              <dbl>          <dbl> <dbl>  <dbl>  <dbl> <dbl>  <dbl>
## 1 IL6     shEWS.Hrs6_down    914.         -0.0935 0.937 -2.48  0.0131 0.658 0.0584
## 2 CCL2    shEWS.Hrs6_up     3119.          0.0489 1.03   0.881 0.378  0.997 0.0561
## 3 CSF2    shEWS.Hrs6_up       38.4         0.123  1.09   0.672 0.501  0.997 0.0317
## # ℹ 2 more variables: ensembl_id <chr>, significant <chr>
```

### `Report` options

The `Report` argument controls what `DE()` returns:

| `Report` | Returns |
|----|----|
| `"Details"` (default) | Complete DE result table |
| `"Features"` | Character vector of feature names per contrast |
| `"contrast_dir"` | Feature names named by `<contrast>_up` / `<contrast>_down` |
| `"Length"` | Number of significant features per contrast |
| `"Ranks"` | Named numeric vector of stats (for GSEA input) |
| `"FC"` | Named numeric vector of log2 fold changes |

``` r
DE(htgs, Report = "Length", padj_cut = 0.1)
## sh_EWS_vs_Ctrl     Hrs_2_vs_0     Hrs_6_vs_0     shEWS.Hrs2     shEWS.Hrs6 
##             52             40             63              5              1
```

``` r
DE(htgs,
   Report   = "Features",
   contrasts = "sh_EWS_vs_Ctrl",
   # limit to top 5
   Topn = 5,
   padj_cut = 0.1) 
## $sh_EWS_vs_Ctrl
## [1] "C1R"   "CCL2"  "C1S"   "MMP3"  "HMGB2"
```

``` r
# Ranks are sorted by the stat column (used as GSEA input)
DE(htgs,
   Report    = "Ranks",
   contrasts  = "sh_EWS_vs_Ctrl",
   
   # limit to top 5
   Topn = 5,
   Rank_name = "Feature",
   padj_cut  = 1) 
## $sh_EWS_vs_Ctrl
##     MMP3    HMGB2     RAC1    HMGB1    MEF2D 
## 8.804869 7.924608 6.377267 6.256207 6.249395
```

``` r
DE(htgs,
   Report    = "FC",
   contrasts  = "sh_EWS_vs_Ctrl",
   # limit to top 5
   Topn = 5,
   Rank_name = "Feature",
   padj_cut  = 1) 
## $sh_EWS_vs_Ctrl
## # A tibble: 5 × 3
##   Feature log2FoldChange ensembl_id     
##   <chr>            <dbl> <chr>          
## 1 C1R             -0.799 ENSG00000159403
## 2 C1R             -0.799 ENSG00000288512
## 3 CCL2            -0.945 ENSG00000108691
## 4 C1S             -0.594 ENSG00000182326
## 5 MMP3             0.578 ENSG00000149968
```

------------------------------------------------------------------------

## 4. `Output_DE_()` — Low-level Access

`Output_DE_()` returns the raw list of DE tables stored in the object,
without further filtering or formatting. Useful when you need to pass
results downstream in list form.

``` r
raw_DE <- Output_DE_(htgs, padj_cut = 1, as_list = TRUE)
names(raw_DE)
## [1] "sh_EWS_vs_Ctrl" "Hrs_2_vs_0"     "Hrs_6_vs_0"     "shEWS.Hrs2"    
## [5] "shEWS.Hrs6"
raw_DE[["sh_EWS_vs_Ctrl"]] |> head(3)
## # A tibble: 3 × 11
##   Feature contrast_dir  baseMean log2FoldChange    FC  stat   pvalue     padj  lfcSE
##   <chr>   <chr>            <dbl>          <dbl> <dbl> <dbl>    <dbl>    <dbl>  <dbl>
## 1 C1R     sh_EWS_vs_Ct…    5967.         -0.799 0.575 -15.1 1.16e-51 1.99e-49 0.0427
## 2 C1R     sh_EWS_vs_Ct…    5967.         -0.799 0.575 -15.1 1.16e-51 1.99e-49 0.0427
## 3 CCL2    sh_EWS_vs_Ct…    3119.         -0.945 0.519 -14.8 1.61e-49 1.37e-47 0.0474
## # ℹ 2 more variables: ensembl_id <chr>, significant <chr>
```

------------------------------------------------------------------------

## 5. `DExps()` — Expression Data Merged with Metadata

`DExps()` joins normalized expression values with sample metadata,
making it easy to build plots or run downstream models without extra
data-wrangling steps.

``` r
DExps(htgs,
      hotList    = c("IL6", "CSF2"),
      coldata_ids = c("Hrs", "sh"),
      Query_set  = TRUE) |> head()
##                  CSF2       IL6 Hrs   sh
## shCON_0hrs_1 4.622126  8.560798   0 Ctrl
## shCON_0hrs_2 4.662570  8.475861   0 Ctrl
## shCON_2hrs_1 4.891036 10.760152   2 Ctrl
## shCON_2hrs_2 4.837244 10.659633   2 Ctrl
## shCON_6hrs_1 5.202213 10.185363   6 Ctrl
## shCON_6hrs_2 5.044985 10.142222   6 Ctrl
```

------------------------------------------------------------------------

## 6. `DECoefs()` — Per-Feature Coefficients

`DECoefs()` extracts the model coefficients (one column per contrast)
for a set of features of interest.

``` r
DECoefs(htgs, hotList = c("CSF2", "IL6", "CCL2"))
##      sh_EWS_vs_Ctrl_log2FoldChange sh_EWS_vs_Ctrl_stat Hrs_2_vs_0_log2FoldChange
## CCL2                    -0.9454428          -14.793754                 0.8075294
## IL6                     -0.4620370           -2.410107                 2.3043826
## CSF2                     0.7835810            1.647501                 0.9656573
##      Hrs_2_vs_0_stat Hrs_6_vs_0_log2FoldChange Hrs_6_vs_0_stat
## CCL2       13.497033                  1.251154       20.809222
## IL6        26.682272                  1.732941       19.975903
## CSF2        2.937682                  2.146082        4.987054
##      shEWS.Hrs2_log2FoldChange shEWS.Hrs2_stat shEWS.Hrs6_log2FoldChange
## CCL2                -0.1065631      -1.1934152                0.04893024
## IL6                 -0.2037340      -3.5671775               -0.09350500
## CSF2                 0.1212042       0.4370389                0.12280155
##      shEWS.Hrs6_stat
## CCL2       0.8812305
## IL6       -2.4798332
## CSF2       0.6723109
```

------------------------------------------------------------------------

## 7. Updating the Mapper

The mapper table can be replaced or augmented at any time:

``` r
current_mapper <- Mapper_(htgs)
head(current_mapper)
## # A tibble: 6 × 2
##   Feature ensembl_id     
##   <chr>   <chr>          
## 1 AGER    ENSG00000204305
## 2 AGER    ENSG00000231268
## 3 AGER    ENSG00000229058
## 4 AGER    ENSG00000237405
## 5 AGER    ENSG00000234729
## 6 AGER    ENSG00000206320

# Add a synthetic annotation column
updated_mapper <- current_mapper |>
  dplyr::mutate(custom_note = paste0("gene:", .data$Feature))

Mapper_(htgs) <- updated_mapper
Mapper_(htgs) |> head()
## # A tibble: 6 × 3
##   Feature ensembl_id      custom_note
##   <chr>   <chr>           <chr>      
## 1 AGER    ENSG00000204305 gene:AGER  
## 2 AGER    ENSG00000231268 gene:AGER  
## 3 AGER    ENSG00000229058 gene:AGER  
## 4 AGER    ENSG00000237405 gene:AGER  
## 5 AGER    ENSG00000234729 gene:AGER  
## 6 AGER    ENSG00000206320 gene:AGER
```

------------------------------------------------------------------------

## 8. Combining Multiple Hotgenes Objects

When you have objects from different experiments or platforms, store
them in a named list. The Shiny app and several reporting functions
accept this list directly.

``` r
# Load the limma-based example as a second object
fit_Hotgenes <- readRDS(
  system.file("extdata", "fit_Hotgenes.RDS",
              package = "Hotgenes",
              mustWork = TRUE)
)

Hotgenes_list <- list(
  DESeq2_experiment = htgs,
  limma_experiment  = fit_Hotgenes
)

# Access individual objects from the list
Hotgenes_list[["limma_experiment"]]
## class: Hotgenes 
## Original class/package:  EList/limma
## 
## Differential expression (default thresholds): 
## |contrast       | total|
## |:--------------|-----:|
## |Hrs_2_vs_0     |    45|
## |Hrs_6_vs_0     |    71|
## |sh_EWS_vs_Ctrl |    51|
## |shEWS.Hrs2     |     8|
## |shEWS.Hrs6     |     1|
## 
## Available feature mapping:  Feature, ensembl_id 
## ExpressionSlots:  logCPM, counts 
## Total auxiliary assays:  0 
## Total samples:  12
```

------------------------------------------------------------------------

## Summary of Key Accessors

| Function              | Purpose                                  |
|-----------------------|------------------------------------------|
| `coldata_()`          | Sample metadata data.frame               |
| `Normalized_Data_()`  | Named list of expression matrices        |
| `ExpressionSlots_()`  | Names of available expression slots      |
| `SampleIDs_()`        | Sample identifiers                       |
| `Features_()`         | Feature (gene/protein) names             |
| `contrasts_()`        | Available contrast names                 |
| `Mapper_()`           | Feature alias mapping table              |
| `designMatrix_()`     | Model design matrix                      |
| `contrastMatrix_()`   | Contrast matrix                          |
| `O_()`                | Original DE object                       |
| `auxiliary_assays_()` | Auxiliary sample-level data              |
| `DE()`                | Query DE results with flexible reporting |
| `Output_DE_()`        | Raw DE list (low-level)                  |
| `DExps()`             | Expression + metadata joined table       |
| `DECoefs()`           | Per-feature model coefficients           |
