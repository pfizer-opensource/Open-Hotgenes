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
##  [1] "shCON_0hrs_1" "shCON_0hrs_2" "shCON_2hrs_1" "shCON_2hrs_2" "shCON_6hrs_1" "shCON_6hrs_2"
##  [7] "shEWS_0hrs_1" "shEWS_0hrs_2" "shEWS_2hrs_1" "shEWS_2hrs_2" "shEWS_6hrs_1" "shEWS_6hrs_2"
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
##    Feature contrast_dir        baseMean log2FoldChange    FC   stat   pvalue     padj  lfcSE ensembl_id  
##    <chr>   <chr>                  <dbl>          <dbl> <dbl>  <dbl>    <dbl>    <dbl>  <dbl> <chr>       
##  1 C1R     sh_EWS_vs_Ctrl_down    5967.         -0.819 0.567 -48.7  1.16e-51 1.99e-49 0.0547 ENSG0000015…
##  2 C1R     sh_EWS_vs_Ctrl_down    5967.         -0.819 0.567 -48.7  1.16e-51 1.99e-49 0.0547 ENSG0000028…
##  3 CCL2    sh_EWS_vs_Ctrl_down    3119.         -0.959 0.514 -46.9  1.61e-49 1.37e-47 0.0657 ENSG0000010…
##  4 C1S     sh_EWS_vs_Ctrl_down    3514.         -0.602 0.659 -29.1  1.36e-31 7.76e-30 0.0525 ENSG0000018…
##  5 MMP3    sh_EWS_vs_Ctrl_up      3986.          0.597 1.51   16.3  1.31e-18 5.60e-17 0.0698 ENSG0000014…
##  6 HMGB2   sh_EWS_vs_Ctrl_up      1118.          0.606 1.52   13.1  2.29e-15 7.83e-14 0.0794 ENSG0000016…
##  7 STAT1   sh_EWS_vs_Ctrl_down    4448.         -0.349 0.785 -10.9  4.52e-13 1.29e-11 0.0493 ENSG0000011…
##  8 CXCL6   sh_EWS_vs_Ctrl_down     130.         -1.77  0.294 -10.3  2.27e-12 5.54e-11 0.267  ENSG0000012…
##  9 IL1R1   sh_EWS_vs_Ctrl_down    1519.         -0.457 0.728  -9.21 2.88e-11 6.16e-10 0.0720 ENSG0000011…
## 10 RAC1    sh_EWS_vs_Ctrl_up      9136.          0.267 1.20    8.47 1.80e-10 3.43e- 9 0.0427 ENSG0000013…
## # ℹ 1 more variable: significant <chr>
## 
## $Hrs_2_vs_0
## # A tibble: 10 × 11
##    Feature contrast_dir  baseMean log2FoldChange    FC  stat    pvalue      padj  lfcSE ensembl_id     
##    <chr>   <chr>            <dbl>          <dbl> <dbl> <dbl>     <dbl>     <dbl>  <dbl> <chr>          
##  1 IL6     Hrs_2_vs_0_up     914.          2.49   5.60 154.  7.56e-157 1.93e-154 0.0936 ENSG00000136244
##  2 TNFAIP3 Hrs_2_vs_0_up     780.          2.12   4.33  98.1 5.75e-101 7.34e- 99 0.0997 ENSG00000118503
##  3 CXCL8   Hrs_2_vs_0_up    1409.          3.49  11.2   81.1 8.35e- 84 7.10e- 82 0.181  ENSG00000169429
##  4 IL11    Hrs_2_vs_0_up     659.          2.00   3.99  71.9 1.93e- 74 1.23e- 72 0.110  ENSG00000095752
##  5 PTGS2   Hrs_2_vs_0_up    1907.          1.12   2.18  51.8 3.41e- 54 1.74e- 52 0.0730 ENSG00000073756
##  6 SMAD7   Hrs_2_vs_0_up    3674.          0.933  1.91  48.1 1.98e- 50 8.40e- 49 0.0629 ENSG00000101665
##  7 CXCL1   Hrs_2_vs_0_up    1533.          2.08   4.24  46.4 1.08e- 48 3.92e- 47 0.144  ENSG00000163739
##  8 CCL2    Hrs_2_vs_0_up    3119.          0.812  1.76  39.3 1.63e- 41 5.19e- 40 0.0607 ENSG00000108691
##  9 JUN     Hrs_2_vs_0_up    1619.          0.767  1.70  19.0 3.84e- 21 1.09e- 19 0.0823 ENSG00000177606
## 10 CXCL3   Hrs_2_vs_0_up     132.          2.00   3.99  14.5 1.37e- 16 3.50e- 15 0.250  ENSG00000163734
## # ℹ 1 more variable: significant <chr>
## 
## $Hrs_6_vs_0
## # A tibble: 10 × 11
##    Feature contrast_dir  baseMean log2FoldChange    FC  stat   pvalue     padj  lfcSE ensembl_id     
##    <chr>   <chr>            <dbl>          <dbl> <dbl> <dbl>    <dbl>    <dbl>  <dbl> <chr>          
##  1 CCL2    Hrs_6_vs_0_up    3119.          1.24   2.36  93.1 3.57e-96 7.14e-94 0.0600 ENSG00000108691
##  2 PTGS2   Hrs_6_vs_0_up    1907.          1.44   2.72  87.9 1.35e-90 1.35e-88 0.0720 ENSG00000073756
##  3 IL6     Hrs_6_vs_0_up     914.          1.89   3.70  86.2 8.93e-89 5.95e-87 0.0952 ENSG00000136244
##  4 CXCL8   Hrs_6_vs_0_up    1409.          3.55  11.7   84.7 4.19e-87 2.09e-85 0.181  ENSG00000169429
##  5 TNFAIP3 Hrs_6_vs_0_up     780.          1.79   3.47  69.6 6.41e-72 2.56e-70 0.101  ENSG00000118503
##  6 CXCL1   Hrs_6_vs_0_up    1533.          2.33   5.03  58.3 1.46e-60 4.86e-59 0.143  ENSG00000163739
##  7 RELB    Hrs_6_vs_0_up     436.          1.90   3.72  55.6 9.66e-58 2.76e-56 0.120  ENSG00000104856
##  8 IL11    Hrs_6_vs_0_up     659.          1.44   2.71  36.5 1.21e-38 3.03e-37 0.112  ENSG00000095752
##  9 BIRC2   Hrs_6_vs_0_up    5074.          0.727  1.66  31.5 1.42e-33 3.15e-32 0.0610 ENSG00000110330
## 10 NFKB1   Hrs_6_vs_0_up     705.          0.985  1.98  24.1 3.94e-26 7.88e-25 0.0949 ENSG00000109320
## # ℹ 1 more variable: significant <chr>
## 
## $shEWS.Hrs2
## # A tibble: 5 × 11
##   Feature contrast_dir  baseMean log2FoldChange    FC  stat  pvalue   padj   lfcSE ensembl_id significant
##   <chr>   <chr>            <dbl>          <dbl> <dbl> <dbl>   <dbl>  <dbl>   <dbl> <chr>      <chr>      
## 1 TNFAIP3 shEWS.Hrs2_d…     780.      -0.495    0.710 -1.93 4.62e-5 0.0118 0.153   ENSG00000… *          
## 2 IL6     shEWS.Hrs2_d…     914.      -0.400    0.758 -1.48 3.61e-4 0.0330 0.156   ENSG00000… *          
## 3 TGFBR1  shEWS.Hrs2_d…    3125.      -0.214    0.862 -1.48 3.88e-4 0.0330 0.0866  ENSG00000… *          
## 4 CXCL3   shEWS.Hrs2_d…     132.      -0.773    0.585 -1.03 1.46e-3 0.0932 0.416   ENSG00000… *          
## 5 CXCL8   shEWS.Hrs2_d…    1409.      -0.000781 0.999 -1.01 1.92e-3 0.0980 0.00815 ENSG00000… *          
## 
## $shEWS.Hrs6
## # A tibble: 1 × 11
##   Feature contrast_dir  baseMean log2FoldChange    FC  stat   pvalue   padj lfcSE ensembl_id  significant
##   <chr>   <chr>            <dbl>          <dbl> <dbl> <dbl>    <dbl>  <dbl> <dbl> <chr>       <chr>      
## 1 FOS     shEWS.Hrs6_up     263.          0.656  1.58  1.04 0.000359 0.0916 0.262 ENSG000001… *
```

### Select specific contrasts

``` r
# List available contrasts
contrasts_(htgs)
## [1] "sh_EWS_vs_Ctrl" "Hrs_2_vs_0"     "Hrs_6_vs_0"     "shEWS.Hrs2"     "shEWS.Hrs6"
```

``` r
DE(htgs,
   contrasts      = "sh_EWS_vs_Ctrl",
   padj_cut       = 0.05,
   .log2FoldChange = 1) 
## $sh_EWS_vs_Ctrl
## # A tibble: 3 × 11
##   Feature contrast_dir        baseMean log2FoldChange    FC   stat   pvalue     padj lfcSE ensembl_id    
##   <chr>   <chr>                  <dbl>          <dbl> <dbl>  <dbl>    <dbl>    <dbl> <dbl> <chr>         
## 1 CXCL6   sh_EWS_vs_Ctrl_down    130.           -1.77 0.294 -10.3  2.27e-12 5.54e-11 0.267 ENSG000001248…
## 2 C3      sh_EWS_vs_Ctrl_down     49.3          -1.62 0.325  -4.78 1.76e- 6 1.67e- 5 0.391 ENSG000001257…
## 3 MAP3K9  sh_EWS_vs_Ctrl_up       26.8           1.59 3.02    3.06 1.18e- 4 8.78e- 4 0.528 ENSG000000064…
## # ℹ 1 more variable: significant <chr>
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
##   Feature contrast_dir        baseMean log2FoldChange    FC    stat   pvalue     padj  lfcSE ensembl_id  
##   <chr>   <chr>                  <dbl>          <dbl> <dbl>   <dbl>    <dbl>    <dbl>  <dbl> <chr>       
## 1 CCL2    sh_EWS_vs_Ctrl_down   3119.          -0.959 0.514 -46.9   1.61e-49 1.37e-47 0.0657 ENSG0000010…
## 2 IL6     sh_EWS_vs_Ctrl_down    914.          -0.218 0.860  -1.23  1.59e- 2 5.83e- 2 0.108  ENSG0000013…
## 3 CSF2    sh_EWS_vs_Ctrl_up       38.4          0.100 1.07    0.624 9.95e- 2 2.37e- 1 0.216  ENSG0000016…
## # ℹ 1 more variable: significant <chr>
## 
## $Hrs_2_vs_0
## # A tibble: 3 × 11
##   Feature contrast_dir  baseMean log2FoldChange    FC   stat    pvalue      padj  lfcSE ensembl_id     
##   <chr>   <chr>            <dbl>          <dbl> <dbl>  <dbl>     <dbl>     <dbl>  <dbl> <chr>          
## 1 IL6     Hrs_2_vs_0_up    914.           2.49   5.60 154.   7.56e-157 1.93e-154 0.0936 ENSG00000136244
## 2 CCL2    Hrs_2_vs_0_up   3119.           0.812  1.76  39.3  1.63e- 41 5.19e- 40 0.0607 ENSG00000108691
## 3 CSF2    Hrs_2_vs_0_up     38.4          1.30   2.46   1.57 3.31e-  3 2.72e-  2 0.661  ENSG00000164400
## # ℹ 1 more variable: significant <chr>
## 
## $Hrs_6_vs_0
## # A tibble: 3 × 11
##   Feature contrast_dir  baseMean log2FoldChange    FC  stat   pvalue     padj  lfcSE ensembl_id     
##   <chr>   <chr>            <dbl>          <dbl> <dbl> <dbl>    <dbl>    <dbl>  <dbl> <chr>          
## 1 CCL2    Hrs_6_vs_0_up   3119.            1.24  2.36 93.1  3.57e-96 7.14e-94 0.0600 ENSG00000108691
## 2 IL6     Hrs_6_vs_0_up    914.            1.89  3.70 86.2  8.93e-89 5.95e-87 0.0952 ENSG00000136244
## 3 CSF2    Hrs_6_vs_0_up     38.4           2.60  6.05  5.39 6.13e- 7 4.09e- 6 0.572  ENSG00000164400
## # ℹ 1 more variable: significant <chr>
## 
## $shEWS.Hrs2
## # A tibble: 3 × 11
##   Feature contrast_dir    baseMean log2FoldChange    FC     stat   pvalue   padj   lfcSE ensembl_id     
##   <chr>   <chr>              <dbl>          <dbl> <dbl>    <dbl>    <dbl>  <dbl>   <dbl> <chr>          
## 1 IL6     shEWS.Hrs2_down    914.      -0.400     0.758 -1.48    0.000361 0.0330 0.156   ENSG00000136244
## 2 CCL2    shEWS.Hrs2_down   3119.      -0.000981  0.999 -0.00746 0.233    0.983  0.00815 ENSG00000108691
## 3 CSF2    shEWS.Hrs2_up       38.4      0.0000373 1.00   0.00746 0.662    0.983  0.00810 ENSG00000164400
## # ℹ 1 more variable: significant <chr>
## 
## $shEWS.Hrs6
## # A tibble: 3 × 11
##   Feature contrast_dir baseMean log2FoldChange    FC     stat pvalue  padj   lfcSE ensembl_id significant
##   <chr>   <chr>           <dbl>          <dbl> <dbl>    <dbl>  <dbl> <dbl>   <dbl> <chr>      <chr>      
## 1 IL6     shEWS.Hrs6_…    914.    -0.0000187   1.000 -0.182   0.0131 0.658 0.00144 ENSG00000… ""         
## 2 CCL2    shEWS.Hrs6_…   3119.    -0.00000324  1.000 -0.00120 0.378  0.997 0.00144 ENSG00000… ""         
## 3 CSF2    shEWS.Hrs6_…     38.4    0.000000956 1.00   0.00120 0.501  0.997 0.00144 ENSG00000… ""
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
# Ranks are sorted (descending order) by the stat column (used as GSEA input)
DE(htgs,
   Report    = "Ranks",
   contrasts  = "sh_EWS_vs_Ctrl",
   
   # limit to top 5
   Topn = 5,
   Rank_name = "Feature",
   padj_cut  = 1) 
## $sh_EWS_vs_Ctrl
##      MMP3     HMGB2      RAC1     HMGB1     MEF2D 
## 16.251774 13.106389  8.465307  8.193452  8.193452
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
## 1 C1R             -0.819 ENSG00000159403
## 2 C1R             -0.819 ENSG00000288512
## 3 CCL2            -0.959 ENSG00000108691
## 4 C1S             -0.602 ENSG00000182326
## 5 MMP3             0.597 ENSG00000149968
```

------------------------------------------------------------------------

## 4. `DExps()` — Expression Data Merged with Metadata

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

## 5. `DECoefs()` — Per-Feature Coefficients

`DECoefs()` extracts the model coefficients (one column per contrast)
for a set of features of interest.

``` r
DECoefs(htgs, hotList = c("CSF2", "IL6", "CCL2"))
##      sh_EWS_vs_Ctrl_log2FoldChange sh_EWS_vs_Ctrl_stat Hrs_2_vs_0_log2FoldChange Hrs_2_vs_0_stat
## CCL2                    -0.9588353         -46.8619299                 0.8115152       39.284892
## IL6                     -0.2175810          -1.2341166                 2.4864430      153.714945
## CSF2                     0.1000838           0.6243459                 1.3006816        1.565419
##      Hrs_6_vs_0_log2FoldChange Hrs_6_vs_0_stat shEWS.Hrs2_log2FoldChange shEWS.Hrs2_stat
## CCL2                  1.240551       93.146191             -9.812312e-04    -0.007460348
## IL6                   1.887527       86.225454             -3.995497e-01    -1.481969498
## CSF2                  2.597729        5.388582              3.725334e-05     0.007460348
##      shEWS.Hrs6_log2FoldChange shEWS.Hrs6_stat
## CCL2             -3.236547e-06    -0.001199919
## IL6              -1.868006e-05    -0.181614602
## CSF2              9.562535e-07     0.001199919
```

------------------------------------------------------------------------

## 6. Updating the Mapper

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

## 7. Combining Multiple Hotgenes Objects

When you have objects from different experiments or platforms, store
them in a named list. The Shiny app can accept this list directly.

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
| `DExps()`             | Expression + metadata joined table       |
| `DECoefs()`           | Per-feature model coefficients           |
