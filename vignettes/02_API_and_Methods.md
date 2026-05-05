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
DE(htgs, padj_cut = 0.05, .log2FoldChange = 1) |> head()
## $sh_EWS_vs_Ctrl
## # A tibble: 2 × 11
##   Feature contrast_dir  baseMean log2FoldChange    FC  stat   pvalue     padj lfcSE ensembl_id significant
##   <chr>   <chr>            <dbl>          <dbl> <dbl> <dbl>    <dbl>    <dbl> <dbl> <chr>      <chr>      
## 1 CXCL6   sh_EWS_vs_Ct…    130.           -1.65 0.318 -7.02 2.27e-12 5.54e-11 0.107 ENSG00000… *          
## 2 C3      sh_EWS_vs_Ct…     49.3          -1.23 0.428 -4.78 1.76e- 6 1.67e- 5 0.138 ENSG00000… *          
## 
## $Hrs_2_vs_0
## # A tibble: 8 × 11
##   Feature contrast_dir  baseMean log2FoldChange    FC  stat    pvalue      padj  lfcSE ensembl_id     
##   <chr>   <chr>            <dbl>          <dbl> <dbl> <dbl>     <dbl>     <dbl>  <dbl> <chr>          
## 1 IL6     Hrs_2_vs_0_up     914.           2.30  4.94 26.7  7.56e-157 1.93e-154 0.0725 ENSG00000136244
## 2 TNFAIP3 Hrs_2_vs_0_up     780.           1.92  3.79 21.3  5.75e-101 7.34e- 99 0.0750 ENSG00000118503
## 3 CXCL8   Hrs_2_vs_0_up    1409.           2.88  7.35 19.4  8.35e- 84 7.10e- 82 0.120  ENSG00000169429
## 4 IL11    Hrs_2_vs_0_up     659.           1.80  3.49 18.3  1.93e- 74 1.23e- 72 0.0821 ENSG00000095752
## 5 PTGS2   Hrs_2_vs_0_up    1907.           1.12  2.17 15.5  3.41e- 54 1.74e- 52 0.0595 ENSG00000073756
## 6 CXCL1   Hrs_2_vs_0_up    1533.           1.87  3.67 14.7  1.08e- 48 3.92e- 47 0.101  ENSG00000163739
## 7 CXCL3   Hrs_2_vs_0_up     132.           1.30  2.47  8.27 1.37e- 16 3.50e- 15 0.143  ENSG00000163734
## 8 CXCL2   Hrs_2_vs_0_up     128.           1.23  2.34  7.40 1.34e- 13 3.11e- 12 0.141  ENSG00000081041
## # ℹ 1 more variable: significant <chr>
## 
## $Hrs_6_vs_0
## # A tibble: 16 × 11
##    Feature contrast_dir    baseMean log2FoldChange    FC  stat   pvalue     padj  lfcSE ensembl_id     
##    <chr>   <chr>              <dbl>          <dbl> <dbl> <dbl>    <dbl>    <dbl>  <dbl> <chr>          
##  1 CCL2    Hrs_6_vs_0_up     3119.            1.25 2.38  20.8  3.57e-96 7.14e-94 0.0507 ENSG00000108691
##  2 PTGS2   Hrs_6_vs_0_up     1907.            1.42 2.68  20.2  1.35e-90 1.35e-88 0.0576 ENSG00000073756
##  3 IL6     Hrs_6_vs_0_up      914.            1.73 3.32  20.0  8.93e-89 5.95e-87 0.0724 ENSG00000136244
##  4 CXCL8   Hrs_6_vs_0_up     1409.            3.04 8.24  19.8  4.19e-87 2.09e-85 0.122  ENSG00000169429
##  5 TNFAIP3 Hrs_6_vs_0_up      780.            1.67 3.18  17.9  6.41e-72 2.56e-70 0.0741 ENSG00000118503
##  6 CXCL1   Hrs_6_vs_0_up     1533.            2.23 4.69  16.4  1.46e-60 4.86e-59 0.101  ENSG00000163739
##  7 RELB    Hrs_6_vs_0_up      436.            1.83 3.55  16.0  9.66e-58 2.76e-56 0.0874 ENSG00000104856
##  8 IL11    Hrs_6_vs_0_up      659.            1.26 2.40  13.0  1.21e-38 3.03e-37 0.0823 ENSG00000095752
##  9 FOS     Hrs_6_vs_0_down    263.           -1.25 0.421 -9.69 3.48e-22 5.80e-21 0.112  ENSG00000170345
## 10 CXCL6   Hrs_6_vs_0_up      130.            1.63 3.10   9.60 7.94e-22 1.22e-20 0.141  ENSG00000124875
## 11 CXCL3   Hrs_6_vs_0_up      132.            1.52 2.87   8.68 3.89e-18 5.55e-17 0.146  ENSG00000163734
## 12 CCL7    Hrs_6_vs_0_up       63.2           1.80 3.48   6.96 3.47e-12 3.66e-11 0.182  ENSG00000108688
## 13 CXCR4   Hrs_6_vs_0_up       88.0           1.27 2.41   6.28 3.37e-10 2.93e- 9 0.161  ENSG00000121966
## 14 CXCL5   Hrs_6_vs_0_up       47.5           1.61 3.06   5.52 3.48e- 8 2.40e- 7 0.207  ENSG00000163735
## 15 CSF2    Hrs_6_vs_0_up       38.4           2.15 4.43   4.99 6.13e- 7 4.09e- 6 0.232  ENSG00000164400
## 16 MAP2K6  Hrs_6_vs_0_down     61.2          -1.05 0.485 -3.39 6.88e- 4 3.13e- 3 0.189  ENSG00000108984
## # ℹ 1 more variable: significant <chr>
## 
## $shEWS.Hrs2
## # A tibble: 0 × 11
## # ℹ 11 variables: Feature <chr>, contrast_dir <chr>, baseMean <dbl>, log2FoldChange <dbl>, FC <dbl>,
## #   stat <dbl>, pvalue <dbl>, padj <dbl>, lfcSE <dbl>, ensembl_id <chr>, significant <chr>
## 
## $shEWS.Hrs6
## # A tibble: 0 × 11
## # ℹ 11 variables: Feature <chr>, contrast_dir <chr>, baseMean <dbl>, log2FoldChange <dbl>, FC <dbl>,
## #   stat <dbl>, pvalue <dbl>, padj <dbl>, lfcSE <dbl>, ensembl_id <chr>, significant <chr>
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
   .log2FoldChange = 1) |> head()
## $sh_EWS_vs_Ctrl
## # A tibble: 2 × 11
##   Feature contrast_dir  baseMean log2FoldChange    FC  stat   pvalue     padj lfcSE ensembl_id significant
##   <chr>   <chr>            <dbl>          <dbl> <dbl> <dbl>    <dbl>    <dbl> <dbl> <chr>      <chr>      
## 1 CXCL6   sh_EWS_vs_Ct…    130.           -1.65 0.318 -7.02 2.27e-12 5.54e-11 0.107 ENSG00000… *          
## 2 C3      sh_EWS_vs_Ct…     49.3          -1.23 0.428 -4.78 1.76e- 6 1.67e- 5 0.138 ENSG00000… *
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
##   Feature contrast_dir        baseMean log2FoldChange    FC   stat   pvalue     padj  lfcSE ensembl_id    
##   <chr>   <chr>                  <dbl>          <dbl> <dbl>  <dbl>    <dbl>    <dbl>  <dbl> <chr>         
## 1 CCL2    sh_EWS_vs_Ctrl_down   3119.          -0.945 0.519 -14.8  1.61e-49 1.37e-47 0.0474 ENSG000001086…
## 2 IL6     sh_EWS_vs_Ctrl_down    914.          -0.462 0.726  -2.41 1.59e- 2 5.83e- 2 0.0596 ENSG000001362…
## 3 CSF2    sh_EWS_vs_Ctrl_up       38.4          0.784 1.72    1.65 9.95e- 2 2.37e- 1 0.147  ENSG000001644…
## # ℹ 1 more variable: significant <chr>
## 
## $Hrs_2_vs_0
## # A tibble: 3 × 11
##   Feature contrast_dir  baseMean log2FoldChange    FC  stat    pvalue      padj  lfcSE ensembl_id     
##   <chr>   <chr>            <dbl>          <dbl> <dbl> <dbl>     <dbl>     <dbl>  <dbl> <chr>          
## 1 IL6     Hrs_2_vs_0_up    914.           2.30   4.94 26.7  7.56e-157 1.93e-154 0.0725 ENSG00000136244
## 2 CCL2    Hrs_2_vs_0_up   3119.           0.808  1.75 13.5  1.63e- 41 5.19e- 40 0.0522 ENSG00000108691
## 3 CSF2    Hrs_2_vs_0_up     38.4          0.966  1.95  2.94 3.31e-  3 2.72e-  2 0.210  ENSG00000164400
## # ℹ 1 more variable: significant <chr>
## 
## $Hrs_6_vs_0
## # A tibble: 3 × 11
##   Feature contrast_dir baseMean log2FoldChange    FC  stat   pvalue     padj  lfcSE ensembl_id significant
##   <chr>   <chr>           <dbl>          <dbl> <dbl> <dbl>    <dbl>    <dbl>  <dbl> <chr>      <chr>      
## 1 CCL2    Hrs_6_vs_0_…   3119.            1.25  2.38 20.8  3.57e-96 7.14e-94 0.0507 ENSG00000… *          
## 2 IL6     Hrs_6_vs_0_…    914.            1.73  3.32 20.0  8.93e-89 5.95e-87 0.0724 ENSG00000… *          
## 3 CSF2    Hrs_6_vs_0_…     38.4           2.15  4.43  4.99 6.13e- 7 4.09e- 6 0.232  ENSG00000… *          
## 
## $shEWS.Hrs2
## # A tibble: 3 × 11
##   Feature contrast_dir   baseMean log2FoldChange    FC   stat  pvalue   padj  lfcSE ensembl_id significant
##   <chr>   <chr>             <dbl>          <dbl> <dbl>  <dbl>   <dbl>  <dbl>  <dbl> <chr>      <chr>      
## 1 IL6     shEWS.Hrs2_do…    914.          -0.204 0.868 -3.57  3.61e-4 0.0330 0.0672 ENSG00000… ""         
## 2 CCL2    shEWS.Hrs2_do…   3119.          -0.107 0.929 -1.19  2.33e-1 0.983  0.0616 ENSG00000… ""         
## 3 CSF2    shEWS.Hrs2_up      38.4          0.121 1.09   0.437 6.62e-1 0.983  0.0421 ENSG00000… ""         
## 
## $shEWS.Hrs6
## # A tibble: 3 × 11
##   Feature contrast_dir    baseMean log2FoldChange    FC   stat pvalue  padj  lfcSE ensembl_id  significant
##   <chr>   <chr>              <dbl>          <dbl> <dbl>  <dbl>  <dbl> <dbl>  <dbl> <chr>       <chr>      
## 1 IL6     shEWS.Hrs6_down    914.         -0.0935 0.937 -2.48  0.0131 0.658 0.0584 ENSG000001… ""         
## 2 CCL2    shEWS.Hrs6_up     3119.          0.0489 1.03   0.881 0.378  0.997 0.0561 ENSG000001… ""         
## 3 CSF2    shEWS.Hrs6_up       38.4         0.123  1.09   0.672 0.501  0.997 0.0317 ENSG000001… ""
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
   padj_cut = 0.1) |> head()
## $sh_EWS_vs_Ctrl
##  [1] "C1R"      "CCL2"     "C1S"      "MMP3"     "HMGB2"    "STAT1"    "CXCL6"    "IL1R1"    "RAC1"    
## [10] "HMGB1"    "MEF2D"    "JUN"      "GNB1"     "HIF1A"    "PTGS1"    "TUBB"     "BCL2L1"   "C3"      
## [19] "HMGN1"    "ROCK2"    "BIRC2"    "CLTC"     "MAP3K9"   "CEBPB"    "CSF1"     "CXCL3"    "MEF2A"   
## [28] "MASP1"    "MAPKAPK2" "PTGER3"   "SMAD7"    "HPRT1"    "RHOA"     "GNAQ"     "TCF4"     "TRAF2"   
## [37] "MX2"      "HRAS"     "RAF1"     "CXCL2"    "RIPK2"    "SHC1"     "MAPK1"    "CD40"     "FLT1"    
## [46] "IL6"      "NR3C1"    "LY96"     "CCL7"     "CFD"      "NFKB1"    "MAPK14"
```

``` r
# Ranks are sorted by the stat column (used as GSEA input)
DE(htgs,
   Report    = "Ranks",
   contrasts  = "sh_EWS_vs_Ctrl",
   Rank_name = "Feature",
   padj_cut  = 1) |> head()
## $sh_EWS_vs_Ctrl
##          MMP3         HMGB2          RAC1         HMGB1         MEF2D           JUN          GNB1 
##   8.804868849   7.924607722   6.377266648   6.256206649   6.249395109   5.840344544   5.749781072 
##         PTGS1          TUBB        BCL2L1         HMGN1         ROCK2         BIRC2          CLTC 
##   5.421186672   4.987660002   4.923589326   4.545913648   4.533313571   4.372898503   4.033724594 
##        MAP3K9         CXCL3         MEF2A        PTGER3      MAPKAPK2         HPRT1         SMAD7 
##   3.850076949   3.582782303   3.571696499   3.200366989   3.189255010   3.187594235   3.182802516 
##          RHOA          GNAQ          TCF4         TRAF2          HRAS          RAF1         CXCL2 
##   3.121686791   3.072143484   3.047615199   2.773197690   2.686051482   2.669045161   2.629158783 
##         RIPK2          SHC1         MAPK1          CD40         NR3C1          FLT1        MAPK14 
##   2.588689220   2.581941646   2.495721897   2.432595037   2.420764606   2.408188384   2.167608036 
##         PRKCA        LTB4R2         HSPB2         MAPK8          BCL6        MAP2K1         PTGS2 
##   2.125741450   2.090346433   2.034842705   2.033232169   2.007788679   1.986748183   1.973559554 
##         PLCB1          IL18          CFL1          GNAS        MAP3K5       TNFAIP3          CSF2 
##   1.966381213   1.873535913   1.776091942   1.762105575   1.691958207   1.671318830   1.647500550 
##          RELA         TGFB1         C3AR1         MEF2C         CDC42          AGER          CCR3 
##   1.644807323   1.506099888   1.412369257   1.345720779   1.326425941   1.313589034   1.255322961 
##          NOD1         CREB1         HDAC4       RPS6KA5          PGK1           MYC         CXCR4 
##   1.151482524   1.150424036   1.134380457   1.112521229   1.107297848   1.054431109   0.975577902 
##        MAP3K7        NFATC3         GAPDH        MAP2K6        TOLLIP        TBXA2R         PDGFA 
##   0.944440546   0.932119682   0.869583654   0.860028124   0.812850737   0.755614819   0.738878451 
##          MAFF         LIMK1          PTK2           CD4            C5        PTGER4          IRF5 
##   0.699018443   0.663112269   0.608641824   0.596682422   0.550849684   0.542046403   0.450286775 
##          GUSB          MAFG          ATF2        PTGER2         LTB4R        PTGER1         KEAP1 
##   0.369993351   0.329969739   0.317138296   0.303307640   0.295656729   0.263652465   0.262764676 
##         HSPB1          ELK1          OASL         FXYD2          AREG         IL12A          IL1B 
##   0.178962936   0.150118893   0.137302381   0.111550319   0.104278358   0.098285814   0.091926422 
##         OXER1         TGFB3         TGFB2         CXCL8        MAP2K4          NOS2          IRF3 
##   0.084396468   0.065136821   0.012037462   0.006409391  -0.022472649  -0.067861697  -0.168294243 
##          DAXX         MYD88        TWIST2        TGFBR1      PPP1R12B          MAFK           FOS 
##  -0.185261760  -0.185671730  -0.198595920  -0.270431958  -0.273608487  -0.274258945  -0.291351601 
##        ALOX15         IL12B          GRB2        IL1RAP           C4A         NLRP3         RIPK1 
##  -0.341430559  -0.345559872  -0.369719241  -0.376122683  -0.382127887  -0.388376547  -0.433518466 
##          IL15         CXCL5         IFIT3          ARG1         DDIT3          NOX1         PTGFR 
##  -0.443176576  -0.447594655  -0.450145072  -0.460057368  -0.501850995  -0.542581776  -0.554283208 
##          IL11        NFE2L2         IFIT1      MAPKAPK5        MAP3K1       PLA2G4A          TLR3 
##  -0.572604125  -0.606611605  -0.620080780  -0.625389854  -0.628078514  -0.656366226  -0.728130723 
##         IL17A         TRADD       RAPGEF2         STAT3  BORCS8-MEF2B          TLR4          IL6R 
##  -0.790217493  -0.800894996  -0.880768554  -0.947527413  -1.030454978  -1.059652572  -1.101442157 
##       TNFSF14          TLR2          TLR1           IL4          CD55         IFIT2           CFB 
##  -1.128538064  -1.136517693  -1.140050128  -1.207157325  -1.210286467  -1.263788109  -1.317486278 
##         PTGIR           MAX         STAT2          IRF7          TLR6        IL10RB         MAPK3 
##  -1.324225420  -1.344429982  -1.365072239  -1.385378120  -1.388137352  -1.412017722  -1.454135184 
##           MX1         CXCL1         MKNK1          OAS2          IRF1          RELB         IFI44 
##  -1.472299599  -1.638167897  -1.652418849  -1.802262738  -1.912021949  -1.925189762  -1.937948286 
##         NFKB1           CFD          CCL7          LY96           IL6           MX2         MASP1 
##  -2.178631940  -2.214359443  -2.254103775  -2.327393367  -2.410107240  -2.750364692  -3.223947475 
##          CSF1         CEBPB            C3         HIF1A         IL1R1         CXCL6         STAT1 
##  -3.756692946  -3.764880703  -4.779164610  -5.614461573  -6.652586493  -7.016910026  -7.238916793 
##           C1S          CCL2           C1R 
## -11.694360410 -14.793754290 -15.121750438
```

``` r
DE(htgs,
   Report    = "FC",
   contrasts  = "sh_EWS_vs_Ctrl",
   Rank_name = "Feature",
   padj_cut  = 1) |> head()
## $sh_EWS_vs_Ctrl
## # A tibble: 208 × 3
##    Feature log2FoldChange ensembl_id     
##    <chr>            <dbl> <chr>          
##  1 C1R             -0.799 ENSG00000159403
##  2 C1R             -0.799 ENSG00000288512
##  3 CCL2            -0.945 ENSG00000108691
##  4 C1S             -0.594 ENSG00000182326
##  5 MMP3             0.578 ENSG00000149968
##  6 HMGB2            0.621 ENSG00000164104
##  7 STAT1           -0.355 ENSG00000115415
##  8 CXCL6           -1.65  ENSG00000124875
##  9 IL1R1           -0.414 ENSG00000115594
## 10 RAC1             0.265 ENSG00000136238
## # ℹ 198 more rows
```

------------------------------------------------------------------------

## 4. `Output_DE_()` — Low-level Access

`Output_DE_()` returns the raw list of DE tables stored in the object,
without further filtering or formatting. Useful when you need to pass
results downstream in list form.

``` r
raw_DE <- Output_DE_(htgs, padj_cut = 1, as_list = TRUE)
names(raw_DE)
## [1] "sh_EWS_vs_Ctrl" "Hrs_2_vs_0"     "Hrs_6_vs_0"     "shEWS.Hrs2"     "shEWS.Hrs6"
raw_DE[["sh_EWS_vs_Ctrl"]] |> head(3)
## # A tibble: 3 × 11
##   Feature contrast_dir baseMean log2FoldChange    FC  stat   pvalue     padj  lfcSE ensembl_id significant
##   <chr>   <chr>           <dbl>          <dbl> <dbl> <dbl>    <dbl>    <dbl>  <dbl> <chr>      <chr>      
## 1 C1R     sh_EWS_vs_C…    5967.         -0.799 0.575 -15.1 1.16e-51 1.99e-49 0.0427 ENSG00000… *          
## 2 C1R     sh_EWS_vs_C…    5967.         -0.799 0.575 -15.1 1.16e-51 1.99e-49 0.0427 ENSG00000… *          
## 3 CCL2    sh_EWS_vs_C…    3119.         -0.945 0.519 -14.8 1.61e-49 1.37e-47 0.0474 ENSG00000… *
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
##      sh_EWS_vs_Ctrl_log2FoldChange sh_EWS_vs_Ctrl_stat Hrs_2_vs_0_log2FoldChange Hrs_2_vs_0_stat
## CCL2                    -0.9454428          -14.793754                 0.8075294       13.497033
## IL6                     -0.4620370           -2.410107                 2.3043826       26.682272
## CSF2                     0.7835810            1.647501                 0.9656573        2.937682
##      Hrs_6_vs_0_log2FoldChange Hrs_6_vs_0_stat shEWS.Hrs2_log2FoldChange shEWS.Hrs2_stat
## CCL2                  1.251154       20.809222                -0.1065631      -1.1934152
## IL6                   1.732941       19.975903                -0.2037340      -3.5671775
## CSF2                  2.146082        4.987054                 0.1212042       0.4370389
##      shEWS.Hrs6_log2FoldChange shEWS.Hrs6_stat
## CCL2                0.04893024       0.8812305
## IL6                -0.09350500      -2.4798332
## CSF2                0.12280155       0.6723109
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
