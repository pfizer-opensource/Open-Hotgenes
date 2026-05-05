Interactive Exploration with Shiny
================

## Overview

`Shiny_Hotgenes()` launches an interactive web application that lets you
explore any Hotgenes object—or a named list of Hotgenes objects—without
writing additional code. The app bundles the most common analysis tasks
into point-and-click panels:

- **DE Summary** — filterable table of differential expression results
- **Volcano plots** — interactive volcano plots per contrast
- **Heatmaps** — top-hit heatmaps with flexible sample/gene selection
- **Expression plots** — per-gene trajectory plots
- **Venn diagrams** — overlap across contrasts
- **GSEA** — gene-set enrichment via `fgsea` with `msigdbr` or custom
  gene sets
- **GSVA** — sample-wise pathway activity scores via GSVA/ssGSEA

``` r
library(Hotgenes)
```

------------------------------------------------------------------------

## 1. Launching the App with a Single Object

The simplest call passes one Hotgenes object:

``` r
# Load a pre-built example
fit_Hotgenes <- readRDS(
  system.file("extdata", "fit_Hotgenes.RDS",
              package = "Hotgenes",
              mustWork = TRUE)
)

Shiny_Hotgenes(fit_Hotgenes)
```

The app opens in your default browser (or the RStudio viewer pane).
Stop it by pressing **Escape** or closing the browser tab.

------------------------------------------------------------------------

## 2. Comparing Multiple Objects Side-by-Side

Pass a **named list** of Hotgenes objects to switch between experiments
in the app. The names become the labels in the dataset selector:

``` r
# Load a second object (DESeq2-based example)
dds_Hotgenes <- readRDS(
  system.file("extdata", "dds_Hotgenes.RDS",
              package = "Hotgenes",
              mustWork = TRUE)
) |> update_object()

# Combine into a named list
Hotgenes_list <- list(
  limma_Ewing  = fit_Hotgenes,
  DESeq2_Ewing = dds_Hotgenes
)

Shiny_Hotgenes(Hotgenes_list)
```

------------------------------------------------------------------------

## 3. Gene-Set Enrichment Analysis (GSEA)

### Built-in msigdbr gene sets

`msigdbr_wrapper()` returns a named list of gene sets sourced from the
MSigDB via the `msigdbr` package. Any collection supported by `msigdbr`
can be used.

``` r
# Retrieve KEGG and Reactome pathways for human
H_paths <- msigdbr_wrapper(
  species  = "human",
  set      = c("CP:KEGG", "CP:REACTOME"),
  gene_col = "gene_symbol"
)

length(H_paths)
H_paths |> names() |> head(5)
```

### Running GSEA with `fgsea_()`

`fgsea_()` is a wrapper around `fgsea::fgsea()` that accepts the ranked
vector returned by `DE(..., Report = "Ranks")`:

``` r
fit_Hotgenes <- readRDS(
  system.file("extdata", "fit_Hotgenes.RDS",
              package = "Hotgenes",
              mustWork = TRUE)
)

# Get ranked statistics for one contrast
InputRanks <- fit_Hotgenes |>
  DE(
    Report    = "Ranks",
    contrasts  = "Hrs_2_vs_0",
    Rank_name = "Feature",
    padj_cut  = 1
  )

head(InputRanks)
```

``` r
# Run GSEA
Out_GSEA <- fgsea_(
  Ranks    = InputRanks,
  pathways = H_paths,
  nproc    = 1,
  minSize  = 5,
  maxSize  = Inf
)
```

### Inspecting GSEA results

``` r
# Tabular summary of significant pathways
Out_GSEA |>
  fgsea_Results(
    contrasts = "Hrs_2_vs_0",
    padj_cut  = 0.2,
    mode      = "D"
  ) |> head()
```

``` r
# Leading-edge genes for one pathway
Out_GSEA |>
  fgsea_Results(
    contrasts = "Hrs_2_vs_0",
    padj_cut  = 0.2,
    mode      = "leadingEdge"
  ) |> head()
```

### Visualizing GSEA results

``` r
Out_GSEA |>
  GSEA_Plots(
    contrasts = "Hrs_2_vs_0",
    padj_cut  = 0.2,
    Topn      = 3,
    width     = 20
  )
```

``` r
# Enrichment plot for a single pathway
sig_paths <- Out_GSEA |>
  fgsea_Results(contrasts = "Hrs_2_vs_0",
                padj_cut  = 0.2,
                mode      = "D")

if (nrow(sig_paths) > 0) {
  plotEnrichment_(Out_GSEA, "Hrs_2_vs_0", sig_paths$pathway[1])
}
```

``` r
# Retrieve the leading-edge gene list for one pathway
if (nrow(sig_paths) > 0) {
  leadingGenes(Out_GSEA, "Hrs_2_vs_0", sig_paths$pathway[1])
}
```

------------------------------------------------------------------------

## 4. Sample-wise Pathway Activity with `HotgeneSets()`

`HotgeneSets()` runs GSVA (or ssGSEA, PLAGE, etc.) on the expression
data to produce per-sample pathway activity scores, then returns a new
Hotgenes object with the pathway scores as the expression matrix.

``` r
# Use the same gene sets as above
HotgeneSets_out <- HotgeneSets(
  Hotgenes = fit_Hotgenes,
  geneSets = H_paths,
  kcdf     = "Gaussian",
  method   = "ssgsea",
  minSize  = 2,
  maxSize  = Inf
)

HotgeneSets_out
```

The result is a Hotgenes object whose expression matrix rows are pathway
names and whose DE slot contains pathway-level differential activity
results. You can pass it directly to `Shiny_Hotgenes()`.

------------------------------------------------------------------------

## 5. Custom Gene Sets in the Shiny App

You can configure the Shiny app to use your own gene-set retrieval
function (or expose multiple databases) via `OntologyMethods()` and
`OntologyFunctions()`:

``` r
# Define custom database functions
Custom_db <- OntologyMethods(
  Ontology_Function  = list("msigdbr" = msigdbr_wrapper),
  InputChoices       = list("msigdbr" = c("CP:REACTOME", "CP:KEGG", "H")),
  gene_col_choices   = list("msigdbr" = c("gene_symbol",
                                           "entrez_gene",
                                           "ensembl_gene")),
  species_choices    = list("msigdbr" = c("human", "mouse", "rat")),
  versions           = list("msigdbr" = packageVersion("msigdbr"))
)

# Pass to Shiny app
Shiny_Hotgenes(fit_Hotgenes, OntologyDB = Custom_db)
```

To supply pre-computed gene sets (bypassing in-app database queries):

``` r
# Use the pathways list returned by msigdbr_wrapper()
Shiny_Hotgenes(fit_Hotgenes, gene_sets = H_paths)
```

------------------------------------------------------------------------

## 6. Tips for Large Datasets

| Tip | Rationale |
|---|---|
| Store objects as `.RDS` files and load them with `readRDS()` | Avoids re-running expensive analyses every session |
| Use `update_object()` when loading saved objects | Ensures slot structure matches the current package version |
| Subset contrasts with `contrasts = c("A_vs_B")` | Reduces memory usage in `DE()` and plot functions |
| Set `padj_cut` and `.log2FoldChange` tightly | Speeds up Venn diagrams and heatmaps for large gene lists |
| Prefer `method = "ssgsea"` over `"gsva"` for sparse data | ssGSEA is more robust when many zeros are present |
| Pre-compute GSEA results and pass `Out_GSEA` to the app | Avoids waiting for `fgsea_` inside the app for large pathway databases |

------------------------------------------------------------------------

## Summary

| Function | Purpose |
|---|---|
| `Shiny_Hotgenes()` | Launch the interactive Shiny app |
| `msigdbr_wrapper()` | Retrieve MSigDB gene sets via msigdbr |
| `fgsea_()` | Run fgsea gene-set enrichment |
| `fgsea_Results()` | Extract GSEA results table |
| `GSEA_Plots()` | Bar plot of top enriched pathways |
| `plotEnrichment_()` | Enrichment curve for one pathway |
| `leadingGenes()` | Leading-edge genes for one pathway |
| `HotgeneSets()` | Compute GSVA/ssGSEA pathway activity |
| `OntologyMethods()` | Define custom gene-set database methods for Shiny |
