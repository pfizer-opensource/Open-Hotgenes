if(interactive()) {
  # load package
  library(Hotgenes)
  
  fit_Hotgenes_dir <- system.file(
    "extdata",
    paste0("fit_Hotgenes", ".RDS"),
    package = "Hotgenes",
    mustWork = TRUE
  )
  
  # from limma
  fit_Hotgenes <- readRDS(fit_Hotgenes_dir)
  
  # object summary on call
  fit_Hotgenes
  
  
  # list available data
  ExpressionSlots_(fit_Hotgenes)
  
  # View plot
  
  BoxPlot(fit_Hotgenes, SampleIDs = SampleIDs_(fit_Hotgenes)[1:6])
  
  # DE related utilities ----------------------------------------------------
  
  # summarize across contrasts
  fit_Hotgenes %>% DEPlot()
  
  # check genes of interest
  genesOfinterest <- c("IL6", "CSF2")
  
  fit_Hotgenes %>%
    DEPlot(hotList = genesOfinterest)
  
  # genes of interest
  fit_Hotgenes %>%
    DE(Report = "Details",
       hotList = genesOfinterest)
  
  # expression data merged with coldata
  DExps(
    fit_Hotgenes,
    hotList = genesOfinterest,
    coldata_ids = c("Hrs", "sh"),
    Query_set = TRUE
  ) %>% head()
  
  # Explore all DEGs --------------------------------------------------------
  
  fit_Hotgenes %>%
    DE(Report = "Details")
  
  # Export Ranks for GSEA
  fit_Hotgenes %>%
    DE(Report = "Ranks")
  
  # or FC
  fit_Hotgenes %>%
    DE(Report = "FC")
  
  # DECoefs ---------------------------------------------
  DECoefs(fit_Hotgenes, hotList = genesOfinterest)
  
}