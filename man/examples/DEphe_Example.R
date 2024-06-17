if(interactive()) {
  # load packages
  library(Hotgenes)
  
  dds_Hotgenes_dir <- system.file(
    "extdata",
    paste0("dds_Hotgenes", ".RDS"),
    package = "Hotgenes",
    mustWork = TRUE
  )
  # from DESeq2
  htgs <- readRDS(dds_Hotgenes_dir)
  
  # DEphe -------------------------------------------------------------------
  ss <- SampleIDs_(htgs)[1:8]
  ss
  
  feat <- Features_(htgs)
   
  
  htgs %>%
    DEphe(
      Topn = 2,
      SampleIDs = ss,
      annotation_colors = coldata_palettes(htgs),
      arrangeby = c("Hrs", "sh"),
      annotations = c("Hrs", "sh"),
      cellheight = 10,
      cellwidth = 8
    )
  
  
  
  
  
  # label by alias in mapper
  
  htgs %>%
    DEphe(
      Topn = 3,
      SampleIDs = ss,
      annotation_colors = coldata_palettes(htgs),
      arrangeby = c("Hrs", "sh"),
      annotations = c("Hrs", "sh"),
      label_by = "ensembl_id",
      cellheight = 10,
      cellwidth = 8
    )
  
 
  
}