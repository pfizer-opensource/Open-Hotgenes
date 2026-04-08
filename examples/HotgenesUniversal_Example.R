if(interactive()){
  # load packages
  library(Hotgenes)
  
  
  # load example data -------------------------------------------------------
  
  dds_Hotgenes_dir <- system.file("extdata",
                                  paste0("dds_Hotgenes", ".RDS"),
                                  package = "Hotgenes",
                                  mustWork = TRUE
  )
  
  htgs <- readRDS(dds_Hotgenes_dir)
  
  # creating new Hotgenes object --------------------------------------------
  
  # Getting example named list of DE statistics
  NewDE <- Output_DE_(htgs, as_list = TRUE, padj_cut = 1)
  NewDE %>% purrr::map(~ head(.))
  
  # Getting example named list of normalized data
  NormlData <- Normalized_Data_(htgs)
  NormlData %>% purrr::map(~ head(.))
  
  # Getting example coldata
  ExpColdata <- coldata_(htgs)
  ExpColdata
  
  # Getting example original data object used for DE analysis
  # This example was generated from DESeq2
  OrigDEObj <- O_(htgs)
  OrigDEObj %>% class()
  
  # Getting example design matrix
  DE_design <- designMatrix_(htgs)
  DE_design
  
  # Getting example mapper
  MapperDF <- Mapper_(htgs)
  MapperDF %>% head()
  
  # Converting example objects to hotgenes
  Hotgenes_Object <- HotgenesUniversal(
    Output_DE = NewDE,
    Normalized_Expression = NormlData,
    coldata = ExpColdata,
    Original_Object = OrigDEObj,
    designMatrix = DE_design,
    Mapper = MapperDF
  )
  
  Hotgenes_Object
  
  
}