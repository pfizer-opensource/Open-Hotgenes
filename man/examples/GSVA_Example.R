if(interactive()){
  # load package
  library(Hotgenes)
  
  dds_Hotgenes_dir <- system.file("extdata",
                                  paste0("dds_Hotgenes", ".RDS"),
                                  package = "Hotgenes",
                                  mustWork = TRUE
  )
  
  # from DESeq2
  htgs <- readRDS(dds_Hotgenes_dir)
  
  # msigdbr_wrapper ---------------------------------------------------------
  choice_set <- "CP:KEGG"
  choice_id <- "gene_symbol"
  
  gsList <- msigdbr_wrapper(
    species = "human",
    set = choice_set,
    gene_col = choice_id
  )
  
  
  names(gsList)
  
  # HotgeneSets -------------------------------------------------------------
  htgs %>% Mapper_()
  
  HotgeneSets_out <- HotgeneSets(
    Hotgenes = htgs,
    geneSets = gsList,
    kcdf = "Gaussian",
    method = "ssgsea",
    minSize = 2,
    maxSize = Inf
  )
  
  
  HotgeneSets_out 
  
}