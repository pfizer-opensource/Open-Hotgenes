if(interactive()){
  # Getting example data ----------------------------------------------------
  
  # load package
  library(Hotgenes)
  
  dds_Hotgenes_dir <- system.file("extdata",
                                  paste0("dds_Hotgenes", ".RDS"),
                                  package = "Hotgenes",
                                  mustWork = TRUE
  )
  # from DESeq2
  htgs <- readRDS(dds_Hotgenes_dir)


# clustering gene ranks ---------------------------------------------------
  InputFeatures <- htgs %>%
    DE(
      Report = "Features",
      padj_cut = 0.1
    )  |> unlist(use.names = FALSE) |> 
    unique()
  
  # Get ranks
  InputRanks <- htgs %>%
    Output_DE_df(
     hotList = InputFeatures
      #padj_cut = 0.1
    ) |> 
    dplyr::ungroup() |> 
    dplyr::select(c("Feature", "group"= "contrast", "stat")) |> 
    dplyr::distinct() 
  

  Input_Frame <-  make_stat_frame(input_stats = InputRanks) 
  

  out_1 <- FactoWrapper_DFs(
    df = Input_Frame)
    
  out_1$TopTibble
  
# with more plotting options
  factoExtra_DFs(PCA_obj = out_1, 
  select.ind = list( contrib = 10),
  plot_type = "biplot",
  repel = TRUE)
  
  
  # clustering genesets -----------------------------------------------------

  # load signatures using msigdbr wrapper
  # This returns a named list of pathways
  
  H_paths <- msigdbr_wrapper(
    species = "human",
    set = c("GO:BP"),
    gene_col = "gene_symbol"
  )
  
  # These genes sets are mapped to gene symbols
  # Verify Feature col contains gene symbols, too
  # In this example the "Feature" column contains gene symbols
  
  # Get ranks
  InputRanks <- htgs %>%
    DE(
      Report = "Ranks",
      Rank_name = "Feature", # see above
      padj_cut = 1
    )
  
  
  
  # fgsea wrapper --------
  Out_GSEA <- fgsea_(
    Ranks = InputRanks,
    pathways = H_paths,
    nproc = 2,
    minSize = 5,
    maxSize = Inf
  )
  
  # Get details for all
  list_out <- Out_GSEA %>%
    fgsea_Results(
      padj_cut = 1,
      mode = "D"
    ) |> 
    purrr::keep(~nrow(.x) > 0) |> 
    purrr::list_rbind(names_to = "group")
  
  Input_Frame_2 <-  make_stat_frame(input_stats = list_out,
                                  feature_var = "pathway",
                                  value_var = "NES") 
  
  head(Input_Frame_2)
  
  # clusters
  out_2 <- FactoWrapper_DFs(
    df = Input_Frame_2
    
  )
  
  out_2$TopTibble
  
  # with more plotting options
  pca_genesets_clusters <- factoExtra_DFs(PCA_obj = out_2, 
                # select.ind = list( contrib = 100),
                label_size = 0,
                point_size = 2,
                col.var = "black",
                alpha.ind  = 0.1,
                ellipse.alpha = 1,
                 #plot_type = "biplot",
                 repel = TRUE)
  
  pca_genesets_clusters
  
}