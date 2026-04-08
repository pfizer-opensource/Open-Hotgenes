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
  
  # load signatures using msigdbr wrapper
  # This returns a named list of pathways
  
  H_paths <- msigdbr_wrapper(
    species = "human",
    set = c("CP:REACTOME", "CP:KEGG"),
    gene_col = "gene_symbol"
  )
  
  H_paths %>% names() %>% head()
  # custom geneset functions can be supplied like this
  # This can be used in the Shiny_Hotgenes() function
  Custom_geneset_function_list <- OntologyMethods(
    Ontology_Function = list("msigdbr" = msigdbr_wrapper),
    InputChoices = list("msigdbr" = c("CP:REACTOME", "CP:KEGG")),
    gene_col_choices = list("msigdbr" = c(
      "gene_symbol", "entrez_gene", "ensembl_gene"
    )),
    species_choices = list("msigdbr" = c("human", "mouse", "rat", "dog")),
    versions = list("msigdbr" = packageVersion("msigdbr"))
  )
  
  
  msigdbr_pthyways <- OntologyFunctions(
    Methods = Custom_geneset_function_list,
    db = "msigdbr",
    species = "human",
    set = c("CP:REACTOME", "CP:KEGG"),
    gene_col = "gene_symbol"
  )
  
  testthat::expect_equal(H_paths, msigdbr_pthyways)
  
  # These genes sets are mapped to gene symbols
  # Verify Feature col contains gene symbols, too
  # In this example the "Feature" column contains gene symbols
  
  htgs %>% Mapper_()
  
  # Get ranks
  InputRanks <- htgs %>%
    DE(
      Report = "Ranks",
      contrasts = "Hrs_2_vs_0",
      Rank_name = "Feature", # see above
      padj_cut = 1
    )
  
  
  # fgsea wrapper --------
  Out_GSEA <- fgsea_(
    Ranks = InputRanks,
    pathways = H_paths,
    nproc = 8,
    minSize = 5,
    maxSize = Inf
  )
  
  # Get details for all
  Out_GSEA %>%
    fgsea_Results(
      padj_cut = 0.2,
      mode = "D"
    )
  
  # Or for one
  Out_GSEA %>%
    fgsea_Results(
      contrasts = "Hrs_2_vs_0",
      padj_cut = 0.2,
      mode = "leadingEdge"
    )
  
  # Generate a summary plot
  Out_GSEA %>%
    GSEA_Plots(
      contrasts = "Hrs_2_vs_0",
      padj_cut = 0.2,
      width = 30,
      Topn = 2
    )
  
  # plotEnrichment_
  plotEnrichment_(
    Out_GSEA, "Hrs_2_vs_0",
    "reactome_interleukin_10_signaling"
  )
  
  
  # leadingGenes
  leadingGenes(
    Out_GSEA, "Hrs_2_vs_0",
    "reactome_interleukin_10_signaling"
  )
  
}