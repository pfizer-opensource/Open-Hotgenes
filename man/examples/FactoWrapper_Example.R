if(interactive()){
  # load packages
  library(Hotgenes)
  library(ggplot2)
  
  
  dds_Hotgenes_dir <- system.file("extdata",
                                  paste0("dds_Hotgenes", ".RDS"),
                                  package = "Hotgenes",
                                  mustWork = TRUE
  )
  # from DESeq2
  htgs <- readRDS(dds_Hotgenes_dir)
  
  
  # Adding auxiliary assay data
  set.seed(12)
  
  max_len <- length(SampleIDs_(htgs))
  
  AssayData  <- auxiliary_assays_default(htgs)  %>% 
    dplyr::mutate(assay1 = rnorm(max_len),
                  assay2 = rnorm(max_len))
  
  auxiliary_assays_(htgs) <- AssayData
  auxiliary_assays_(htgs)
  
  # FactoWrapper --------------------------------------
  
  # run PCA
  FactoOutput <- FactoWrapper(htgs,
                              contrasts = "Hrs_6_vs_0",
                              coldata_ids = c("text_Hrs", "Hrs", "Time"),
                              aux_features = "assay1",
                              padj_cut = 0.1,
                              .log2FoldChange = 0
  )
  
  # getting HCPC details
  FactoOutput$TopTibble # Feature
  FactoOutput$TopTibble_sup # Coldata Quanti.variable
  FactoOutput$TopGroups # Coldata Factors
  
  # Feature Ranks
  FactoOutput$Ranks 
  
  # Getting FactoMiner PCA object
  res.pca <- FactoOutput$res
  
  # loadings
  sweep(res.pca$var$coord, 2, sqrt(res.pca$eig[1:ncol(res.pca$var$coord), 1]), FUN = "/")
  
  # heatmap
  FactoOutput$TopTibble %>% 
    plyr::dlply("Cluster", identity) %>% 
    purrr::imap(function(x, y){
      x %>% 
        # tibble::enframe() %>% 
        dplyr::select(dplyr::any_of(c("Cluster", #"Interpretation", 
                                      "Feature"))) %>% 
        dplyr::slice(1:5) %>% 
        tibble::column_to_rownames("Feature")
      # dplyr::mutate()
    }) 
  

  annotation_row = FactoOutput$TopTibble %>% 
    plyr::dlply("Cluster", identity) %>% 
    purrr::imap(function(x, y){
      x %>% 
        # tibble::enframe() %>% 
        dplyr::select(dplyr::any_of(c("clust"="Cluster", #"Interpretation", 
                                      "Feature"))) %>% 
        dplyr::slice(1:5) %>% 
        tibble::column_to_rownames("Feature")
      # dplyr::mutate()
    })  %>% 
    purrr::list_rbind()
    
 df<- coldata_palettes(annotation_row)
 
 
  
  htgs %>%
    DEphe(
     # Topn = 2,
      show_colnames = FALSE,
      cutree_rows = 3,
      cluster_rows = TRUE,
      hotList = rownames(annotation_row),
      annotation_row = NULL ,
      annotation_colors = coldata_palettes(htgs) %>% append(df),
      arrangeby = c("Hrs", "sh"),
      annotations = c("Hrs", "sh"),
      cellheight = 10,
      cellwidth = 8
    )
  
  
  # update plot with fviz_pca_ind -------------------------------------------
  
  factoextra::fviz_pca_ind(FactoOutput$res,
                           axes = c(1, 2),
                           repel = FALSE, label = "none",
                           habillage = FactoOutput$res.hcpc$data.clust[, c("clust")],
                           col.quanti.sup = "black",
                           col.var = c("black"),
                           pointsize = 3,
                           labelsize = 4,
                           select.ind = list(contrib = 50),
                           col.ind.sup = "black",
                           ellipse.alpha = 0,
                           # title = title,
                           legend.title = "clust",
                           addEllipses = TRUE, ellipse.level = 0.5
  ) +
    theme_classic() +
    scale_shape_manual(values = rep(20, 100))
  
  
  
  # FactoWrapper_DFs --------------------------------------------------------
  
  
  ExpressionDat <- Normalized_Data_(htgs) %>%
    purrr::chuck(1) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Features") %>%
    tidyr::pivot_longer(cols = -"Features", names_to = "sample_id") %>%
    tidyr::pivot_wider(
      id_cols = "sample_id",
      names_from = "Features"
    )
  
  supp_data <- coldata_(htgs
  ) %>%
    tibble::rownames_to_column("sample_id") %>% 
    dplyr::left_join(auxiliary_assays_(htgs)%>%
                       tibble::rownames_to_column("sample_id"))
  
  
  PC <- FactoWrapper_DFs(
    ExpressionDat = ExpressionDat,
    supp_data = supp_data,
    sampleID_col = "sample_id"
  )
  
  PC$res_PPI_pa_1
  
  # with more plotting options
  factoExtra_DFs(PCA_obj = PC, repel = TRUE)
  
}