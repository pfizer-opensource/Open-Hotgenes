if(interactive()){
  # load packages
  library(Hotgenes)
  
  dds_Hotgenes_dir <- system.file("extdata",
                                  paste0("dds_Hotgenes", ".RDS"),
                                  package = "Hotgenes",
                                  mustWork = TRUE
  )
  # from DESeq2
  htgs <- readRDS(dds_Hotgenes_dir)
  
  # VPlot -------------------------------------------------------------------
  
  # see all contrasts in object summary
  htgs %>% Mapper_()
  
  # select by name
 htgs %>%
    VPlot(
      contrasts = "sh_EWS_vs_Ctrl",
      .log2FoldChange = 0.5,
      hotList = "C1S",
      Hide_labels = FALSE,
      point_label_size = 5,
      base_size = 20,
      padj_cut = 0.1
    )

  # select by name
 
  htgs %>%
    VPlot(
      contrasts = "sh_EWS_vs_Ctrl",
      .log2FoldChange = 0.5,
      interactive = TRUE,
      hotList = "CCL2",
      Hide_labels = FALSE,
      pointSize = 2,
      max.overlaps = 50,
      legend_pointSize = 50,
      point_label_size = 5,
      padj_cut = 0.1
    ) 
  
}