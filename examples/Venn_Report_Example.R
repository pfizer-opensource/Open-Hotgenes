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
  
  # Venn_Report -----------------------------------------
  # only four lists at a time
  idc <- htgs %>% contrasts_()
  idc[1:4]
  
 
  # set Report = "Features" for overlap (up to 4 contrasts)
  # set Report = "contrast_dir" to include direction
  ## limit 2 contrast for "contrast_dir"
  V <- htgs %>%
    DE(
     # Report = "Features",
      Report = "contrast_dir",
      contrasts = idc[1:2]
    ) %>%
    Venn_Report(set_name_size = 4,
                stroke_size = 0.5,
                text_size = 4)
  
  # view
  V$vennD
  
  # names of lists that intersect
  V$Names
  
  # all lists
  V$Intsect
  
  # Label venn diagram groups by query gene ---------------------------------
  
  Venn_hotList <- c("CCL2", "CXCL1")
  
  htgs %>%
    DE(contrasts = idc[1:4]) %>%
    purrr::map(function(x) {
      x %>%
        dplyr::filter(.data$Feature %in% Venn_hotList) %>%
        dplyr::pull("Feature")
    }) %>%
    Venn_Report(
      label_sep = "\n",
      show_elements = TRUE
    )
  
}