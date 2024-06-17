
# VPlot -------------------------------------------------------------------

# No error here
htgs %>%
  VPlot(
    contrasts = "sh_EWS_vs_Ctrl",
    .log2FoldChange = 0.5,
    hotList = "C1R",
    Hide_labels = FALSE,
    point_label_size = 5,
    base_size = 20,
    padj_cut = 0.1
  )


# No error here
htgs %>%
  VPlot(
    contrasts = "sh_EWS_vs_Ctrl",
    .log2FoldChange = 0.5,
    Hide_labels = FALSE,
    point_label_size = 5,
    base_size = 20,
    padj_cut = 0.1,
    repel_labels = FALSE
  )

# no error here
  VPlot(htgs)

# contrast doesn't exit
   testthat::expect_error( VPlot(htgs,
  contrasts = "DoesNotExist"))
   
   
   
 
  
  