# load package
library(Hotgenes)

fit_Hotgenes_dir <- system.file("extdata",
                                paste0("fit_Hotgenes", ".RDS"),
                                package = "Hotgenes",
                                mustWork = TRUE
)


# from limma
fit_Hotgenes <- readRDS(fit_Hotgenes_dir) %>% 
  update_object()


enlist_ <-c("ENSG00000170345", 
            "ENSG00000159403",
            'FOS',
            "ENSG00000011600")


hotList_mapper(fit_Hotgenes,
             hotList = NULL)

mapped_f<-hotList_Feature(fit_Hotgenes,
                hotList = enlist_)
mapped_f

# with utilities
DE(fit_Hotgenes, hotList = mapped_f)

fit_Hotgenes %>%
  VPlot(
    contrasts = "shEWS.Hrs6",
    .log2FoldChange = 0.5,
    hotList = mapped_f,
    Hide_labels = FALSE,
    point_label_size = 5,
    base_size = 10,
    padj_cut = 0.1
  )


fit_Hotgenes %>%
  DEPlot()

fit_Hotgenes %>%
  DEPlot(
    hotList = mapped_f)


fit_Hotgenes %>%
  DEphe(
    hotList = mapped_f)


