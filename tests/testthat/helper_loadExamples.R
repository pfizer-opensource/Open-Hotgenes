library(Hotgenes)

fit_Hotgenes_dir <- system.file("extdata",
                                paste0("fit_Hotgenes", ".RDS"),
                                package = "Hotgenes",
                                mustWork = TRUE
)

# from limma
fit_Hotgenes <- readRDS(fit_Hotgenes_dir) 

dds_Hotgenes_dir <- system.file("extdata",
                                paste0("dds_Hotgenes", ".RDS"),
                                package = "Hotgenes",
                                mustWork = TRUE)
# from DESeq2
htgs <- readRDS(dds_Hotgenes_dir) 

# load dds_con.Rdata
dds_con_dir <- system.file("extdata",
                           "dds_con.Rdata",
                           package = "Hotgenes",
                           mustWork = TRUE
)
load(dds_con_dir)
