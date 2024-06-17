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

# call the object to get summary
fit_Hotgenes


# auxiliary_assays_ -------------------------------------------------------
fit_Hotgenes %>% 
auxiliary_assays_()

# ExpressionSlots_ ------------------------------------
fit_Hotgenes %>% ExpressionSlots_()

# Normalized_Data_ --------------------------------------------------------
# should return object class data.frame when slot specified
fit_Hotgenes %>%
  Normalized_Data_(slot = "logCPM") %>% head()

# O_ --------------------------------------------------
# get original object
fit_Hotgenes %>%
  O_() %>%
  class()

# designMatrix_ ---------------------------------------
fit_Hotgenes %>%
  designMatrix_() %>%
  head()

# coldata_ --------------------------------------------
fit_Hotgenes %>% 
  coldata_()

# SampleIDs_ ------------------------------------------
fit_Hotgenes %>% SampleIDs_()

# contrasts_ ------------------------------------------
fit_Hotgenes %>% contrasts_()

# Mapper_ ---------------------------------------------
fit_Hotgenes %>%
  Mapper_() %>%
  head()

# Features_ -------------------------------------------
fit_Hotgenes %>% Features_()

# Background_ -----------------------------------------
# this returns Features mapped to Mapper aliases
fit_Hotgenes %>% 
  Background_(Col = "ensembl_id")


# updating object ---------------------------------------------------------

# Adding column of strings
coldata_(fit_Hotgenes)$New <- "value"
coldata_(fit_Hotgenes)


# update mapper slot

# in case you wanted to include aliases for your genes
# requires a "Feature" column that contains gene names in expression matrix
dbCon <- org.Hs.eg.db::org.Hs.eg_dbconn()

sqlQuery <- "SELECT * FROM ENSEMBL, gene_info WHERE ENSEMBL._id == gene_info._id;"

ensembl_Symbol <- DBI::dbGetQuery(dbCon, sqlQuery) %>%
  dplyr::select(dplyr::any_of(c("Feature" = "symbol", 
                                "gene_name", 
                                "ensembl_id"))) %>%
  tibble::as_tibble()

Test_Hotgenes <- fit_Hotgenes

Mapper_(Test_Hotgenes) <- ensembl_Symbol
Mapper_(Test_Hotgenes)


# auxiliary_assays_ --------------------------------------------------------
fit_Hotgenes_2 <- fit_Hotgenes

coldata_(fit_Hotgenes_2)$Samples <- as.character(SampleIDs_(fit_Hotgenes_2))

df_assays_slot <- auxiliary_assays_default(fit_Hotgenes_2) 
df_assays_slot

# making example data
set.seed(12)

max_len <- length(SampleIDs_(fit_Hotgenes_2))

AssayData  <- df_assays_slot %>% 
  dplyr::mutate(assay1 = rnorm(max_len),
                assay2 = rnorm(max_len))

AssayData

auxiliary_assays_(fit_Hotgenes_2) <- AssayData
auxiliary_assays_(fit_Hotgenes_2) 

