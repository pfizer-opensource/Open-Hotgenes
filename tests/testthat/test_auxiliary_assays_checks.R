

# auxiliary_assays_ --------------------------------------------------------
dds_Hotgenes_2 <- htgs

df_assays_slot <- auxiliary_assays_default(dds_Hotgenes_2, coldata_ids = NULL) 

# making example data
set.seed(12)

max_len <- length(SampleIDs_(dds_Hotgenes_2))

AssayData  <- df_assays_slot %>% 
  dplyr::mutate(assay1 = rnorm(max_len),
                assay2 = rnorm(max_len))

AssayData

# testing missing SampleIDs column
MissingSampleCol <- AssayData %>% 
  dplyr::select(-"SampleIDs")

testthat::expect_error(auxiliary_assays_(dds_Hotgenes_2) <- MissingSampleCol,
                       regexp = "SampleIDs")

# test typo
Typo_SampleCol <- AssayData %>% 
  dplyr::rename("sampleIDs" = "SampleIDs")

testthat::expect_error(auxiliary_assays_(dds_Hotgenes_2) <- Typo_SampleCol,
                       regexp = "SampleIDs")

# test no matching with metadata
NoMatch_SampleCol <- AssayData %>% 
  dplyr::mutate(SampleIDs = stringr::str_remove(.data$SampleIDs,
                                                      pattern = "sh"))

testthat::expect_error(auxiliary_assays_(dds_Hotgenes_2) <- NoMatch_SampleCol,
                       regexp = "Sample")

# partial overlap ok
PartialMatch_SampleCol <- AssayData %>% 
  dplyr::slice_head(n = 5)

testthat::expect_no_error(auxiliary_assays_(dds_Hotgenes_2) <- PartialMatch_SampleCol)


# Checking for overlap across slots ---------------------------------------

# making example data match with Expression

set.seed(12)

max_len <- length(SampleIDs_(dds_Hotgenes_2))

AssayDataExp  <- df_assays_slot %>% 
  dplyr::mutate(CXCL8 = rnorm(max_len),
                IL6 = rnorm(max_len))

testthat::expect_message(auxiliary_assays_(dds_Hotgenes_2) <- AssayDataExp,
                       regexp = "already in")

# making example data match with coldata

set.seed(12)

max_len <- length(SampleIDs_(dds_Hotgenes_2))

AssayColdata  <- df_assays_slot %>% 
  dplyr::mutate(Hrs = rnorm(max_len),
                assay2 = rnorm(max_len))


testthat::expect_message(auxiliary_assays_(dds_Hotgenes_2) <- AssayColdata,
                       regexp = "already in")
