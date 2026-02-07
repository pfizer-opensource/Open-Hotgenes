

# Ranks are not nameless --------------------------------------------------

df_missing_Prot <- data.frame(
  ProteinID = c("IL-8", "GM-CSF", ""),
  Feature = c("CXCL8", "CSF2", "IL6")
)

New_dds_Hotgenes <- htgs
Mapper_(New_dds_Hotgenes) <- df_missing_Prot


DE_1 <- New_dds_Hotgenes %>%
  DE(
    Report = "Features",
    contrasts = "shEWS.Hrs2",
    padj_cut = 1
  ) %>%
  purrr::chuck(1)


DE_2 <- New_dds_Hotgenes %>%
  DE(
    Report = "Ranks",
    contrasts = "shEWS.Hrs2",
    padj_cut = 1
  ) %>%
  purrr::chuck(1)

testthat::expect_equal(length(DE_1), length(names(DE_2)))

DE_2_bad_mapping <- New_dds_Hotgenes %>%
  DE(
    Report = "Ranks",
    contrasts = "shEWS.Hrs2",
    padj_cut = 1,
    Rank_name = "ProteinID",
    mapFeatures = TRUE
  ) %>%
  purrr::chuck(1)


testthat::expect_equal(length(DE_2_bad_mapping), length(DE_2)) 


# Testing mapFeatures -----------------------------------------------------
# mapFeatures should not impact "Length" | "Features" Report

# how many DEGs
Lengths_mapped <- htgs %>%
  DE(Report = "Length", mapFeatures = TRUE)

# Lengths_mapped

Lengths_Not_mapped <- htgs %>%
  DE(Report = "Length", mapFeatures = FALSE)

Lengths_Not_mapped

# verify that mapping does not change values
testthat::expect_identical(
  Lengths_mapped,
  Lengths_Not_mapped
)

# getting DEGs
DEGs_mapped <- htgs %>%
  DE(Report = "Features", mapFeatures = TRUE)

DEGs_mapped

DEGs_Not_mapped <- htgs %>%
  DE(Report = "Features", mapFeatures = FALSE)

DEGs_Not_mapped

# verify that mapping does not change values
testthat::expect_identical(
  DEGs_mapped,
  DEGs_Not_mapped
)


# check minimum -----------------------------------------------------------
min_map <- Mapper_(htgs) %>% 
  dplyr::mutate(species = "hooman")

htgs2 <- htgs
Mapper_(htgs2) <- min_map

# identical when min = 0
testthat::expect_true(
  identical(names(Mapper_(htgs2, min_n = 0)), names(min_map)))

# not identical with defaults
testthat::expect_false(
  identical(names(Mapper_(htgs2)), names(min_map)))



