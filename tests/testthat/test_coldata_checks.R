# Getting example data ----------------------------------------------------

# load package
library(Hotgenes)

# min level set
htgs_mod <- htgs
coldata_(htgs_mod)$New_character <- "value"
coldata_(htgs_mod)$New_factor <- factor("value")


df_coldata_min_all <- coldata_(htgs_mod, mode = "all")
df_coldata_min_factor <- coldata_(htgs_mod, mode = "quali")

testthat::expect_true(all(c("New_character", "New_factor") %in% 
                            names(df_coldata_min_all)) )

testthat::expect_true(
  all(c("New_factor") %in% names(df_coldata_min_factor),
                          
!"New_character" %in% names(df_coldata_min_factor)) )

# set min = 2 should only be "Hrs"
df_coldata_min_2 <- coldata_(htgs, mode = "quali", min_col_levels = 2)
testthat::expect_true(names(df_coldata_min_2) == "Hrs" )

# valid_factors -----------------------------------------------------------
# removes invalid factors
df_coldata <- coldata_(htgs, mode = "quali")
df_coldata

# single level or all NAs removed
df_coldata_w_invalid <- df_coldata
df_coldata_w_invalid$inv1 <- NA
df_coldata_w_invalid$inv2 <- c("Human")
df_coldata_w_invalid$inv3 <- 4
df_coldata_w_invalid

testthat::expect_equal(
  df_coldata,
  valid_factors(df_coldata_w_invalid)
)

# SampleIDs can't be added 
testthat::expect_error(
  coldata_(htgs)$SampleIDs <- as.character(SampleIDs_(htgs)),
  regexp = "SampleIDs can't"
)

# UpdateLevelsbyList updates levels ---------------------------------------

df_coldata$Bio_Rep %>% levels()
df_coldata$Hrs %>% levels()

# change these levels
named_refs <- list(
  Bio_Rep = c("DoesntExist", "2"),
  Hrs = c("6")
)

New_df_coldata <- df_coldata %>%
  UpdateLevelsbyList(named_refs)

testthat::expect_equal(levels(New_df_coldata$Bio_Rep)[1], "2")
testthat::expect_equal(levels(New_df_coldata$Hrs)[1], "6")

# coldata_palettes works with Hotgenes or dfs -----------------------------

NewCOlFromHotgenes <- coldata_palettes(htgs, coldata_ids = coldata_names(htgs))
# NewCOlFromHotgenes

# or from data.frame
DF_coldata <- coldata_(htgs)
# DF_coldata

NewCOlFromData.frame <- coldata_palettes(DF_coldata)
# NewCOlFromData.frame

# verify equal results
testthat::expect_equal(NewCOlFromHotgenes, NewCOlFromData.frame)
