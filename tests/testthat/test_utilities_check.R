

# no error here
BoxPlot(fit_Hotgenes,
  SampleIDs = SampleIDs_(fit_Hotgenes)[0]
)


# summary plot
fit_Hotgenes %>% DEPlot()

# check genes of interest
fit_Hotgenes %>%
  DEPlot(hotList = c("IL6", "CSF2"))

# top 5 with details
fit_Hotgenes %>%
  DE(
    Report = "Details",
    Topn = 5
  )

# genes of interest
fit_Hotgenes %>%
  DE(
    Report = "Details",
    hotList = c("CSF2", "CXCL8", "TRAF2")
  )

# Export Ranks for GSEA
fit_Hotgenes %>%
  DE(Report = "Ranks")

# or FC
fit_Hotgenes %>%
  DE(Report = "FC")

# Expression Plots
fit_Hotgenes %>%
  ExpsPlot(
    yVar = "IL6",
    xVar = "sh", 
    group = "Hrs"
  ) +

  ggplot2::facet_wrap("Hrs",
    labeller = ggplot2::label_both
  )



# all hotList items must be truthy ----------------------------------------

example_1 <-fit_Hotgenes %>%
  Output_DE_(
    
    hotList = character(0)
    
  )


example_2 <-fit_Hotgenes %>%
  Output_DE_(
    
    hotList = c("IL6", "CSF2", "")
    
  )


example_3 <-fit_Hotgenes %>%
  Output_DE_(
    
    hotList = c( "")
    
  )

testthat::expect_identical(example_1, example_2)

testthat::expect_identical(example_2, example_3)


exp_4_q <- c("IL6", "CSF2")

example_4 <-fit_Hotgenes %>%
  Output_DE_(
    
    hotList = exp_4_q
    
  )

testthat::expect_true(all(example_4$Feature %in% exp_4_q))
# DECoefs ---------------------------------------------
DECoefs(fit_Hotgenes) %>% head()

# DExps -----------------------------------------------
# list available data
ExpressionSlots_(fit_Hotgenes)

# only matched and unique values returned
DExps(fit_Hotgenes,
  hotList = c("Doesn't Exist", "TUBB", "TUBB"),
  coldata_ids = c("Hrs", "sh"),
  Query_set = TRUE
) %>% head()

# all data merged as data.frame
DExps(fit_Hotgenes,
  Query_set = TRUE, Q = TRUE
) %>% head()


# BoxPlot -----------------------------------------------------------------
# View plot
BoxPlot(fit_Hotgenes)

# list available data
ExpressionSlots_(fit_Hotgenes)




# ExpsPlot ----------------------------------------------------------------

yvar <- c("CSF2", "CXCL8", "IL6")
xvar <- "Hrs"

# name_col param must be available in Mapper_()
M_Plot_aliases <- ExpsPlot(fit_Hotgenes,
                           xVar = xvar,
                           color = "Hrs",
                           yVar = yvar,
                           scales = "fixed", 
                           name_col = "ensembl_id"
)

testthat::expect_true("ensembl_id" %in% names(M_Plot_aliases$data))




# Test keep_na_padj parameter --------------------------------------------

testthat::test_that("keep_na_padj parameter works correctly", {
  
  # Create a modified fit_Hotgenes with some NA padj values
  fit_Hotgenes_na <- fit_Hotgenes
  
  # Get the first contrast and introduce some NA padj values
  contrast_name <- names(fit_Hotgenes_na@Output_DE)[1]
  modified_de <- fit_Hotgenes_na@Output_DE[[contrast_name]]
  
  # Set padj to NA for the first 5 features
  modified_de$padj[1:5] <- NA
  
  # Put it back
  fit_Hotgenes_na@Output_DE[[contrast_name]] <- modified_de
  
  # Test 1: keep_na_padj = FALSE (default) should exclude NA padj values
  result_exclude_na <- Output_DE_(
    fit_Hotgenes_na,
    contrasts = contrast_name,
    padj_cut = 0.1,
    keep_na_padj = FALSE,
    mapFeatures = FALSE
  )
  
  testthat::expect_true(
    all(!is.na(result_exclude_na$padj)),
    info = "keep_na_padj = FALSE should exclude all NA padj values"
  )
  
  # Test 2: keep_na_padj = TRUE should include NA padj values
  result_include_na <- Output_DE_(
    fit_Hotgenes_na,
    contrasts = contrast_name,
    padj_cut = 0.1,
    keep_na_padj = TRUE,
    mapFeatures = FALSE
  )
  
  testthat::expect_true(
    any(is.na(result_include_na$padj)),
    info = "keep_na_padj = TRUE should include NA padj values"
  )
  
  # Test 3: Verify that NA features are actually included
  na_feature_count <- sum(is.na(result_include_na$padj))
  testthat::expect_gt(
    na_feature_count, 0,
    label = "Should have at least one feature with NA padj when keep_na_padj = TRUE"
  )
  
  # Test 4: keep_na_padj = FALSE should have fewer rows than TRUE
  testthat::expect_lt(
    nrow(result_exclude_na),
    nrow(result_include_na),
    label = "Excluding NAs should result in fewer rows"
  )
  
  # Test 5: Test Output_DE_df with keep_na_padj
  result_df_exclude <- Output_DE_df(
    fit_Hotgenes_na,
    contrasts = contrast_name,
    padj_cut = 0.1,
    keep_na_padj = FALSE
  )
  
  result_df_include <- Output_DE_df(
    fit_Hotgenes_na,
    contrasts = contrast_name,
    padj_cut = 0.1,
    keep_na_padj = TRUE
  )
  
  testthat::expect_true(
    all(!is.na(result_df_exclude$padj)),
    info = "Output_DE_df with keep_na_padj = FALSE should exclude NAs"
  )
  
  testthat::expect_true(
    any(is.na(result_df_include$padj)),
    info = "Output_DE_df with keep_na_padj = TRUE should include NAs"
  )
  
  # Test 6: Verify consistency between Output_DE_ and Output_DE_df
  testthat::expect_equal(
    nrow(result_exclude_na),
    nrow(result_df_exclude),
    info = "Output_DE_ and Output_DE_df should return same number of rows"
  )
  
})

# Test Ranks report with NA padj handling --------------------------------

testthat::test_that("Ranks report handles infinite values correctly", {
  
  # Create a modified fit_Hotgenes with some NA padj values
  fit_Hotgenes_imputing <- fit_Hotgenes
  
  # Get the first contrast and introduce some NA padj values
  contrast_name <- names(fit_Hotgenes_imputing@Output_DE)[1]
  modified_de <- fit_Hotgenes_imputing@Output_DE[[contrast_name]]
  
  # 
  testthat::expect_true(
    all(is.finite(modified_de$stat)),
    label = "These should all be finite as the downstream code will add infinites"
  )
  
  index_feature_mods <- nrow(modified_de)
  
  mid_feature <- floor(index_feature_mods/2)
  
  fea_index <- seq(mid_feature, mid_feature+5)
  
  # Set padj to NA for the first 5 features
  modified_de$padj[fea_index] <- 0
  
  mod_these <- modified_de$Feature[fea_index] 
  
  modified_de <- modified_de %>% 
  dplyr::mutate(stat = (-log10(.data$padj) * sign(.data$log2FoldChange)))
  
  # Put it back
  fit_Hotgenes_imputing@Output_DE[[contrast_name]] <- modified_de
  
  # # Test 1: Ranks report should exclude features with NA in rank column
  
  
  Imputed_ranks <- testthat::expect_warning({
    DE(
      fit_Hotgenes_imputing,
      Report = "Ranks",
      contrasts = contrast_name,
      Rank_name = "Feature",
      padj_cut = 1,
      impute_finite_ranks = TRUE  # Even with TRUE, Ranks should filter NAs
    )
    
  },
    regexp = "Inputing max"
  )
  
  
  
  Not_imputed <- testthat::expect_warning({
    DE(
      fit_Hotgenes_imputing,
      Report = "Ranks",
      contrasts = contrast_name,
      Rank_name = "Feature",
      padj_cut = 1,
      impute_finite_ranks = FALSE  # Even with TRUE, Ranks should filter NAs
    )
    
  },
  regexp = "No imputation"
  )
  
  
 
  # Extract the ranks vector for the contrast
  
  
  testthat::expect_true(
    any(is.infinite(Not_imputed[[contrast_name]])),
    info = "some should be infinite"
  )
  
  
  testthat::expect_true(
    all(is.finite(Imputed_ranks[[contrast_name]])),
    info = "all should be finite"
  )
  
  # Test 2: Verify that NA padj features are filtered out from ranks
  # (They should be removed by the na_check in the Ranks code)
  ranks_features <- names(Not_imputed[[contrast_name]][is.infinite(Not_imputed[[contrast_name]])])
  
  

  # Check if any of the NA padj features made it into ranks
  # (they shouldn't if stat values are also problematic)
  testthat::expect_true(
    all(ranks_features %in% mod_these),
    info = "Ranks should return a numeric vector with length > 0"
  )
  
  
  
  
  
  
})

# Test keep_na_padj with hotList -----------------------------------------

testthat::test_that("keep_na_padj works correctly with hotList", {
  
  # Create a modified fit_Hotgenes with some NA padj values
  fit_Hotgenes_na <- fit_Hotgenes
  
  # Get the first contrast
  contrast_name <- names(fit_Hotgenes_na@Output_DE)[1]
  modified_de <- fit_Hotgenes_na@Output_DE[[contrast_name]]
  
  # Set padj to NA for the first 5 features and keep their names
  na_features <- modified_de$Feature[1:5]
  modified_de$padj[1:5] <- NA
  
  # Put it back
  fit_Hotgenes_na@Output_DE[[contrast_name]] <- modified_de
  
  # Test 1: hotList with keep_na_padj = TRUE should include NA features
  result_with_na <- Output_DE_(
    fit_Hotgenes_na,
    contrasts = contrast_name,
    hotList = na_features,
    padj_cut = 0.1,
    keep_na_padj = TRUE,
    annotateSig = TRUE,
    mapFeatures = FALSE
  )
  
  testthat::expect_equal(
    nrow(result_with_na), 5,
    info = "hotList should return all requested features when keep_na_padj = TRUE"
  )
  
  testthat::expect_true(
    any(is.na(result_with_na$padj)),
    info = "hotList with keep_na_padj = TRUE should include NA padj values"
  )
  
  # Test 2: hotList with keep_na_padj = FALSE should exclude NA features
  result_without_na <- Output_DE_(
    fit_Hotgenes_na,
    contrasts = contrast_name,
    hotList = na_features,
    padj_cut = 0.1,
    keep_na_padj = FALSE,
    annotateSig = TRUE,
    mapFeatures = FALSE
  )
  
  testthat::expect_equal(
    nrow(result_without_na), 5,
    info = "hotList should include features with NA padj even when keep_na_padj = FALSE"
  )
  
  # Test 3: Test significance annotation with NA padj values
  testthat::expect_true(
    "significant" %in% colnames(result_with_na),
    info = "annotateSig should create 'significant' column"
  )
  
  # Features with NA padj should be marked as significant
  na_rows <- is.na(result_with_na$padj)
  testthat::expect_true(
    all(result_with_na$significant[na_rows] == "*"),
    info = "Features with NA padj should be marked as significant"
  )
  
})
