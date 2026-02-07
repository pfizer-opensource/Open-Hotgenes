
# getting a mapper 
alt_Mapper <- Mapper_(htgs) %>% 
  dplyr::slice(1, .by = "Feature") %>% 
  dplyr::rename("egid" = "ensembl_id")


# preparing to randomly remove some rows

index_seq <- seq(1,nrow(alt_Mapper))
set.seed(12)

index_replace <- sample(index_seq, size = 50)

# creating a gnm column
alt_Mapper <- alt_Mapper %>% 
  dplyr::mutate(gnm = .data$Feature)
  
# replacing rows with NA
alt_Mapper[index_replace,"gnm"] <- NA

alt_Mapper[index_replace,"egid"] <- NA
  
# verifying NAs
testthat::expect_true(all(is.na(alt_Mapper$gnm[index_replace])))
testthat::expect_true(all(is.na(alt_Mapper$egid[index_replace])))

# verifying identical NA missing rows
testthat::expect_identical(alt_Mapper$gnm[index_replace],
                           alt_Mapper$egid[index_replace])



# Filling in NAs from Feature
alt_Mapper_updated <- alt_Mapper %>% 
  set_missing_alias_default(preferred_alias = c("egid", "gnm"))

# these should be identical because they should all match Feature
testthat::expect_identical(alt_Mapper_updated$gnm[index_replace],
                           alt_Mapper_updated$egid[index_replace])

# verifying that they match Feature
testthat::expect_identical(alt_Mapper_updated$Feature[index_replace],
                           alt_Mapper_updated$egid[index_replace])

# verifying that they match Feature
testthat::expect_identical(alt_Mapper_updated$gnm[index_replace],
                           alt_Mapper_updated$Feature[index_replace])

# Should not have NA
testthat::expect_false(all(is.na(alt_Mapper_updated$gnm[index_replace])))
testthat::expect_false(all(is.na(alt_Mapper_updated$egid[index_replace])))

testthat::expect_false(all(is.na(alt_Mapper_updated$gnm)))
testthat::expect_false(all(is.na(alt_Mapper_updated$egid)))


# verifying behavior

# no change if preferred_alias is NULL
alt_Mapper_updated_2 <- alt_Mapper %>% 
  set_missing_alias_default(preferred_alias = NULL)

testthat::expect_identical(alt_Mapper_updated_2,
                           alt_Mapper)

# no change if default_alias is NULL
alt_Mapper_updated_3 <- alt_Mapper %>% 
  set_missing_alias_default(preferred_alias = c("egid", "gnm"),
                            default_alias = NULL)

testthat::expect_identical(alt_Mapper_updated_3,
                           alt_Mapper)

# expect error if no match returned by match.arg
testthat::expect_error({
  
  alt_Mapper_updated_4 <- alt_Mapper %>% 
    set_missing_alias_default(preferred_alias = c("gene_symbol"),
                              default_alias = "Feature")
}, regexp = "should be one of")



# testing  ----------------------------------------------------------------

main_Mapper <- Mapper_(htgs) 

Out_df_1 <- Output_DE_(htgs, mapFeatures = FALSE, padj_cut = 1)

Out_df_2 <- Output_DE_(htgs, mapFeatures = TRUE, padj_cut = 1)

Out_df_3 <-  df_append_mapperDF(
  df = Out_df_1,
  mapperDF = main_Mapper, 
  by = "Feature",
  relationship = "many-to-many")

testthat::expect_identical(Out_df_3, Out_df_2)
