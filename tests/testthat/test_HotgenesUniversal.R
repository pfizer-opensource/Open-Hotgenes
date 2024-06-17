# creating new Hotgenes object --------------------------------------------

# Getting example named list of DE statistics
NewDE <- Output_DE_(htgs, as_list = TRUE, padj_cut = 1)
NewDE %>% purrr::map(~ head(.))

# Getting example named list of normalized data
NormlData <- Normalized_Data_(htgs)
NormlData %>% purrr::map(~ head(.))

# Getting example coldata
ExpColdata <- coldata_(htgs)
ExpColdata

# Getting example original data object used for DE analysis
# This example was generated from DESeq2
OrigDEObj <- O_(htgs)
OrigDEObj %>% class()

# Getting example design matrix
DE_design <- designMatrix_(htgs)
DE_design

# Getting example mapper
MapperDF <- Mapper_(htgs)
MapperDF %>% head()

# Converting example objects to hotgenes
Hotgenes_Object <- HotgenesUniversal(
  Output_DE = NewDE,
  Normalized_Expression = NormlData,
  coldata = ExpColdata,
  Original_Object = OrigDEObj,
  designMatrix = DE_design,
  Mapper = MapperDF
)



# Verifying tests for hotgenes object generation --------------------------
# This checks if rownames have blanks
NormlData2 <- NormlData

# Replace name with blank
rownames(NormlData2[[1]])
rownames(NormlData2[[1]])[255] <- " "
rownames(NormlData2[[1]])

# expect an error for empty rowname
testthat::expect_error(
  HotgenesUniversal(
    Output_DE = NewDE,
    Normalized_Expression = NormlData2,
    coldata = ExpColdata,
    Original_Object = NULL,
    designMatrix = DE_design,
    Mapper = MapperDF
  ),
  "Empty name detected! Check your expression matrix rownames"
)

# This verifies that list of normalized data is named
NormlData3 <- NormlData

# remove list names
names(NormlData3) <- NULL

# expect an error for empty names
testthat::expect_error(
  HotgenesUniversal(
    Output_DE = NewDE,
    Normalized_Expression = NormlData3,
    coldata = ExpColdata,
    Original_Object = NULL,
    designMatrix = DE_design,
    Mapper = MapperDF
  ),
  "Normalized_Expression elements must be named"
)

# This verifies that coldata matches normalized expression data -----------
# coldata is reordered to not match expression data
ExpColdata2 <- ExpColdata %>%
  dplyr::arrange(.data$Time)

# expect an error for resorted coldata
testthat::expect_error(
  HotgenesUniversal(
    Output_DE = NewDE,
    Normalized_Expression = NormlData,
    coldata = ExpColdata2,
    Original_Object = NULL,
    designMatrix = DE_design,
    Mapper = MapperDF
  ),
  "Sample IDs from coldata and expression data aren't equal"
)

# This verifies that at least one column contains factors
ExpColdata_characters <- ExpColdata %>%
  dplyr::mutate_all(as.character)

# expect an error for not having at least one factor
testthat::expect_error(
  HotgenesUniversal(
    Output_DE = NewDE,
    Normalized_Expression = NormlData,
    coldata = ExpColdata_characters,
    Original_Object = NULL,
    designMatrix = DE_design,
    Mapper = MapperDF
  ),
  "coldata must have at least one factor"
)


# Verifies mapper ---------------------------------------------------------
MapperDF2 <- MapperDF %>%
  dplyr::select(-Feature)

# expect an error for missing Feature col
testthat::expect_error(
  Mapper_(Hotgenes_Object) <- MapperDF2,
  "Problem with `Feature`"
)


# Verifies output_DE ------------------------------------------------------
NewDE2 <- NewDE %>%
  purrr::map(~ .x %>%
               dplyr::select(-padj))

# expect an error for missing padj column
testthat::expect_error(
  HotgenesUniversal(
    Output_DE = NewDE2,
    Normalized_Expression = NormlData,
    coldata = ExpColdata,
    Original_Object = NULL,
    designMatrix = DE_design,
    Mapper = MapperDF
  ),
  regexp = "padj"
)

# expect an error for "Output_DE elements must be named"
NewDE3 <- NewDE
names(NewDE3) <- NULL

testthat::expect_error(
  HotgenesUniversal(
    Output_DE = NewDE3,
    Normalized_Expression = NormlData,
    coldata = ExpColdata,
    Original_Object = NULL,
    designMatrix = DE_design,
    Mapper = MapperDF
  ),
  "Output_DE elements must be named"
)


# This checks if Feature have blanks
NewDE_error <- NewDE

# Replace name with blank
NewDE_error$sh_EWS_vs_Ctrl
NewDE_error$sh_EWS_vs_Ctrl$Feature[5] <- " "
NewDE_error$sh_EWS_vs_Ctrl

# expect an error for empty Feature
testthat::expect_error(
  HotgenesUniversal(
    Output_DE = NewDE_error,
    Normalized_Expression = NormlData,
    coldata = ExpColdata,
    Original_Object = NULL,
    designMatrix = DE_design,
    Mapper = MapperDF
  ),
  "Empty Feature detected! Check your Output_DE slot"
)
