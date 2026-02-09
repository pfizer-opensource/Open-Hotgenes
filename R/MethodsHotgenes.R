# Hotgenes class  ---------------------------------------------------------
#' Hotgenes class
#'
#' @importFrom methods is new
#' @slot Output_DE Named list of data.frames
#' Each data.frame contains atleast: "Feature", "baseMean","log2FoldChange",
#' "stat", "pvalue", "padj".
#' Additional values dependent on original DE analysis.
#' Accessible via Output_DE_.
#' @slot Normalized_Expression named list of data.frames available
#' expression data. Names of available slots accessible via ExpressionSlots_.
#' Expression data can be accessed via Normalized_Data_.
#' @slot coldata The coldata, as a data.frame with Sample Metadata.
#' Accessible via coldata_.
#' @slot designMatrix The model matrix used for DE analysis.
#' Accessible via designMatrix_.
#' @slot contrastMatrix The contrast matrix used for DE analysis.
#' Accessible via contrastMatrix_.
#' @slot Original_Object This is the original object imported into Hotgenes.
#' Accessible via O_.
#' @slot auxiliary_assays A data.frame with auxiliary assays 
#' @slot Mapper Slot for a data.frame that can be used
#' to map aliases to original assay names. Must be a
#' data.frame containing a "Feature" column,
#' which will be used for mapping to results. This object
#' can be updated or viewed using Mapper_
#' @example man/examples/HotgenesClass_Example.R

setClass("Hotgenes", slots = list(
  Output_DE = "list",
  Normalized_Expression = "list",
  coldata = "data.frame",
  auxiliary_assays = "data.frame",
  designMatrix = "matrix",
  contrastMatrix = "ANY",
  Original_Object = "ANY",
  Mapper = "data.frame"
))


# object setValidity ------------------------------------------------------


setValidity("Hotgenes", function(object) {
  
  # Check Output_DE
  # verifies the presence of required columns
  # verifies that requires statistical data columns are numeric
  # verifies that no Features are empty
  object %>% Output_DE_check()

  # check Features
  # Get rownames from all expression data
  # verifies that rownames do not have empty names
  object %>% Features_check()

  # check ExpressionSlots_
  # check if Normalized_Expression slot is a named list
  object %>% expression_data_check()
  
  # check auxiliary_assays
  # must be data.frame and have SampleIDs column
  object %>% auxiliary_assays_check()
  
  # check designMatrix 
  # check designMatrix col length > 0 
  # verify designMatrix has rownames and
  # that they match with exps data
  object %>% designMatrix_check()
  
  
  # check coldata
  # check coldata row lengths
  # verify coldata has rownames and
  # that they match with exps data
  object %>% coldata_check()

  # check Mapper slot
  # verifies that it has "Feature" column
  object %>% Mapper_()
  TRUE
})


#' update_object updates object class
#' @importFrom BiocGenerics updateObject getObjectSlots
#' @param Hotgenes Hotgenes object.
#' @export

update_object <- function(Hotgenes = NULL){
  
 
  
  Hotgenes<- HotgenesUniversal(Output_DE = Hotgenes@Output_DE,
                    Normalized_Expression = Hotgenes@Normalized_Expression,
                    auxiliary_assays = Hotgenes@auxiliary_assays,
                    coldata = Hotgenes@coldata,
                    Original_Object = Hotgenes@Original_Object,
                    designMatrix = Hotgenes@designMatrix,
                    contrastMatrix  = Hotgenes@contrastMatrix,
                    Mapper = Hotgenes@Mapper)
    


  return(Hotgenes)
  
}


# Internal functions ------------------------------------------------------
#'
#'
auxiliary_assays_check <- function(Hotgenes = NULL) {
  
  aux_df <- Hotgenes@auxiliary_assays
  
  if(isFALSE(is.data.frame(aux_df))){
    stop("auxiliary_assays slot must contain a data.frame")
  }
  
  if(isFALSE("SampleIDs" %in% colnames(aux_df))){
    stop("data.frame must contain a SampleIDs column")
  }
  
  # verify at least some overlap
  sample_check <- aux_df %>% 
    dplyr::filter(.data$SampleIDs %in% SampleIDs_(Hotgenes))
  
  
  # verify the SampleIDs column
  if (nrow(sample_check) == 0) {
    warning("no overlap with metadata")
  }
  
  # checks for duplicate names with coldata slots
  aux_df_features <- aux_df %>% 
    dplyr::select(-dplyr::any_of(c("SampleIDs"))) %>% 
    names()
  
  overl_names <- aux_df_features[aux_df_features %in% names(Hotgenes@coldata)]
  
  if(length(overl_names) > 0){

    message(
        glue::glue(
          "{glue::glue_collapse(overl_names, sep = ', ')} 
           not added, since already in coldata"
          ))
  }
  
  # checks for duplicate names with expression slots
  overl_features <- aux_df_features[aux_df_features %in% Features_(Hotgenes)]
  
  if(length(overl_features) > 0){
    message(
      glue::glue(
        "{glue::glue_collapse(
        overl_features, sep = ', ')
        }  not added, since already in Expression"))
  }
  

  
}

# checks normalized data
expression_data_check <- function(Hotgenes = NULL){
  
  
  # Check object
  stopifnot(is(Hotgenes, "Hotgenes"))
  
  
  # check if Normalized_Expression has names
  if (is.null(names(Hotgenes@Normalized_Expression))) {
    stop("Normalized_Expression elements must be named")
  }
  
  
  
}

#' internal function returns Mapper_() contents as named list
#' @noRd
Mapper_list <- function(Hotgenes = NULL){
  
  Mapper_(Hotgenes) %>% 
    purrr::imap(~unique(.x))
}

#'
#'
designMatrix_check  <- function(Hotgenes = NULL) {
  
 # Hotgenes_Object@designMatrix
  
  # check coldata row lengths
  if (nrow(Hotgenes@designMatrix) == 0) {
    stop("designMatrix slot has no rows")
  }
  
  # verify designMatrix has rownames
  
  if (!tibble::has_rownames(
    base::as.data.frame(Hotgenes@designMatrix) 
    )) {
    stop("designMatrix slot has no row names")
  }
  
  # verify that all output_DE names are in 
  # designMatrix has atleast one column
  # of factors
  check_cols_names <- Hotgenes@designMatrix %>%
    colnames() 
  
  length_cols <- check_cols_names %>%
    length()
  
  if (length_cols == 0) {
    stop("designMatrix must have at least one column")
  }
  
  
  # verify match with exps data
  Normalized_Data_(Hotgenes) %>%
    purrr::map(function(x) {
      Exps_SampleNames <- colnames(x)
      
      coldata_SampleNames <- rownames(Hotgenes@designMatrix)
      if (!all(coldata_SampleNames == Exps_SampleNames)) {
        stop("Sample IDs from designMatrix and expression data aren't equal")
      } else if (all(coldata_SampleNames == Exps_SampleNames)) {
        return("Sample IDs from designMatrix and expression data are equal")
      }
    })
}

#'
#'
coldata_check <- function(Hotgenes = NULL) {
  # check coldata row lengths
  if (nrow(Hotgenes@coldata) == 0) {
    stop("coldata slot has no rows")
  }
  
  # verify coldata has rownames
  if (!has_rownames(Hotgenes@coldata)) {
    stop("coldata slot has no row names")
  }
  
  # can't have SampleIDs in coldata
  if ("SampleIDs" %in% base::colnames(Hotgenes@coldata)) {
    stop("SampleIDs can't be in coldata")
  }
  
  # verify coldata has atleast one column
  # of factors
  check_factor <- Hotgenes@coldata %>%
    dplyr::select_if(is.factor) %>%
    names() %>%
    length()
  
  if (check_factor == 0) {
    stop("coldata must have at least one factor")
  }
  
  # verify match with exps data
  Normalized_Data_(Hotgenes) %>%
    purrr::map(function(x) {
      Exps_SampleNames <- colnames(x)
      
      coldata_SampleNames <- rownames(Hotgenes@coldata)
      if (!all(coldata_SampleNames == Exps_SampleNames)) {
        stop("Sample IDs from coldata and expression data aren't equal")
      } else if (all(coldata_SampleNames == Exps_SampleNames)) {
        return("Sample IDs from coldata and expression data are equal")
      }
    })
}


#'
#'
Features_check <- function(Hotgenes = NULL) {
  Normalized_Data_(Hotgenes) %>%
    purrr::map(function(x) {
      Features_names <- rownames(x) %>%
        sort()
      
      # Check for empty spaces
      if (any_missingness(Features_names)) {
        stop("Empty name detected! Check your expression matrix rownames")
      } else {
        return(Features_names)
      }
    })
}


#' internal required DE table column names
#' 
required_DE_cols <- function(){
  DE_Reqs <- c(
    "Feature", "contrast_dir",
    "baseMean", "log2FoldChange", "FC",
    "stat", "pvalue", "padj", "contrast"
  )
  
  return(DE_Reqs)
}

#' internal required DE table column names
#' 
required_DE_numeric_cols <- function(){
  
  
  numeric_Reqs <- c(
    
    "baseMean", "log2FoldChange",
    "stat", "pvalue", "padj"
  )
 
  return(numeric_Reqs)
}
#'
#'
Output_DE_check <- function(Hotgenes = NULL) {
  # check if Output_DE has names
  if (is.null(names(Hotgenes@Output_DE))) {
    stop("Output_DE elements must be named")
  }
  
  # verifies that require col are present
  # sets required cols
  DE_Reqs <- required_DE_cols()
  
  # Check Output_DE cols
  DE_Check <- Hotgenes@Output_DE %>%
    purrr::map(function(x) {
      all(DE_Reqs %in% names(x))
    }) %>%
    unlist(use.names = FALSE) %>%
    unique()
  
  if (isFALSE(all(DE_Check))) {
    stop(paste0(
      "All data.frames in Output_DE must contain: ",
      paste0(DE_Reqs, collapse = ", ")
    ))
  }
  
  # verify numeric
  # Check Output_DE cols
  
  numeric_Reqs <- required_DE_numeric_cols()
 
  
  numeric_Check <- Hotgenes@Output_DE %>%
    purrr::map(function(x) {
      
      dat_x <- x %>% 
        dplyr::select_if(is.numeric)
      
      
      all(numeric_Reqs %in% names(dat_x))
    }) %>%
    unlist(use.names = FALSE) %>%
    unique()
  
  # numeric_Check
  if (isFALSE(all(numeric_Check))) {
    stop(
      "Missing numeric columns in Output_DE"
    )
  }
  
  # check for empties
  DE_Check_blanks <- Hotgenes@Output_DE %>%
    purrr::map_dfr(function(x) {
      x %>%
        dplyr::filter(detected_missingness(.data$Feature))
    }) %>%
    nrow()
  
  if (DE_Check_blanks > 0) {
    stop("Empty Feature detected! Check your Output_DE slot")
  }
}

# show methods ------------------------------------------------------------
setMethod(
  "show", "Hotgenes",
  function(object) {
    cat("class:", class(object), "\n")

    cat("Original class/package: ", stringr::str_c(
      class(O_(object)), "/",
      attr(class(O_(object)), "package")
    ))

    cat("\n\n")
    
    cat("Differential expression (default thresholds): ",
      knitr::kable(
        tibble::enframe(DE(object, Report = "Length")),
        format = "pipe",
        col.names = c("contrast", "total")
      ),
      sep = "\n"
    )

    cat("\n")

   
    cat("Available feature mapping: ",
      glue::glue("{names(Mapper_(object))  }") %>%
      glue::glue_collapse(sep = ", "), "\n")

    
    cat("ExpressionSlots: ", 
      glue::glue("{ExpressionSlots_(object)}") %>%
      glue::glue_collapse(sep = ", "), "\n")
    
    
    cat("Total auxiliary assays: ", 
        glue::glue("{length(auxiliary_assays_features(object))}") %>%
          glue::glue_collapse(sep = ", "), "\n")

    
    cat("Total samples: ", 
        glue::glue("{length(SampleIDs_(object))}"), "\n")


  }
)



# support auxiliary_assays_ -----------------------------------------------
# Internal
# removes cols from new_auxiliary_assays that are already in other slots
# also removes empty cols
auxiliary_assays_filter <- function(Hotgenes = NULL, 
                                    new_auxiliary_assays = NULL) {
  
  valid_samples <- SampleIDs_(Hotgenes)
  
  aux_df <-new_auxiliary_assays 
  
  if(isFALSE(is.data.frame(aux_df))){
    stop("auxiliary_assays slot must contain a data.frame")
  }
  
  # exclude_these
  check_for_these <-c( names(Hotgenes@coldata),
                       Features_(Hotgenes)) %>% unique()
  exclude_these <- names(aux_df)[names(aux_df) %in% check_for_these]
  
  
  if(length(exclude_these) > 0){
    
    message(
      glue::glue(
        "{glue::glue_collapse(exclude_these, sep = ', ')} not added, 
        since already in coldata or Features slot"
      ))
  }
  
  # checks for duplicate names with coldata slots
  aux_df_features <- aux_df %>% 
    dplyr::select(-dplyr::any_of(c( exclude_these))) %>% 
    dplyr::filter(.data$SampleIDs %in% valid_samples) %>% 
    janitor::remove_empty(which = "cols")
  
  
  return(aux_df_features)
}
# Creates the initial aux assay dataframe
auxiliary_assays_frame <-function(SampleIDs = NULL){
  aux_initial <- data.frame(SampleIDs = SampleIDs) %>%
  tibble::as_tibble()
  
  return(aux_initial)
}


#' Returns features from auxiliary_assays slot
#' @rdname Hotgenes-class
#' @importFrom janitor remove_empty
#' @export

auxiliary_assays_features <- function(Hotgenes = NULL) {
  assays_names <- Hotgenes@auxiliary_assays %>%
    dplyr::select(-"SampleIDs") %>% 
    colnames()
  
  return(assays_names)
}
#' @rdname Hotgenes-class
#' @export
#' @param Hotgenes Hotgenes object.
#' @return auxiliary_assays_default() returns
#' a data.frame with the minimum requirements for a
#' Hotgenes object. The coldata_ids parameter can be used to
#' define metadata columns to include in the data.frame, to help
#' align data. 
#' 
auxiliary_assays_default <- function(
    Hotgenes = NULL,
    coldata_ids = NULL) {
  
  # Check object
  stopifnot(is(Hotgenes, "Hotgenes"))
  
  # making clean auxiliary_assays_slot data
  auxiliary_assays_empty <- data.frame(
    SampleIDs = SampleIDs_(Hotgenes)) %>%
    tibble::as_tibble() 
  
  auxiliary_assays_empty <- coldata_(Hotgenes = Hotgenes,
                         coldata_ids = coldata_ids) %>% 
    tibble::rownames_to_column("SampleIDs") %>% 
    dplyr::left_join(auxiliary_assays_empty, by = "SampleIDs")
  
  
  
  
  return(auxiliary_assays_empty)
  
  
}

#' @rdname Hotgenes-class
#' @export
#' @param Hotgenes Hotgenes object.
#' @param features vector of feature to return. If NULL (default)
#' all data will be returned. No need to specifiy "SampleIDs".
#' @return auxiliary_assays_ returns a data.frame
#' containing auxiliary assays to complement expression data
auxiliary_assays_ <- function(Hotgenes = NULL, features = NULL,
                              mode = "all") {
  
  
  
  # pulling Mapper slot
  # new merges to Feature
  auxiliary_assays_Slot <- auxiliary_assays_default(Hotgenes = Hotgenes) %>%
    dplyr::left_join(Hotgenes@auxiliary_assays, by = "SampleIDs") %>% 
    tibble::column_to_rownames("SampleIDs")
  
  if(!is.null(features)){
    #complete_features <- c("SampleIDs", features) %>% unique()
    complete_features <- features %>% unique()
    
    auxiliary_assays_Slot <- auxiliary_assays_Slot %>% 
      dplyr::select(dplyr::any_of(c(complete_features)))
  }
  
  # check mode
  selected_mode <- match.arg(mode, c("all", "quanti", "quali"),
                             several.ok = FALSE
  )
  
  if (selected_mode == "quanti") {
    auxiliary_assays_Slot <- auxiliary_assays_Slot %>%
      dplyr::select_if(~ dplyr::n_distinct(na.omit(.x)) > 1) %>%
      dplyr::select_if(is.numeric)
  } else if (selected_mode == "quali") {
    auxiliary_assays_Slot <- auxiliary_assays_Slot %>%
      dplyr::select_if(~ dplyr::n_distinct(na.omit(.x)) > 1) %>%
      dplyr::select_if(is.factor)
  }
  
  return(auxiliary_assays_Slot)
  
  
}

#' @rdname Hotgenes-class
#' @export
#' @importFrom methods validObject
#' @param object Hotgenes object.
#' @param value a data.frame containing a "Feature"
#' for alias mapping.
"auxiliary_assays_<-" <- function(object, value) {
  object
}
setGeneric("auxiliary_assays_<-")

#' @rdname Hotgenes-class
#' @export
setMethod(
  "auxiliary_assays_<-", signature(object = "Hotgenes"),
  function(object, value) {
    
    new_obj <-auxiliary_assays_filter(Hotgenes = object, 
                                      new_auxiliary_assays = value)
    
    object@auxiliary_assays <- new_obj
    
    validObject(object)
    
    object
  }
)


# Returns original import object ------------------------------------------



#' @rdname Hotgenes-class
#' @export
#' @param Hotgenes Hotgenes object.
#' @return O_ Returns the original R object imported into Hotgenes.
#' For HotgenesDESeq2 this a DESeqDataSet class object.
#' For HotgenesMonocle3 this is a cell_data_set class object.
#' For Hotgeneslimma this is an EList class object (Expression data),
#' with the limmafit (MArrayLM class object) add in in the "limmafit"
#' slot.
O_ <- function(Hotgenes = NULL) {
  # Check object
  stopifnot(is(Hotgenes, "Hotgenes"))
  return(Hotgenes@Original_Object)
}


# Returns design/DE setup objects -----------------------------------------

#' @rdname Hotgenes-class
#' @export
#' @md
#' @param Hotgenes Hotgenes object.
#' @param min_col_levels numeric value used to exclude
#' columns that do not contain unique values greater than this.
#' @param coldata_ids Variable(s) stored in coldata slot can be selected
#' by column name, using a vector strings. Default is [coldata_names()]. If
#' set to "", an empty data.frame is returned.
#' @param mode string indicating the type of coldata to call.
#' Choices include: "all", "quanti", "quali". If "all", the
#' default, all columns will be returned. If "quanti", only
#' numeric coldata will be returned. If "quali", only
#' qualitative coldata will be returned.
#' @return coldata_ Returns Coldata (Experimental Variables).
#' @details Except for when mode = "all", single level columns will be
#' dropped.
coldata_ <- function(Hotgenes = NULL,
                     coldata_ids = coldata_names(Hotgenes),
                     min_col_levels = 0,
                     mode = "all") {
  # Check object
  stopifnot(is(Hotgenes, "Hotgenes"))
  
  
  
  coldata_Out <- Hotgenes@coldata %>%
    dplyr::select(dplyr::any_of(coldata_ids))
  
  # check mode
  selected_mode <- match.arg(mode, c("all", "quanti", "quali"),
                             several.ok = FALSE
  )
  
  
  if (selected_mode == "quanti") {
    coldata_Out <- coldata_Out %>%
      dplyr::select_if(~ dplyr::n_distinct(.x,
                                           na.rm = TRUE
      ) > min_col_levels) %>%
      dplyr::select_if(is.numeric)
  } else if (selected_mode == "quali") {
    coldata_Out <- coldata_Out %>%
      dplyr::select_if(~ dplyr::n_distinct(.x,
                                           na.rm = TRUE
      ) > min_col_levels) %>%
      dplyr::select_if(is.factor)
  }
  
  return(coldata_Out)
}


#' Returns colnames for coldata slot
#' @rdname Hotgenes-class
#' @export

coldata_names <- function(Hotgenes = NULL) {
  coldata_names <- Hotgenes@coldata %>%
    colnames()
  
  return(coldata_names)
}

#' @rdname Hotgenes-class
#' @export
#' @importFrom methods validObject
#' @param object Hotgenes object.
#' @param value a data.frame containing a "Feature"
#' for alias mapping.
"coldata_<-" <- function(object, value) {
  object
}
setGeneric("coldata_<-")

#' @rdname Hotgenes-class
#' @export
setMethod(
  "coldata_<-", signature(object = "Hotgenes"),
  function(object, value) {
    object@coldata <- value
    
    validObject(object)
    
    object
  }
)


#' @rdname Hotgenes-class
#' @export
#' @param Hotgenes Hotgenes object.
#' @return designMatrix_ Returns the model matrix used
#' for the original DE analysis (if available).
designMatrix_ <- function(Hotgenes = NULL) {
  # Check object
  stopifnot(is(Hotgenes, "Hotgenes"))
  return(Hotgenes@designMatrix)
}


#' @rdname Hotgenes-class
#' @export
#' @param Hotgenes Hotgenes object.
#' @return contrastMatrix_ Returns the contrast matrix used
#' for the original DE analysis (if available).
contrastMatrix_ <- function(Hotgenes = NULL) {
  # Check object
  stopifnot(is(Hotgenes, "Hotgenes"))
  return(Hotgenes@contrastMatrix)
}


# Returns expression data slot and related details ------------------------



#' @rdname Hotgenes-class
#' @export
#' @param Hotgenes Hotgenes object.
#' @param slot Single name of expression data.frame in
#' Normalized_Expression slot to be returned.
#' @return Named list from Normalized_Expression slot
#' if "slot" is NULL.
#' If slot argument is provided, a data.frame is returned.
Normalized_Data_ <- function(Hotgenes = NULL,
                             slot = NULL) {
  # Check object
  stopifnot(is(Hotgenes, "Hotgenes"))
  
  # pulling Normalized_Expression slot
  NormDat <- Hotgenes@Normalized_Expression
  
  # get slot names
  slotNames <- names(NormDat)
  
  if (is.null(slot)) {
    NormDat
  } else if (!is.null(slot)) {
    slotID <- match.arg(slot, slotNames)
    
    NormDat %>%
      purrr::chuck(slotID) %>%
      base::as.data.frame()
  }
}

#' @rdname Hotgenes-class
#' @export
#' @param Hotgenes Hotgenes object.
#' @return SampleIDs_ Returns the colnames of the expression data in the
#' first ExpressionSlot.
SampleIDs_ <- function(Hotgenes = NULL) {
  # Check object
  stopifnot(is(Hotgenes, "Hotgenes"))
  
  allsamples <- coldata_(Hotgenes) %>% 
    base::rownames()
  
  return(allsamples)
}

#' @rdname Hotgenes-class
#' @export
#' @param Hotgenes Hotgenes object.
#' @return ExpressionSlots_ Returns the names for all the available
#' expression data. This depends on the original
#' R package used for DE analysis. For example,
#' DESeq2 data can be normalized by "rld" or "VST".
ExpressionSlots_ <- function(Hotgenes = NULL) {
  # Check object
  stopifnot(is(Hotgenes, "Hotgenes"))
  
  
  # check if Normalized_Expression has names
  if (is.null(names(Hotgenes@Normalized_Expression))) {
    stop("Normalized_Expression elements must be named")
  }

  return(names(Hotgenes@Normalized_Expression))
}

#' @rdname Hotgenes-class
#' @export
#' @param Hotgenes Hotgenes object.
#' @param slot Single name of expression data.frame in
#' Normalized_Expression slot to be returned. Default is
#' first string returned by ExpressionSlots_()
#' @return returns a data.frame
ExpressionData_ <- function(Hotgenes = NULL,
  slot = ExpressionSlots_(Hotgenes)[1]) {
  
  # Check object
  stopifnot(is(Hotgenes, "Hotgenes"))
  
  # pulling Normalized_Expression slot
  NormDat_List <- Normalized_Data_(Hotgenes)
  
  # get slot names
  slotNames <- ExpressionSlots_(Hotgenes)
  
  slotID <- match.arg(slot, slotNames)
  
  NormDat <- NormDat_List %>%
    purrr::chuck(slotID) %>%
    base::as.data.frame()
  
  return(NormDat)
  
}


# Returns DE stats results slot -------------------------------------------



#' @rdname Hotgenes-class
#' @export
#' @importFrom purrr list_rbind
#' @param Hotgenes Hotgenes object.
#' @param contrasts vector for DE data.frames in Output_DE slot to return.
#' See ?contrasts_. If NULL (default), all contrasts used.
#' @param mapFeatures If FALSE, no changes will be made. If TRUE (default),
#' features will be mapped to information provided in the mapper slot.
#' See details below and  \code{\link[Hotgenes]{Mapper_}}.
#' @param hotList vector of Features to select. This will override
#' all cut offs.
#' @param padj_cut numeric value for padj limit (0.1 is Default).
#' @param .log2FoldChange absolute log2FoldChange for filtering.
#' Default is 0.
#' @param as_list logical, default is FALSE. If true, results returned
#' as a named list other object with be a dataframe grouped by "contrast".
#' @param annotateSig logical, if TRUE (default), a column is added
#' to show significance for any hotList feature (based on cut offs).
#' @param keep_na_padj logical, if TRUE
#' any features with NA padj will be included with features that
#' meet padj_cut. If FALSE (default), only features that meet the padj_cut
#' threshold will be returned.
#' @return Output_DE_ Returns complete unfiltered Output_DE slot object.
Output_DE_ <- function(Hotgenes = NULL,
                       contrasts = NULL,
                       mapFeatures = TRUE,
                       hotList = NULL,
                       padj_cut = 0.1,
                       .log2FoldChange = 0,
                       annotateSig = TRUE,
                       as_list = FALSE,
                       keep_na_padj = FALSE) {
  # Check object
  stopifnot(is(Hotgenes, "Hotgenes"))
  
  
  # Build padj filter condition based on keep_na_padj
  padj_condition <- internal_padj_handling(padj_cut = padj_cut,
                                           keep_na_padj = keep_na_padj)
  
 

  # sets required cols
  DE_Reqs <- required_DE_cols()

  # getting available contrasts, excluding Intercept
  available_contrasts <- contrasts_(Hotgenes)

  if (is.null(contrasts)) {
    all_contrasts <- available_contrasts
  } else if (!is.null(contrasts)) {
    all_contrasts <- match.arg(contrasts, 
    available_contrasts, several.ok = TRUE)
  }



  # pulling slot
  Output_DE <- Hotgenes@Output_DE[all_contrasts] %>% 
    purrr::list_rbind(names_to = "contrast") 

 
  # check mapFeatures
  
  if (isTRUE(mapFeatures)) {
    
    mapped_check <- Mapper_aliases(Hotgenes)
    
    if(length(mapped_check) > 0 ) {
      
    
    mappingDF <- hotList_mapper(Hotgenes,
                                hotList = hotList)
    
    Output_DE <-  df_append_mapperDF(
        df = Output_DE,
        mapperDF = mappingDF, 
        by = "Feature",
        relationship = "many-to-many")
    }
    
  }

  Output_DE <- Output_DE %>%
    dplyr::relocate(dplyr::all_of(DE_Reqs), .before = 1) %>%
    dplyr::mutate(contrast = forcats::fct_relevel(
      .data$contrast,
      all_contrasts
    )) %>%
    dplyr::group_by(.data$contrast)
  
  
  if (isTRUE(annotateSig)) {
    Output_DE <- Output_DE %>%
      dplyr::mutate(significant = dplyr::case_when(
        eval(padj_condition) &
          abs(.data$log2FoldChange) >= .log2FoldChange ~ "*",
        TRUE ~ ""
      ))
  }


  if (all_truthy(hotList)) {
    #updateDetails <- TRUE

    # DF version


    Output_DE <- Output_DE %>%
       dplyr::filter(.data$Feature %in% hotList, .preserve = TRUE) 
    # tibble::as_tibble()

    # if (isTRUE(annotateSig)) {
    #   Output_DE <- Output_DE %>%
    #     dplyr::mutate(significant = dplyr::case_when(
    #        eval(padj_condition) &
    #         abs(.data$log2FoldChange) >= .log2FoldChange ~ "*",
    #       TRUE ~ ""
    #     ))
    # }
   
  } else {
   # updateDetails <- FALSE


    Output_DE <- Output_DE %>%
      dplyr::filter(abs(.data$log2FoldChange) >= .log2FoldChange,
        .preserve = TRUE
      ) %>%
      dplyr::filter( eval(padj_condition), .preserve = TRUE) 

  }

  if (isTRUE(as_list)) {
    Output_DE <- Output_DE %>%

      dplyr::group_map(~.x) %>%
      purrr::set_names(levels(Output_DE[[dplyr::group_vars(Output_DE)]])) %>%
      purrr::map(function(x) {
        x %>%
          tibble::as_tibble() %>%
          dplyr::select(-any_of("contrast"))
      })
  }

  return(Output_DE)
}

  
# internal
# returns DE results with minimal modification for faster performance
Output_DE_df <- function(Hotgenes = NULL,
                       contrasts = NULL,
                       hotList = NULL,
                       padj_cut = 0.1,
                       .log2FoldChange = 0,
                       keep_na_padj = FALSE) {
  # Check object
  stopifnot(is(Hotgenes, "Hotgenes"))
  
  
  # Build padj filter condition based on keep_na_padj
  padj_condition <- internal_padj_handling(padj_cut = padj_cut,
                                           keep_na_padj = keep_na_padj)

  
  # getting available contrasts, excluding Intercept
  available_contrasts <- contrasts_(Hotgenes)
  
  if (is.null(contrasts)) {
    all_contrasts <- available_contrasts
  } else if (!is.null(contrasts)) {
    all_contrasts <- match.arg(contrasts,
    available_contrasts, several.ok = TRUE)
  }
  
  
  
  # pulling slot
  Output_DE <- Hotgenes@Output_DE[all_contrasts] %>%
    # purrr::imap( ~.x) %>% 
    purrr::list_rbind(names_to = "contrast") 
  
  if(length(all_contrasts) > 1){
    
    Output_DE <- Output_DE   %>%
      dplyr::group_by(.data$contrast)
    
  }
 
  
  if (all_truthy(hotList)) {
    
    Output_DE <- Output_DE %>%
      dplyr::filter(.data$Feature %in% hotList, .preserve = TRUE)  
    
  } else {
    
    Output_DE <- Output_DE %>%
      dplyr::filter(abs(.data$log2FoldChange) >= .log2FoldChange,
                    .preserve = TRUE
      ) %>%
      dplyr::filter( eval(padj_condition), .preserve = TRUE) 
    
    
  }
  
  return(Output_DE)
}

#' @rdname Hotgenes-class
#' @export
#' @param Hotgenes Hotgenes object.
#' @return contrasts_ Returns exported contrasts.
contrasts_ <- function(Hotgenes = NULL) {
  # Check object
  stopifnot(is(Hotgenes, "Hotgenes"))

  # check if Output_DE has names
  if (is.null(names(Hotgenes@Output_DE))) {
    stop("Output_DE elements must be named")
  }
  return(names(Hotgenes@Output_DE))
}

#' @noRd
internal_padj_handling <- function(padj_cut = 0.1,
                                   keep_na_padj = FALSE){
  # Build padj filter condition based on keep_na_padj
  if (keep_na_padj) {
    padj_condition <- quote(
      detected_missingness(.data$padj) | .data$padj < padj_cut)
  } else {
    padj_condition <- quote( .data$padj < padj_cut)
  }
  return(padj_condition)
}


# Calls/updates Features/mapper -------------------------------------------


#' @rdname Hotgenes-class
#' @export
#' @param Hotgenes Hotgenes object.
#' @return Features_ Returns the rownames of the expression data in the
#' first ExpressionSlot. It will also check for empty rownames and
#' will return an error if any are detected. This will be used as
#' a last check before returning a hotgenes object during generation.
Features_ <- function(Hotgenes = NULL) {
  # Check object
  stopifnot(is(Hotgenes, "Hotgenes"))

  # Get rownames from all expression data
  Normalized_Data_(Hotgenes) %>%
    purrr::map(~ rownames(.x)) %>%
    unlist(use.names = FALSE) %>%
    unique() %>%
    sort()
}

#' @rdname Hotgenes-class
#' @export
# @importFrom BiocGenerics getObjectSlots updateObject
#' @param Hotgenes Hotgenes object.
#' @param min_n numeric, sets the minimum number of unique values
#' a column can have, default is 1. 
#' @param use_default_aliases logical, if TRUE, 
#' uses [set_missing_alias_default(]). If FALSE (Default), 
#' Mapper returened without modification.
#' @return Mapper_ Returns Mapper containing
#' aliases mapped to Features. Only aliases in 
#' expression data are returned
#'
Mapper_ <- function(Hotgenes = NULL, min_n = 1, 
                    use_default_aliases = FALSE) {
  # Check object
  stopifnot(is(Hotgenes, "Hotgenes"))

  
  MapperSlot <- data.frame(Feature = Features_(Hotgenes)) %>%
    tibble::as_tibble() %>%
    dplyr::left_join(Hotgenes@Mapper, by = "Feature") %>% 
  
  dplyr::select_if(~ dplyr::n_distinct(.x, na.rm = TRUE) > min_n) 

  if(use_default_aliases) {
    choices_aliases <- Mapper_aliases(Hotgenes)
    
    MapperSlot <- MapperSlot  %>% 
      set_missing_alias_default(
        default_alias = "Feature",
        preferred_alias = choices_aliases)
  }
    
   
  
  return(MapperSlot)
 
}

#' @rdname Hotgenes-class
#' @export
#' @return Mapper_aliases Returns a vector of Mapper column names other
#' than "Feature"
Mapper_aliases <- function(Hotgenes = NULL, min_n = 1) {
  # Check object
  stopifnot(is(Hotgenes, "Hotgenes"))
  
  Mapperfeatures <-  Mapper_(Hotgenes, min_n = min_n) %>% 
    dplyr::select(-dplyr::any_of(c("Feature"))) %>% 
    colnames()
  
  
  return(Mapperfeatures)
  
}



#' @rdname Hotgenes-class
#' @export
#' @param object Hotgenes object.
#' @param value a data.frame containing a "Feature"
#' for alias mapping.
"Mapper_<-" <- function(object, value) {
  object
}
setGeneric("Mapper_<-")

#' @rdname Hotgenes-class
#' @export
setMethod(
  "Mapper_<-", signature(object = "Hotgenes"),
  function(object, value) {
    
    
    object@Mapper <- value %>% 
      dplyr::mutate_all(as.character)

    validObject(object)

    object
  }
)





#' @rdname Hotgenes-class
#' @export
# @importFrom BiocGenerics getObjectSlots updateObject
#' @param Hotgenes Hotgenes object.
#' @param Col String matching the column id from Mapper
#' slot to pull. Default is "Feature".
#' @return Background_ Returns a vector of
#' all Features in the expression data.
Background_ <- function(Hotgenes = NULL,
                        Col = "Feature") {
  # Check object
  stopifnot(is(Hotgenes, "Hotgenes"))

  # Check for Mapper slot
  # Mapper_() automatically filters for expression
  Mapped_df <- Hotgenes %>% Mapper_()
  #AllFeatures <- Hotgenes %>% Features_()

  Mapped_df %>%
    dplyr::pull(Col) %>%
    unique()
}



# missing feature names check ---------------------------------------------


#' checks a vector for empty strings, NA, and "NA"
#' @export
#' @importFrom stringr str_replace_na str_detect
#' @inheritParams stringr::str_detect
#' @details
#' The string is first passed to the
#' \code{\link[stringr]{str_replace_na}} function to
#' convert NA into "".
#'
#' Then \code{\link[stringr]{str_detect}} is used to
#' detect spaces and "NA" with: "^[[:space:]]*$|^NA$"
#' @examples
#' this_string <- c("a", NA, "", "b", "NA", "NANA", "na")
#'
#' na_check(this_string)
na_check <- function(string = NULL) {
  string %>%
    stringr::str_replace_na(replacement = "NA") %>%
    trimws() %>% 
    stringr::str_detect("^[[:space:]]*$|^NA$")
}



# missingness_checks ------------------------------------------------------


#' Vector-aware missingness detection
#' @name missingness_checks
#' @title Check for missing or empty values in vectors
#' @description
#' These functions provide vector-aware detection of missing or empty values.
#' Unlike \code{na_check()}, these functions treat the string "NA" as valid data.
#' 
NULL

#' [detected_missingness()] checks for non-whitespace values in specified vector
#' @description
#' Use [detected_missingness()] if you want logical values for each element in
#' a vector.
#' 
#' Use [any_missingness()] if you want a single logical value returned for 
#' a vector. This is the equivalent to ```any(detected_missingness(x))```.
#' 
#' @param x vector to check
#' @return [detected_missingness()] returns vector of logical values
#' indicating the presence of empty values TRUE 
#' or non-empty values FALSE.
#' @export
#' @rdname missingness_checks
#' @md
detected_missingness <- function(x = NULL) {
  
  # Handle NULL input - return single TRUE
  if (is.null(x)) {
    return(TRUE)
  }
  
  # Handle empty vector - return single TRUE
  if (length(x) == 0) {
    return(TRUE)
  }
  
  # grepl returns FALSE for NA values (they don't match the pattern)
  # !grepl returns TRUE for empty strings, whitespace-only, and NA
  vector_out <- !grepl(x = x, pattern = "\\S") 
  
  return(vector_out)
}


#' [any_missingness()] checks if any element has missingness
#' 
#' @param x vector to check
#' @return [any_missingness()] returns a single logical value indicating 
#' the presence of any empty values TRUE, otherwise FALSE
#' @export
#' @rdname missingness_checks
#' @md
any_missingness <- function(x = NULL) {
  
  vector_out <- any(detected_missingness(x = x))
  
  return(vector_out)
}


#' [all_truthy()] checks if all elements in a vector are truthy
#' @description
#' Vector-aware equivalent to checking if all elements would pass 
#' shiny::isTruthy(). Returns FALSE if any element is NA, NULL, "", 
#' or whitespace-only.
#' 
#' @param x vector to check
#' @return Single logical value: TRUE if all elements are truthy, FALSE otherwise
#' @export
#' @rdname missingness_checks
#' @md
all_truthy <- function(x = NULL) {
  
  # NULL or empty vector is not truthy
  if (is.null(x) || length(x) == 0) {
    return(FALSE)
  }
  
  # Check if any elements are missing/empty
  !any_missingness(x)
}


#' [any_truthy()] checks if any element in a vector is truthy
#' @description
#' Vector-aware check for at least one truthy element.
#' 
#' @param x vector to check
#' @return Single logical value: TRUE if any element is truthy, FALSE otherwise
#' @export
#' @rdname missingness_checks
#' @md
any_truthy <- function(x = NULL) {
  
  # NULL or empty vector has no truthy values
  if (is.null(x) || length(x) == 0) {
    return(FALSE)
  }
  
  # Inverse of all elements being missing
  !all(detected_missingness(x))
}