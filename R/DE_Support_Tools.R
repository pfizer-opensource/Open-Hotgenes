#' Tools to support DE analysis
#' @name DE_Support_Tools
NULL


#' @rdname DE_Support_Tools
#' @export
#' @param raw_levels string of factor levels in data to be sorted.
#' @param ordered_levels string of factor levels in the desired order
#' @md
#' @details [level_matcher()] returns overlapping strings between
#' raw_levels and ordered_levels, with the ordering of the later.
#' No partial string matching. This is useful when applying a generic
#' factor level order to different metadata that do not share
#' all the same levels.
level_matcher <- function(raw_levels = NULL,
                          ordered_levels = NULL) {
  # making as character
  raw_levels <- as.character(raw_levels)
  ordered_levels <- as.character(ordered_levels)

  set_order <- seq(1, length(ordered_levels)) %>%
    purrr::set_names(ordered_levels)


  matched_levels <- names(sort(set_order[raw_levels]))

  return(matched_levels)
}


#' @export
#' @rdname DE_Support_Tools
#' @importFrom purrr imap_dfc
#' @importFrom forcats fct_relevel
#' @param metadata A data.frame containing factors
#' @param named_levels a named list containing the
#' levels to set as reference via
#' \code{\link[forcats]{fct_relevel}}
#' @details
#' The names from named_levels map to the
#' columns in metadata to update. List slots
#' contain strings passed to the ... argument of
#' \code{\link[forcats]{fct_relevel}}
#'
#' This is all done with \code{\link[purrr]{imap_dfc}}
UpdateLevelsbyList <- function(
    metadata = NULL,
    named_levels = NULL) {

   # forcats::fct_relevel requires vars be character of factor
  metadata <- metadata %>%
    dplyr::mutate(dplyr::across(names(named_levels),as.character ))
  
  metadata[names(named_levels)] <- metadata[names(named_levels)] %>%
    purrr::imap_dfc(function(x, y) {
      named_levels[[y]] <- level_matcher(
        raw_levels = unique(x),
        ordered_levels = named_levels[[y]]
      )

      forcats::fct_relevel(x, named_levels[[y]])
    }) 

  return(metadata)
}


#' @export
#' @rdname DE_Support_Tools
#' @importFrom dplyr select_if n_distinct
#' @details drops single level factors with
#'  metadata %>%
#'  select_if(~n_distinct(na.omit(.x))>1)
#'
#' @return data.frame with factors
valid_factors <- function(metadata = NULL) {
  metadata %>%
    dplyr::select_if(~ dplyr::n_distinct(na.omit(.x)) > 1) %>%
    dplyr::select_if(is.factor)
}

#' Makes string model from variables
#' @export
#' @importFrom purrr imap_chr set_names
#' @importFrom glue glue
#' @importFrom limma makeContrasts
#' @param ... pairs of character vector for specifying contrasts. For
#' each pair, position 1 is query and position 2 is reference. 
#' @param modelMatrix output from model.matrix() function.
#' @returns
#' Value object returned by limma::makeContrasts(), named by contrasts
#' @example man/examples/Hotgeneslimma_Example.R
#' 
limma_paired_contrasts <- function(..., modelMatrix = NULL){
  
  
  contrast_input <- list(...)
  
  # generating names
  contrast_names <- contrast_input %>% 
    purrr::imap_chr(function(x,y){
      
      
      c_names <- glue::glue("{x[1]}_v_{x[2]}")
      
      return(c_names)
      
    })
  
  # generating values
  contrast_values <- contrast_input %>% 
    purrr::set_names(contrast_names) %>% 
    purrr::imap_chr(function(x,y){
      
      val <- glue::glue("{x[1]} - {x[2]}")
      return(val)
      
    }) 
  
  
  # name making contrasts
  contrast.matrix <- limma::makeContrasts(
    
    contrasts = contrast_values,
    levels=modelMatrix)
  
  # contrast.matrix
  
  # adding names
  colnames(contrast.matrix) <-  names(contrast_values)
  
  return(contrast.matrix) 
  
}


