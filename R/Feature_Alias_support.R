
#' Returns a named list showing where the features are found
#' @export
#' @inheritParams DE
#' @importFrom rlang is_empty
#' @param features strings to check
#' @rdname utilities

parse_features <- function(Hotgenes = NULL,
                           features = NULL){
  
  # auxiliary_assays_features
  
  aux_features <-  features[features %in% auxiliary_assays_features(Hotgenes) ]
  
  #exps_features <-  features[features %in% Features_(Hotgenes) ]
  
  exps_features <- hotList_Feature(Hotgenes = Hotgenes,
                  hotList = features)
  
  coldata_features <-  features[features %in% coldata_names(Hotgenes) ]
  
  
  output <- list(aux_features = aux_features,
                 exps_features = exps_features,
                 coldata_features = coldata_features) %>% 
    purrr::imap(function(x,y){
      if(rlang::is_empty(x)){
        return("")
      } else {
        return(unique(x))
      }
    })
  return(output)
  
}


#' Returns a named list showing where the aliases are found
#' @export
#' @inheritParams DE
#' @param features strings to check

parse_mapper_aliases <- function(Hotgenes = NULL,
                                 features = NULL){
  
  

 initial_list <- Mapper_(Hotgenes) %>% 
    purrr::imap(function(x,y){
      
      features[features %in% x]
    })
  

  output <- initial_list %>% 
    purrr::imap(function(x,y){
      if(rlang::is_empty(x)){
        return("")
      } else {
        return(x)
      }
    })
  return(output)
  
}

#' filters Mapper across all columns with hotList
#' @export
#' @inheritParams DE
#' @importFrom shiny isTruthy
#' @importFrom dplyr filter any_of if_any
#' @example man/examples/Parse_Aliases_Example.R
hotList_mapper <- function(Hotgenes = NULL,
                          hotList = NULL){
  mapper_df <- Mapper_(Hotgenes)
  if(all_truthy(hotList)){
    
    mapper_col<- names(mapper_df)
    
    out<- mapper_df %>% 
      dplyr::filter(
        dplyr::if_any(
          .cols = dplyr::any_of(mapper_col),
          .fns = ~ .x %in% hotList
        )
      ) 
  } else {
    out <- mapper_df
  }
  
  
  
  return(out)
  
  
  
}



#' returns mapped Features for any hotList
#' @export
#' @rdname hotList_mapper
hotList_Feature <- function(Hotgenes = NULL,
    hotList = NULL){
  
  if(all_truthy(hotList)){
    out <- hotList_mapper(Hotgenes = Hotgenes,
                          hotList = hotList)$Feature
    
  } else {
    out <- NULL
  }
  
  
  
  return(out)
  
  
  
}

