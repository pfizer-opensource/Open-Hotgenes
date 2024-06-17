#' Generates VennDiagrams from lists
#' @export
#' @importFrom ggvenn ggvenn
#' @param FeatureList Named list with up to five vectors
#' @param ... Additional paramters for ggvenn \code{\link[ggvenn]{ggvenn}}
#' @inheritParams ggvenn::ggvenn
#' @inheritParams ggplot2::theme
#' @return List containing: venn diagram, corresponding lists,
#' and names of lists with overlapping features.
#' @details Up to five vectors are compared.
#' @seealso gplots::venn
#' @seealso ggvenn::ggvenn
#' @example man/examples/Venn_Report_Example.R

Venn_Report <- function(FeatureList = NULL,
                        show_percentage = FALSE,
                        fill_alpha = 0.1,
                        set_name_size = 6,
                        text_size = 4,
                        ...) {
  # check lengths
  if(length(FeatureList) > 4){
    stop("no more that 4 groups")
  }
  
  
  if (isFALSE(all(lengths(FeatureList) == 0))) {
    Venn_FeatureList <- FeatureList %>%
      ggvenn::ggvenn(
        show_percentage = show_percentage,
        fill_alpha = fill_alpha,
        set_name_size = set_name_size,
        text_size = text_size,
        ...
      ) +
    #  ggplot2::theme(plot.margin = ggplot2::unit(c(1,1,1,1), "cm"),
     #                aspect.ratio = aspect.ratio) +
      ggplot2::coord_fixed(clip = "off" ) +
      ggplot2::theme(plot.margin = ggplot2::unit(c(top = 1,
                                                   right = 2,
                                                   bottom = 1,
                                                   left = 2), "in"))

    Int_FeatureList <- FeatureList %>%
      find_intersections()
    
    Int_Names <-
      names(Int_FeatureList)[grepl(":", names(Int_FeatureList))]
  } else {
    Venn_FeatureList <- empty_plot()
    
    Int_FeatureList <- FeatureList
    
    Int_Names <- NULL
  }
  
  ReportList <- list(vennD = Venn_FeatureList,
                     Intsect = Int_FeatureList,
                     Names = Int_Names)
  return(ReportList)
}



#' gets overlapping values
#' @noRd
find_intersections <- function(named_list = NULL) {
  obj_sorted <- named_list %>%
    purrr::imap( ~ data.frame(value = .x)) %>% 
    purrr::list_rbind(names_to = "contrast" )  %>%
    
    dplyr::summarise(col = stringr::str_c(.data$contrast, collapse = ":"),
                     .by = "value")
  
  obj_sorted <- split(obj_sorted$value, obj_sorted$col)
  
  ne <- names(obj_sorted)
  
  ne <- ne[order(stringr::str_count(ne, pattern = ":"))]
  return(obj_sorted[ne])
}
