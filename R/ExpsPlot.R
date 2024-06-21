#' Makes expression plots via ggplot2
#' @inheritParams DExps
#' @inheritParams ggplot2::facet_wrap
#' @importFrom ggplot2 coord_flip element_blank element_text
#' geom_boxplot geom_hline
#' geom_line geom_segment ggplot label_both labs scale_fill_manual theme
#' theme_bw theme_minimal unit
#' @param group String indicating the grouping variable.
#' @param xVar String indicating the x variable. If NULL (default),
#' the first column name reported by coldata_ function will be used.
#' @param yVar String indicating the y variable
#' @param basemean Logical, if TRUE a horizontal line
#' indicating the mean value for all samples will be shown. If FALSE (default)
#' no line will be added
#' @param fill String variable to set fill, default is NULL.
#' @param color String variable to set color, default is NULL.
#' @param linevar String variable to set line, default is NULL.
#' @param boxplot logical default is FALSE. If TRUE a boxplot will
#' be generated
#' @param pointplot logical if TRUE (default), point will be generated.
#' @param name_col string for naming column with feature names,
#' default is "Feature". Options include columns in mapper slot.
#' @param value_col string for naming column with expression
#' values, default is "value".
#' @param filter_eval to be passed to
#'  \code{\link[dplyr]{filter}}.
#' @inheritParams UpdateLevelsbyList
#' @importFrom rlang enexprs !!!
#' @param ... additional arguments for
#' \code{\link[ggplot2]{facet_wrap}}
#' @inheritParams DEphe
#' @export
#' @return ggplot object
#' @details xVar and group must be a column name accessible by coldata_
#' function. String for yVar must be accessible by parse_features() function.
#' @example man/examples/ExpsPlot_Example.R
#'
ExpsPlot <- function(
    Hotgenes = NULL,
    ExpressionSlots = NULL,
    SampleIDs = NULL,
    group = NULL,
    fill = NULL,
    color = NULL,
    xVar = NULL,
    yVar = NULL,
    linevar = NULL,
    boxplot = FALSE,
    pointplot = TRUE,
    basemean = FALSE,
    name_col = "Feature",
    value_col = "value",
    scales = "fixed",
    facets = NULL,
    filter_eval = NULL,
    named_levels = NULL,
    ...) {
  # Check object
  stopifnot(is(Hotgenes, "Hotgenes"))

  # color pal

  if (!shiny::isTruthy(fill)) {
    fill_pal <- NULL
  } else {
    fill_pal <- Hotgenes %>%
      coldata_palettes(coldata_ids = fill, SampleIDs = SampleIDs) %>%
      purrr::chuck(fill)
  }

  if (!shiny::isTruthy(color)) {
    color_pal <- NULL
  } else {
    color_pal <- Hotgenes %>%
      coldata_palettes(coldata_ids = color, SampleIDs = SampleIDs) %>%
      purrr::chuck(color)
  }



  # Prep to building long expression data
  Matched_ExpSel <- match.arg(
    ExpressionSlots,
    ExpressionSlots_(Hotgenes)
  )

 # parsing variables
 
  parsed_features <- parse_features(Hotgenes,
                                features = c(xVar, yVar, "SampleIDs",
                                             color, linevar,
                                             group, fill, facets))
 
  
  yVar_parsed <- parsed_features[c("aux_features",
                                   "exps_features")] %>% 
    purrr::imap(function(x,y){
      x %>% stringr::str_subset( ".+")
      
    }) %>% purrr::compact() %>% 
    unlist(use.names = FALSE)
  
  if(!all(yVar_parsed %in% parsed_features$exps_features)){
    Update_Y_lab <- "Relative values\n(mixed assays)"
  } else {
    Update_Y_lab <- Matched_ExpSel 
  }
  
  # yVar_parsed can be multiple Features, but not NULL
  if (is.null(yVar_parsed) | is.null(xVar)) {
    ExpDat <- NULL
  } else if (!is.null(yVar_parsed) && !is.null(xVar)) {
    # expression data
    ExpDat <- Hotgenes %>%
      DExps(
        Query_set = TRUE,
        hotList = parsed_features$exps_features,
        aux_features = parsed_features$aux_features,
        coldata_ids = parsed_features$coldata_features,
        ExpressionSlots = Matched_ExpSel,
        SampleIDs = SampleIDs
      ) %>%
      dplyr::arrange(dplyr::across(dplyr::any_of(xVar)))# %>%
      
    
    if(name_col != "Feature"){
      # label_by
      label_aliases <- hotList_mapper(Hotgenes,
                                      hotList = parsed_features$exps_features)%>%
        dplyr::select(dplyr::any_of(c("Feature", name_col))) 
      
      ExpDat <- ExpDat %>% 
        tidyr::pivot_longer(
          cols = dplyr::any_of(c(yVar_parsed)),
          names_to = "Feature",
          values_to = value_col
        ) %>% 
        dplyr::left_join(label_aliases, by = "Feature")
    } else {
      
      
      ExpDat <- ExpDat %>% 
      # converting to long format
      tidyr::pivot_longer(
        cols = dplyr::any_of(c(yVar_parsed)),
        names_to = name_col,
        values_to = value_col
      ) 
    }
    
    
    
  }

  # setting blank
  if (is.null(ExpDat) || nrow(ExpDat) == 0) {
    gplot <- empty_plot()
   
  } else {
    # plot --------------------------------------------------------------------

    # setting mapping
    mapping <- ggplot2::aes(
      x = .data[[xVar]],
      y = .data[[value_col]],
      fill = .data[[fill]],
      color = .data[[color]],
      group = .data[[group]]
    )



    # Preparing to drop empty mappings
    list_mappings <- list(
      x = xVar,
      y = yVar_parsed,
      fill = fill,
      colour = color,
      group = group
    ) %>%
      purrr::keep(~ all(!shiny::isTruthy(.x)))

    # removing empties
    mapping[names(list_mappings)] <- NULL

    # updating expression data
    # can't use is.null
    if (!missing(filter_eval)) {
      args <- rlang::enquos(filter_eval)

      # enexprs this works
      # args <- rlang::enexprs(filter_eval)
      ExpDat <- ExpDat %>%
        dplyr::filter(!!!args)
    }

    # updating levels
    if (!is.null(named_levels)) {
      ExpDat <- ExpDat %>%
        UpdateLevelsbyList(named_levels = named_levels)
    }

    # making plot
    gplot <- ExpDat %>%
      ggplot2::ggplot(mapping) +
      list(
        # geom_point(size=2),
        # scale_fill_manual(values = alpha(fill_pal,0.3)) ,
        # scale_color_manual(values = alpha(color_pal,1)),
        ggplot2::scale_fill_manual(values = fill_pal),
        ggplot2::scale_color_manual(values = color_pal),
        ggplot2::labs(y = Update_Y_lab),
        ggplot2::theme_classic(),
        ggplot2::theme(
          aspect.ratio = 1,
          axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
          title = ggplot2::element_text(hjust = 0.5),
          plot.margin = ggplot2::unit(c(1, 1, 1, 1), "cm")
        ),
        ggplot2::facet_wrap(c(name_col, facets), # labeller = labeller ,
          scales = scales, ...
        )
      )

    # setting baseMean
    if (isTRUE(basemean)) {
      # by Feature
      baseMean <- ExpDat %>%
        dplyr::group_by(.data[[name_col]]) %>%
        dplyr::summarise(baseMean = mean(.data[[value_col]]))

      gplot <- gplot +
        ggplot2::geom_hline(data = baseMean, aes(yintercept = .data$baseMean))
    }
    # line check
    if (shiny::isTruthy(linevar)) {
      gplot <- gplot +
        ggplot2::geom_line(ggplot2::aes( # color = .data[[linevar]],
          group = .data[[linevar]]
        ))
    }

    # point check
    if (isTRUE(pointplot)) {
      gplot <- gplot + ggplot2::geom_point()
    }

    # boxplot check
    if (isTRUE(boxplot)) {
      gplot <- gplot + ggplot2::geom_boxplot()
    }
  }

  # final return

  return(gplot)
}

#' Converts ExpsPlot data to matrix format
#' @inheritParams ExpsPlot
#' @importFrom tidyr pivot_longer pivot_wider
#' @param ggplot_obj ggplot object
#' @export
ExpsPlot_Table <- function(
    ggplot_obj = NULL,
    xVar = NULL,
    facets = NULL,
    name_col = "Feature",
    value_col = "value") {
  data <- ggplot_obj$data
  # preparing to reformat data
  xx_vars <- c(facets, xVar)

  D2 <- data %>%
    dplyr::select(dplyr::any_of(c(
      xx_vars,
      name_col,
      value_col
    )))


  # adding column name as prefix
  D2[xx_vars] <- D2[xx_vars] %>%
    purrr::imap(function(x, y) {
      paste(y, x, sep = "")
    })

  # pivot table
  D_table_Out <- D2 %>%
    dplyr::group_by(dplyr::across(dplyr::any_of(xx_vars))) %>%
    dplyr::mutate(Uniq_id = row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::across(dplyr::any_of(xx_vars))) %>%
    tidyr::pivot_wider(
      id_cols = c(name_col, "Uniq_id"),
      names_from = dplyr::any_of(xx_vars),
      values_from = dplyr::any_of(value_col)
    ) %>%
    dplyr::select(-dplyr::any_of("Uniq_id")) %>%
    dplyr::rename(yVar = .data[[name_col]]) %>%
    dplyr::relocate(dplyr::any_of(c("yVar")), .before = 1) %>%
    dplyr::arrange(.data$yVar)

  return(D_table_Out)
}

