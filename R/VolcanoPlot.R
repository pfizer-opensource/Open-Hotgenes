#' Internal VolcanoPlot
#' @noRd
#' @importFrom ggrepel geom_text_repel
#' @inheritParams ggrepel::geom_text_repel
#' 
ggplot2_Volcano_prep <- function(data = NULL,
                       hotList = NULL,
                       FeatureCol = "Feature",
                       padj_cut = 0.1,
                       mapper_df = NULL,
                       .log2FoldChange = 0,
                       pointSize = 2,
                       max.overlaps = 50,
                       legend_pointSize = 5,
                       point_label_size = 10,
                       base_size = 15,
                       Hide_labels = FALSE,
                       repel_labels = TRUE,
                       col = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C"),
                       colAlpha = 1,
                       matched_contrast = NULL) {
  # set up legend
  output_levels <- c(
    "NS", "log2FC",
    "padj", "log2FC_padj"
  )

  col_set <- col %>%
    purrr::set_names(output_levels)

  # specify labeling
  if (isFALSE(Hide_labels)) {
    if (!is.null(hotList)) {
      label_these_features <- hotList
    } else if (is.null(hotList)) {
      label_these_features <- data %>%
        dplyr::pull(FeatureCol) %>%
        unique()
    }
  } else if (isTRUE(Hide_labels)) {
    label_these_features <- NA_character_
  }


  # annotating data
  Updated_data <- data %>%
    dplyr::mutate(FDR = -log10(.data$padj)) %>%
    dplyr::mutate(
      Expression = dplyr::case_when(

        # output_levels[4] "log2FC & padj"
        abs(.data$log2FoldChange) >= .log2FoldChange &
          .data$padj < padj_cut ~ output_levels[4],

        #  output_levels[3] "padj"
        abs(.data$log2FoldChange) < .log2FoldChange &
          .data$padj < padj_cut ~ output_levels[3],

        # output_levels[2] "log2FC"
        abs(.data$log2FoldChange) >= .log2FoldChange &
          .data$padj >= padj_cut ~ output_levels[2],

        # output_levels[1] NS
        TRUE ~ output_levels[1]
      )
    ) %>%
    dplyr::mutate(
      Expression = factor(.data$Expression, output_levels)
    )

  # adding showGenes col
  Updated_data <- Updated_data %>%
    dplyr::mutate(showGenes = dplyr::case_when(
      .data$Expression == output_levels[4] &
        .data$Feature %in% label_these_features ~ .data$Feature
    ))



  total_above_threshold <- Updated_data %>%
    dplyr::count(.data$Expression, .drop = FALSE) %>%
    unique() %>%
    tibble::deframe()
 

  This_q <- total_above_threshold %>%
    
    purrr::imap(function(x, y) {
      
      glue::glue("bquote(~{y}^{x})") %>%
        str2lang()
    })
  
  # mapper join
  p2 <- Updated_data %>% 
    dplyr::left_join(mapper_df, by = "Feature") %>% 
  
   ggplot2::ggplot(
  #  Updated_data,
    ggplot2::aes(.data$log2FoldChange, .data$FDR)
  ) +
    ggplot2::geom_point(ggplot2::aes(color = .data$Expression),
      alpha = colAlpha,
      size = pointSize
    ) +
    ggplot2::xlab(expression("log"[2] * "FC")) +
    ggplot2::ylab(expression("-log"[10] * "padj")) +
    ggplot2::scale_color_manual(
      values = col_set, drop = FALSE,
      labels = purrr::map(This_q, ~ eval(.x))
    ) +
    ggplot2::labs(
      title = glue::glue("{matched_contrast}")
     
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      title = glue::glue(""),
      override.aes = list(size = legend_pointSize)
    )) +

    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.position = "top",
      title = ggplot2::element_text(hjust = 0.5),
      plot.margin = ggplot2::unit(c(1, 1, 1, 1), "cm")
    ) +
    ggplot2::geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
    ggplot2::geom_vline(xintercept = -.log2FoldChange, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = .log2FoldChange, linetype = "dashed")

  if(isFALSE(Hide_labels)){
    
  
  if (isTRUE(repel_labels)) {
    p2 <- p2 +
      ggrepel::geom_text_repel(
        mapping = aes(label = .data$showGenes),
        size = point_label_size,
        max.overlaps = max.overlaps,
        na.rm = TRUE
      )
  } else {
    p2 <- p2 +
      ggplot2::geom_text(
        mapping = aes(label = .data$showGenes),
        size = point_label_size,
        na.rm = TRUE
      )
  }

  }

  # p2
  return(p2)
}

#' Hotgenes volcano plot
#' @export
#' @inheritParams DE
#' @importFrom ggrepel geom_text_repel
#' @importFrom plotly layout ggplotly
#' @inheritParams ggrepel::geom_text_repel
#' @param Hide_labels Logical, if TRUE (Default), labels will not be shown.
#' @param ... additional internal parameters.
#' @param col vector of colors to use for plot. Default is
#' c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C") and
#' maps with c("NS", "log2FC", "padj", "log2FC & padj").
#' @param repel_labels logical, if TRUE (default), labels will
#' be generated using ggrepel::geom_text_repel. If FALSE,
#' ggplot2::geom_text will be used.
#' @param interactive logical, if TRUE, plotly object returned.
#' If FALSE (default), ggplot object will be returned.
#' @param shinymode logical, if TRUE a named list of 
#' parameters will be exported to processing by internal plotting functions.
#' If FALSE, a volcano plot is returned. 
#' @example man/examples/VPlot_Example.R

VPlot <- function(
    Hotgenes = NULL,
    contrasts = NULL,
    hotList = NULL,
    padj_cut = 0.1,
    .log2FoldChange = 0,
    Hide_labels = TRUE,
    max.overlaps = 50,
    repel_labels = TRUE,
    interactive = FALSE,
    col = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C"),
    shinymode = FALSE,
    ...) {
  # checks if is hotgenes
  stopifnot(is(Hotgenes, "Hotgenes"))

  # check contrasts


  # only one allowed for volcano plots
  matched_contrast <- match.arg(contrasts, contrasts_(Hotgenes),
    several.ok = FALSE
  )
  
  
  # DE Frame
  DE_Frame <- Output_DE_df(Hotgenes,
    padj_cut = 1,
    .log2FoldChange = 0,
    contrasts = matched_contrast)

  mapper_df <- Mapper_(Hotgenes)
  
 
  
  
  ###
  if(isTRUE(shinymode)){
    
    # sets up mapper
    
    outVplot <- list(
      data = DE_Frame,
      hotList = hotList,
      FeatureCol = "Feature",
      mapper_df = mapper_df,
      .log2FoldChange = .log2FoldChange,
      padj_cut = padj_cut,
      Hide_labels = Hide_labels,
      col = col,
      repel_labels = repel_labels,
      max.overlaps = max.overlaps,
      matched_contrast = matched_contrast, ...
    ) %>% purrr::compact()
    
    return(outVplot)
  } else {
    
  

  if (isFALSE(interactive)) {
    outVplot <- ggplot2_Volcano_prep(
      data = DE_Frame,
      hotList = hotList,
      FeatureCol = "Feature",
      mapper_df = mapper_df,
      .log2FoldChange = .log2FoldChange,
      padj_cut = padj_cut,
      Hide_labels = Hide_labels,
      col = col,
      repel_labels = repel_labels,
      max.overlaps = max.overlaps,
      matched_contrast = matched_contrast, ...
    )
  } else {
    outVplot <- plotly_Volcano_prep(
      data = DE_Frame,
      hotList = hotList,
      mapper_df = mapper_df,
      FeatureCol = "Feature",
      .log2FoldChange = .log2FoldChange,
      padj_cut = padj_cut,
      Hide_labels = Hide_labels,
      col = col,
      repel_labels = repel_labels,
      max.overlaps = max.overlaps,
      matched_contrast = matched_contrast, ...
    ) 
  }



  return(outVplot)
  }
}




#' Internal VolcanoPlot
#' @noRd
plotly_Volcano_prep <- function(data = NULL,
                                hotList = NULL,
                                FeatureCol = "Feature",
                                mapper_df = NULL,
                                padj_cut = 0.1,
                                .log2FoldChange = 0,
                                pointSize = 2,
                                max.overlaps = 50,
                                legend_pointSize = 5,
                                point_label_size = 10,
                                base_size = 15,
                                Hide_labels = FALSE,
                                repel_labels = TRUE,
                                col = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C"),
                                colAlpha = 1,
                                matched_contrast = NULL) {
  # set up legend
  output_levels <- c(
    "NS", "log2FC",
    "padj", "log2FC_padj"
  )
  
  # specify labeling
  if (isFALSE(Hide_labels)) {
    if (!is.null(hotList)) {
      label_these_features <- hotList
    } else if (is.null(hotList)) {
      label_these_features <- data %>%
        dplyr::pull(FeatureCol) %>%
        unique()
    }
  } else if (isTRUE(Hide_labels)) {
    label_these_features <- NA_character_
  }


  # annotating data
Updated_data <- data %>% 
  
    dplyr::mutate(FDR = -log10(.data$padj)) %>%
    dplyr::mutate(
      Expression = dplyr::case_when(

        # output_levels[4] "log2FC & padj"
        abs(.data$log2FoldChange) >= .log2FoldChange &
          .data$padj < padj_cut ~ output_levels[4],

        #  output_levels[3] "padj"
        abs(.data$log2FoldChange) < .log2FoldChange &
          .data$padj < padj_cut ~ output_levels[3],

        # output_levels[2] "log2FC"
        abs(.data$log2FoldChange) >= .log2FoldChange &
          .data$padj >= padj_cut ~ output_levels[2],

        # output_levels[1] NS
        TRUE ~ output_levels[1]
      )
    ) %>%
    dplyr::mutate(
      Expression = factor(.data$Expression, levels = output_levels)
    ) 


# checks for missing output_levels and then add empty row to tibble
# otherwise plotly will drop color 
Updated_data <- output_levels %>% 
  purrr::set_names(~.x) %>% 
  purrr::imap(function(x,y){
    
    Tbale <-  Updated_data %>% 
      dplyr::filter(.data$Expression == y)
    
    count_ <- nrow(Tbale)
    
    if(nrow(Tbale) == 0){
      Tbale <- Tbale%>% 
        tibble::add_row()%>% 
        dplyr::mutate(Expression = y)
    }
    Count_ <- glue::glue("{y}<sup>{count_}</sup>") 
    
    Tbale <- Tbale %>% 
      dplyr::mutate(count = Count_,
      g_name = y, .before = 1)

    return( Tbale)
  }) %>% 
  purrr::list_rbind()%>%
  dplyr::mutate(
    Expression = factor(.data$Expression, levels = output_levels)
  ) 

  
This_q <- Updated_data %>%
    dplyr::select(c("count", "g_name")) %>% 
    unique() %>%
    tibble::deframe()
  
# set color scheme
  col_set <- col %>%
    purrr::set_names(names(This_q))
  
# updates names
    Updated_data <- Updated_data %>%
      dplyr::mutate(showGenes = dplyr::case_when(
        .data$Expression == output_levels[4] &
          .data$Feature %in% label_these_features ~ .data$Feature 
        
      )) %>%
      
     dplyr::mutate(
       Expression = forcats::fct_recode(.data$Expression, !!!This_q)
     ) 
      

 if("PointIDs" %in% names(Updated_data)){
    
   Updated_data <- Updated_data %>% 
     dplyr::mutate(!!FeatureCol := .data$PointIDs)
   
     
 }
   
    
   # Updated_data$Expression %>% levels() %>% print()
new_col_id <- "{matched_contrast}: "
    
  p2_test <- Updated_data %>% 
     ggplot2::ggplot(
   # Updated_data,
    ggplot2::aes(x = .data$log2FoldChange, y = .data$FDR)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
         color = .data[[glue::glue("Expression")]],
        group = .data[[FeatureCol]]  
      ),
      alpha = colAlpha,
      size = pointSize
    ) +
    ggplot2::ylab("-log<sub>10</sub>padj") +
    ggplot2::xlab("-log<sub>2</sub>FC") +
    
    ggplot2::scale_color_manual(values = col_set, drop = FALSE ) +
                                
    
ggplot2::ggtitle(label = glue::glue("{matched_contrast}")) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
    ggplot2::geom_vline(xintercept = -.log2FoldChange, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = .log2FoldChange, linetype = "dashed") +
    ggplot2::theme(plot.margin = ggplot2::unit(c(3,1,1,1), "cm"))

  if (isFALSE(Hide_labels)) {
    p2_test <- p2_test +
      ggplot2::geom_text(
        mapping = ggplot2::aes(label = .data$showGenes),
        check_overlap = TRUE,
        size = point_label_size,
        na.rm = TRUE
      )
  }


  p2_test <- p2_test %>%
    plotly::ggplotly() %>%
    plotly::layout(legend = list(
      title = "",
      
      xanchor = "center",
      yanchor = "top",
      x = 0.5, y = 1.1,
      orientation = "h"
      
    ))
  
  return(p2_test)
  
}
