#' Heatmap generator
#' @export
#' @importFrom pheatmap pheatmap
#' @inheritParams DExps
#' @inheritParams coldata_
#' @inheritParams pheatmap::pheatmap
#' @inheritParams UpdateLevelsbyList
#' @param ... Additional parameters for
#' \code{\link[pheatmap]{pheatmap}}.
#' @param arrangeby String mapping to a column name returned by
#' coldata_. Samples in the heatmap will be sorted by this variable.
#' @param annotations, string indicating the coldata_ids to select.
#' If NULL (default), none are returned (unlike all with coldata_).
#' @param label_by String matching a column name returned by
#' Mapper_ function, default is "Feature".
#' @example man/examples/DEphe_Example.R
#' @details If Feature maps to multiple names, only one will be
#' used.

DEphe <- function(
    Hotgenes = NULL,
    ExpressionSlots = NULL,
    SampleIDs = NULL,
    hotList = NULL,
    contrasts = NULL,
    padj_cut = 0.1,
    keep_na_padj = FALSE,
    .log2FoldChange = 0,
    Topn = Inf,
    annotations = NULL,
    mode = "all",
    border_color = "white",
    scale = "row",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    arrangeby = NULL,
    label_by = "Feature",
    named_levels = NULL,
    ...) {
  # Check object
  stopifnot(is(Hotgenes, "Hotgenes"))

  # heatmap
  DM <- DExps(
    Hotgenes = Hotgenes,
    hotList = hotList,
    SampleIDs = SampleIDs,
    Topn = Topn,
    ExpressionSlots = ExpressionSlots,
    contrasts = contrasts,
    padj_cut = padj_cut,
    keep_na_padj = keep_na_padj,
    .log2FoldChange = .log2FoldChange
  )

  
  
 
  
  Top_Features <- rank_hotList(Hotgenes = Hotgenes, 
                               contrasts = contrasts, 
                              hotList = rownames(DM), 
                              Topn = Topn) %>% 
    dplyr::pull("Feature")
  
  DM <- DM[Top_Features, ]

  if (nrow(DM) == 0 | ncol(DM) == 0) {
    
    empty_plot()
    
  } else {
    # check annotations
    # needs to return null, otherwise coldata_ returns all by default
    if (is.null(annotations)) {
      annot_input <- NULL
    } else if (!is.null(annotations)) {
      annot_input <- Hotgenes %>%
        coldata_(coldata_ids = annotations, mode = mode) %>%
        dplyr::arrange_at(arrangeby)
      
      # updating levels
      if (!is.null(named_levels)) {
        annot_input <- annot_input %>%
          UpdateLevelsbyList(named_levels = named_levels ) %>%
          dplyr::arrange_at(arrangeby)
      }
    }

    # checking if DM needs to be sorted
    if (is.null(annot_input)) {
      # print("no change to DM")
    } else if (!is.null(annot_input)) {
      DM_s <- rownames(annot_input)[rownames(annot_input) %in% colnames(DM)]
      DM <- DM[, DM_s]
    }

    # checking Features
    mappingDF <- Mapper_(Hotgenes, use_default_aliases = TRUE)
    Feature_Ids <- colnames(mappingDF)

    Matched_FeatureID <- match.arg(label_by,
      choices = Feature_Ids
    )

    if (Matched_FeatureID != "Feature") {
      de_annot <- mappingDF %>%
        dplyr::select(dplyr::any_of(c("Feature", Matched_FeatureID))) %>%
        tibble::deframe()

      labels_row <- de_annot[rownames(DM)]
    } else if (Matched_FeatureID == "Feature") {
      labels_row <- NULL
    }


   
    
    # pheatmap
    pheat_p <- pheatmap::pheatmap(DM,
      annotation_col = annot_input,
      scale = scale,
      labels_row = labels_row,
      border_color = border_color,
      cluster_rows = cluster_rows,
      cluster_cols = cluster_cols,
      # main = paste0("Expression: ", ExpressionSlots),
      ...
    )


    pheat_p$topHits <- rownames(DM)
    return(pheat_p)
  }
}


#' Converts a coldata data.frame to a color palette.
#' @export
#' @param Input_Object A data.frame containing factors or a hotgenes object.
#' @param brewer.pals A vector of strings indicating the desired RColorBrewer
#' palettes to use. See \link[RColorBrewer]{brewer.pal}.
#' Default is set to "Dark2", "Set2". Multiple maybe used.
#' @param BinaryColorScheme A vector of strings to for binary levels.
#' The default is set to "lightgrey","black". If NULL, two colors will
#' be assigned.
#' @inheritParams coldata_
#' @inheritParams DExps
#' @inheritParams UpdateLevelsbyList
#' @example man/examples/coldata_palettes_Example.R
#'


coldata_palettes <- function(
    Input_Object = NULL,
    brewer.pals = c("Dark2", "Set2"),
    BinaryColorScheme = c("lightgrey", "black"),
    coldata_ids = NULL,
    mode = "all",
    named_levels = NULL,
    SampleIDs = NULL) {
  if (is(Input_Object, "Hotgenes")) {
    # Get SampleIDs
    if (is.null(SampleIDs)) {
      Out_SampleIDs <- SampleIDs_(Input_Object)
    } else if (!is.null(SampleIDs)) {
      Out_SampleIDs <- SampleIDs
    }
    
    
    #if(!)
    coldata_df <- Input_Object %>%
      coldata_( # coldata_ids=coldata_ids,
        mode = mode
      )
    
    
    coldata_df <- coldata_df %>%
      tibble::rownames_to_column("rs") %>%
      dplyr::filter(.data$rs %in% Out_SampleIDs) %>%
      droplevels() %>%
      tibble::column_to_rownames("rs") %>%
      valid_factors()
  } else if (is(Input_Object, "data.frame")) {
    coldata_df <- Input_Object %>%
      tibble::as_tibble() %>%
      valid_factors()
  }
  
  
  
  coldata_df %>%
    dplyr::select_all(factor) %>%
    UpdateLevelsbyList(named_levels = named_levels) %>% 
    purrr::imap(function(x, y) {
      # getting levels
      Factor_levels <- x %>%
        levels()
      
      # Getting colors
      Start_Colors <- RColorBrewer::brewer.pal.info %>%
        tibble::rownames_to_column("Name") %>%
        dplyr::filter(.data$Name %in% c(brewer.pals)) %>%
        plyr::dlply("Name", identity) %>%
        purrr::imap(function(xPal, yPal) {
          RColorBrewer::brewer.pal(xPal$maxcolors, yPal)
        }) %>%
        unlist(use.names = FALSE)
      
      
      # making sure that we have enough colors
      if (length(Start_Colors) < length(Factor_levels)) {
        All_Colors <- colorRampPalette(
          Start_Colors
        )(length(Factor_levels))
      } else if (length(Start_Colors) >= length(Factor_levels)) {
        All_Colors <- Start_Colors
      }
      
      
      # checking if BinaryColorScheme is supplied
      if (is.null(BinaryColorScheme)) {
        # Assigning colors to levels
        GroupColors <- All_Colors[seq(1, length(Factor_levels))] %>%
          purrr::set_names(Factor_levels)
        
        GroupColors
      } else if (!is.null(BinaryColorScheme)) {
        # checking BinaryColorScheme
        if (length(BinaryColorScheme) > 2) {
          stop("BinaryColorScheme length > 2")
        }
        
        # Binary color
        if (length(Factor_levels) <= 2) {
          # verifies that BinaryColorScheme length
          # matches length(Factor_levels)
          GroupColors <- BinaryColorScheme[c(seq(1, length(Factor_levels)))] %>%
            purrr::set_names(Factor_levels)
          
          GroupColors
        } else if (length(Factor_levels) > 2) {
          GroupColors <- All_Colors[seq(1, length(Factor_levels))] %>%
            purrr::set_names(Factor_levels)
          
          GroupColors
        }
      }
    }) %>% # list end
    
    purrr::compact() # dropping NAs
}

# internal function
# sorts all features by abs(stat)
rank_hotList <- function(
    Hotgenes = NULL,
    contrasts = NULL,
    hotList = NULL, 
    Topn = Inf){
  
  # getting available contrasts, excluding Intercept
  available_contrasts <- contrasts_(Hotgenes)

  if (is.null(contrasts)) {
    all_contrasts <- available_contrasts
  } else if (!is.null(contrasts)) {
    all_contrasts <- match.arg(contrasts, available_contrasts, several.ok = TRUE)
  }
  
  Out <-Hotgenes@Output_DE[all_contrasts] %>%
    purrr::list_rbind() 
  
  if(!is.null(hotList)){
    Out <- Out %>%
      dplyr::filter(.data$Feature %in%  hotList)

  }
  Out <- Out %>% 
    dplyr::select(all_of(c("Feature", "stat"))) %>%
    # tibble::as_tibble() %>%
    # Getting absolutate ranks
    dplyr::mutate(stat = abs(.data$stat)) %>%
    #dplyr::group_by(.data$Feature) %>%
    dplyr::slice_max(.data$stat, by = "Feature", n = 1) %>%
    #dplyr::ungroup() %>% 
    dplyr::slice_max(n = Topn, order_by = .data$stat) %>%
    dplyr::arrange(dplyr::desc(.data$stat)) %>% 
    tibble::as_tibble()
  
  return(Out)
  
}