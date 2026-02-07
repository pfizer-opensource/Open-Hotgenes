#' FactoWrapper for PCA/HCPC of Hotgenes objects
#' @importFrom FactoMineR summary.PCA HCPC PCA
#' @importFrom factoextra fviz_pca_biplot fviz_pca_ind
#' @importFrom purrr list_rbind
#' @importFrom dplyr group_split
#' @inheritParams DExps
#' @inheritParams DE
#' @inheritParams FactoMineR::HCPC
#' @param biplot Boolean value. If TRUE (Default), both individuals and variables
#' will be shown on the plot. If FALSE, only individuals will be plotted.
#' @param Top_var a numeric value indicating the number of top contributing
#' genes to show on the biplot.
#' @param ellipse.level numeric value for ellipse size, from 0 to 1.0
#' @param ellipse.alpha numeric value for ellipse transparency, from 0 to 1.0.
#' @param habillage_selection string indicating the grouping to be shown.
#' Default is "clust", which represents the HCPC ward clusters.
#' @param label_sel String indicating if individuals "ind" and variables
#' "var" should be labeled on the biplot.
#' @param pointsize numeric value for point size.
#' @param labelsize numeric value for label size.
#' @param Hotgenes Hotgenes R object
#' @export
#' @return FactoMiner plot showing HCPC of genes, samples, and conditions
#' @details Character coldata are removed prior to PCA. Please convert to
#' factors if these columns are needed.
#' @example man/examples/FactoWrapper_Example.R

FactoWrapper <- function(Hotgenes = NULL,
                         biplot = TRUE,
                         min = 3,
                         max = 5,
                         ExpressionSlots = NULL,
                         aux_features = "",
                         coldata_ids = coldata_names(Hotgenes),
                         Top_var = 10,
                         ellipse.level = 0.5,
                         ellipse.alpha = 0.5,
                         habillage_selection = "clust",
                         label_sel = c("all"),
                         pointsize = 1,
                         labelsize = 1,
                         contrasts = NULL,
                         hotList = NULL,
                         padj_cut = 0.1,
                         keep_na_padj = FALSE,
                         signif_ = NULL,
                         .log2FoldChange = 0,
                         SampleIDs = NULL) {
  # Check object
  stopifnot(is(Hotgenes, "Hotgenes"))
  
  title <- paste0("Expression: ", ExpressionSlots)
  

# making supplementary data -----------------------------------------------

  # get aux_assays
  aux_data_wide <- Hotgenes %>% 
    auxiliary_assays_(features = aux_features) %>% 
    tibble::rownames_to_column("SampleIDs")
  
  
  # merging with coldata
  supp_data <- coldata_(Hotgenes, coldata_ids = coldata_ids) %>% 
    tibble::rownames_to_column("SampleIDs") %>% 
    dplyr::left_join(aux_data_wide, by = "SampleIDs") %>%
    dplyr::select_if(~ dplyr::n_distinct(na.omit(.x)) > 1) 
  
  # Normalized expression
  # making expression data with coldata data ---------------------------------
  
  ExpressionDat <- DExps(
    Hotgenes = Hotgenes,
    ExpressionSlots = ExpressionSlots,
    SampleIDs = SampleIDs,
    hotList = hotList,
    aux_features = "",
    contrasts = contrasts,
    padj_cut = padj_cut,
    keep_na_padj = keep_na_padj,
    .log2FoldChange = .log2FoldChange,
    Topn = Inf,
    Query_set = TRUE,
    coldata_ids = ""
  ) %>%
    dplyr::select_if( ~ !is.character(.x)) %>%
    tibble::rownames_to_column("SampleIDs")
  
  output_PCA <- FactoWrapper_DFs(
    ExpressionDat = ExpressionDat,
    supp_data = supp_data,
    sampleID_col = "SampleIDs",
    Top_var = Top_var,
    biplot = biplot,
    min = min,
    max = max,
    ellipse.level = ellipse.level,
    ellipse.alpha = ellipse.alpha,
    habillage_selection = habillage_selection,
    label_sel = label_sel,
    pointsize = pointsize,
    labelsize = labelsize
  )
  
  #output_PCA$TopTibble
  
  # Setting decimal points --------------------------------------------------
  summary_ids <- c("TopTibble", "TopGroups", "TopTibble_sup")
  
  available <- Mapper_aliases(Hotgenes)
  
  if(all(length(available) > 0,
     nrow(output_PCA$TopTibble)> 0)) {
    
    
    cli::cli_inform(
      "Appending TopTibble with available aliases: {available}")
  
  
  output_PCA["TopTibble"] <- output_PCA["TopTibble"] %>%
    purrr::imap(function(x,y) {
      
    
      x <-  df_append_mapperDF(
        df = x,
        mapperDF = Mapper_(Hotgenes), 
        by = "Feature",
        relationship = "many-to-many")
      
      return(x)
      
    })
  
  } 
  
  if (!is.null(signif_)) {
    output_PCA[summary_ids] <- output_PCA[summary_ids] %>%
      purrr::imap(function(x,y) {
        
 
  
  out_x <- x %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                ~ signif(.x, signif_)))
  
        
      })
  }
  
  
  
  return(output_PCA)
}


#' FactoWrapper_DFs for PCA/HCPC of expression data
#' @rdname FactoWrapper
#' @inheritParams factoextra::fviz_pca_biplot
#' @param ExpressionDat wide expression data, rows are samples
#' @param supp_data wide format data.frame to be used as supplementary data
#' @param sampleID_col string column name for merging ExpressionDat with
#' supp_data
#' @export
#' @return FactoMiner plot showing HCPC of genes, samples, and conditions

FactoWrapper_DFs <- function(ExpressionDat = NULL,
                             supp_data = NULL,
                             sampleID_col = NULL,
                             biplot = TRUE,
                             min = 3,
                             max = 5,
                             Top_var = 10,
                             ellipse.level = 0.5,
                             ellipse.alpha = 0.5,
                             habillage_selection = "clust",
                             label_sel = c("all"),
                             pointsize = 1,
                             labelsize = 1,
                             title = NULL,
                             signif_ = NULL) {
  # Normalized expression
  # making expression data with coldata data ---------------------------------
  
  quanti_ <- supp_data %>% 
    dplyr::select(-dplyr::any_of(sampleID_col)) %>% 
    dplyr::select_if(is.numeric) %>% names()
  
  quali_ <- supp_data %>% 
    dplyr::select(-dplyr::any_of(sampleID_col)) %>% 
    dplyr::select_if(~!is.numeric(.x)) %>% names()
  
  # filtered
  filtered_coldata_ids <- c(quanti_, quali_)
  
  # ensuring only filtered_coldata_ids in
  # supp_data get passed along
  final_col <- supp_data %>%
    dplyr::select(dplyr::any_of(c(sampleID_col, filtered_coldata_ids)))
  
  
  dm_pheno <- ExpressionDat %>%
    dplyr::left_join(final_col, by = sampleID_col) %>%
    dplyr::relocate(dplyr::any_of(filtered_coldata_ids), .before = 1) %>%
    base::as.data.frame() %>%
    tibble::column_to_rownames(sampleID_col)
  
  
  # Getting supplemental variables
  # Converting quanti_ to quanti_sup
  
  if (length(quanti_) == 0) {
    quanti_sup <- NULL
  } else if (length(quanti_) > 0) {
    quanti_sup <- match(quanti_, names(dm_pheno))
  }
  
  
  # Converting quali_ to quali_sup
  if (length(quali_) == 0) {
    quali_sup <- NULL
  } else if (length(quali_) > 0) {
    quali_sup <- match(quali_, names(dm_pheno))
  }
  
  # PCA start ---------------------------------------------------------------
  Feature_check <- dm_pheno %>%
    dplyr::select(-dplyr::any_of(filtered_coldata_ids)) %>%
    names() %>%
    length()
  
  # adding check
  if (nrow(dm_pheno) < 2 | Feature_check <= 2) {
    res_PPI_pa_1 <- empty_plot()
    
    summary_PCA <- data.frame() %>% as.matrix()
    
    res <- list()
    
    Object_slots <- c("data.clust", "desc.var",
                      "desc.axes", "desc.ind", "call")
    
    res.hcpc <- Object_slots %>%
      as.list() %>%
      purrr::set_names(Object_slots)
  } else {
    res <- FactoMineR::PCA(
      dm_pheno,
      scale.unit = TRUE,
      ncp = 5,
      quali.sup = quali_sup,
      quanti.sup = quanti_sup,
      graph = FALSE
    )
    
    # Clustering
    res.hcpc <- FactoMineR::HCPC(
      res,
      nb.clust = -1,
      consol = FALSE,
      min = min,
      max = max,
      graph = FALSE
    )
    
    # summary dat
    summary_PCA <- FactoMineR::summary.PCA(
      res,
      nb.dec = 3,
      nbelements = 100,
      nbind = 100,
      ncp = 3,
      file = nullfile()
    )
    
    
    
    # Functions
    if (biplot == TRUE) {
      fviz_ <- function(...)
        factoextra::fviz_pca_biplot(...)
    } else if (biplot == FALSE) {
      fviz_ <- function(...)
        factoextra::fviz_pca_ind(...)
    }
    
    
    res_PPI_pa_1 <- fviz_(
      res,
      axes = c(1, 2),
      repel = TRUE,
      label = label_sel,
      habillage = res.hcpc$data.clust[, habillage_selection],
      col.quanti.sup = "red",
      col.var = c("black"),
      pointsize = pointsize,
      labelsize = labelsize,
      select.var = list(contrib = Top_var),
      col.ind.sup = "black",
      ellipse.alpha = ellipse.alpha,
      title = title,
      legend.title = habillage_selection,
      addEllipses = TRUE,
      ellipse.level = ellipse.level
    ) +
      ggplot2::theme_classic() +
      ggplot2::scale_shape_manual(values = rep(20,
  length(res.hcpc$data.clust[, habillage_selection])))
  }
  
  
  # set for NULL
  TopTibble_All <- tibble::as_tibble(NULL)
  TopTibble <- tibble::as_tibble(NULL)
  TopTibble_sup <- tibble::as_tibble(NULL)
  TopGroups <- tibble::as_tibble(NULL)
  Ranks <- NULL
  
  # saving Feature table as tibble
  if (all(is.list(res.hcpc$desc.var),
          "quanti" %in% names(res.hcpc$desc.var))) {
    TopTibble_All <- res.hcpc$desc.var$quanti %>%
      purrr::map( ~ .x %>%
                    base::as.data.frame() %>%
                    tibble::rownames_to_column(var = "Feature")) %>%
      purrr::list_rbind(names_to = "Cluster") %>%
      tibble::as_tibble()
  }
  
  # verify that v.test exists
  if ("v.test" %in% names(TopTibble_All)) {
    TopTibble_sup <- TopTibble_All %>%
      dplyr::filter(.data$Feature %in% quanti_) %>%
      dplyr::rename(dplyr::any_of(c(Quanti.variable = "Feature"))) %>%
      dplyr::mutate(
        "Interpretation" = dplyr::case_when(
          .data$v.test > 0 ~ "Cluster has elevated",
          .data$v.test < 0 ~ "Cluster has reduced"
        ),
        .after = "Cluster"
      ) %>%
      dplyr::mutate_at(c("Cluster", "Interpretation"), as.factor)
    
    TopTibble <- TopTibble_All %>%
      dplyr::filter(!.data$Feature %in% quanti_) %>%
      dplyr::mutate(
        "Interpretation" = dplyr::case_when(
          .data$v.test > 0 ~ "Cluster has elevated",
          .data$v.test < 0 ~ "Cluster has reduced"
        ),
        .after = "Cluster"
      ) %>%
      dplyr::mutate_at(c("Cluster", "Interpretation"), as.factor)
    
    # Making Ranks for Feature
    Ranks <- plyr::dlply(TopTibble, "Cluster", identity) %>%
      purrr::map(
        ~ .x %>%
          dplyr::arrange(dplyr::desc(.data$v.test)) %>%
          dplyr::select(dplyr::any_of(c(
            "Feature", "v.test"
          ))) %>%
          tibble::deframe()
      ) %>% 
      purrr::set_names(~glue::glue("clust_{.x}"))
  }
  
  
  # Getting category
  if (all(is.list(res.hcpc$desc.var),
          "category" %in% names(res.hcpc$desc.var))) {
    TopGroups <- res.hcpc$desc.var$category %>%
      purrr::map( ~ .x %>%
                    base::as.data.frame() %>%
                    tibble::rownames_to_column(var = "Category")) %>%
      purrr::list_rbind(names_to = "Cluster") %>%
      tibble::as_tibble() %>%
      dplyr::mutate(
        "Interpretation" = dplyr::case_when(
          .data$v.test > 0 ~ "Cluster has elevated",
          .data$v.test < 0 ~ "Cluster has reduced"
        ),
        .after = "Cluster"
      ) %>%
      dplyr::mutate_at(c("Cluster", "Interpretation"), as.factor)
  }
  
  # preparing to export
  output_PCA <- list(
    res = res,
    summary_PCA = summary_PCA,
    res.hcpc = res.hcpc,
    res_PPI_pa_1 = res_PPI_pa_1,
    TopTibble = TopTibble,
    TopTibble_sup = TopTibble_sup,
    TopGroups = TopGroups,
    Ranks = Ranks,
    quali_sup = quali_,
    quanti_sup = quanti_,
    coldata_dat = supp_data
  )
  
  
  # Setting decimal points --------------------------------------------------
  summary_ids <- c("TopTibble", "TopGroups", "TopTibble_sup")
  if (!is.null(signif_)) {
    output_PCA[summary_ids] <- output_PCA[summary_ids] %>%
      purrr::map(function(x) {
        x %>%
          dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                      ~ signif(.x, signif_)))
      })
  }
  
  
  return(output_PCA)
}


#' factoExtra_DFs a factoextra wrapper for FactoWrapper_DFs
#' @md
#' @importFrom ggplot2 theme_classic geom_point aes stat_ellipse
#' scale_color_brewer scale_shape guides guide_legend
#' @importFrom factoextra fviz_pca_ind
#' @param PCA_obj output of [FactoWrapper_DFs()]
#' @param ellipse.alpha 0 to 1 numeric value to control
#' transparency of ellipse.
#' @param ellipse.level 0 to 1 numeric value to control
#' The level at which to draw an ellipse.
#' @param habillage_id string with the qualitative variable
#' to use for ellipses. Default is "clust".
#' @param habillage_shape_id string with the qualitative variable
#' to use for individual shapes. Default is "clust".
#' @param point_size numeric value. Default is 3.
#' @param label_size numeric value. Default is 3.
#' @param legend_text_size numeric value. Default is 12.
#' @param ... additional parameter to pass to [factoextra::fviz_pca_ind()]
#' @rdname FactoWrapper
#' @export

factoExtra_DFs <- function(PCA_obj = NULL,
                           ellipse.alpha = 1,
                           ellipse.level = 0.5,
                           habillage_id = "clust",
                           habillage_shape_id = "clust",
                           point_size = 3,
                           label_size = 3,
                           legend_text_size = 12,
                           ...) {
  # getting all quali IDs + clust
  quali_IDs <- c(colnames(PCA_obj$coldata_dat), "clust")
  
  # building data.frame of these values
  quali_frame <- PCA_obj$res.hcpc$data.clust %>%
    tibble::rownames_to_column("name") %>%
    dplyr::select(dplyr::any_of(c(quali_IDs, "name")))
  
  if(label_size == 0){
    geom <- c("point")
    
  } else {
    geom <- c("point", "text")
  }
  
  # label_size
  PCA_OUT <- factoextra::fviz_pca_ind(
    PCA_obj$res,
    axes = c(1, 2),
    geom  = geom ,
    # repel = TRUE,
    # label=c("all"),
    # col.var = c("black"),
    pointsize = NA,
    addEllipses = FALSE,
    labelsize = label_size,
    # scol.ind.sup = "black",
    ellipse.alpha = ellipse.alpha,
    legend.title = habillage_id,
    ...
  )
  
  # To allows users to select var for point shape,
  # ellipse must be added after fviz_PCA_ind
  
  # adding data to update plot
  New_data <- PCA_OUT$data %>%
    dplyr::select(dplyr::all_of(c("x", "y", "name"))) %>%
    dplyr::left_join(quali_frame, by = "name")
  
  # new plot
  Final_PCA_OUT <- PCA_OUT +
    ggplot2::theme_classic() +
    
    # Adding shape and ellipse
    ggplot2::geom_point(
      data = New_data,
      na.rm = TRUE,
      mapping = ggplot2::aes(
        shape = .data[[habillage_shape_id]],
        group = .data[[habillage_id]],
        color = .data[[habillage_id]]
      ),
      size = point_size
    ) +
    
    ggplot2::stat_ellipse(
      data = New_data,
      na.rm = TRUE,
      mapping = ggplot2::aes(color = .data[[habillage_id]]),
      type = "norm",
      level = ellipse.level,
      alpha = ellipse.alpha
    ) +
    
    ggplot2::scale_color_brewer(palette = "Set2") +
    
    ggplot2::scale_shape(solid = TRUE) +
    
    ggplot2::guides(
      shape = ggplot2::guide_legend(title = habillage_shape_id,
                                    title.position = "top"),
      color = ggplot2::guide_legend(title = habillage_id,
                                    title.position = "top")
    ) +
    ggplot2::theme(
      legend.text = ggplot2::element_text(size = legend_text_size)
    )
  
  
  return(Final_PCA_OUT)
}
