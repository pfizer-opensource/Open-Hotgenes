#' Hotgenes object utilities
#' @name utilities
#' @title Hotgenes object utilities
#' @description
#' `DE()`, `DEPlot()`, `DECoefs()`, and `DExps()` are used
#' to review DE details, generate summary plots,
#' get DE coefficients, or get expression data.
#'
#' `BoxPlot()` generates a box plot for expression data
#'
NULL

#' Returns DE details from Hotgenes Object
#' @export
#' @rdname utilities
#' @importFrom utils head
#' @importFrom dplyr pull where group_map group_vars
# @importFrom tidyselect vars_select_helpers
#' @inheritParams Output_DE_
#' @param Hotgenes Hotgenes object.
#' @param Topn numeric value for the maximum number of DEGs to show.
#' Arranged by padj. Default is Inf for all available DEGs.
#' @param Report vector indicating the details to return.
#' Options include,
#' Details (Complete DE details per contrast),
#' Features (Feature column of DE details as a vector),
#' contrast_dir (Feature column of DE details as a vector) named by contrast_dir,
#' Length (Number Feature per contrast),
#' Ranks (named vector with feature-level statistics for GSEA),
#' FC (Features with log2FoldChange)
#' @param signif_ integer indicating the number of significant digits (signif)
#' to be used. If NULL (NULL), original values returned.
#' @param Rank_name String matching a column name returned by
#' Mapper_ function, default is "Feature". This parameter is only valid when
#' Report = "Ranks" or "FC" and mapFeatures = TRUE.
#' @details By default, DE results are sorted by padj (increasing). Filtering by padj
#' or log2FoldChange only applies when hotList = NULL.
#'
#' If a hotList is provided, details for all matched Features will be returned.
#' An additional column named "significant" is added to highlight features that
#' show significant differences for the given condition. Threshold for this can be
#' set by the "padj_cut" and ".log2FoldChange" arguments.
#'
#' When calling "Ranks" features are sorted (descending) by the "stats"
#' column (see Hotgenes class). When mapFeatures = TRUE and Rank_name
#' is not set to "Feature", due to mapping issues, non-unique names
#' may be produced. In this case, the rank with the highest absolute
#' value is returned via \code{\link[dplyr]{slice_max}}, with Rank_name
#' set as group.
#'
#' The dataframe supplied to mapFeatures must contain a "Feature" column,
#' which will be used for mapping to results. When Report is set to "FC"
#' the first column will be set by Rank_name argument (default is "Feature"),
#' followed by "log2FoldChange", and then the remaining columns of returned
#' by Mapper_ function.
#'
#'
#' @example man/examples/Utilities_Example.R

DE <- function(
    Hotgenes = NULL,
    contrasts = NULL,
    mapFeatures = TRUE,
    hotList = NULL,
    Report = "Details",
    padj_cut = 0.1,
    .log2FoldChange = 0,
    Topn = Inf,
    signif_ = NULL,
    Rank_name = "Feature",
    annotateSig = TRUE) {
  # Check object
  stopifnot(is(Hotgenes, "Hotgenes"))


  # Pulling Output_DE -------------------------------------------------------

  if (Report != "Details") {
    annotateSig <- FALSE
  }


  Output_DE <- Output_DE_(Hotgenes,
    contrasts = contrasts,
    mapFeatures = mapFeatures,
    hotList = hotList,
    .log2FoldChange = .log2FoldChange,
    padj_cut = padj_cut,
    annotateSig = annotateSig
  )


  # Setting decimal points --------------------------------------------------

  if (!is.null(signif_)) {
    
    # DF version
    Output_DE <- Output_DE %>%
      dplyr::mutate(dplyr::across(
        dplyr::where(is.numeric),
        ~ signif(.x, signif_)
      ))

  
  }

  # mapping Features --------------------------------------------------------
  if (isTRUE(mapFeatures)) {
    mappingDF <- Hotgenes %>% Mapper_()
    Feature_Ids <- colnames(mappingDF)
  } else if (isFALSE(mapFeatures)) {

    Feature_Ids <- c("Feature")
  }


  # checking Rank_name ------------------------------------------------------

  if (is.null(Rank_name)) {
    Rank_name_Final <- c("Feature")
  } else if (!is.null(Rank_name)) {
    Rank_name_Final <- match.arg(Rank_name,
      Feature_Ids,
      several.ok = FALSE
    )
  }


  # match print Report ------
  ReportID <- match.arg(Report, c(
    "Details", "Features", "contrast_dir",
    "Length", "Ranks", "FC"
  ))

  
  if (any(ReportID == "Details")) {
   
    ReportOut <- Output_DE %>%
      
      dplyr::group_map(~.x) %>%
      purrr::set_names(levels(Output_DE[[dplyr::group_vars(Output_DE)]])) %>%
      purrr::map(function(x) {
        x %>%
          dplyr::select(-any_of("contrast")) %>%
          tibble::as_tibble() %>%
          head(n = Topn)
      })
    return(ReportOut)
  }

  if (any(ReportID == "Length")) {
    ReportOut <- Output_DE %>%
      
      dplyr::group_map(~.x) %>%
      purrr::set_names(levels(Output_DE[[dplyr::group_vars(Output_DE)]])) %>%
      purrr::map_int(function(x) {
        x %>%
          head(n = Topn) %>%
          dplyr::pull(.data$Feature) %>%
          # to bypass issue with multiple mappings due to mapFeatures = TRUE
          unique() %>%
          length()
      })
    return(ReportOut)
  }
  
  if (any(ReportID == "Features")) {
    ReportOut <- Output_DE %>%
      
      dplyr::group_map(~.x) %>%
      purrr::set_names(levels(Output_DE[[dplyr::group_vars(Output_DE)]])) %>%
      purrr::map(function(x) {
        x %>%
          head(n = Topn) %>%
          dplyr::pull(.data$Feature) %>%
          
          unique()
      })
    return(ReportOut)
  }
  
  if (any(ReportID == "contrast_dir")) {
    ReportOut <- Output_DE %>%
      
      plyr::dlply("contrast_dir", identity) %>% 
      purrr::map(function(x) {
        x %>%
          head(n = Topn) %>%
          dplyr::pull(.data$Feature) %>%
          unique()
      })
    return(ReportOut)
  }
  
  
  if (any(ReportID == "Ranks")) {
    ReportOut <- Output_DE %>%
      
      dplyr::group_map(~.x) %>%
      purrr::set_names(levels(Output_DE[[dplyr::group_vars(Output_DE)]])) %>%
      purrr::map(function(x) {
        x %>%
          head(n = Topn) %>%
         
          dplyr::select(dplyr::all_of(c(Rank_name_Final, "stat"))) %>%
          
          dplyr::filter(
            dplyr::if_any(
              .cols = dplyr::any_of(Rank_name_Final),
              .fns = ~ !na_check(.x)
            )
          ) %>%
          dplyr::group_by(across(dplyr::any_of(Rank_name_Final))) %>%
          dplyr::slice_max(abs(.data$stat), n = 1, with_ties = FALSE) %>%
          dplyr::ungroup() %>%
          dplyr::arrange(dplyr::desc(.data$stat)) %>%
          unique() %>%
          tibble::deframe()
      })
    return(ReportOut)
  }
  
  if (any(ReportID == "FC")) {
    ReportOut <- Output_DE %>%
      
      dplyr::group_map(~.x) %>%
      purrr::set_names(levels(Output_DE[[dplyr::group_vars(Output_DE)]])) %>%
      purrr::map(function(x) {
        x %>%
          tibble::as_tibble() %>%
          head(n = Topn) %>%
          dplyr::select(dplyr::any_of(c(Feature_Ids, "log2FoldChange"))) %>%
          dplyr::relocate(
            dplyr::any_of(c(Rank_name_Final, "log2FoldChange")),
            .before = 1
          )
      })
    return(ReportOut)
  }
}

#' Returns DE summary plot from Hotgenes Object
#' @export
#' @importFrom forcats fct_relevel
#' @inheritParams Output_DE_
#' @inheritParams DE
#' @param ncol integer to specific the number of columns to pass to
#' \code{\link[ggplot2]{facet_wrap}}. Valid only when hotList != NULL.
#' @details If hotList = NULL, total DE counts will be presented in
#' a plot. However, if hotList != NULL, a plot showing the
#' -log10(padj) value for corresponding Feature will be shown. A horizontal
#' line will represent the selected -log10(padj_cut).
#' @rdname utilities

DEPlot <- function(
    Hotgenes = NULL,
    contrasts = NULL,
    padj_cut = 0.1,
    .log2FoldChange = 0,
    ncol = NULL,
    hotList = NULL) {
  # Check hotList

  if (!shiny::isTruthy(hotList)) {
    # make df
    Summary_df <- Hotgenes %>%
      Output_DE_df( 
       
        contrasts = contrasts,
        padj_cut = padj_cut,
        .log2FoldChange = .log2FoldChange
      ) %>%

      # xContrast %>%
      dplyr::mutate(DEGs = dplyr::case_when(
        .data$log2FoldChange > 0 ~ "Increased",
        TRUE ~ "Decreased"
      ), .before = 1) %>%
      dplyr::group_by(.data$contrast, .data$DEGs) %>%
      dplyr::summarise(n = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(contrast = as.character(.data$contrast)) %>%
      dplyr::select(dplyr::all_of(c("DEGs", "n", "contrast")))

    

    # Generates summary plot
    Summary_df %>%
      ggplot2::ggplot(ggplot2::aes(
        y = .data$n,
        group = .data$DEGs,
        x = .data$contrast
      )) +
      list(
        ggplot2::geom_bar(
          stat = "identity", ggplot2::aes(fill = .data$DEGs),
          width = .5, position = "dodge"
        ),
        ggplot2::xlab(""),
        ggplot2::ylab("Count"),
        ggplot2::coord_trans(y = "sqrt"),
        ggplot2::ggtitle( # "Total DE counts by contrast",
          label = glue::glue(
            "Total DE by contrast"
          ),
          subtitle = glue::glue(
            "(padj = {padj_cut}; log2FC = {.log2FoldChange})"
          )
        ),

        # scale_fill_brewer(palette = "Set2", breaks = c("Decreased", "Increased")) ,
        ggplot2::scale_fill_manual(values = c(
          "Increased" = "#FC8D62",
          "Decreased" = "#66C2A5"
        )),
        ggplot2::theme_classic(base_size = 12),
        ggplot2::theme(
          aspect.ratio = 1 / 4,
          axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
          title = ggplot2::element_text(hjust = 0.5),
          plot.margin = ggplot2::unit(c(1, 1, 1, 1), "cm")
        )
      )
  } else {
    Summary_df <- Hotgenes %>%
      Output_DE_df(
        contrasts = contrasts,
        hotList = hotList
      ) %>%
      dplyr::ungroup() %>%

      # xContrast  %>%
      dplyr::mutate(`-log10(padj)` = -log10(.data$padj)) %>%
      dplyr::select(dplyr::all_of(c("Feature", "-log10(padj)", "contrast"))) %>%
      dplyr::mutate(contrast = as.character(.data$contrast)) %>%
      dplyr::mutate(AboveThreshold = case_when(
        .data$`-log10(padj)` > -log10(padj_cut) ~ "Y",
        TRUE ~ "N"
      ) %>% as.factor())

    # }) #%>%
    # mutate(AboveThreshold  = forcats::fct_relevel(.data$AboveThreshold, c("N", "Y")))

    # Generates summary plot
    Summary_df %>%
      ggplot2::ggplot(ggplot2::aes(
        y = .data$`-log10(padj)`,
        x = .data$contrast
      )) +
      list(
        ggplot2::geom_bar(stat = "identity", 
          ggplot2::aes(fill = .data$AboveThreshold)),
        ggplot2::xlab(""),
        # ylab("stat"),
        ggplot2::geom_hline(
          yintercept = -log10(padj_cut),
          linetype = "dashed"
        ),
        ggplot2::ggtitle("Query DE by contrast"),
        # scale_fill_manual(values = c("N"="grey11","Y"= "darkorange")),
        ggplot2::scale_fill_manual(values = c("Y" = "#A6D854", "N" = "#B3B3B3")),
        ggplot2::theme_classic(base_size = 14),
        ggplot2::theme(
          aspect.ratio = 1 / 3,
          axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
          title = ggplot2::element_text(hjust = 0.5),
          # legend.position = "none",
          plot.margin = ggplot2::unit(c(1, 1, 1, 1), "cm")
        ),
        ggplot2::facet_wrap("Feature", ncol = ncol)
      )
  }
}

#' Provides a table of DE coefficients to be used for clustering
#' @export
#' @inheritParams DE
#' @param coef_ids vector of coefficients to be returned per contrast.
#' Options include, "log2FoldChange", "stat", "baseMean".
#' @param split_names Logical, if FALSE (default), contrast and coefficient
#' label will be separated by "_". If TRUE, strings will be split onto two
#' rows.
#' @param .baseMean string indicating the name of the baseMean column.
#' If "NULL" the default, the baseMean will not be returned. Options
#' include "baseMean" or "AveExpr". "AveExpr" is for older Hotgenes objects.
#' @details Features may be selected by contrast, log2FoldChange,
#' and/or adjusted palue. Alternatively, features may be selected by name
#' via the "hotList" parameter, which overrides other parameters.
#'
#' With the exception of "baseMean", coefficents will be returned for each
#' contrast. The parameter "baseMean", is the average expression of a given feature
#' across samples.
#'
#' @rdname utilities

DECoefs <- function(
    Hotgenes = NULL,
    hotList = NULL,
    contrasts = NULL,
    coef_ids = c("log2FoldChange", "stat"),
    .baseMean = NULL,
    padj_cut = 0.1,
    .log2FoldChange = 0,
    Topn = Inf,
    split_names = FALSE) {
  # Parameters are only used for feature selection.
  DE_TopIds <- Hotgenes %>%
    DE(
      hotList = hotList,
      contrasts = contrasts,
      padj_cut = padj_cut,
      .log2FoldChange = .log2FoldChange,
      Topn = Topn,
      Report = "Features"
    ) %>%
    unlist(use.names = FALSE) %>%
    unique()

  # selecting coef_ids
  Contrast_coefs <- match.arg(coef_ids, c("log2FoldChange", "stat"), several.ok = TRUE)

  Coefs_Out <- Hotgenes %>%
    DE(hotList = DE_TopIds, Report = "Details", mapFeatures = FALSE) %>%
    purrr::imap(function(x, y) {
      dat_coefs <- x %>%
        base::as.data.frame() %>%
        tibble::column_to_rownames("Feature") %>%
        dplyr::select(dplyr::all_of(Contrast_coefs))

      # splitting names
      if (isTRUE(split_names)) {
        colnames(dat_coefs) <- paste0(y, "\n", colnames(dat_coefs))
      } else if (isFALSE(split_names)) {
        colnames(dat_coefs) <- paste0(y, "_", colnames(dat_coefs))
      }

      dat_coefs %>%
        tibble::rownames_to_column("rwnames")
    }) %>%
    purrr::reduce(dplyr::full_join,
      by = "rwnames"
    ) %>%
    tibble::column_to_rownames("rwnames")

  # Getting baseMean

  if (!is.null(.baseMean)) {
    Sel_baseMean <- match.arg(.baseMean, c("baseMean", "AveExpr"))

    baseMean_dat <- Hotgenes %>%
      DE(
        hotList = DE_TopIds,
        contrasts = contrasts_(Hotgenes)[1],
        Report = "Details"
      ) %>%
      purrr::map_df(function(x) {
        x %>%
          dplyr::select(dplyr::any_of(c(Sel_baseMean, "Feature")))
      })

    Final_Coefs <- Coefs_Out %>%
      tibble::rownames_to_column("Feature") %>%
      dplyr::left_join(baseMean_dat) %>%
      tibble::column_to_rownames("Feature")
  } else if (is.null(.baseMean)) {
    Final_Coefs <- Coefs_Out
  }

  return(Final_Coefs)
}


#' Provides expression data for selected, contrasts, and samples
#' @export
#' @inheritParams DE
#' @inheritParams coldata_
#' @param SampleIDs vector of sample ids to select
#' @param ExpressionSlots name of normalized data to use.
#' If NULL (default), the first name returned by ExpressionSlots_
#' will be used.
#' See HotgenesObj.
#' @param Query_set if FALSE (Default), expression data
#' returned, with rows are genes and columns as Samples. If TRUE,
#' the data are transposed and merged with coldata data. This is
#' suitable for ggplots or linear regression analysis.
#' @param Q Logical, if TRUE complete expression data returned.
#' Overrides contrast selection. Default is FALSE.
#' @param aux_features vector of strings. Options include
#' any value returned by auxiliary_assays_features() function. If provided,
#' corresponding auxiliary_assays values will be returned along with
#' expression data. If "" (default) no auxiliary_assays data will be included.
#' Query_set must be TRUE.
#' @rdname utilities

DExps <- function(
    Hotgenes = NULL,
    ExpressionSlots = NULL,
    SampleIDs = NULL,
    hotList = NULL,
    aux_features = "",
    contrasts = NULL,
    padj_cut = 0.1,
    .log2FoldChange = 0,
    Topn = Inf,
    Query_set = FALSE,
    Q = FALSE,
    coldata_ids = coldata_names(Hotgenes)) {
  # Check object
  stopifnot(is(Hotgenes, "Hotgenes"))

  # Get normalized data ----
  ExprOptions <- Hotgenes %>%
    ExpressionSlots_()

  Matched_ExpSel <- match.arg(ExpressionSlots, ExprOptions)

  NormalizedData <- Normalized_Data_(Hotgenes,
    slot = Matched_ExpSel
  )

  # Get SampleIDs
  if (is.null(SampleIDs)) {
    Out_SampleIDs <- SampleIDs_(Hotgenes)
  } else if (!is.null(SampleIDs)) {
    Out_SampleIDs <- SampleIDs
  }

  # gene_ids

  if (is.null(hotList)) {
    FeaturesOut <- Output_DE_df(
      Hotgenes = Hotgenes,
      hotList = hotList,
      contrasts = contrasts,
     # Report = "Features",
      #mapFeatures = FALSE,
      padj_cut = padj_cut,
      .log2FoldChange = .log2FoldChange
    ) %>% 
      dplyr::slice_head(n=Topn) %>% 
      dplyr::pull("Feature") %>% unique()
    gene_ids <- unique(unlist(FeaturesOut, use.names = FALSE))
  } else {
    gene_ids <- unique(hotList)
  }

  if (isFALSE(Q)) {
    gene_ids <- gene_ids
  } else if (isTRUE(Q)) {
    gene_ids <- rownames(NormalizedData)
  }

  # Return expression matrix
  if (isFALSE(Query_set)) {
    dm <- NormalizedData[row.names(NormalizedData) %in% gene_ids,
      colnames(NormalizedData) %in% Out_SampleIDs,
      drop = FALSE
    ]
  } else if (isTRUE(Query_set)) {
    Exps_Dat <- NormalizedData[row.names(NormalizedData) %in% gene_ids,
      colnames(NormalizedData) %in% Out_SampleIDs,
      drop = FALSE
    ]

    tExps_Dat <- data.frame(t(Exps_Dat),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    # get aux_assays
      aux_data_wide <- Hotgenes %>% 
        auxiliary_assays_(features = aux_features) %>% 
        tibble::rownames_to_column("SampleIDs_rw")
        #dplyr::mutate(SampleIDs_rw = .data$SampleIDs)
    
    # merging with phenodata
    coldata_Out <- coldata_(Hotgenes, coldata_ids = coldata_ids) %>% 
      tibble::rownames_to_column("SampleIDs_rw") %>% 
      dplyr::left_join(aux_data_wide, by = "SampleIDs_rw") %>% 
     # dplyr::mutate(rw = .data$SampleIDs) %>% 
      tibble::column_to_rownames("SampleIDs_rw")

    dm <- merge(tExps_Dat, coldata_Out,
      by = "row.names"
    )

    rownames(dm) <- dm$Row.names
    dm$Row.names <- NULL
  }
  
  # dm <- tibble::as_tibble(dm) %>% 
  #   dplyr::relocate("SampleIDs", .before = 1)

  return(dm)
}






#' @rdname utilities
#' @inheritParams grDevices::boxplot.stats
#' @importFrom grDevices boxplot.stats
#' @param fillby string indicating the variable for setting
#' fill, default is NULL. Options include any variable returned by
#' the coldata_ function.
#' @export

BoxPlot <- function(
    Hotgenes = NULL,
    ExpressionSlots = NULL,
    SampleIDs = NULL,
    fillby = NULL,
    coef = 1.5) {
  # Check object
  stopifnot(is(Hotgenes, "Hotgenes"))

  # Get normalized data ----
  ExprOptions <- ExpressionSlots_(Hotgenes)
  Matched_ExpSel <- match.arg(ExpressionSlots, ExprOptions)


  InitialMeta <- coldata_(Hotgenes)

  # calc stats
  Initial_stats <- Normalized_Data_(
    Hotgenes = Hotgenes,
    slot = Matched_ExpSel
  )

  # Get SampleIDs
  if (!shiny::isTruthy(SampleIDs)) {
    Out_SampleIDs <- names(Initial_stats)
  } else if (shiny::isTruthy(SampleIDs)) {
    Initial_stats <- Initial_stats %>%
      dplyr::select(dplyr::any_of(SampleIDs))

    # returns matched names
    Out_SampleIDs <- names(Initial_stats)
  }




  # building initial stats
  Initial_stats <- Initial_stats %>%
    purrr::imap(
      function(.x, .y) {
        int_stats <- grDevices::boxplot.stats(.x,
          coef = coef, do.conf = FALSE,
          do.out = TRUE
        )

        return(int_stats)
      }
    )



  all_stats <- Initial_stats %>%
    purrr::imap(~ .x$stats %>%
      purrr::set_names(c(
        "lw", "lh", "median",
        "uh", "uw"
      )) %>%
      rbind() %>%
      tibble::as_tibble()) %>%
    purrr::list_rbind(names_to = "Samples") %>%
    dplyr::left_join(
      InitialMeta %>%
        tibble::rownames_to_column("Samples"),
      by = "Samples"
    )

  # getting outliers
  all_outliers <- Initial_stats %>%
    purrr::imap(
      ~ tibble::as_tibble(.x$out)
    ) %>%
    purrr::list_rbind(names_to = "Samples")

  # empty plot
  if (nrow(all_stats) == 0 | ncol(all_stats) == 0) {
    empty_plot()
  } else {
    P <- all_stats %>%
      ggplot(aes(
        x = .data$Samples, ymin = .data$lw, lower = .data$lh, middle = .data$median,
        upper = .data$uh, ymax = .data$uw
      )) +
      geom_point(aes(y = .data$value, x = .data$Samples),
        data = all_outliers,
        inherit.aes = FALSE
      )




    if (shiny::isTruthy(fillby)) {
      fillby <- match.arg(fillby,
        choices = names(InitialMeta)
      )


      c_breaks <- coldata_palettes(InitialMeta,
        coldata_ids = fillby
      )[[fillby]]

      # updating plot
      P <- P +
        list(
          geom_boxplot(aes(fill = .data[[fillby]]),
            outlier.shape = NA,
            stat = "identity",
            orientation = "x"
          ),
          scale_fill_manual(values = c_breaks),
          theme_classic(),
          labs(
            title = paste0(length(Out_SampleIDs), " Samples"),
            y = Matched_ExpSel
          ),
          theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            plot.margin = unit(c(1, 1, 1, 1), "cm"),
            title = element_text(hjust = 0.5)
          )
        )

      return(P)
    } else if (!shiny::isTruthy(fillby)) {
      P <- P +
        list(
          geom_boxplot(
            outlier.shape = NA,
            stat = "identity",
            orientation = "x"
          ),
          theme_classic(),
          labs(
            title = paste0(length(Out_SampleIDs), " Samples"),
            y = Matched_ExpSel
          ),
          theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            plot.margin = unit(c(1, 1, 1, 1), "cm"),
            title = element_text(hjust = 0.5)
          )
        )

      return(P)
    }






    return(P)
  }
}

