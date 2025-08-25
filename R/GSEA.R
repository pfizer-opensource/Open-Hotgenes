#' Support functions for GSEA
#' @name GSEA_Support
NULL

#' Wrapper for fgsea
#' @rdname GSEA_Support
#' @importFrom fgsea fgsea
#' @inheritParams fgsea::fgseaMultilevel
#' @param Ranks A named list containing gene stats in
#' descending order.
#' @param ... Additional parameters for \link[fgsea]{fgsea}.
#' @export
#' @details `fgsea_()` detects for ranks that are all positive or
#' all negative and then switches to "pos" or "neg", depending on the outcome.
#' Otherwise scoreType is "std"
#' @example man/examples/GSEA_tools_Example.R
fgsea_ <- function(
    Ranks = NULL,
    pathways = NULL,
    minSize = 1,
    maxSize = Inf,
    nproc = 1,
    ...) {
  # Check if Ranks is a list
  if (!is(Ranks, "list")) {
    stop("Ranks is not a list")
  }

  # Check if all elements in Ranks are numeric
  Ranks %>%
    purrr::imap(function(x, y) {
      if (!is(x, "numeric")) {
        stop(paste0(y, " ranks are not numeric"))
      }
    })

  # Check if any ranks have missing names
  Ranks %>%
    purrr::imap(function(x, y) {
      empty_check <- names(x) %>%
        na_check() %>%
        any()

      if (empty_check) {
        stop(paste0(y, " has empty names"))
      }
    })



  # fsgea for list
  List_fgsea_Out <- Ranks %>%
    purrr::imap(function(x, y) {
      print(y)


      allSigns <- sign(x)# %>%
        #unique() %>%
        #sum()


      if (all(allSigns == 1)) {
        "allSigns == 1, switching to scoreType <-'pos' " %>%
          message()

        scoreType <- "pos"
      } else if (all(allSigns == -1)) {
        "allSigns == -1, switching to scoreType <-'neg' " %>%
          message()

        scoreType <- "neg"
      } else {
        scoreType <- "std"
      }


      fgsea::fgsea(
        pathways = pathways,
        stats = x,
        minSize = minSize,
        maxSize = maxSize,
        nproc = nproc,
        scoreType = scoreType,
        ...
      )
    })



  GSEA_Out <- list(
    Results = List_fgsea_Out,
    Ranks = Ranks,
    pathways = pathways
  )

  return(GSEA_Out)
}




#' Support functions for GSEA
#' @export
#' @inheritParams DE
#' @rdname GSEA_Support
#' @param fgseaRes A named list containing the output from
#' fgsea_.
#' @param Topn number indicating the top number of pathways to show.
#' Results grouped by sign(NES) and then Topn max
#' absolute NES selected.
#' @param mode vector indicating the details to return.
#' Options include,
#' Details (Complete pathway details per contrast),
#' leadingEdge (vector for of the column with the same name
#' Length (Number pathways per contrast),
#'
fgsea_Results <- function(
    fgseaRes = NULL,
    Topn = Inf,
    padj_cut = 0.1,
    contrasts = NULL,
    mode = "Details") {
  # Check if fgseaRes is a list
  if (!is(fgseaRes, "list")) {
    stop("fgseaRes is not a list")
  }

  # Getting Results
  Output_Results <- fgseaRes %>%
    purrr::pluck("Results")

  # Check if fgseaRes is a list
  if (!is(Output_Results, "list")) {
    stop("Output_Results is not a list")
  }

  # getting available contrasts
  available_contrasts <- names(Output_Results)

  if (is.null(contrasts)) {
    all_contrasts <- available_contrasts
  } else if (!is.null(contrasts)) {
    all_contrasts <- match.arg(contrasts, available_contrasts, several.ok = TRUE)
  }


  # match print Report
  modeID <- match.arg(mode, c("Details", "leadingEdge", "Length"))

  # Filtering Results

  GSEA_Report <- Output_Results[all_contrasts] %>%
    purrr::imap(function(x, y) {
      # Filtering by padj
      x %>%
        dplyr::as_tibble() %>%
        dplyr::filter(.data$padj < padj_cut) %>%
        # slicing by absolute NES
        dplyr::mutate(sign_NES = sign(.data$NES)) %>%
        dplyr::group_by(.data$sign_NES) %>%
        dplyr::slice_max(abs(.data$NES), n = Topn) %>%
        dplyr::ungroup() %>%
        # Then sorting
        dplyr::arrange(desc(.data$NES))
    })


  # Details
  if ("Details" %in% modeID) {
    ReportOut <- GSEA_Report %>%
      purrr::map(function(x) {
        # Filtering by padj
        x
      }) %>%
      purrr::compact()
  }

  # Length
  if ("Length" %in% modeID) {
    ReportOut <- GSEA_Report %>%
      purrr::map_chr(function(x) {
        x %>%
          dplyr::pull("pathway") %>%
          length()
      })
  }


  # leadingEdge
  if ("leadingEdge" %in% modeID) {
    ReportOut <- GSEA_Report %>%
      purrr::map(function(x) {
        x %>%
          dplyr::select(dplyr::all_of(c("pathway", "leadingEdge"))) %>%
          tibble::deframe()
      }) %>%
      purrr::compact()
  }


  return(ReportOut)
}


#' `wrap_fixed_names()` supports wrapping gene set names that
#' have '_' instead of ' '
#' @export
#' @rdname GSEA_Support
#' @importFrom stringr str_replace_all str_wrap
#' @inheritParams stringr::str_wrap
#' @inheritParams stringr::str_replace_all
wrap_fixed_names <- function(string,
                             width = 80, 
                             indent = 0,
                             exdent = 0,
                             pattern = "_", 
                             replacement = " "){
  
   out_string <- string %>% 
    stringr::str_replace_all(pattern = pattern,
                             replacement = replacement) %>% 
  stringr::str_wrap( width = width,
                     indent = indent,
                     exdent = exdent, 
                     whitespace_only = FALSE ) %>% 
    
    stringr::str_replace_all(replacement = pattern,
                             pattern = replacement) 
  
 
  return(out_string)
  
}


#' Shows top pathways by contrast
#' @export
#' @rdname GSEA_Support
#' @importFrom forcats fct_rev
#' @importFrom stats na.omit reorder
#' @importFrom stringr str_wrap
#' @inheritParams fgsea_Results
#' @inheritParams DE
#' @inheritParams stringr::str_wrap
#'

GSEA_Plots <- function(
    fgseaRes = NULL,
    Topn = 10,
    padj_cut = 0.1,
    contrasts = NULL,
    width = 50) {
  # Setting new colors
  red <- "#E7298A"
  green <- "#1F78B4"

  # Getting Results

  OutResults <- fgsea_Results(
    fgseaRes = fgseaRes,
    Topn = Topn,
    padj_cut = padj_cut,
    contrasts = contrasts,
    mode = "Details"
  )

  OutResults %>%
    purrr::imap(function(x, y) {
      filtRes <- x %>%
        dplyr::mutate(Enrichment = ifelse(.data$NES > 0,
          "Up-regulated", "Down-regulated"
        )) %>%
        dplyr::mutate(Enrichment = as.factor(.data$Enrichment)) %>%
        dplyr::mutate(Enrichment = forcats::fct_rev(.data$Enrichment)) %>%
        dplyr::mutate(pathway_breaks = wrap_fixed_names(.data$pathway,
          width = width
        )) %>%
        dplyr::as_tibble()

      if (nrow(filtRes) > 0) {
        # plot 
        ggplot2::ggplot(filtRes, 
                        ggplot2::aes(reorder(.data$pathway_breaks, .data$NES), .data$NES)) +
          ggplot2::geom_segment(ggplot2::aes(reorder(.data$pathway_breaks, .data$NES),
            xend = .data$pathway_breaks, y = 0, yend = .data$NES
          )) +
          ggplot2::geom_point(
            size = 3, ggplot2::aes(fill = .data$Enrichment),
            shape = 21, stroke = 0.6
          ) +
          ggplot2::scale_fill_manual(values = c(
            "Up-regulated" = red,
            "Down-regulated" = green
          )) +
          ggplot2::coord_flip() +
          ggplot2::labs(title = y, x = "Pathway", y = "Normalized Enrichment Score") +
          ggplot2::theme_minimal(base_size = 8) +
          ggplot2::theme(
            axis.text = ggplot2::element_text(size = 8),
            plot.margin = ggplot2::unit(c(1, 1, 1, 1), "cm")
          )
      } else if (nrow(filtRes) == 0) {
        
        empty_plot()
        
      }
    })
}


#' Returns leadingEdge genes for specified contrast and pathway
#' @export
#' @rdname GSEA_Support
#' @inheritParams fgsea_Results
#' @inheritParams DE
#' @param genesetName String indicating the pathway to return.
#' @param contrast String indicating the contrast to select.
#'

leadingGenes <- function(
    fgseaRes = NULL,
    contrast = NULL,
    genesetName = NULL) {
  fgsea_Results(
    fgseaRes = fgseaRes,
    Topn = Inf,
    padj_cut = Inf,
    contrasts = contrast,
    mode = "leadingEdge"
  ) %>%
    purrr::chuck(contrast) %>%
    purrr::chuck(genesetName)
}

#' Returns gene signature plot for specified contrast and pathway
#' @export
#' @rdname GSEA_Support
#' @importFrom rlang .env
#' @inheritParams fgsea_Results
#' @param genesetName String indicating the pathway to return.
#' @param contrast String indicating the contrast to select.
#'
plotEnrichment_ <- function(
    fgseaRes = NULL,
    contrast = NULL,
    genesetName = NULL) {
  
 

  # Getting top pathway
  TopPathway <- fgsea_Results(
    fgseaRes = fgseaRes,
    Topn = Inf,
    padj_cut = Inf,
    contrasts = contrast,
    mode = "Details"
  ) %>%
    purrr::chuck(contrast) %>%
    dplyr::filter(.data$pathway == .env$genesetName)

  # getting leadingEdge
  leadingGenes_msigdbr <- TopPathway %>%
    dplyr::pull("leadingEdge") %>%
    purrr::chuck(1)

  # getting padj
  TopPathway_padj <- TopPathway %>%
    dplyr::mutate(padj = signif(.data$padj, 3)) %>%
    dplyr::pull("padj")

  # getting NEWS
  TopPathway_NES <- TopPathway %>%
    dplyr::mutate(NES = signif(.data$NES, 3)) %>%
    dplyr::pull("NES")

  # message to show genes
  message(paste0("leading edge genes for ", genesetName))
  leadingGenes_msigdbr %>% dput()

  # Making plot
  fgPlot(
    pathway = fgseaRes$pathways[[genesetName]],
    stats = fgseaRes$Ranks[[contrast]]
  ) +

    ggplot2::labs(
      title = genesetName,
      subtitle = paste0(
        contrast,
        " (padj: ",
        TopPathway_padj,
        "; NES: ",
        TopPathway_NES,
        ")"
      )
    )
}


#' wrapper inspired by fgsea::plotEnrichment()
#' @noRd
fgPlot <- function(
    pathway = NULL,
    stats = NULL,
    gseaParam = 1,
    ticksSize = 0.2) {
  
  rnk <- rank(-stats)
  ord <- order(rnk)
  
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)

  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  pathway <- unique(pathway)
  
  gseaRes <- fgsea::calcGseaStat(statsAdj,
                                 selectedStats = pathway,
                                 returnAllExtremes = TRUE)
  
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms)) / 8
  x <- y <- NULL
  
  # Setting new colors
  red <- "#E7298A"
  green <- "#1F78B4"
  
  g <- ggplot2::ggplot(toPlot, 
                       ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(
      color = green,
      size = 0.1
    ) +
    ggplot2::geom_hline(
      yintercept = max(tops), colour = red,
      linetype = "dashed"
    ) +
    ggplot2::geom_hline(
      yintercept = min(bottoms),
      colour = red, linetype = "dashed"
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      colour = "black"
    ) +
    ggplot2::geom_line(color = green) +
    ggplot2::theme_bw() +
    ggplot2::geom_segment(data = data.frame(x = pathway), 
                          mapping = ggplot2::aes(
                            x = x,
                            y = -diff / 2, xend = x, yend = diff / 2
                          ), size = ticksSize) +
    ggplot2::theme(panel.border = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::labs(x = "rank", y = "enrichment score")
  
  g +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(plot.margin = ggplot2::unit(c(1, 1, 1, 1), "cm"))
}


#' Support functions for GSEA
#' @rdname GSEA_Support
#' @export
#' @description
#'
#' `OntologyMethods()` and `OntologyFunctions()` work together
#' to supply geneset library functions for the `Shiny_Hotgenes()` app.
#'
#' The `msigdbr_wrapper()` function is a wrapper that
#' is used supply the default Ontology_Function argument for
#' `OntologyMethods()`, although it can be used as a standalone
#' function for supplying genesets.
#'
#' @param Ontology_Function a named list containing
#' functions for returning a named list of geneset
#' for GSVA/GSEA analysis.
#' @param InputChoices a named list containing a
#' vector with string choices corresponding to ontology choices
#' for supplied ontology functions.
#' @param gene_col_choices a named list containing a
#' vector with string choices corresponding to gene id choices
#' for supplied ontology functions.
#' @param species_choices a named list containing a
#' vector with string choices corresponding to species choices
#' for supplied ontology functions.
#' @param versions a named list containing
#' Ontology database version number as a reference.

OntologyMethods <- function(
    Ontology_Function = list("msigdbr" = msigdbr_wrapper),
    InputChoices = list("msigdbr" = msigdbr_wrapper_choices()$set),
    gene_col_choices = list("msigdbr" = c(
      "gene_symbol",
      "entrez_gene", "ensembl_gene"
    )),
    species_choices = list("msigdbr" = c("human", "mouse", "rat", "dog")),
    versions = list("msigdbr" = packageVersion("msigdbr"))) {
  # check Ontology_Function
  
  if (!is(Ontology_Function, "list")) {
    stop("Ontology_Function is not a list")
  }

  # check InputChoices
  if (!is(InputChoices, "list")) {
    stop("InputChoices is not a list")
  }

  # check gene_col_choices
  if (!is(gene_col_choices, "list")) {
    stop("gene_col_choices is not a list")
  }

  # check species_choices
  if (!is(species_choices, "list")) {
    stop("species_choices is not a list")
  }

  # verify names
  if (isFALSE(all(names(Ontology_Function) == names(InputChoices)))) {
    stop("names don't match")
  }

  # Final methods
  Final_methods <- Ontology_Function %>%
    purrr::imap(function(x, y) {
      # check if function
      
      if (!is(Ontology_Function[[y]], "function")) {
        stop("function not valid")
      }

      # check if character
      if (!is(InputChoices[[y]], "character")) {
        stop("choices not valid")
      }

      # check if character gene_col_choices
      if (!is(gene_col_choices[[y]], "character")) {
        stop("choices not valid")
      }

      # check if character species_choices
      if (!is(species_choices[[y]], "character")) {
        stop("choices not valid")
      }

      # check version
      if (is.null(versions[[y]])) {
        stop("versions not valid")
      }

      list(
        Function = Ontology_Function[[y]],
        Choices = InputChoices[[y]],
        gene_cols = gene_col_choices[[y]],
        species_choices = species_choices[[y]],
        Version = versions[[y]]
      )
    })


  return(Final_methods)
}

#' Support functions for GSVA
#' @export
#' @rdname GSEA_Support
#' @param Methods output of OntologyMethods function
#' @param db string indicating the ontology database to use.
#' @param species string
#' @param set string
#' @param gene_col string
#' @param ... additional args passed to supplied methods

OntologyFunctions <- function(Methods = NULL,
                              db = NULL,
                              species = NULL,
                              set = NULL,
                              gene_col = NULL, ...) {
  
  Final_set <- match.arg(arg = set, 
                         choices = Methods[[db]]$Choices,
                         several.ok = TRUE)
  
  Final_db <- match.arg(arg = db, choices = names(Methods))
  
  Final_gene_col <- match.arg(arg = gene_col, 
                              choices = Methods[[db]]$gene_cols)
  Final_species_choices <- match.arg(arg = species,
                                     choices = Methods[[db]]$species)
  
  arg_list <- list(
    species = Final_species_choices,
    set = Final_set,
    gene_col = Final_gene_col,
    ...
  )

  out <- base::do.call(what = Methods[[db]]$Function, args = arg_list)
  
  return(out)
}

#' msigdbr_wrapper_choices() returns table of set choices.
#' @export
#' @rdname GSEA_Support
#' @importFrom msigdbr msigdbr_collections
msigdbr_wrapper_choices <- function(){
  # set up choices
  df_choices <- msigdbr::msigdbr_collections() %>% 
    # dplyr::rename(
    #   dplyr::any_of(
    #     c("gs_collection" ="gs_cat",
    #       "gs_subcollection" = "gs_subcat"
    #       )
    #   )
    # ) %>% 
 
    dplyr::mutate(set = dplyr::case_when(
      .data$gs_subcollection != "" ~ .data$gs_subcollection,
      TRUE ~ .data$gs_collection
    ))
  
  return(df_choices)
}

#' Support functions for GSVA
#' @export
#' @rdname GSEA_Support
#' @importFrom msigdbr msigdbr
#' @param species maps to the species argument of \link[msigdbr]{msigdbr}.
#' @param set maps to the gs_subcollection argument of \link[msigdbr]{msigdbr}, 
#' when not an empty string. If gs_subcollection is an empty string, uses gs_collection. 
#' Multiple okay.
#' @param gene_col string specifying id to use.
#' Choices include "gene_symbol" (default),
#' @param clean_names Logical, if TRUE (default) geneset names are passed to
#' \link[janitor]{make_clean_names}
#' @details `msigdbr_wrapper()` returns a named list containing genesets provided
#' by msigdbr. List object names are supplied by the "gs_name" column
#' and are populated by the column specify by the gene_col argument.
#'
#'

msigdbr_wrapper <- function(species = "human",
                            set = NULL,
                            gene_col = NULL,
                            clean_names = TRUE) {
  # check cols
  Final_gene_col <- match.arg(gene_col,
    choices = c(
      "gene_symbol", "entrez_gene",
      "ensembl_gene"
    ),
    several.ok = FALSE
  )

  # set up choices
  df_choices <- msigdbr_wrapper_choices()
  
  matched_set <- match.arg(
    arg = set,
    choices = df_choices$set,
    several.ok = TRUE
  )
  
  df_subcats <- df_choices %>% 
    dplyr::filter(.data$set %in% matched_set)
  
  
arg_list <- list(collection = unique(df_subcats$gs_collection ),
                 species = species)

final_sigs <- do.call(msigdbr::msigdbr, arg_list) %>% 
  dplyr::filter(.data$gs_subcollection %in% unique(df_subcats$gs_subcollection)) %>% 
 
    dplyr::select(dplyr::all_of(c(Final_gene_col, "gs_name"))) %>%
    plyr::dlply("gs_name", identity) %>%
    purrr::map(function(x) x %>%
                 dplyr::pull(Final_gene_col) %>% 
                  unique())
 

if(isTRUE(clean_names)){
  final_sigs <- final_sigs %>% 
    purrr::set_names(~janitor::make_clean_names(.x ))
  
 
}
 return(final_sigs)
}



#' Support functions for GSVA
#' @export
#' @rdname GSEA_Support
#' @importFrom msigdbr msigdbr
#' @param library_name string to call this. 
#' @param source_genesets a long format dataframe with the following 
#' columns: 'species', 'set', 'geneset_names', 'aliase_category', 'aliases', 
#' corresponding to species, geneset category, 
#' gene set names, gene alias category, and alias values (respectively).
#' @param source_requirements a named list of
#' input requirements. Must include:
#' 'species_choices', 'InputChoices', and 'gene_col_choices', 
#' popluated with vector choices for species, geneset categories, 
#' and gene alias categories (respectively). 
#' @param versions stringr for version to append
#' to database name. This will be concantenated to name
#' in shiny app.
#'


make_custom_geneset_library <- function(library_name = "custom", 
                                        source_genesets = NULL, 
                                        source_requirements = NULL,
                                        versions = NULL){
  
  
  custom_geneset_function <- function(
    species = NULL, 
    set = NULL, 
    gene_col = NULL, 
    
    source_genesets = NULL,
    source_requirements = NULL,
    
    clean_names = TRUE) {
    
   
    Final_gene_col <- match.arg(arg = gene_col, 
                                choices = unique(source_requirements$gene_col_choices),
                                several.ok = FALSE)
    
    matched_set <- match.arg(arg = set, 
                             choices = unique(source_requirements$InputChoices), 
                             several.ok = TRUE)
    
    matched_species <- match.arg(arg = species, 
                                 choices = unique(source_requirements$species_choices), 
                                 several.ok = FALSE)
    
    
    filtered_sigs <- source_genesets %>% 
      dplyr::filter(dplyr::if_any("species", ~.x == .env$matched_species)) %>% 
      dplyr::filter(dplyr::if_any("set", ~.x %in% .env$matched_set)) %>%
      dplyr::filter(dplyr::if_any("aliase_category", ~.x == .env$Final_gene_col)) %>% 
      dplyr::select(dplyr::all_of(c("geneset_names", "aliases"))) %>% 
      dplyr::distinct()
    
    
    final_sigs_list <- filtered_sigs %>%
      plyr::dlply("geneset_names", identity) %>% 
      purrr::map(function(x) x %>% 
                   dplyr::pull("aliases"))
    
    if (isTRUE(clean_names)) {
      final_sigs_list <- final_sigs_list %>%
        purrr::set_names(~janitor::make_clean_names(.x))
    }
    return(final_sigs_list)
  }
  
  
  
  OntologyMethods_new <- source_requirements %>% 
    append(list(versions = versions)) %>% 
    purrr::imap(~  list(.x) %>% 
                  purrr::set_names(library_name) )
  
  
  base::formals(custom_geneset_function)$source_requirements <- source_requirements
  base::formals(custom_geneset_function)$source_genesets <- source_genesets
  
  
  OntologyMethods_new$Ontology_Function <-  list(custom_geneset_function)%>% 
    purrr::set_names(library_name)
  
  out<-base::do.call(what = OntologyMethods, args = OntologyMethods_new)
  return(out) 
}

