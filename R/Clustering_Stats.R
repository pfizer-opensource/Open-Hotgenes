
#' Build a wide dataframe from long stat to use with FactoWrapper_DFs()
#'
#' @param input_stats dataframe in long format. 
#' @param group_var string for column containing grouping var. Usually
#' contrast names from DE analysis. 
#' @param feature_var string for column containing features. Can be
#' features from DE analysis, or geneset names. 
#' @param value_var string for column containing values to cluster. Ideally,
#' with negative and positive values, such normalized enrichment scores (GSEA),
#' log2FC, or feature ranks calculated via -log10(pvalue) * sign(log2FC).
#' @importFrom stats complete.cases
#' @returns a dataframe in wide format with features as
#' rownames and column names generated using ["group_var__value_var"].
#' @param strict logical, if TRUE, incomplete observations will be dropped.
#' Default is FALSE.
#' @example examples/Clustering_Stats_Example.R
make_stat_frame <- function(
  input_stats = NULL,
  feature_var = "Feature",
  value_var = "stat",
  group_var = "group",
  strict = FALSE) {
  
  
  de <- input_stats |>
    dplyr::ungroup()
  
  output <- de |>
    dplyr::select(
      dplyr::all_of(c(feature_var, group_var)),
      dplyr::any_of(value_var)
    ) |>
    tidyr::pivot_wider(
      names_from  = .env$group_var,
      values_from = dplyr::any_of(value_var),
      names_glue  = "{.name}__{.value}",
      values_fill = NA_real_
    ) |>
    tibble::column_to_rownames(var = feature_var)
  
  if(strict){
    
    output_complete <- stats::complete.cases(output)
    output <- output[output_complete, , drop = FALSE]
    
    cli::cli_inform(
  "{sum(output_complete)} features with complete stats \\
  ({sum(!output_complete)} dropped)"
  )
    
  }
  
  return(output)
}
