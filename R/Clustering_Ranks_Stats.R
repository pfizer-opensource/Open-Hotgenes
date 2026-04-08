
#' Build a wide stat / log2FC / neg_log10_padj matrix from a Hotgenes object
#'
#' @param obj       Hotgenes object
#' @param contrasts character vector of contrast names to extract
#' @param padj_cut  significance threshold (default 0.05)
#' @param value_vars one or more of "stat", "log2FoldChange", "neg_log10_padj"
make_stat_wide <- function(obj, contrasts, padj_cut = 0.05,
                           value_vars = "stat") {
  de_init <- Hotgenes::Output_DE_(
    obj,
    padj_cut  = padj_cut,
    contrasts = contrasts
  ) |>
    dplyr::ungroup()
  
  de <- Hotgenes::Output_DE_(
    obj,
    padj_cut  = padj_cut,
    contrasts = contrasts,
    hotList   = unique(de_init$Feature)
  ) |>
    dplyr::ungroup()
  
  if ("neg_log10_padj" %in% value_vars) {
    de <- de |>
      dplyr::mutate(neg_log10_padj = -log10(.data$padj))
  }
  
  de |>
    dplyr::select(
      "Feature", "contrast",
      dplyr::any_of(value_vars)
    ) |>
    tidyr::pivot_wider(
      names_from  = "contrast",
      values_from = dplyr::any_of(value_vars),
      names_glue  = "{contrast}__{.value}",
      values_fill = NA_real_
    ) |>
    tibble::column_to_rownames("Feature")
}

make_col_names <- function(contrasts, value_vars) {
  tidyr::expand_grid(contrast = contrasts, value = value_vars) |>
    glue::glue_data("{contrast}__{value}") |>
    as.character()
}

run_pca_hcpc <- function(stat_wide,
                         active_cols,
                         quanti_sup_cols = NULL,
                         title,
                         ncp           = 5,
                         min           = 2,
                         max           = 5,
                         Top_var       = 20L,
                         ellipse.level = 0.5,
                         ellipse.alpha = 0.3) {
  
  active_complete <- stats::complete.cases(stat_wide[, active_cols, drop = FALSE])
  stat_wide       <- stat_wide[active_complete, , drop = FALSE]
  
  cli::cli_inform(
    "{sum(active_complete)} proteins with complete active-column stats \\
     ({sum(!active_complete)} dropped)"
  )
  
  quanti_sup_idx <- if (!is.null(quanti_sup_cols)) {
    match(quanti_sup_cols, colnames(stat_wide))
  } else {
    NULL
  }
  
  res_pca <- FactoMineR::PCA(
    stat_wide,
    scale.unit = TRUE,
    ncp        = ncp,
    quanti.sup = quanti_sup_idx,
    graph      = FALSE
  )
  
  res_hcpc <- FactoMineR::HCPC(
    res_pca,
    nb.clust = -1,
    consol   = FALSE,
    min      = min,
    max      = max,
    graph    = FALSE
  )
  
  biplot <- factoextra::fviz_pca_biplot(
    res_pca,
    axes           = c(1, 2),
    repel          = TRUE,
    label          = "all",
    habillage      = res_hcpc$data.clust[, "clust"],
    col.var        = "black",
    col.quanti.sup = "firebrick",
    select.ind     = list(contrib = Top_var),
    addEllipses    = TRUE,
    ellipse.level  = ellipse.level,
    ellipse.alpha  = ellipse.alpha
  ) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle(title)
  
  varplot <- factoextra::fviz_pca_var(
    res_pca,
    axes           = c(1, 2),
    repel          = TRUE,
    col.var        = "black",
    col.quanti.sup = "firebrick"
  ) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle(paste0(title, " [variables]"))
  
  ranks <- res_hcpc$data.clust |>
    tibble::rownames_to_column("Feature") |>
    tibble::as_tibble() |>
    dplyr::mutate(
      mean_active_stat = rowMeans(
        dplyr::pick(dplyr::all_of(active_cols)),
        na.rm = TRUE
      )
    ) |>
    dplyr::arrange(.data$clust, dplyr::desc(.data$mean_active_stat))
  
  list(
    res_pca  = res_pca,
    res_hcpc = res_hcpc,
    biplot   = biplot,
    varplot  = varplot,
    ranks    = ranks
  )
}
