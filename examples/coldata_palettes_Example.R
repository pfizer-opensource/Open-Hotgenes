require(Hotgenes)

df_mtcars <- mtcars %>%
  dplyr::mutate(dplyr::across(c(
    "gear", "carb",
    "vs", "cyl", "am"
  ), as.factor))
df_mtcars

# get factors
df_anno <- dplyr::select_if(df_mtcars, is.factor)
df_anno

# make palette
color_pal <- coldata_palettes(df_mtcars)

df_mtcars %>%
  dplyr::select_if(is.numeric) %>%
  pheatmap::pheatmap(
    scale = "column",
    annotation_colors = color_pal,
    annotation_row = df_anno,
    border_color = "white"
  )
