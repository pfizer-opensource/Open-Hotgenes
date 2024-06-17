

# DEphe -------------------------------------------------------------------

ss <- SampleIDs_(htgs)[1:8]
ss

htgs %>%
  DEphe(
    Topn = 3,
    SampleIDs = ss,
    annotation_colors = coldata_palettes(htgs),
    arrangeby = c("Hrs", "sh"),
    annotations = c("Hrs", "sh"),
    cellheight = 10,
    cellwidth = 8
  )

# label by alias in mapper
htgs %>%
  DEphe(
    Topn = 3,
    SampleIDs = ss,
    annotation_colors = coldata_palettes(htgs),
    arrangeby = c("Hrs", "sh"),
    annotations = c("Hrs", "sh"),
    label_by = "ensembl_id",
    cellheight = 10,
    cellwidth = 8
  )

# when too stringent, no error
htgs %>%
  DEphe(
    Topn = 0,
    SampleIDs = ss,
    annotation_colors = coldata_palettes(htgs),
    arrangeby = c("Hrs", "sh"),
    annotations = c("Hrs", "sh"),
    cellheight = 10,
    cellwidth = 8
  )


htgs %>%
  DEphe(
    Topn = 10,
    SampleIDs = "ss_noexist",
    annotation_colors = coldata_palettes(htgs),
    arrangeby = c("Hrs", "sh"),
    annotations = c("Hrs", "sh"),
    cellheight = 10,
    cellwidth = 8
  )

