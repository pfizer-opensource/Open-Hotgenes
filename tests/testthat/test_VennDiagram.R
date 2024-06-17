# load packages

# Venn_Report -----------------------------------------
# only four lists at a time
idc <- htgs %>% contrasts_()
idc[1:4]

# even when too stringent no error
htgs %>%
  DE(
    Report = "Features",
    padj_cut = 0,
    contrasts = idc[1:4]
  ) %>%
  Venn_Report()


# Report = Features
V <- htgs %>%
  DE(
    Report = "Features",
    contrasts = idc[1:4]
  ) %>%
  Venn_Report()

# view
V$vennD

# names of lists that intersect
V$Names

# all lists
V$Intsect

# Label venn diagram groups by query gene ---------------------------------

Venn_hotList <- c("CCL2", "CXCL1")

htgs %>%
  DE(contrasts = idc[1:4]) %>%
  purrr::map(function(x) {
    x %>%
      dplyr::filter(.data$Feature %in% Venn_hotList) %>%
      dplyr::pull("Feature")
  }) %>%
  Venn_Report(
    label_sep = "\n",
    show_elements = TRUE
  )
