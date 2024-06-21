

# no error here
BoxPlot(fit_Hotgenes,
  SampleIDs = SampleIDs_(fit_Hotgenes)[0]
)


# summary plot
fit_Hotgenes %>% DEPlot()

# check genes of interest
fit_Hotgenes %>%
  DEPlot(hotList = c("IL6", "CSF2"))

# top 5 with details
fit_Hotgenes %>%
  DE(
    Report = "Details",
    Topn = 5
  )

# genes of interest
fit_Hotgenes %>%
  DE(
    Report = "Details",
    hotList = c("CSF2", "CXCL8", "TRAF2")
  )

# Export Ranks for GSEA
fit_Hotgenes %>%
  DE(Report = "Ranks")

# or FC
fit_Hotgenes %>%
  DE(Report = "FC")

# Expression Plots
fit_Hotgenes %>%
  ExpsPlot(
    yVar = "IL6",
    xVar = "sh"
  ) +
  ggplot2::facet_wrap("Hrs",
    labeller = ggplot2::label_both
  )

# DECoefs ---------------------------------------------
DECoefs(fit_Hotgenes) %>% head()

# DExps -----------------------------------------------
# list available data
ExpressionSlots_(fit_Hotgenes)

# only matched and unique values returned
DExps(fit_Hotgenes,
  hotList = c("Doesn't Exist", "TUBB", "TUBB"),
  coldata_ids = c("Hrs", "sh"),
  Query_set = TRUE
) %>% head()

# all data merged as data.frame
DExps(fit_Hotgenes,
  Query_set = TRUE, Q = TRUE
) %>% head()


# BoxPlot -----------------------------------------------------------------
# View plot
BoxPlot(fit_Hotgenes)

# list available data
ExpressionSlots_(fit_Hotgenes)




# ExpsPlot ----------------------------------------------------------------

yvar <- c("CSF2", "CXCL8", "IL6")
xvar <- "Hrs"

# name_col param must be available in Mapper_()
M_Plot_aliases <- ExpsPlot(fit_Hotgenes,
                           xVar = xvar,
                           color = "Hrs",
                           yVar = yvar,
                           scales = "fixed", 
                           name_col = "ensembl_id"
)

testthat::expect_true("ensembl_id" %in% names(M_Plot_aliases$data))

