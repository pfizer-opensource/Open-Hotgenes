# load packages
library(Hotgenes)
library(ggplot2)

dds_Hotgenes_dir <- system.file("extdata",
  paste0("dds_Hotgenes", ".RDS"),
  package = "Hotgenes",
  mustWork = TRUE
)
# from DESeq2
htgs <- readRDS(dds_Hotgenes_dir)

# ExpsPlot_Multi ----------------------------------------------------------

yvar <- c("CSF2", "CXCL8", "IL6")
xvar <- "Hrs"

# Multiple Features
M_Plot <- ExpsPlot(htgs,
  xVar = xvar,
  color = "Hrs",
  yVar = yvar,
  scales = "fixed"
)

M_Plot

M_Plot %>%
  ExpsPlot_Table(xVar = xvar)


# ExpsPlot with aliases ---------------------------------------------------

# name_col param must be available in Mapper
M_Plot_aliases <- ExpsPlot(htgs,
                   xVar = xvar,
                   color = "Hrs",
                   yVar = yvar,
                   scales = "fixed", 
                   name_col = "ensembl_id"
)

M_Plot_aliases

M_Plot_aliases %>%
  ExpsPlot_Table(xVar = xvar,  
                 name_col = "ensembl_id" )


# with auxiliary assays ----------------------------------------------------


htgs <- htgs %>% update_object()

# making example data
set.seed(12)

max_len <- length(SampleIDs_(htgs))

AssayData  <- auxiliary_assays_default(htgs)  %>% 
  dplyr::mutate(assay1 = rnorm(max_len),
                assay2 = rnorm(max_len))

auxiliary_assays_(htgs) <- AssayData
auxiliary_assays_(htgs)

# add colors
ExpsPlot(htgs,
         xVar = "Hrs",
         color = "Hrs",
         yVar = c("CSF2", "ENSG00000169429", "assay2")
)

# add colors
ExpsPlot(htgs,
         xVar = "Hrs",
         color = "Hrs",
         yVar = c("CSF2", "CXCL8")
)

# add colors
ExpsPlot(htgs,
         xVar = "sh",
         color = "sh",
         yVar = c("CSF2", "CXCL8")
)

# add colors
ExpsPlot(htgs,
         xVar = "sh",
         color = "sh",
         yVar = c("assay2")
)


# ExpsPlot --------------------------------------------

# yvar<-"CSF2"
xvar <- "Hrs"

# add colors
ExpsPlot(htgs,
  xVar = xvar,
  color = "Hrs",
  yVar = yvar
)

# add lines
ExpsPlot(htgs,
  xVar = xvar,
  linevar = "Bio_Rep",
  yVar = yvar
)

# set new level factor level
ExpsPlot(htgs,
  xVar = "Hrs",
  color = "sh",
  named_levels = list(Hrs = c("6")),
  yVar = yvar
)

# set subset data
ExpsPlot(htgs,
  xVar = "Hrs",
  color = "sh",
  filter_eval = sh == "EWS",
  yVar = yvar
)

# add facets
facets <- c("sh")

htgs %>%
  coldata_(coldata_ids = facets)

P <- ExpsPlot(htgs,
  xVar = xvar,
  yVar = yvar,
  # color = "sh",
  facets = facets,
  scales = "free",
  nrow = 3, labeller = label_both
) +

  theme_classic(base_size = 11) +
  theme(aspect.ratio = 1 / 2)

P

# Reformat data for graphpad ----------------------------------------------

P %>%
  ExpsPlot_Table(xVar = xvar, facets = facets)

