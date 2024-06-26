% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExpsPlot.R
\name{ExpsPlot}
\alias{ExpsPlot}
\title{Makes expression plots via ggplot2}
\usage{
ExpsPlot(
  Hotgenes = NULL,
  ExpressionSlots = NULL,
  SampleIDs = NULL,
  group = NULL,
  fill = NULL,
  color = NULL,
  xVar = NULL,
  yVar = NULL,
  linevar = NULL,
  boxplot = FALSE,
  pointplot = TRUE,
  basemean = FALSE,
  name_col = "Feature",
  value_col = "value",
  scales = "fixed",
  facets = NULL,
  filter_eval = NULL,
  named_levels = NULL,
  ...
)
}
\arguments{
\item{Hotgenes}{Hotgenes object.}

\item{ExpressionSlots}{name of normalized data to use.
If NULL (default), the first name returned by ExpressionSlots_
will be used.
See HotgenesObj.}

\item{SampleIDs}{vector of sample ids to select}

\item{group}{String indicating the grouping variable.}

\item{fill}{String variable to set fill, default is NULL.}

\item{color}{String variable to set color, default is NULL.}

\item{xVar}{String indicating the x variable. If NULL (default),
the first column name reported by coldata_ function will be used.}

\item{yVar}{String indicating the y variable}

\item{linevar}{String variable to set line, default is NULL.}

\item{boxplot}{logical default is FALSE. If TRUE a boxplot will
be generated}

\item{pointplot}{logical if TRUE (default), point will be generated.}

\item{basemean}{Logical, if TRUE a horizontal line
indicating the mean value for all samples will be shown. If FALSE (default)
no line will be added}

\item{name_col}{string for naming column with feature names,
default is "Feature". Options include columns in mapper slot.}

\item{value_col}{string for naming column with expression
values, default is "value".}

\item{scales}{Should scales be fixed (\code{"fixed"}, the default),
free (\code{"free"}), or free in one dimension (\code{"free_x"},
\code{"free_y"})?}

\item{facets}{A set of variables or expressions quoted by \code{\link[ggplot2:vars]{vars()}}
and defining faceting groups on the rows or columns dimension.
The variables can be named (the names are passed to \code{labeller}).

For compatibility with the classic interface, can also be a
formula or character vector. Use either a one sided formula, \code{~a + b},
or a character vector, \code{c("a", "b")}.}

\item{filter_eval}{to be passed to
\code{\link[dplyr]{filter}}.}

\item{named_levels}{a named list containing the
levels to set as reference via
\code{\link[forcats]{fct_relevel}}}

\item{...}{additional arguments for
\code{\link[ggplot2]{facet_wrap}}}
}
\value{
ggplot object
}
\description{
Makes expression plots via ggplot2
}
\details{
xVar and group must be a column name accessible by coldata_
function. String for yVar must be accessible by parse_features() function.
}
\examples{
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

M_Plot \%>\%
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

M_Plot_aliases \%>\%
  ExpsPlot_Table(xVar = xvar,  
                 name_col = "ensembl_id" )


# with auxiliary assays ----------------------------------------------------


htgs <- htgs \%>\% update_object()

# making example data
set.seed(12)

max_len <- length(SampleIDs_(htgs))

AssayData  <- auxiliary_assays_default(htgs)  \%>\% 
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

htgs \%>\%
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

P \%>\%
  ExpsPlot_Table(xVar = xvar, facets = facets)

}
