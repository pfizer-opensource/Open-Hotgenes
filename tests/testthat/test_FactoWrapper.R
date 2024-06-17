# load packages
library(Hotgenes)
library(ggplot2)



# FactoWrapper --------------------------------------
# using FactoWrapper


# which DEGs do you want to use?
contrasts_(htgs) # shows contrasts in hotgenes object

# columns containing characters will be shown in PCA
coldata_(htgs)
coldata_(htgs)$text_Hrs <- as.character(coldata_(htgs)$Hrs)
coldata_(htgs)


# run PCA
FactoOutput <- FactoWrapper(htgs,
  contrasts = "Hrs_6_vs_0",
  coldata_ids = c("text_Hrs", "Hrs", "Time"),
  padj_cut = 0.1,
  .log2FoldChange = 0
)

# getting HCPC details
FactoOutput$TopTibble # Feature
FactoOutput$TopTibble_sup # Coldata Quanti.variable
FactoOutput$TopGroups # Coldata Factors

# verifies Hrs, but not text_Hrs is detected
dat_names <- names(FactoOutput$res.hcpc$data.clust)
testthat::expect_true("Hrs" %in% dat_names)
testthat::expect_true("text_Hrs" %in% dat_names)

# Feature Ranks
FactoOutput$Ranks

# Getting FactoMiner PCA object
res.pca <- FactoOutput$res

# loadings
sweep(res.pca$var$coord, 2, sqrt(res.pca$eig[1:ncol(res.pca$var$coord), 1]), FUN = "/")

# Should be empty
FactoOutput$TopTibble %>%
  dplyr::filter(.data$Feature == "Time")

FactoOutput$TopGroups %>%
  dplyr::slice(1:10)

FactoOutput$TopTibble_sup


# update plot
factoextra::fviz_pca_ind(FactoOutput$res,
  axes = c(1, 2),
  repel = FALSE, label = "none",
  habillage = FactoOutput$res.hcpc$data.clust[, c("clust")],
  col.quanti.sup = "black",
  col.var = c("black"),
  pointsize = 3,
  labelsize = 4,
  select.ind = list(contrib = 50),
  col.ind.sup = "black",
  ellipse.alpha = 0,
  # title = title,
  legend.title = "clust",
  addEllipses = TRUE, ellipse.level = 0.5
) +
  theme_classic() +
  scale_shape_manual(values = rep(20, 100))

# sup
FactoOutput$TopTibble_sup # cluster stats

# Get clusts
FactoOutput$TopGroups # cluster stats
FactoOutput$TopTibble %>%
  dplyr::filter(.data$Cluster == "4")



# unit testing minimum features -------------------------------------------

FactoWrapper(htgs,
  hotList = Features_(htgs)[0],
  coldata_ids = "",
  aux_features = ""
)

FactoWrapper(htgs,
  hotList = Features_(htgs)[1:1],
  coldata_ids = "",
  aux_features = ""
)

FactoWrapper(htgs,
  hotList = Features_(htgs)[1:2],
  coldata_ids = "",
  aux_features = ""
)

FactoWrapper(htgs,
  SampleIDs = "sd",
  coldata_ids = "",
  aux_features = ""
)


# testing  ----------------------------------------------------------------

if (FALSE) {
  require(FactoMineR)
  require(factoextra)
  require(tidyverse)

  x <- as.list(formals(FactoWrapper))
  list2env(x, .GlobalEnv)
  Hotgenes <- htgs
  contrasts <- "Hrs_6_vs_0"
  coldata_ids <- c("text_Hrs", "Hrs", "Time")
  padj_cut <- 0.1
  .log2FoldChange <- 0
  # biplot<-FALSE
  label_sel <- "all"
}
