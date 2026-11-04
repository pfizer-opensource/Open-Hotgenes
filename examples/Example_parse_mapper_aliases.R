library(package = Hotgenes)

hotgenes <- readRDS(file = system.file(
    "extdata", "fit_Hotgenes.RDS",
    package = "Hotgenes", mustWork = TRUE
))
features <- head(x = Mapper_(Hotgenes = hotgenes)$Feature, n = 2)
parse_mapper_aliases(Hotgenes = hotgenes, features = features)
