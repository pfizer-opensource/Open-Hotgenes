library(package = Hotgenes)

hotgenes <- readRDS(file = system.file(
    "extdata", "fit_Hotgenes.RDS",
    package = "Hotgenes", mustWork = TRUE
))
updated <- update_object(Hotgenes = hotgenes)
methods::validObject(object = updated)
