library(package = Hotgenes)

values <- c("gene", "", " ", NA_character_, "NA")
detected_missingness(x = values)
any_missingness(x = values)
all_truthy(x = c("gene1", "gene2"))
any_truthy(x = c("", "gene"))
