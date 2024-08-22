

# msigdbr_wrapper ---------------------------------------------------------
choice_set <- "CP:KEGG"
choice_id <- "gene_symbol"

gsList <- msigdbr_wrapper(
  species = "human",
  set = choice_set,
  gene_col = choice_id
)

class(gsList)
names(gsList)

# HotgeneSets -------------------------------------------------------------

HotgeneSets_out <- HotgeneSets(
  Hotgenes = htgs,
  geneSets = gsList,
  kcdf = "Gaussian",
  method = "gsva",
  minSize = 5,
  maxSize = Inf,
  BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
)


# check names -------------------------------------------------------------

gsList_not_clean <- gsList %>% 
  purrr::set_names(~.x %>% stringr::str_replace_all("_", " "))
  
names_check <-checkisValidAndUnreserved(names(gsList_not_clean))

testthat::expect_false(all(names_check))

# these are valid
names_check_valid <-checkisValidAndUnreserved(names(gsList))

testthat::expect_true(all(names_check_valid))


testthat::expect_warning(
  HotgeneSets(
    Hotgenes = htgs,
    geneSets = gsList_not_clean,
    kcdf = "Gaussian",
    method = "gsva",
    minSize = 5,
    maxSize = Inf,
    BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
  ),
  
  regexp = "correcting invalid geneSets names")



