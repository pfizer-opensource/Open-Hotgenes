

# msigdbr_wrapper ---------------------------------------------------------
#choice_set <- "CP:KEGG"
choice_set <- c("CP:REACTOME", "CP:KEGG", "CP:WIKIPATHWAYS")

choice_id <- "gene_symbol"

gsList <- msigdbr_wrapper(
  species = "human",
  set = choice_set,
  gene_col = choice_id
)

class(gsList)
names(gsList)

# test wrap_fixed_names() -------------------------------------------------

tst <- "kegg_nod_like_receptor_signaling_pathway"

tst2 <- tst %>% 
  wrap_fixed_names(width = 8) 

testthat::expect_identical(
 tst2,
"kegg\nnod_like\nreceptor\nsignaling\npathway")

# test Hotgsva _-----------------------------------------------------------

Hotgsva_O2 <- Hotgsva(
  Hotgenes = htgs,
  geneSets = gsList,
  kcdf = "Gaussian",
  method = "gsva",
  minSize = 5,
  maxSize = Inf,
  BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
)

# Hotgsva_O2 %>% names()
# Hotgsva_O2$Expression

# HotgeneSets -------------------------------------------------------------

HotgeneSets(
  Hotgenes = htgs,
  geneSets = gsList,
  kcdf = "Gaussian",
  method = "gsva",
  minSize = 5,
  maxSize = Inf,
  use_weights = TRUE,
  BPPARAM = BiocParallel::SnowParam(workers = 1)
)


HotgeneSets_out_1 <- HotgeneSets(
  Hotgenes = htgs,
  geneSets = gsList,
  kcdf = "Gaussian",
  method = "gsva",
  minSize = 5,
  maxSize = Inf,
  use_weights = TRUE,
  #use_vooma = FALSE,
  BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
)


HotgeneSets_out_2 <- HotgeneSets(
  Hotgenes = htgs,
  geneSets = gsList,
  kcdf = "Gaussian",
  method = "gsva",
  minSize = 5,
  maxSize = Inf,
  use_weights = FALSE,
 # use_vooma = TRUE,
  BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
)

# HotgeneSets_out_1
# HotgeneSets_out_2

testthat::expect_false(identical(HotgeneSets_out_1,
                                 HotgeneSets_out_2))

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



