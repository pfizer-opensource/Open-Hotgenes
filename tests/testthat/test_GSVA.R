

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

HotgeneSets_out



