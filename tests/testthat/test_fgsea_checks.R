# Getting example data ----------------------------------------------------


# GO is C5 with BP

H_paths <- msigdbr::msigdbr(
  species = "Homo sapiens",
  # category = "C5", subcategory = "BP"
  category = "H"
) %>%
  # options for ids include: gene_symbol, entrez_gene, ensembl_gene

  dplyr::select(c("gene_symbol", "gs_name")) %>%
  plyr::dlply("gs_name", identity) %>%
  purrr::map(function(x) x %>% dplyr::pull("gene_symbol")) %>% 
  purrr::set_names(~ janitor::make_clean_names(.x ))

# These genes sets are mapped to gene symbols
# Verify Feature col contains gene symbols, too
# In this example the "Feature" column contains gene symbols

htgs %>% Mapper_()

# Get ranks
InputRanks <- htgs %>%
  DE(
    Report = "Ranks",
    contrasts = "sh_EWS_vs_Ctrl",
    Rank_name = "Feature", # see above
    padj_cut = 1, Topn = 5000
  )


# check if a rank has a missing name
InputRanks$missingGenes <- InputRanks$sh_EWS_vs_Ctrl %>%
  tibble::enframe() %>%
  dplyr::mutate(name = dplyr::case_when(
    name %in% c("JUN") ~ NA_character_,
    name %in% c("MMP3") ~ "",
    TRUE ~ name
  )) %>%
  tibble::deframe()

InputRanks$sh_EWS_vs_Ctrl[c(1:6)]
InputRanks$missingGenes[c(1:6)]

# fgsea wrapper --------
testthat::expect_error(
  fgsea_(
    Ranks = InputRanks,
    pathways = H_paths,
    minSize = 5,
    maxSize = Inf
  ),
  regexp = "has empty names"
)

InputRanks$missingGenes <- NULL

# This simulates all positive DEGs
InputRanks$TestPos <- sort(abs(InputRanks$sh_EWS_vs_Ctrl), 
                           decreasing = TRUE)
InputRanks$TestNeg <- sort(abs(InputRanks$sh_EWS_vs_Ctrl) * -1, 
                           decreasing = TRUE)

InputRanks$TestPos %>%
  sign() %>%
  unique() %>%
  sum()

InputRanks$TestNeg %>%
  sign() %>%
  unique() %>%
  sum()

InputRanks$sh_EWS_vs_Ctrl %>%
  sign() %>%
  unique() %>%
  sum()

# fgsea wrapper --------
testthat::expect_message(
  fgsea_(
    Ranks = InputRanks,
    pathways = H_paths,
    minSize = 5,
    maxSize = Inf
  ),
  regexp = "switching to scoreType"
)


H_paths <- msigdbr::msigdbr(
  species = "Homo sapiens",
  # category = "C5", subcategory = "BP"
  subcategory = "CP:KEGG"
) %>%
  # options for ids include: gene_symbol, entrez_gene, ensembl_gene

  dplyr::select(gene_symbol, gs_name) %>%
  plyr::dlply("gs_name", identity) %>%
  purrr::map(function(x) x %>% dplyr::pull("gene_symbol")) %>% 
  purrr::set_names(~ janitor::make_clean_names(.x ))


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

# OntologyMethods ---------------------------------------------------------

OutMe <- OntologyMethods(
  Ontology_Function = list(
    "msigdbr" = msigdbr_wrapper),
  InputChoices = list(
    "msigdbr" = c("CP:REACTOME", "CP:KEGG")),
  gene_col_choices = list(
    "msigdbr" = c("gene_symbol", "entrez_gene", "ensembl_gene")),
  
  species_choices = list(
    "msigdbr" = c("human", "mouse", "rat", "dog")),
  versions = list("msigdbr" = packageVersion("msigdbr"))
)

# OntologyFunctions -------------------------------------------------------


msigdbr_pthyways <- OntologyFunctions(
  Methods = OutMe,
  db = "msigdbr",
  species = "human",
  set = choice_set,
  gene_col = choice_id
)

testthat::expect_equal(gsList, msigdbr_pthyways)
testthat::expect_equal(gsList, H_paths)
