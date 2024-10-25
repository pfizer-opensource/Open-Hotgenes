library(Hotgenes)

# load example data
dds_Hotgenes_dir <- system.file("extdata",
  paste0("dds_Hotgenes", ".RDS"),
  package = "Hotgenes",
  mustWork = TRUE
)

# from DESeq2
HotgenesObj <- readRDS(dds_Hotgenes_dir) %>% Hotgenes::update_object()

fit_Hotgenes_dir <- system.file("extdata",
                                paste0("fit_Hotgenes", ".RDS"),
                                package = "Hotgenes",
                                mustWork = TRUE
)

# from limma
fit_Hotgenes <- readRDS(fit_Hotgenes_dir) 

# converting to list of Hotgenes
HotgenesList <- list(
  DESeq2_Hotgenes = HotgenesObj,
  limma_Hotgenes = fit_Hotgenes
) 

# app with single Hotgenes object
if(FALSE){
  
  Shiny_Hotgenes(fit_Hotgenes)
}

# app with multiple Hotgenes objects
if(FALSE){
 
  Shiny_Hotgenes(HotgenesList)
  
}


# or gsva -----------------------------------------------------------------
if(FALSE){
# msigdbr_wrapper ---------------------------------------------------------
Hotgenes::msigdbr_wrapper_choices() 

choice_set <- c("CP:KEGG", "CP:REACTOME", "CP:WIKIPATHWAYS")
choice_id <- "gene_symbol"


gsList <- msigdbr_wrapper(
  species = "human",
  set = choice_set,
  gene_col = choice_id
)

# HotgeneSets -------------------------------------------------------------
HotgenesObj %>% Mapper_()

HotgenesObj %>% contrastMatrix_()

# ?HotgeneSets
HotgeneSets_out <- HotgeneSets(
  Hotgenes = HotgenesObj,
  gset.idx.list = gsList,
  kcdf = "Gaussian",
  method = "gsva",
  min.sz = 5,
  max.sz = Inf,
  parallel.sz = 16L
)

HotgeneSets_out

# converting to list of Hotgenes
Hotgenes_sets_List <- list(
  DESeq2_Hotgenes = HotgenesObj,
  HotgeneSets = HotgeneSets_out
) 



# app with multiple Hotgenes objects

  
  Shiny_Hotgenes(Hotgenes_sets_List,
                 theme = "shinythemes/css/united.min.css")
  
}


# adding custom signatures ------------------------------------------------
# custom library ----------------------------------------------------------
  
if(FALSE){
  
  example_sigs <- msigdbr::msigdbr( species = "human", category = "H") %>% 
    dplyr::mutate_at( c("gene_symbol",  "entrez_gene", "ensembl_gene"),
                      as.character) #%>% 
  
  # create custom genesets in long format
  source_genesets <- example_sigs %>% 
    dplyr::mutate(species = "human",
                  set = .data$gs_cat,
                  geneset_names = .data$gs_name) %>% 
    dplyr::select(c("species", "set", "geneset_names", 
                    "gene_symbol",  "entrez_gene", "ensembl_gene"
    )) %>% 
    tidyr::pivot_longer(cols = c("gene_symbol",  "entrez_gene", "ensembl_gene"),
                        names_to = "aliase_category", values_to = "aliases") %>% 
    unique()
  
  # must have these column names:
  c("species", "set", "geneset_names", "aliase_category", "aliases")
  source_genesets 
  
  # create source requirements
  # this will feed into the shiny UI
  source_requirements <- list(
    species_choices = unique(source_genesets$species), 
    InputChoices = unique(source_genesets$set), 
    gene_col_choices = unique(source_genesets$aliase_category))
  
  
  # building library --------------------------------------------------------
  
  # give it a good name
  my_fun <-make_custom_geneset_library(library_name = "my_custom",
                                       version = "1", 
                                       source_genesets = source_genesets, 
                                       source_requirements = source_requirements)
  
  
  # append your custom library with default OntologyMethods 
  OntologyMethods_default <- OntologyMethods() %>% 
    append(my_fun)
  
  # now yor custom library is available in the app
  Shiny_Hotgenes(HotgenesObj,
                 OntologyMethods = OntologyMethods_default)
  
  
}

