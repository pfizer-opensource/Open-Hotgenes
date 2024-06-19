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

