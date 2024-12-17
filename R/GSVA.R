#' Support functions for GSVA
#' @importFrom GSVA gsva
#' @inheritParams GSVA::gsva
#' @inheritParams GSVA::gsvaParam
#' @inheritParams DExps
#' @inheritParams Hotgeneslimma
#' @inheritParams HotgenesUniversal
#' @param MapperCol String matching a column name returned by
#' Mapper_ function, default is "Feature".
#' @param use_weights logical, if TRUE (default) gene set sizes applied to 
#' ebayes trend. If FALSE, trend set to FALSE
#' @param use_vooma logical, if TRUE gsva matrix
#' passed to vooma before lmFit. If FALSE (default),
#' gsva matrix converted to EList object
#' and directly passed lmFit.
#' @param voomaGroup string indicating the column name to be
#' passed to voomaByGroup(). Default is NULL. Ignored when
#' use_vooma = FALSE.
#' @param ... Additional parameters for \code{\link[GSVA]{gsva}}.
#' @param contrast_matrix an optional matrix with contrasts to be passed to
#' \code{\link[limma]{contrasts.fit}}.  
#' @param method string for method to use for GSVA package. Options
#' include: c("ssgsea","gsva", "zscore", "plage").
#' @export
#' @example man/examples/Hotgeneslimma_Example.R

HotgeneSets <- function(Hotgenes = NULL,
ExpressionSlots = NULL,
geneSets = NULL,
method = c("ssgsea"),
kcdf = c("Gaussian"),
minSize = 2,
maxSize = Inf,
MapperCol = "Feature",
voomaGroup = NULL,
contrast_matrix = NULL,
use_weights = TRUE,
use_vooma = FALSE,
...) {


gsva_outList <- Hotgsva(Hotgenes = Hotgenes,
ExpressionSlots = ExpressionSlots,
geneSets = geneSets,
method = method,
kcdf = kcdf,
minSize = minSize,
maxSize = maxSize,
MapperCol = MapperCol, 
...)


# converting
gsva_es <-  gsva_outList$Expression
Mapper <- gsva_outList$Mapper
geneset_weights <- gsva_outList$geneset_weights


if(use_vooma){

glue::glue("using vooma") %>% message()

# limma method
if(is.null(voomaGroup)){
vm_exp <- limma::vooma(y = gsva_es,
design = designMatrix_(Hotgenes),
plot = TRUE) 
} else {

stopifnot(voomaGroup %in% coldata_names(Hotgenes) )

vm_exp <- limma::voomaByGroup(y = gsva_es,
group = as.character(coldata_(Hotgenes)[,voomaGroup]),
design = designMatrix_(Hotgenes),
plot = TRUE
) }

vooma_plot <- recordPlot()
dev.off() ## clean up device

}  else {
vm_exp <- new("EList", list(E = as.matrix(gsva_es),
design = designMatrix_(Hotgenes)))

vooma_plot <- NULL

}

if(use_weights){
glue::glue("using geneset weights") %>% message()


trend <- geneset_weights

} else {
trend <- FALSE
}



# make fit
fit <- limma::lmFit(vm_exp)

# using contrast matrix
if(!is.null(contrast_matrix)){
glue::glue("using contrast matrix") %>% message()

fit_final <- limma::contrasts.fit(fit = fit,
contrasts = contrast_matrix) 

fit_final <- limma::eBayes(fit = fit_final, 
trend = trend)
} else {

fit_final <- limma::eBayes(fit = fit, 
trend = trend)
}



# vm_exp
vm_exp$vooma_plot <- vooma_plot

fit_Hotgenes_base <- Hotgeneslimma(
limmafit = fit_final,
coldata = coldata_(Hotgenes),
auxiliary_assays = Hotgenes@auxiliary_assays,
Expression = vm_exp,
Expression_name = method,
Mapper = Mapper
)



return(fit_Hotgenes_base)
}

# internal function for preparing HotgeneSets
#' @noRd
Hotgsva <- function(Hotgenes = NULL,
ExpressionSlots = NULL,
geneSets = NULL,
method = c("ssgsea"),
kcdf = c("Gaussian"),
minSize = 2,
maxSize = Inf,
MapperCol = "Feature",
...){


GSVA_check <- utils::packageVersion("GSVA")

if(GSVA_check < '1.50.0'){
stop("This function requires GSVA (>= 1.50.0).")
}

# Preparing expression data
ExprOptions <- Hotgenes %>%
ExpressionSlots_()

Matched_ExpSel <- match.arg(ExpressionSlots, ExprOptions)

# prepare to remap
if (MapperCol == "Feature") {
print("using Feature col")

NormalizedData <- Normalized_Data_(
Hotgenes,
slot = Matched_ExpSel
) %>%
as.matrix()
} else if (MapperCol != "Feature") {
paste0("using ", MapperCol, " col") %>%
print()

# This converts ids from expression data to
# ids supplied in original mapper slot

mapperD <- Hotgenes %>%
Mapper_() %>%
dplyr::select(dplyr::all_of(c("Feature", MapperCol)))

NormalizedData <- Normalized_Data_(Hotgenes,
slot = Matched_ExpSel
) %>%
tibble::rownames_to_column("Feature") %>%
dplyr::left_join(mapperD, by = "Feature") %>%
# dropping columns that don't have a mapped id
dplyr::filter(!is.na(.data[[MapperCol]])) %>%
dplyr::select(-dplyr::all_of(c("Feature"))) %>%
# unique %>%
dplyr::group_by(.data[[MapperCol]]) %>%
dplyr::summarise_all(mean) %>%
tibble::column_to_rownames(MapperCol) %>%
as.matrix()
}


# build a new mapper for geneset names

print("building mapper")

Featur_s <- rownames(NormalizedData)

# check names
are_names_valid <-checkisValidAndUnreserved(names(geneSets))

if(!all(are_names_valid)){
warning("correcting invalid geneSets names")

geneSets <- geneSets %>% 
purrr::set_names(~.x %>% janitor::make_clean_names())
}


# msigdbr_Mapper
Final_Mapper <- geneSets %>%
purrr::imap(function(x, y) {
filt_x <- x[x %in% Featur_s] %>% unique()
gene_x <- paste0(filt_x, collapse = "|")

gene_x %>%
tibble::as_tibble() %>%
dplyr::rename("original_features" = "value" ) %>%
dplyr::mutate(Feature = y, 
size = length(filt_x)) %>%
dplyr::relocate("Feature", .before = 1)
}) %>% 
purrr::list_rbind()

# GSVA building
method_sel <- match.arg(method, 
choices = c("gsva", "ssgsea", "zscore", "plage"))

method_list <- list(
"gsva" = GSVA::gsvaParam, 
"ssgsea" = GSVA::ssgseaParam,
"zscore" = GSVA::zscoreParam, 
"plage" = GSVA::plageParam)[[method_sel]]

method_formalArgs <- methods::formalArgs(method_list)


# filter matching params
param_list <-list(
exprData = NormalizedData,
geneSets = geneSets,
minSize = minSize,
maxSize = maxSize,
kcdf = kcdf,
...)

obj_names <- base::intersect(
names(param_list), method_formalArgs)

# creating method object
gsva_es_method_obj <- base::do.call(what = method_list, 
args = param_list[obj_names])

# preparing gsva parameter list
gsva_params <- list(gsva_es_method_obj) %>% 
  
  # to allow for backwards compatibility
  purrr::set_names(methods::formalArgs(GSVA::gsva)[1]) %>% 
  append(param_list)

# filtering parameter
gsva_params_names <- base::intersect(
names(gsva_params), methods::formalArgs(GSVA::gsva))

# converting
gsva_params_list <- gsva_params[gsva_params_names] %>% 
  append(list(...))

gsva_es <-  base::do.call(what = GSVA::gsva, args = gsva_params_list)

# making weights
size_weights <-geneSets[rownames(gsva_es)] %>% 
lengths() %>% as.vector() %>% 
purrr::set_names(rownames(gsva_es))


gen_gsva_Out <- list(
coldata = coldata_(Hotgenes),
auxiliary_assays = Hotgenes@auxiliary_assays,
Expression = gsva_es,
Expression_name = method,
Mapper = Final_Mapper,
geneset_weights = size_weights
)

return(gen_gsva_Out)

}

# internal names check function
#' @noRd
checkisValidAndUnreserved <- function(string, ...) {
make.names(string, ...) == string
}


