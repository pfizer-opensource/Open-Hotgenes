require(Hotgenes)

# incase you wanted to include aliases for your genes
# requires a "Feature" column that contains gene names in expression matrix
dbCon <- org.Hs.eg.db::org.Hs.eg_dbconn()
sqlQuery <- "SELECT * FROM ENSEMBL, gene_info WHERE ENSEMBL._id == gene_info._id;"

ensembl_Symbol <- DBI::dbGetQuery(dbCon, sqlQuery) %>%
  dplyr::select(c("Feature" = "symbol", "ensembl_id"))


ensembl_Symbol %>% head()

# Hotgeneslimma -----------------------------------------------------------
require(DESeq2)
#require(limma)
#require(edgeR)


# Example Expression data and coldata
cts <- counts(dds_con) %>% as.data.frame()
Design <- colData(dds_con) %>%
  as.data.frame() %>%
  dplyr::select_if(is.factor) %>%
  dplyr::mutate(Time = as.numeric(levels(.data$Hrs))[.data$Hrs])

# Create DGEList object
# and calculate normalization factors
d0 <- edgeR::DGEList(cts)
d0 <- edgeR::calcNormFactors(d0)

# Filter low-expressed genes
# disabled in this example
if (FALSE) {
  cutoff <- 1
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  d <- d0[-drop,]
  dim(d) # number of genes lef
}

d <- d0

# make a model.matrix
model_Matrix <- model.matrix( ~ sh * Hrs,
                              data = Design)
# voom
vm_exp <- limma::voom(d, model_Matrix, plot = TRUE)

# make fit
fit <- limma::lmFit(vm_exp, model_Matrix)
fit <- limma::eBayes(fit)

# Get alternative exps
alt_Exp <- list(counts = data.matrix(d0))


# Convert to Hotgenes Object
fit_Hotgenes_noMapper <- Hotgeneslimma(
  limmafit = fit,
  coldata = Design,
  Expression = vm_exp,
  Expression_name = "logCPM",
  Exps_list = alt_Exp
)


# Convert to Hotgenes Object
df_aux_assays <- auxiliary_assays_default(fit_Hotgenes_noMapper) 

# making example data
set.seed(12)

max_len <- length(SampleIDs_(fit_Hotgenes_noMapper))

AssayDataWithExport  <- df_aux_assays %>% 
  dplyr::mutate(assay1 = rnorm(max_len),
                assay2 = rnorm(max_len))

fit_HotgenesOneContrast <- Hotgeneslimma(
  limmafit = fit,
  coldata = Design,
  contrasts = "shEWS",
  Expression = vm_exp,
  Expression_name = "logCPM",
  Exps_list = alt_Exp,
  auxiliary_assays = AssayDataWithExport,
  Mapper = ensembl_Symbol
)

testthat::expect_equivalent(
  tibble::column_to_rownames(AssayDataWithExport, "SampleIDs"),
  auxiliary_assays_(fit_HotgenesOneContrast))


# Convert to Hotgenes Object
fit_HotgenesWithAuxAssays <- Hotgeneslimma(
  limmafit = fit,
  coldata = Design,
  contrasts = "shEWS",
  Expression = vm_exp,
  Expression_name = "logCPM",
  Exps_list = alt_Exp,
  Mapper = ensembl_Symbol
)



# Convert to Hotgenes Object
fit_Hotgenes <- Hotgeneslimma(
  limmafit = fit,
  coldata = Design,
  Expression = vm_exp,
  Expression_name = "logCPM",
  Exps_list = alt_Exp,
  Mapper = ensembl_Symbol
)

class(vm_exp) == class(O_(fit_Hotgenes))
class(fit_Hotgenes) == "Hotgenes"
Mapper_(fit_Hotgenes)

fit_Hotgenes %>% DE()

# saving example
if (FALSE) {
  saveRDS(fit_Hotgenes,
          file = file.path(
            getwd(),
            "inst", "extdata",
            paste0("fit_Hotgenes", ".RDS")
          )
  )
}
