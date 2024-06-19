
# Download example data -------------------------------
#library(readr)

Web_NanoString_dir <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE110102&format=file"

tempDir <- tempdir()
temp.NanoString.file <- file.path(tempDir, "NanoString.csv.gz")
download.file(Web_NanoString_dir, temp.NanoString.file, mode = "wb")

NanoString_allFiles <- untar(temp.NanoString.file, list = TRUE)
NanoString_allFiles

untar(temp.NanoString.file, exdir = tempDir)

NanoString_allFiles_paths <- file.path(tempDir, NanoString_allFiles)
NanoString_allFiles_paths

Crude_LTbRNanoString <- NanoString_allFiles_paths %>%
  purrr::map(~ read_csv(.x)) %>%
  purrr::reduce(full_join) %>%
  as.data.frame() %>%
  dplyr::rename("ID" = 1) %>%
  # remove markers
  dplyr::slice(1:255)

head(Crude_LTbRNanoString)

# Getting annotations
annot_wb <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=GPL19963&id=19167&db=GeoDb_blob125"

# accession table
acc <- readr::read_table(annot_wb) %>%
  dplyr::select(1, 3) %>%
  dplyr::slice(26:280) %>%
  dplyr::rename("ID" = 1) %>%
  dplyr::rename("accession" = 2) %>%
  dplyr::mutate(accession = gsub('["]', "", .data$accession)) %>%
  dplyr::mutate(accession = gsub(".*>(.+)..<.*", "\\1", .data$accession))

head(acc)

# updating names
dbCon <- org.Hs.eg.db::org.Hs.eg_dbconn()

# write your SQL query
sqlQuery <- "SELECT * FROM RefSeq, gene_info WHERE RefSeq._id == gene_info._id;"

# execute the query on the database
aliasSymbol <- DBI::dbGetQuery(dbCon, sqlQuery)
aliasSymbol %>% head()

# merging
result <- acc %>% dplyr::left_join(aliasSymbol %>% dplyr::select(accession, symbol))
result %>% head()

# final output
LTbRNanoString <- dplyr::left_join(Crude_LTbRNanoString, result) %>%
  dplyr::mutate(rw = dplyr::case_when(
    is.na(.data$symbol) ~ .data$ID,
    !is.na(.data$symbol) ~ .data$symbol
  )) %>%
  dplyr::select(-c(ID, accession, symbol)) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("rw")

LTbRNanoString %>% head()

# make coldata ----------------------------------------

coldata <- colnames(LTbRNanoString) %>%
  tibble::as_tibble() %>%
  dplyr::rename("SampleIDs" = "value") %>%
  dplyr::mutate(sh = dplyr::case_when(
    grepl("EWS", .data$SampleIDs) ~ "EWS",
    grepl("CON", .data$SampleIDs) ~ "Ctrl"
  )) %>%
  dplyr::mutate(Bio_Rep = dplyr::case_when(
    grepl("_1$", .data$SampleIDs) ~ "1",
    grepl("_2$", .data$SampleIDs) ~ "2"
  )) %>%
  dplyr::mutate(Hrs = dplyr::case_when(
    grepl("0hrs", .data$SampleIDs) ~ "0",
    grepl("2hrs", .data$SampleIDs) ~ "2",
    grepl("6hrs", .data$SampleIDs) ~ "6"
  )) %>%
  dplyr::mutate_at(dplyr::vars(c("sh", "Bio_Rep", "Hrs")), as.factor) %>%
  tibble::column_to_rownames("SampleIDs") %>%
  as.data.frame()

coldata

# incase you wanted to include aliases for your genes
# requires a "Feature" column that contains gene names in expression matrix
dbCon <- org.Hs.eg.db::org.Hs.eg_dbconn()
sqlQuery <- "SELECT * FROM ENSEMBL, gene_info WHERE ENSEMBL._id == gene_info._id;"

ensembl_Symbol <- DBI::dbGetQuery(dbCon, sqlQuery) %>%
  dplyr::select(dplyr::any_of(c("symbol", "ensembl_id"))) %>%
  dplyr::rename("Feature" = "symbol")

ensembl_Symbol %>% head()
# setup DESeq2 ----------------------------------------

#library(DESeq2)
library(Hotgenes)

model <- eval(~ sh * Hrs)

dds_con <- DESeq2::DESeqDataSetFromMatrix(
  countData = LTbRNanoString[c(1:255), ], # DESeq2 Object
  colData = coldata,
  design = model
)

modeMatrix <- cleanModelMatrix(
  model, # model matrix
  coldata
)

dds_con <- DESeq2::DESeq(dds_con, full = modeMatrix)
# DESeq2::resultsNames(dds_con) # lists the coefficients
# DESeq2::plotDispEsts(dds_con)

# Convert to Hotgenes Object
htgs <- Hotgenes::HotgenesDEseq2(DEseq2_object = dds_con, 
                       Mapper = ensembl_Symbol)


temp_coldata <- coldata_(htgs)

temp_coldata <- temp_coldata %>% 
  dplyr::mutate(Time = as.numeric(as.character(.data$Hrs)))

coldata_(htgs) <- temp_coldata

# Saving data ----------------------------------------------

# Create directory
dir.create(file.path(
  getwd(),
  "inst", "extdata"
), recursive = TRUE)

save(dds_con,
  file = file.path(
    getwd(),
    "inst", "extdata", "dds_con.Rdata"
  )
)

saveRDS(htgs,
  file = file.path(
    getwd(),
    "inst", "extdata", paste0("dds_Hotgenes", ".RDS")
  )
)

list.files(file.path(getwd(), "inst", "extdata"))


# Hotgeneslimma -----------------------------------------------------------
#require(DESeq2)
#require(limma)
#require(edgeR)



# Example Expression data and coldata
cts <- DESeq2::counts(dds_con) %>% as.data.frame()
Design <- SummarizedExperiment::colData(dds_con) %>%
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

# make a cleanModelMatrix
model_Matrix <- cleanModelMatrix( ~ sh * Hrs,
                              data = Design)
# voom
vm_exp <- limma::voom(d, model_Matrix, plot = TRUE)

# make fit
fit <- limma::lmFit(vm_exp, model_Matrix)
fit <- limma::eBayes(fit)

# Get alternative exps
alt_Exp <- list(counts = data.matrix(d0))


# Convert to Hotgenes Object
fit_Hotgenes <- Hotgenes::Hotgeneslimma(
  limmafit = fit,
  coldata = Design,
  Expression = vm_exp,
  Expression_name = "logCPM",
  Exps_list = alt_Exp,
  Mapper = ensembl_Symbol
)

class(vm_exp) == class(O_(fit_Hotgenes))
class(fit_Hotgenes) == "Hotgenes"

# saving example
#if (FALSE) {
  saveRDS(fit_Hotgenes,
          file = file.path(
            getwd(),
            "inst", "extdata",
            paste0("fit_Hotgenes", ".RDS")
          )
  )
#}
rm(list = ls())
