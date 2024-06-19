require(Hotgenes)

# incase you wanted to include aliases for your genes
# requires a "Feature" column that contains gene names in expression matrix
dbCon <- org.Hs.eg.db::org.Hs.eg_dbconn()
sqlQuery <- "SELECT * FROM ENSEMBL, gene_info WHERE ENSEMBL._id == gene_info._id;"

ensembl_Symbol <- DBI::dbGetQuery(dbCon, sqlQuery) %>%
  dplyr::select(c("Feature" = "symbol", "ensembl_id"))


ensembl_Symbol %>% head()

# HotgenesDEseq2 ----------------------------------------------------------
require(DESeq2)

# Example Expression data and coldata
cts <- counts(dds_con) %>% as.data.frame()
Design <- colData(dds_con) %>%
  base::as.data.frame() %>%
  dplyr::select_if(is.factor) %>%
  dplyr::mutate(Time = as.numeric(levels(.data$Hrs))[.data$Hrs])

# set up DESeq2
model_DESeq <- eval(~ sh * Hrs)

dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = Design,
  design = model_DESeq
)

dds <- DESeq(dds)

# Convert to Hotgenes Object
# ?HotgenesDEseq2



# no mapper
dds_Hotgenes_noMapper <- HotgenesDEseq2(
  DEseq2_object = dds,
  lfcShrink_type = "apeglm",
  contrasts = "sh_EWS_vs_Ctrl", # set this for speed
  ExpressionData = "log2", # set this for speed
  Mapper = NULL
)

# no_contrast
dds_Hotgenes_no_contrast <- HotgenesDEseq2(
  DEseq2_object = dds,
  lfcShrink_type = "apeglm",
  Mapper = NULL
)

# one contrast
dds_Hotgenes <- HotgenesDEseq2(
  DEseq2_object = dds,
  lfcShrink_type = "apeglm",
  contrasts = "sh_EWS_vs_Ctrl", # set this for speed
  ExpressionData = "log2", # set this for speed
  Mapper = ensembl_Symbol
)

dds_Hotgenes

# check class
class(dds_Hotgenes) == "Hotgenes"
class(dds) == class(O_(dds_Hotgenes))
Mapper_(dds_Hotgenes)

dds_Hotgenes %>% DE()

# comparing output

dds_Hotgenes_dir <- system.file("extdata",
  paste0("dds_Hotgenes", ".RDS"),
  package = "Hotgenes",
  mustWork = TRUE
)

# from DESeq2
htgs <- readRDS(dds_Hotgenes_dir)

testthat::expect_equivalent(dds_Hotgenes, htgs)


