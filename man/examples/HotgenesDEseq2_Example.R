if (interactive()) {
  require(Hotgenes)

  # incase you wanted to include aliases for your genes
  # requires a "Feature" column that contains gene names in expression matrix
  dbCon <- org.Hs.eg.db::org.Hs.eg_dbconn()
  sqlQuery <- "SELECT * FROM ENSEMBL, gene_info WHERE ENSEMBL._id == gene_info._id;"

  ensembl_Symbol <- DBI::dbGetQuery(dbCon, sqlQuery) %>%
    dplyr::select(c("symbol", "ensembl_id")) %>%
    dplyr::rename("Feature" = .data$symbol)

  ensembl_Symbol %>% head()

  # HotgenesDEseq2 ----------------------------------------------------------
  require(DESeq2)


  dds_con_dir <- system.file("extdata",
    "dds_con.Rdata",
    package = "Hotgenes",
    mustWork = TRUE
  )
  load(dds_con_dir)

  # Example Expression data and coldata
  cts <- DESeq2::counts(dds_con) %>% as.data.frame()
  Design <- SummarizedExperiment::colData(dds_con) %>%
    base::as.data.frame() %>%
    dplyr::select_if(is.factor) %>%
    dplyr::mutate(Time = as.numeric(levels(.data$Hrs))[.data$Hrs])

  # set up DESeq2
  model_DESeq <- eval(~ sh * Hrs)

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = cts,
    colData = Design,
    design = model_DESeq
  )

  dds <- DESeq2::DESeq(dds)

  # Convert to Hotgenes Object

  dds_Hotgenes <- HotgenesDEseq2(
    DEseq2_object = dds,
    lfcShrink_type = "apeglm",
    contrasts = "sh_EWS_vs_Ctrl", # set this for speed
    ExpressionData = "log2", # set this for speed
    Mapper = ensembl_Symbol
  )

  dds_Hotgenes
  
  
}
