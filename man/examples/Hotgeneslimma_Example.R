if(interactive()) {
  require(Hotgenes)
  
  # incase you wanted to include aliases for your genes
  # requires a "Feature" column that contains gene names in expression matrix
  dbCon <- org.Hs.eg.db::org.Hs.eg_dbconn()
  sqlQuery <-
    "SELECT * FROM ENSEMBL, gene_info WHERE ENSEMBL._id == gene_info._id;"
  
  ensembl_Symbol <- DBI::dbGetQuery(dbCon, sqlQuery) %>%
    dplyr::select(c("symbol", "ensembl_id")) %>%
    dplyr::rename("Feature" = .data$symbol)
  
  ensembl_Symbol %>% head()
  
  # Hotgeneslimma -----------------------------------------------------------
  require(DESeq2)
  #require(limma)
  #require(edgeR)
  
  dds_con_dir <- system.file("extdata",
                             "dds_con.Rdata",
                             package = "Hotgenes",
                             mustWork = TRUE)
  load(dds_con_dir)
  
  # Example Expression data and coldata
  cts <- counts(dds_con) %>% as.data.frame()
  Design <- colData(dds_con) %>%
    base::as.data.frame() %>%
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
  fit_Hotgenes <- Hotgeneslimma(
    limmafit = fit,
    coldata = Design,
    Expression = vm_exp,
    Expression_name = "logCPM",
    Exps_list = alt_Exp,
    Mapper = ensembl_Symbol
  )
  
  fit_Hotgenes
  
  # with contrasts ----------------------------------------------------------
  
  # make a model.matrix
  model_Matrix_2 <- model.matrix( ~ 0 + sh,
                                data = Design)
  # voom
  contrast_m <- limma_paired_contrasts(c("shEWS", "shCtrl"),
                                       modelMatrix= model_Matrix_2)
  head(contrast_m)
  
  # voom
  vm_exp <- limma::voom(d, model_Matrix_2, plot = TRUE)
  
  # make fit
  fit <- limma::lmFit(vm_exp, model_Matrix_2)
  fit2 <- limma::contrasts.fit(fit = fit, contrasts = contrast_m)
  
  fit2 <- limma::eBayes(fit2)
  
  # Get alternative exps
  alt_Exp <- list(counts = data.matrix(d0))
  
  # Convert to Hotgenes Object
  fit_Hotgenes_contrast <- Hotgeneslimma(
    limmafit = fit2,
    coldata = Design,
    Expression = vm_exp,
    Expression_name = "logCPM",
    Exps_list = alt_Exp,
    Mapper = ensembl_Symbol
  )
  
  fit_Hotgenes_contrast
  
  
  
  
}