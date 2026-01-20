
library(DRomics)
library(Hotgenes)

# example with RNAseqdata -------------------------------------------------
datafilename <- system.file("extdata", "RNAseq_sample.txt", package="DRomics")
o <- DRomics::RNAseqdata(datafilename, check = TRUE, transfo.method = "vst")
s_quad <- DRomics::itemselect(o, select.method = "quadratic", FDR = 0.05)

f <- DRomics::drcfit(itemselect = s_quad, parallel = "no")

bmdcalc_out <- DRomics::bmdcalc(f)

# optional mapper
dbCon <- org.Mm.eg.db::org.Mm.eg_dbconn()

# Your query - joining ENSEMBL and gene_info
sqlQuery <- "SELECT * FROM ensembl, refseq, gene_info 
             WHERE ensembl._id == gene_info._id 
             AND refseq._id == gene_info._id;"

mapper <- DBI::dbGetQuery(dbCon, sqlQuery) %>%
  
  dplyr::select(c("Feature" = "accession", "ensembl_id", "gene_name", "symbol")) %>% 
  tibble::tibble()


hotDR_RNAseqdata <- HotgenesDRomics(
    bmdcalc = bmdcalc_out,
    Mapper = mapper)


# example with continuousomicdata -----------------------------------------
datafilename2 <- system.file("extdata", "metabolo_sample.txt", package = "DRomics")

con_o <- DRomics::continuousomicdata(datafilename2)

s_quad_2 <- DRomics::itemselect(con_o, select.method = "quadratic", FDR = 0.05)
f_2 <- DRomics::drcfit(itemselect = s_quad_2, parallel = "no")
bmdcalc_f_2 <- DRomics::bmdcalc(f_2)

hotDR_continuous <- HotgenesDRomics(
  bmdcalc = bmdcalc_f_2)

hotDR_continuous
