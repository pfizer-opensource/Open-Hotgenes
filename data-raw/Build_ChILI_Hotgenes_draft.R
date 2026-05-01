# Build ChILI QuantSeq Hotgenes Object ---------------------------------------
# Source: Lanz & Sheehan (2025), Cancer Immunol Immunother
#         https://doi.org/10.1007/s00262-025-04033-z
# Data:   GEO GSE287540 (public deposit)
#
# Study design:
#   QuantSeq bulk 3' RNA-seq of PBMCs from cancer patients treated with ICI.
#   Two groups: chili (developed checkpoint inhibitor-induced liver injury)
#               ctrl  (did not develop ChILI)
#   Timepoints:
#     BC  = Baseline (pre-ICI treatment)              [both groups]
#     TOL = Time Of Liver injury onset                [chili only]
#     AC  = On-treatment sample (post-ICI, no ChILI)  [ctrl only]
#     FU1 = Follow-up 1                               [chili only]
#     FU2 = Follow-up 2                               [chili only]
#
# See the original publication for full analysis details.


# 1. Libraries ---------------------------------------------------------------

# library(Hotgenes)
# library(GEOquery)
# library(DESeq2)
# library(dplyr)
# library(tibble)
# library(purrr)
# library(readr)


# 2. Download count matrix from GEO ------------------------------------------

counts_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE287nnn/",
  "GSE287540/suppl/GSE287540_ChILI_bulk_readcount.tsv.gz"
)

counts_file <- file.path(tempdir(), "GSE287540_ChILI_bulk_readcount.tsv.gz")

if (!file.exists(counts_file)) {
  download.file(counts_url, destfile = counts_file, mode = "wb")
}

raw_counts <- readr::read_tsv(counts_file, show_col_types = FALSE) |>
  tibble::column_to_rownames(var = "...1")


# 3. Fetch sample metadata from GEO -----------------------------------------

gse <- GEOquery::getGEO("GSE287540", GSEMatrix = FALSE)

gsm_list <- GEOquery::GSMList(gse)

pheno <- gsm_list |>
  purrr::imap_dfr(function(gsm, gsm_id) {
    
    meta  <- GEOquery::Meta(gsm)
    chars <- meta[["characteristics_ch1"]]
    
    extract_char <- function(key) {
      val <- chars[grepl(paste0("^", key, ":"), chars)]
      ifelse(length(val) == 0, NA_character_,
             gsub(paste0("^", key, ": *"), "", val))
    }
    
    data.frame(
      GSM       = gsm_id,
      Title     = meta[["title"]],
      iTox      = gsub("Library name: ", "", meta[["description"]]),
      Patient   = extract_char("subject"),
      Group     = extract_char("group"),
      Timepoint = extract_char("time"),
      Sex       = extract_char("Sex"),
      Age       = as.integer(extract_char("age")),
      Treatment = extract_char("treatment"),
      Steroid   = extract_char("steroid"),
      Batch     = as.integer(extract_char("batch")),
      stringsAsFactors = FALSE
    )
  }) |>
  dplyr::mutate(Group = dplyr::case_when(
    .data$Group == "control" ~ "ctrl",
    TRUE                     ~ .data$Group
  ))


# 4. Build coldata -----------------------------------------------------------

coldata <- data.frame(iTox = colnames(raw_counts)) |>
  dplyr::left_join(pheno, by = "iTox") |>
  tibble::column_to_rownames("iTox") |>
  dplyr::mutate(
    GroupTime = paste(.data$Group, .data$Timepoint, sep = "_"),
    Patient   = factor(.data$Patient),
    Group     = factor(.data$Group,     levels = c("ctrl", "chili")),
    Timepoint = factor(.data$Timepoint, levels = c("BC", "AC", "TOL",
                                                   "FU1", "FU2")),
    Steroid   = dplyr::case_when(
      .data$Steroid == "y"  ~ "yes",
      .data$Steroid == "NA" ~ "no",
      TRUE                  ~ NA_character_
    ) |> factor(levels = c("no", "yes")),
    Batch     = factor(.data$Batch),
    GroupTime = factor(.data$GroupTime,
                       levels = c("ctrl_BC",   "ctrl_AC",
                                  "chili_BC",  "chili_TOL",
                                  "chili_FU1", "chili_FU2"))
  ) |>
  dplyr::mutate(GroupTime = droplevels(.data$GroupTime))

stopifnot(all(rownames(coldata) == colnames(raw_counts)))


# 5. Build DESeqDataSet ------------------------------------------------------

use_parallel <- TRUE

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = as.matrix(raw_counts),
  colData   = coldata,
  design    = ~ GroupTime + Batch
)


# 6. Pre-filter low-count genes ----------------------------------------------

min_samples <- ceiling(ncol(dds) * 0.10)
keep_these  <- rowSums(DESeq2::counts(dds) >= 10) >= min_samples
dds         <- dds[keep_these, ]

cli::cli_inform("Genes retained after pre-filtering: {sum(keep_these)} / {length(keep_these)}")


# 7. Run DESeq2 --------------------------------------------------------------

dds <- DESeq2::DESeq(dds, parallel = use_parallel)


# 8. Gene ID mapper ----------------------------------------------------------

dbCon    <- org.Hs.eg.db::org.Hs.eg_dbconn()
sqlQuery <- "SELECT * FROM ENSEMBL, gene_info WHERE ENSEMBL._id == gene_info._id;"

mapper_df <- DBI::dbGetQuery(dbCon, sqlQuery) |>
  dplyr::select(c("Feature" = "ensembl_id", "symbol")) |>
  tibble::tibble()


# 9. Create Hotgenes object --------------------------------------------------

chili_rna <- HotgenesDEseq2(
  DEseq2_object  = dds,
  ExpressionData = c("vsd", "log2"),
  Mapper         = mapper_df,
  parallel       = use_parallel,
  contrast_vectors = list(
    ChILI_TOL_v_Ctrl_AC = c("GroupTime", "chili_TOL", "ctrl_AC"),
    ChILI_TOL_v_Ctrl_BC = c("GroupTime", "chili_TOL", "ctrl_BC")
  )
)

chili_rna
# 
# saveRDS(
#   chili_rna,
#   file = file.path(getwd(), "inst", "extdata", "chili_rna.RDS")
# )

# rm(list = ls())