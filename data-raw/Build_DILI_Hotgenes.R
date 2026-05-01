# Build DILI Discovery Proteomics Hotgenes Object -------------------------
# Source: Federspiel et al. (2023), J Hepatol
# Data:   MassIVE MSV000089782 (public deposit)
#
# Differences from the original analysis:
#   The original analysis used two variables from private metadata
#   (1219_p1-p5_protein_KEY.xlsx) that are not available in the public
#   MassIVE deposit:
#
#   (1) subject: patient ID — used with duplicateCorrelation() to control
#       for repeated measures (DO and DF samples from the same patient).
#       Without subject IDs this block cannot be reproduced and is omitted.
#
#   (2) channel: exact TMT channel label (e.g. 126C, 127N) within each
#       pool — used with voomaByGroup() to compensate for channel-specific
#       variance. Pool (P1-P5) is used as a proxy here, as it captures
#       the dominant source of TMT batch variance.
#
#   All other analytical steps — VSN normalization, voomaByGroup, robust
#   lmFit, and BH-adjusted p-value threshold of 0.1 — reproduce those
#   described in Federspiel et al.


# 1. Load data ------------------------------------------------------------
# dili_raw.RDS is a pre-parsed snapshot of the public MassIVE deposit
# (MSV000089782), saved to inst/extdata to avoid a runtime download.
# To refresh it, run the commented block below.

raw_df <- readRDS(
  system.file("extdata", "dili_raw.RDS",
              package = "Hotgenes",
              mustWork = TRUE)
)

# To re-download and refresh dili_raw.RDS from MassIVE:
if (FALSE) {
  url <- paste0(
    "https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?",
    "file=f.MSV000089782%2Fupdates%2F2023-02-28_jfederspiel_1bc96582",
    "%2Fother%2FDILI_discovery_data.xlsx&forceDownload=true"
  )
  
  tmp    <- tempfile(fileext = ".xlsx")
  download.file(url, tmp, mode = "wb")
  
  raw_df <- readxl::read_excel(tmp) |>
    janitor::clean_names()
  
  saveRDS(
    raw_df,
    file = file.path(getwd(), "inst", "extdata", "dili_raw.RDS")
  )
}


# 2. Sample columns -------------------------------------------------------

sample_cols <- colnames(raw_df)[
  grepl("^(do|df|hv|nafld|ndo|ndf)_", colnames(raw_df))
]


# 3. Filter proteins ------------------------------------------------------
# - remove contaminants
# - keep reviewed UniProt entries (^sp|)
# - keep min_peps > 0
# - generate unique gene symbols as Feature (required for GSEA)

filtered_exps <- raw_df |>
  dplyr::filter(!grepl("contaminant", .data$protein)) |>
  dplyr::filter(grepl("^sp[|]", .data$protein)) |>
  dplyr::filter(.data$min_peps > 0) |>
  dplyr::mutate(Feature = make.names(.data$gene_symbol, unique = TRUE))


# 4. Expression matrix ----------------------------------------------------

expr_matrix <- filtered_exps |>
  dplyr::select("Feature", dplyr::any_of(sample_cols)) |>
  tibble::column_to_rownames("Feature") |>
  as.matrix()


# 5. Protein ID mapper ----------------------------------------------------

mapper_df <- filtered_exps |>
  dplyr::select(
    "Feature",
    "Gene"        = "gene_symbol",
    "Protein"     = "protein",
    "Description" = "description"
  )


# 6. Sample metadata (coldata) --------------------------------------------

coldata <- data.frame(
  Sample    = sample_cols,
  Condition = toupper(sub("_p[0-9]+_[0-9]+$", "", sample_cols)),
  Pool      = toupper(sub(".*_(p[0-9]+)_.*", "\\1", sample_cols)),
  row.names = sample_cols,
  stringsAsFactors = TRUE
)

coldata[["Condition"]] <- factor(
  coldata[["Condition"]],
  levels = c("HV", "DO", "DF", "NDO", "NDF", "NAFLD")
)


# 7. VSN normalization ----------------------------------------------------

expr_matrix[is.na(expr_matrix)] <- 0
expr_matrix_vsn <- limma::normalizeVSN(expr_matrix)


# 8. Design matrix --------------------------------------------------------

design <- model.matrix(~ 0 + Condition, data = coldata)
colnames(design) <- gsub("Condition", "", colnames(design))


# 9. voomaByGroup ---------------------------------------------------------

vm_exp <- limma::voomaByGroup(
  y      = expr_matrix_vsn,
  group  = coldata[["Pool"]],
  design = design,
  plot   = FALSE
)


# 10. Robust lmFit --------------------------------------------------------
# Note: duplicateCorrelation() is omitted — subject IDs required to model
# repeated measures (DO/DF pairing) are not available in the public data.

fit <- limma::lmFit(
  vm_exp,
  design,
  method = "robust"
)


# 11. Contrasts -----------------------------------------------------------

contrasts_mat <- limma::makeContrasts(
  DO_vs_HV    = DO    - HV,
  DF_vs_HV    = DF    - HV,
  NDO_vs_HV   = NDO   - HV,
  NDF_vs_HV   = NDF   - HV,
  NAFLD_vs_HV = NAFLD - HV,
  DF_vs_DO    = DF    - DO,
  NDO_vs_DO   = NDO   - DO,
  levels = design
)

fit2 <- limma::contrasts.fit(fit, contrasts_mat)
fit2 <- limma::eBayes(fit2)


# 12. Hotgenes object -----------------------------------------------------

dili_hotgenes <- Hotgenes::Hotgeneslimma(
  limmafit        = fit2,
  coldata         = coldata,
  Expression      = vm_exp,
  Expression_name = "VSN",
  Exps_list       = list(log2 = log2(expr_matrix + 1)),
  Mapper          = mapper_df
)

dili_hotgenes

#rm(list = ls())