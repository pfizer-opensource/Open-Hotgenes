library(Hotgenes)
library(shinytest2)
library(magick)

# ── 1. Load object ────────────────────────────────────────────────────────────
dds_Hotgenes_dir <- system.file("extdata",
                                paste0("dds_Hotgenes", ".RDS"),
                                package = "Hotgenes",
                                mustWork = TRUE
)
HotgenesObj <- readRDS(dds_Hotgenes_dir) %>% Hotgenes::update_object()

# ── 2. Dynamically extract values ─────────────────────────────────────────────
all_contrasts    <- contrasts_(HotgenesObj)
all_expr_slots   <- ExpressionSlots_(HotgenesObj)
all_coldata_cols <- coldata_names(HotgenesObj)
all_features     <- Background_(HotgenesObj, Col = "Feature")
all_mapper_cols  <- names(Mapper_(HotgenesObj))

default_contrast    <- all_contrasts[1]
default_expr_slot   <- all_expr_slots[1]
default_coldata_col <- all_coldata_cols[1]
default_feature     <- all_features[1]
venn_contrasts      <- all_contrasts[seq_len(min(2, length(all_contrasts)))]
pca_contrasts       <- all_contrasts[seq_len(min(3, length(all_contrasts)))]

cat("=== Dynamic values extracted from HotgenesObj ===\n")
cat("Contrasts:          ", paste(all_contrasts,    collapse = ", "), "\n")
cat("Expr slots:         ", paste(all_expr_slots,    collapse = ", "), "\n")
cat("Coldata cols:       ", paste(all_coldata_cols,  collapse = ", "), "\n")
cat("Mapper cols:        ", paste(all_mapper_cols,   collapse = ", "), "\n")
cat("Default contrast:   ", default_contrast,        "\n")
cat("Default expr slot:  ", default_expr_slot,       "\n")
cat("Default coldata col:", default_coldata_col,     "\n")
cat("Venn contrasts:     ", paste(venn_contrasts, collapse = ", "), "\n")
cat("PCA contrasts:      ", paste(pca_contrasts,  collapse = ", "), "\n")

# ── 3. Helpers ────────────────────────────────────────────────────────────────
ns <- function(tab, input_id) {
  paste("Hotgenes_A", tab, input_id, sep = "-")
}

inp <- function(tab, ...) {
  args <- list(...)
  setNames(args, paste("Hotgenes_A", tab, names(args), sep = "-"))
}

take_shot <- function(tab_id, filename, inputs = list(), btn = NULL, sleep = 1.5) {
  app$set_inputs(tabs = tab_id, wait_ = FALSE)
  Sys.sleep(sleep)
  
  if (length(inputs) > 0) {
    do.call(app$set_inputs, c(inputs, list(wait_ = FALSE)))
    Sys.sleep(sleep)
  }
  
  if (!is.null(btn)) {
    app$click(btn)
    Sys.sleep(sleep * 2)
  }
  
  path <- file.path(fig_dir, filename)
  app$get_screenshot(path)
  message("Saved screenshot: ", filename)
}

# annotate a screenshot with a title label for the manuscript
annotate_fig <- function(filename,
                         label,
                         output_filename = NULL,
                         label_size      = 40,
                         bg_color        = "white",
                         text_color      = "black") {
  
  if (is.null(output_filename)) output_filename <- filename
  
  input_path  <- file.path(fig_dir, filename)
  output_path <- file.path(fig_dir, output_filename)
  
  img <- image_read(input_path)
  
  # get dimensions
  info  <- image_info(img)
  width <- info$width
  
  # create a label banner to prepend above the image
  banner <- image_blank(width = width, height = 80, color = bg_color) |>
    image_annotate(
      text     = label,
      gravity  = "West",
      location = "+20+0",
      size     = label_size,
      color    = text_color,
      font     = "Helvetica",
      weight   = 700
    )
  
  # stack banner on top of screenshot
  image_append(c(banner, img), stack = TRUE) |>
    image_write(output_path, format = "png")
  
  message("Annotated: ", output_filename)
}

# ── 4. Launch app ─────────────────────────────────────────────────────────────
app <- AppDriver$new(
  app    = Shiny_Hotgenes(HotgenesObj),
  height = 900,
  width  = 1400
)

# 2x device pixel ratio for retina-quality screenshots
app$get_chromote_session()$Emulation$setDeviceMetricsOverride(
  width             = 1400,
  height            = 900,
  deviceScaleFactor = 2,
  mobile            = FALSE
)

# ── 5. Output dirs ────────────────────────────────────────────────────────────
fig_dir <- file.path(tempdir(), "hotgenes_figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── 6. BoxPlot ────────────────────────────────────────────────────────────────
take_shot(
  tab_id   = "Hotgenes_A-BoxPlot",
  filename = "fig1_boxplot_raw.png",
  inputs   = inp("BoxPlot",
                 NormSlot = default_expr_slot
  )
)
annotate_fig(
  filename        = "fig1_boxplot_raw.png",
  label           = "A  Expression Data (BoxPlot) — normalised per-sample distributions",
  output_filename = "fig1_boxplot.png"
)

# ── 7. DEstats ────────────────────────────────────────────────────────────────
take_shot(
  tab_id   = "Hotgenes_A-DEstats",
  filename = "fig2_destats_raw.png",
  inputs   = inp("DEstats",
                 DE_Contrasts       = default_contrast,
                 padj_cut_DE_tables = 0.1,
                 lfc_Statistics     = 0
  )
)
annotate_fig(
  filename        = "fig2_destats_raw.png",
  label           = paste0("B  Differential Expression — contrast: ", default_contrast),
  output_filename = "fig2_destats.png"
)

# ── 8. ExpsPlot ───────────────────────────────────────────────────────────────
take_shot(
  tab_id   = "Hotgenes_A-ExpsPlot",
  filename = "fig3_expsplot_raw.png",
  inputs   = inp("ExpsPlot",
                 yVar = default_feature,
                 xVar = default_coldata_col
  ),
  btn = ns("ExpsPlot", "makePlot")
)
annotate_fig(
  filename        = "fig3_expsplot_raw.png",
  label           = paste0("C  Expression Plot — feature: ", default_feature),
  output_filename = "fig3_expsplot.png"
)

# ── 9. PCA ────────────────────────────────────────────────────────────────────
take_shot(
  tab_id   = "Hotgenes_A-PCA",
  filename = "fig4_pca_raw.png",
  inputs   = inp("PCA",
                 PCA_contrasts = pca_contrasts
  ),
  btn = ns("PCA", "goButton2")
)
annotate_fig(
  filename        = "fig4_pca_raw.png",
  label           = "D  Principal Component Analysis",
  output_filename = "fig4_pca.png"
)

# ── 10. VennDiag ──────────────────────────────────────────────────────────────
take_shot(
  tab_id   = "Hotgenes_A-VennDiag",
  filename = "fig5_venn_raw.png",
  inputs   = inp("VennDiag",
                 Venn_Contrasts = venn_contrasts
  ),
  btn = ns("VennDiag", "execute_btn")
)
annotate_fig(
  filename        = "fig5_venn_raw.png",
  label           = paste0("E  Venn Diagram — ", paste(venn_contrasts, collapse = " vs ")),
  output_filename = "fig5_venn.png"
)

# ── 11. GSEA ──────────────────────────────────────────────────────────────────
take_shot(
  tab_id   = "Hotgenes_A-GSEA",
  filename = "fig6_gsea_raw.png",
  inputs   = inp("GSEA",
                 fgsea_Contrasts = default_contrast,
                 input_MapperCol = all_mapper_cols[1]
  )
)
annotate_fig(
  filename        = "fig6_gsea_raw.png",
  label           = "F  Gene Set Enrichment Analysis (GSEA)",
  output_filename = "fig6_gsea.png"
)

app$stop()

# ── 12. Composite figure: all tabs in one panel ───────────────────────────────
# useful for a single manuscript figure showing the full UI scope
annotated_pngs <- file.path(fig_dir, c(
  "fig1_boxplot.png",
  "fig2_destats.png",
  "fig3_expsplot.png",
  "fig4_pca.png",
  "fig5_venn.png",
  "fig6_gsea.png"
))

# check all exist before compositing
existing <- annotated_pngs[file.exists(annotated_pngs)]

if (length(existing) >= 2) {
  imgs <- lapply(existing, image_read)
  
  # resize all to the same width for clean stacking
  target_width <- 1400 * 2  # matches deviceScaleFactor = 2
  imgs_resized <- lapply(imgs, function(img) {
    image_scale(img, as.character(target_width))
  })
  
  composite <- image_append(do.call(c, imgs_resized), stack = TRUE)
  
  composite_path <- file.path(fig_dir, "fig_composite_all_tabs.png")
  image_write(composite, composite_path, format = "png")
  message("Saved composite: fig_composite_all_tabs.png")
}

# ── 13. Clean up raw (unannotated) files ──────────────────────────────────────
raw_files <- list.files(fig_dir, pattern = "_raw\\.png$", full.names = TRUE)
file.remove(raw_files)

# ── 14. Summary ───────────────────────────────────────────────────────────────
all_files <- list.files(fig_dir, full.names = FALSE)
message("\nAll done! Files saved to: ", fig_dir)
message(paste0("  ", all_files, collapse = "\n"))

browseURL(fig_dir)  # opens in Finder
# unlink(fig_dir, recursive = TRUE)
