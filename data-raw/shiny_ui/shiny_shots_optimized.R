#library(Hotgenes)
library(shinytest2)
#library(magick)
#library(cli)

if(FALSE){
  # Debug block for interactive development
  # Load HotgenesObj and test parameters
  dds_Hotgenes_dir <- system.file("extdata",
                                  paste0("dds_Hotgenes", ".RDS"),
                                  package = "Hotgenes",
                                  mustWork = TRUE)
  HotgenesObj <- readRDS(dds_Hotgenes_dir)
  fig_dir <- file.path(tempdir(), "hotgenes_figures")
}

# load_hotgenes_object -----------------------------------------------------------

dds_Hotgenes_dir <- system.file("extdata",
                                paste0("dds_Hotgenes", ".RDS"),
                                package = "Hotgenes",
                                mustWork = TRUE)

HotgenesObj <- readRDS(dds_Hotgenes_dir) %>%
  Hotgenes::update_object()

cli::cli_h1("Loaded Hotgenes object")
cli::cli_inform(c("i" = "Object class: {.cls {class(HotgenesObj)}}"))

# extract_dynamic_values -----------------------------------------------------------

all_contrasts <- contrasts_(HotgenesObj)
all_expr_slots <- ExpressionSlots_(HotgenesObj)
all_coldata_cols <- coldata_names(HotgenesObj)
all_features <- Background_(HotgenesObj, Col = "Feature")
all_mapper_cols <- names(Mapper_(HotgenesObj))

default_contrast <- all_contrasts[1]
default_expr_slot <- all_expr_slots[1]
default_coldata_col <- all_coldata_cols[1]
default_feature <- all_features[1]
venn_contrasts <- all_contrasts[seq_len(min(2, length(all_contrasts)))]
pca_contrasts <- all_contrasts[seq_len(min(3, length(all_contrasts)))]

cli::cli_h2("Extracted dynamic values")
cli::cli_bullets(c(
  "Contrasts:" = paste(all_contrasts, collapse = ", "),
  "Expression slots:" = paste(all_expr_slots, collapse = ", "),
  "Coldata cols:" = paste(all_coldata_cols, collapse = ", "),
  "Mapper cols:" = paste(all_mapper_cols, collapse = ", "),
  "Default contrast:" = default_contrast,
  "Default expr slot:" = default_expr_slot,
  "Default coldata col:" = default_coldata_col
))

# helper_ns_id -----------------------------------------------------------

#' Construct namespaced input ID
#'
#' @param tab Character. Tab name.
#' @param input_id Character. Input identifier.
#'
#' @return Character string with namespace pattern.
#'
#' @keywords internal
ns_id <- function(tab = NULL, input_id = NULL) {
  if (is.null(tab) || is.null(input_id)) {
    cli::cli_abort(c(
      "Arguments {.arg tab} and {.arg input_id} are required.",
      "x" = "Both must be non-null character strings."
    ))
  }
  return(paste("Hotgenes_A", tab, input_id, sep = "-"))
}

# helper_construct_inputs -----------------------------------------------------------

#' Construct input list for app$set_inputs
#'
#' @param tab Character. Tab name.
#' @param ... Named arguments for inputs.
#'
#' @return Named list suitable for app$set_inputs.
#'
#' @keywords internal
construct_inputs <- function(tab = NULL, ...) {
  if (is.null(tab)) {
    cli::cli_abort("Argument {.arg tab} is required and cannot be NULL.")
  }

  args <- list(...)
  if (length(args) == 0) {
    return(list())
  }

  return(setNames(
    args,
    paste("Hotgenes_A", tab, names(args), sep = "-")
  ))
}

# helper_wait_for_render -----------------------------------------------------------

#' Wait for Shiny UI to stabilize and render
#'
#' @param app AppDriver object from shinytest2.
#' @param wait_base Numeric. Base wait time in seconds (default: 1.5).
#' @param stability_threshold Numeric. Wait time to confirm stability (default: 1.0).
#' @param max_timeout Numeric. Maximum total wait time in seconds (default: 15).
#'
#' @return Invisible NULL. Waits until rendering is complete.
#'
#' @keywords internal
wait_for_render <- function(app = NULL,
                             wait_base = 1.5,
                             stability_threshold = 1.0,
                             max_timeout = 15) {
  if (is.null(app)) {
    cli::cli_abort("Argument {.arg app} is required and cannot be NULL.")
  }

  if (!inherits(app, "AppDriver")) {
    cli::cli_abort("Argument {.arg app} must be an {.cls AppDriver} object.")
  }

  start_time <- Sys.time()

  Sys.sleep(wait_base)

  elapsed <- difftime(Sys.time(), start_time, units = "secs")
  if (elapsed < max_timeout) {
    Sys.sleep(stability_threshold)
  }

  return(invisible(NULL))
}

# helper_take_screenshot -----------------------------------------------------------

#' Capture a single screenshot from the Shiny app
#'
#' @param app AppDriver object from shinytest2.
#' @param tab_id Character. Tab identifier to activate.
#' @param filename Character. Output PNG filename.
#' @param inputs List. Named inputs to set via app$set_inputs (default: empty list).
#' @param button_id Character. Optional button ID to click after setting inputs.
#' @param wait_time Numeric. Wait time between actions in seconds (default: 1.5).
#' @param output_dir Character. Directory to save screenshot (default: tempdir()).
#' @param max_retries Integer. Number of retry attempts on failure (default: 3).
#'
#' @return Invisible logical. TRUE if successful, FALSE otherwise.
#'
#' @keywords internal
take_screenshot <- function(app = NULL,
                             tab_id = NULL,
                             filename = NULL,
                             inputs = list(),
                             button_id = NULL,
                             wait_time = 1.5,
                             output_dir = NULL,
                             max_retries = 3) {
  if (is.null(app) || is.null(tab_id) || is.null(filename)) {
    cli::cli_abort(c(
      "Arguments {.arg app}, {.arg tab_id}, and {.arg filename} are required.",
      "x" = "All must be non-null."
    ))
  }

  if (is.null(output_dir)) {
    output_dir <- tempdir()
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }

  for (attempt in seq_len(max_retries)) {
    tryCatch({
      app$set_inputs(tabs = tab_id, wait_ = FALSE)
      wait_for_render(app, wait_base = wait_time)

      if (length(inputs) > 0) {
        do.call(app$set_inputs, c(inputs, list(wait_ = FALSE)))
        wait_for_render(app, wait_base = wait_time)
      }

      if (!is.null(button_id)) {
        app$click(button_id)
        wait_for_render(app, wait_base = wait_time * 2)
      }

      output_path <- file.path(output_dir, filename)
      app$get_screenshot(output_path)

      if (!file.exists(output_path)) {
        cli::cli_warn(c(
          "Screenshot file not created at {.file {output_path}}.",
          "i" = "Attempt {attempt}/{max_retries}"
        ))
        return(invisible(FALSE))
      }

      file_size <- file.size(output_path) / 1024
      cli::cli_alert_success(c(
        "Saved {.file {filename}} ({format(file_size, digits = 1)} KB)"
      ))

      return(invisible(TRUE))
    }, error = function(e) {
      cli::cli_alert_warning(c(
        "Attempt {attempt}/{max_retries} failed for {.file {filename}}.",
        "i" = conditionMessage(e)
      ))
      return(invisible(FALSE))
    })
  }

  cli::cli_warn("Failed to capture {.file {filename}} after {max_retries} retries.")
  return(invisible(FALSE))
}

# helper_annotate_screenshot -----------------------------------------------------------

#' Add annotation banner to screenshot
#'
#' @param input_path Character. Path to input PNG file.
#' @param label Character. Label text to prepend.
#' @param output_path Character. Path for output annotated PNG.
#' @param label_size Numeric. Font size in points (default: 40).
#' @param bg_color Character. Banner background color (default: "white").
#' @param text_color Character. Text color (default: "black").
#'
#' @return Invisible logical. TRUE if successful, FALSE otherwise.
#'
#' @keywords internal
annotate_screenshot <- function(input_path = NULL,
                                 label = NULL,
                                 output_path = NULL,
                                 label_size = 40,
                                 bg_color = "white",
                                 text_color = "black") {
  if (is.null(input_path) || is.null(label) || is.null(output_path)) {
    cli::cli_abort(c(
      "Arguments {.arg input_path}, {.arg label}, and {.arg output_path} are required.",
      "x" = "All must be non-null."
    ))
  }

  if (!file.exists(input_path)) {
    cli::cli_warn("Input file does not exist: {.file {input_path}}")
    return(invisible(FALSE))
  }

  tryCatch({
    img <- magick::image_read(input_path)
    info <- magick::image_info(img)
    width <- info$width

    banner <- magick::image_blank(width = width, height = 80, color = bg_color) %>%
      magick::image_annotate(
        text = label,
        gravity = "West",
        location = "+20+0",
        size = label_size,
        color = text_color,
        font = "Helvetica",
        weight = 700
      )

    magick::image_append(c(banner, img), stack = TRUE) %>%
      magick::image_write(output_path, format = "png")

    cli::cli_alert_success("Annotated: {.file {basename(output_path)}}")
    return(invisible(TRUE))
  }, error = function(e) {
    cli::cli_warn(c(
      "Failed to annotate {.file {basename(input_path)}}.",
      "i" = conditionMessage(e)
    ))
    return(invisible(FALSE))
  })
}

# helper_create_composite -----------------------------------------------------------

#' Combine multiple annotated screenshots into composite figure
#'
#' @param png_files Character vector. Paths to PNG files to combine.
#' @param output_path Character. Path for output composite PNG.
#' @param target_width Numeric. Target width for scaling (default: NULL, no scaling).
#'
#' @return Invisible logical. TRUE if successful, FALSE otherwise.
#'
#' @keywords internal
create_composite <- function(png_files = NULL,
                              output_path = NULL,
                              target_width = NULL) {
  if (is.null(png_files) || is.null(output_path)) {
    cli::cli_abort(c(
      "Arguments {.arg png_files} and {.arg output_path} are required.",
      "x" = "Both must be non-null."
    ))
  }

  existing_files <- png_files[file.exists(png_files)]

  if (length(existing_files) < 2) {
    cli::cli_warn(c(
      "Need at least 2 PNG files to create composite.",
      "i" = "Found {length(existing_files)} file(s)."
    ))
    return(invisible(FALSE))
  }

  tryCatch({
    imgs <- lapply(existing_files, magick::image_read)

    if (!is.null(target_width)) {
      imgs <- lapply(imgs, function(img) {
        magick::image_scale(img, as.character(target_width))
      })
    }

    composite <- magick::image_append(do.call(c, imgs), stack = TRUE)
    magick::image_write(composite, output_path, format = "png")

    composite_size <- file.size(output_path) / (1024^2)
    cli::cli_alert_success(c(
      "Saved composite: {.file {basename(output_path)}} ({format(composite_size, digits = 2)} MB)"
    ))

    return(invisible(TRUE))
  }, error = function(e) {
    cli::cli_warn(c(
      "Failed to create composite figure.",
      "i" = conditionMessage(e)
    ))
    return(invisible(FALSE))
  })
}

# setup_output_directory -----------------------------------------------------------

fig_dir <- file.path(tempdir(), "hotgenes_figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cli::cli_h2("Output directory")
cli::cli_inform("Figures will be saved to: {.file {fig_dir}}")

# launch_app -----------------------------------------------------------

cli::cli_h2("Launching Shiny app")

app <- AppDriver$new(
  app = Shiny_Hotgenes(HotgenesObj),
  height = 1000,
  width = 1600
)

cli::cli_alert_success("App launched successfully")

# set_retina_quality -----------------------------------------------------------

cli::cli_h2("Configuring device pixel ratio (retina 2x)")

app$get_chromote_session()$Emulation$setDeviceMetricsOverride(
  width = 1600,
  height = 1000,
  deviceScaleFactor = 2,
  mobile = FALSE
)

cli::cli_alert_success("Device pixel ratio set to 2x (retina quality)")

# capture_boxplot_tab -----------------------------------------------------------

cli::cli_h2("Capturing BoxPlot tab")

take_screenshot(
  app = app,
  tab_id = "Hotgenes_A-BoxPlot",
  filename = "shiny-01-boxplot_raw.png",
  inputs = construct_inputs("BoxPlot", NormSlot = default_expr_slot),
  output_dir = fig_dir,
  wait_time = 1.5
)

annotate_screenshot(
  input_path = file.path(fig_dir, "shiny-01-boxplot_raw.png"),
  label = "A  Expression Data (BoxPlot) — normalised per-sample distributions",
  output_path = file.path(fig_dir, "shiny-01-boxplot.png")
)

# capture_destats_tab -----------------------------------------------------------

cli::cli_h2("Capturing DEstats tab")

take_screenshot(
  app = app,
  tab_id = "Hotgenes_A-DEstats",
  filename = "shiny-02-destats_raw.png",
  inputs = construct_inputs("DEstats",
                             DE_Contrasts = default_contrast,
                             padj_cut_DE_tables = 0.1,
                             lfc_Statistics = 0),
  output_dir = fig_dir,
  wait_time = 1.5
)

annotate_screenshot(
  input_path = file.path(fig_dir, "shiny-02-destats_raw.png"),
  label = paste0("B  Differential Expression — contrast: ", default_contrast),
  output_path = file.path(fig_dir, "shiny-02-destats.png")
)

# capture_expsplot_tab -----------------------------------------------------------

cli::cli_h2("Capturing ExpsPlot tab")

take_screenshot(
  app = app,
  tab_id = "Hotgenes_A-ExpsPlot",
  filename = "shiny-03-expsplot_raw.png",
  inputs = construct_inputs("ExpsPlot",
                             yVar = default_feature,
                             xVar = default_coldata_col),
  button_id = ns_id("ExpsPlot", "makePlot"),
  output_dir = fig_dir,
  wait_time = 1.5
)

annotate_screenshot(
  input_path = file.path(fig_dir, "shiny-03-expsplot_raw.png"),
  label = paste0("C  Expression Plot — feature: ", default_feature),
  output_path = file.path(fig_dir, "shiny-03-expsplot.png")
)

# capture_pca_tab -----------------------------------------------------------

cli::cli_h2("Capturing PCA tab")

take_screenshot(
  app = app,
  tab_id = "Hotgenes_A-PCA",
  filename = "shiny-04-pca_raw.png",
  inputs = construct_inputs("PCA", PCA_contrasts = pca_contrasts),
  button_id = ns_id("PCA", "goButton2"),
  output_dir = fig_dir,
  wait_time = 1.5
)

annotate_screenshot(
  input_path = file.path(fig_dir, "shiny-04-pca_raw.png"),
  label = "D  Principal Component Analysis",
  output_path = file.path(fig_dir, "shiny-04-pca.png")
)

# capture_venn_tab -----------------------------------------------------------

cli::cli_h2("Capturing VennDiag tab")

take_screenshot(
  app = app,
  tab_id = "Hotgenes_A-VennDiag",
  filename = "shiny-05-venn_raw.png",
  inputs = construct_inputs("VennDiag", Venn_Contrasts = venn_contrasts),
  button_id = ns_id("VennDiag", "execute_btn"),
  output_dir = fig_dir,
  wait_time = 1.5
)

annotate_screenshot(
  input_path = file.path(fig_dir, "shiny-05-venn_raw.png"),
  label = paste0("E  Venn Diagram — ", paste(venn_contrasts, collapse = " vs ")),
  output_path = file.path(fig_dir, "shiny-05-venn.png")
)

# capture_gsea_tab -----------------------------------------------------------

cli::cli_h2("Capturing GSEA tab")

take_screenshot(
  app = app,
  tab_id = "Hotgenes_A-GSEA",
  filename = "shiny-06-gsea_raw.png",
  inputs = construct_inputs("GSEA",
                             fgsea_Contrasts = default_contrast,
                             input_MapperCol = all_mapper_cols[1]),
  output_dir = fig_dir,
  wait_time = 1.5
)

annotate_screenshot(
  input_path = file.path(fig_dir, "shiny-06-gsea_raw.png"),
  label = "F  Gene Set Enrichment Analysis (GSEA)",
  output_path = file.path(fig_dir, "shiny-06-gsea.png")
)

# stop_app -----------------------------------------------------------

app$stop()
cli::cli_alert_success("App stopped successfully")

# create_composite_figure -----------------------------------------------------------

cli::cli_h2("Creating composite figure")

annotated_pngs <- file.path(fig_dir, c(
  "shiny-01-boxplot.png",
  "shiny-02-destats.png",
  "shiny-03-expsplot.png",
  "shiny-04-pca.png",
  "shiny-05-venn.png",
  "shiny-06-gsea.png"
))

create_composite(
  png_files = annotated_pngs,
  output_path = file.path(fig_dir, "shiny_composite_all_tabs.png"),
  target_width = 1600 * 2
)

# cleanup_raw_files -----------------------------------------------------------

cli::cli_h2("Cleaning up raw (unannotated) files")

raw_files <- list.files(fig_dir, pattern = "_raw\\.png$", full.names = TRUE)
n_removed <- length(raw_files)

if (n_removed > 0) {
  file.remove(raw_files)
  cli::cli_alert_success("Removed {n_removed} raw screenshot file(s)")
} else {
  cli::cli_inform("No raw files found to clean up")
}

# summary_and_output -----------------------------------------------------------

cli::cli_h1("All done! Vignette figures ready")

all_files <- list.files(fig_dir, full.names = FALSE)
n_files <- length(all_files)

cli::cli_h2("Generated files ({n_files} total)")

file_info <- lapply(all_files, function(f) {
  size_kb <- file.size(file.path(fig_dir, f)) / 1024
  paste0("  • {.file ", f, "} ({format(size_kb, digits = 0)} KB)")
})

cli::cli_bullets(unlist(file_info))

cli::cli_h2("Output directory")
cli::cli_inform("Path: {.file {fig_dir}}")

browseURL(fig_dir)
