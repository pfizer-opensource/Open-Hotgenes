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
#' @param wait_element a check for a specific UI element to be rendered
#' @param output_dir Character. Directory to save screenshot (default: tempdir()).
#' @param output_name name of shiny object to wait for.
#' @param timeout numeric for time in seconds for time out of output_name
#' @param max_retries Integer. Number of retry attempts on failure (default: 3).
#' @param overwrite logical if TRUE file will be overwritten
#' @param button_timeout numeric for time to wait for buttons.
#' @return Invisible logical. TRUE if successful, FALSE otherwise.
#'
#' @keywords internal
take_screenshot <- function(app = NULL,
                            tab_id = NULL,
                            filename = NULL,
                            inputs = list(),
                            button_id = NULL,
                            button_timeout = 20000,
                            output_name = NULL,  # NEW
                            timeout = NULL,
                            wait_element = NULL,  # NEW
                            wait_time = 1.5,
                            output_dir = NULL,
                            overwrite = TRUE,
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
  # attempt <- 1
  for (attempt in seq_len(max_retries)) {
    tryCatch({
      
      app$set_inputs(tabs = tab_id, wait_ = FALSE)
     
      #  #wait_for_render(app, wait_base = wait_time)
      # if(!is.null(value_check)) {
      #   
      #   app$wait_for_value()
      #   
      #   app$get_value(output = "total")
      #   
      # }
      
      # Wait for specific output if provided
      if (!is.null(output_name)) {
        app$wait_for_js(
          glue::glue("Shiny.shinyapp.$values.output.{output_name} !== null"),
          timeout = timeout
        )
      }
      
      # Wait for specific element instead of fixed time
      if (!is.null(wait_element)) {
        app$wait_for_js(glue::glue("$('{wait_element}').is(':visible')"))
      } else {
        wait_for_render(app, wait_base = wait_time)
      }
      
      if (length(inputs) > 0) {
        do.call(app$set_inputs, c(inputs, list(wait_ = FALSE)))
        wait_for_render(app, wait_base = wait_time)
      }
      
      if (!is.null(button_id)) {
        app$click(button_id, timeout_ = button_timeout)
        wait_for_render(app, wait_base = wait_time * 2)
      }
      
      output_path <- file.path(output_dir, filename)
      
      if(overwrite) {
        
        if(file.exists(output_path)) file.remove(output_path)
        
      }
      
      app$get_screenshot(output_path)
      
      if (!file.exists(output_path)) {
        cli::cli_warn(c(
          "Screenshot file not created at {.file {output_path}}.",
          "i" = "Attempt {attempt}/{max_retries}"
        ))
       #return(invisible(FALSE))
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




