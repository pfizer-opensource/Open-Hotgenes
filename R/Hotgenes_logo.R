#' Hex logo for Hotgenes
#' @name Hotgenes_logo
NULL

#' path to package logo
#' @export
#' @rdname Hotgenes_logo
#' @example examples/Hotgenes_logo_Example.R
Hotgenes_logo_path <-function(){
  
  logo_file_path <- internal_system.file(
  file.path("logo",
  "HotgenesLogo_hexSticker.png"),
  package = "Hotgenes")
  
  return(logo_file_path)
}



#' Embed the Hotgenes logo in a Shiny UI
#'
#' @description Creates a Shiny UI element containing the Hotgenes logo
#' and a clickable version number that opens the changelog modal.
#' Use in combination with \code{\link{logo_news}} in the server function.
#'
#' @param logo_dir Character. Path to the logo directory.
#'   Defaults to \code{\link{Hotgenes_logo_path}()}.
#' @param align Character. Horizontal alignment of the logo.
#'   One of \code{"center"}, \code{"left"}, \code{"right"}. Default is \code{"center"}.
#' @param height Character. Height of the logo image as a CSS string (e.g. \code{"100px"}).
#'   Default is \code{"100px"}.
#' @param style Character. CSS style applied to the wrapping \code{div}.
#'   Default is \code{"text-align: center;margin-bottom:10px;"}.
#'   
#' @param id string passed to actionLink()
#'
#' @return A Shiny \code{tags$div} UI element containing the logo image
#'   and a versioned action link.
#'
#' @examples
#' \dontrun{
#' # In your Shiny UI, e.g. inside sidebarMenu():
#' embed_Hotgenes_logo(height = "100px")
#' }
#'
#' @seealso \code{\link{logo_news}}
#'
#' @export
embed_Hotgenes_logo <- function (
  logo_dir = Hotgenes_logo_path(), 
  id = "show_news",
  align = "center", 
  height = "100px", 
  style = "text-align: center;margin-bottom:10px;") {
  
  finaL <- logo_dir %>% basename()
  shiny::addResourcePath(finaL, logo_dir)
  
  em_logo <- shiny::tags$div(
    shiny::img(
      src = finaL,
      align = align, 
      height = height
    ),
    shiny::actionLink(
      inputId = id,
      label = shiny::h5(glue::glue("v {utils::packageVersion('Hotgenes')} news"),
                        style = "color: black;")
    ),
    style = style
  )
  
  
  return(em_logo)
}

#' Internal helper to locate package files during development and installation
#'
#' @description A wrapper around \code{base::system.file} that supports
#' locating files in both installed packages and packages loaded via
#' \code{devtools::load_all()}. When a package is loaded with
#' \code{devtools::load_all()}, files may reside under \code{inst/} rather
#' than at the top level, and this function checks both locations.
#'
#' @param ... Character vectors of path components passed to
#'   \code{base::system.file}.
#' @param package Character. The name of the package to search. Default is
#'   \code{"base"}.
#' @param lib.loc Character vector or \code{NULL}. Location of the package
#'   library. Passed to \code{base::system.file}. Default is \code{NULL}.
#' @param mustWork Logical. If \code{TRUE}, an error is thrown if the file
#'   cannot be found. Default is \code{FALSE}.
#'
#' @return A character string of the normalized file path if found, or
#'   \code{""} if not found and \code{mustWork = FALSE}.
#'
#' @examples
#' \dontrun{
#' internal_system.file("logo", "logo.png", package = "Hotgenes")
#' }
#'
#' @importFrom cli cli_abort
#'
#' @noRd
internal_system.file <- function (..., 
                                  package = "base", 
                                  lib.loc = NULL, 
                                  mustWork = FALSE) {
  
  if (any(!requireNamespace("devtools", quietly = TRUE))) {
    return(base::system.file(..., package = package, lib.loc = lib.loc, 
                             mustWork = mustWork))
  } 
  
  if (!(package %in% devtools::dev_packages())) {
    return(base::system.file(..., package = package, lib.loc = lib.loc, 
                             mustWork = mustWork))
  }
  
  if (rlang::dots_n(...) && rlang::is_string(..1)) {
    if (rlang::is_string("inst", ..1) || grepl("^inst/", ..1)) {
      cli::cli_abort(c("Paths can't start with `inst`", 
                       i = "Files in `inst` are installed at top-level."))
    }
  }
  pkg_path <- find.package(package)
  files_inst <- file.path(pkg_path, "inst", ...)
  present_inst <- file.exists(files_inst)
  files_top <- file.path(pkg_path, ...)
  present_top <- file.exists(files_top)
  files <- files_top
  files[present_inst] <- files_inst[present_inst]
  files <- files[present_inst | present_top]
  if (length(files) > 0) {
    normalizePath(files, winslash = "/")
  } else {
    if (mustWork) {
      cli::cli_abort("Can't find package file.", call = NULL)
    }
    ""
  }
}

#' Server-side handler for the Hotgenes changelog modal
#'
#' @description Registers a Shiny \code{observeEvent} that listens for clicks
#' on the version link created by \code{\link{embed_Hotgenes_logo}} and displays
#' the package changelog in a modal dialog. Supports both plain text
#' \code{NEWS} files and markdown \code{NEWS.md} files — the appropriate
#' rendering function is selected automatically based on the file extension.
#'
#' @param input The Shiny \code{input} object from the server function.
#' @param session The Shiny \code{session} object. Defaults to the current
#'   reactive domain via \code{shiny::getDefaultReactiveDomain()}.
#' @param news_file file name for news file.
#'
#' @return None. Called for its side effects (registers an \code{observeEvent}).
#'
#' @details
#' The function first looks for a \code{NEWS.md} file in the package root.
#' If not found, it falls back to a plain text \code{NEWS} file.
#' \code{NEWS.md} is rendered with \code{shiny::includeMarkdown()} and
#' plain \code{NEWS} is rendered with \code{shiny::includeText()}.
#'
#' @examples
#' \dontrun{
#' # In your Shiny server function:
#' Hotgenes_core_server <- function(input, output, session) {
#'   logo_news(input, session)
#' }
#' }
#'
#' @seealso \code{\link{embed_Hotgenes_logo}}
#'
#' @export
logo_news <- function(input, 
                      session = shiny::getDefaultReactiveDomain(),
                      news_file = c("NEWS", "NEWS.md")) {
  
  news_file <- match.arg(arg = news_file, 
                         choices = c("NEWS", "NEWS.md"),
                         several.ok = FALSE)
  
  shiny::observeEvent(input$show_news, {
    
    news_path <- internal_system.file(news_file, package = "Hotgenes")
    
    news_content <- if (!nzchar(news_path)) {
      shiny::p("No changelog found.")
    } else if (grepl("\\.md$", news_file)) {
      shiny::includeMarkdown(news_path)
    } else {
      shiny::tags$pre(
        style = "white-space: pre-wrap; word-wrap: break-word; font-size: 13px;",
        paste(readLines(news_path), collapse = "\n")
      )
    }
    
    shiny::showModal(shiny::modalDialog(
  title = glue::glue("Hotgenes v{utils::packageVersion('Hotgenes')} - Changelog"),
      news_content,
      easyClose = TRUE,
      size = "l",
      footer = shiny::modalButton("Close")
    ))
  })
}
