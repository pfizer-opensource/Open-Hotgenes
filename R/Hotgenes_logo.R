#' Hex logo for Hotgenes
#' @name Hotgenes_logo
NULL

#' path to package logo
#' @export
#' @rdname Hotgenes_logo
#' @example examples/Hotgenes_logo_Example.R
Hotgenes_logo_path <-function(){
  
  logo_file_path <- system.file(
  file.path("logo",
  "HotgenesLogo_hexSticker.png"),
  package = "Hotgenes", mustWork = TRUE)
  
  return(logo_file_path)
}



#' Embed the Hotgenes logo in a Shiny UI
#'
#' @description Creates a Shiny UI element containing the Hotgenes logo
#' and a clickable version number that opens the changelog modal.
#' Use in combination with \code{\link{logo_news}} in the server function.
#'
#' @param logo_dir Character. Path to the logo image file.
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
  
  logo_name <- basename(logo_dir)
  shiny::addResourcePath(logo_name, dirname(logo_dir))
  
  em_logo <- shiny::tags$div(
    shiny::img(
      src = paste(logo_name, logo_name, sep = "/"),
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
    
    news_path <- system.file(news_file, package = "Hotgenes")
    
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
