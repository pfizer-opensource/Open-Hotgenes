#' Hex logo for Hotgenes
#' @name Hotgenes_logo
NULL

#' path to package logo
#' @export
#' @rdname Hotgenes_logo
#' @example man/examples/Hotgenes_logo_Example.R
Hotgenes_logo_path <-function(){
  
  logo_file_path <- system.file(file.path("logo", "HotgenesLogo_hexSticker.png"),
                                package = "Hotgenes")
  
  #shiny::addResourcePath("www",logo_file_path)
  
  return(logo_file_path)
}



#' Adds logo file to shiny app directory
#' @export
#' @rdname Hotgenes_logo
#' @param logo_dir path to www dir in shiny app working directory.
#' Default is [Hotgenes_logo_path()].
#' @param align string html tag
#' @param height string html tag
#' @param style string html tag
#' @param align string html tag
#' @importFrom shiny h5 tags img
#' @md
embed_Hotgenes_logo <-function(logo_dir = Hotgenes_logo_path(),
                               align = "center", height = "100px",
                               style="text-align: center;margin-bottom:10px;"){
  finaL <- logo_dir %>% basename()
  
  shiny::addResourcePath(finaL, logo_dir)
  
  
  em_logo<- shiny::tags$div(
    shiny::img( src=finaL,
    shiny::h5(glue::glue("v {utils::packageVersion('Hotgenes')}")),
    align = align, height = height),
    style= style)
  
  return(em_logo)
}