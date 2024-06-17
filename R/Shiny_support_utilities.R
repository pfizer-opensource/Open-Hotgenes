#' Provides the condition for when a specific tab is active
#' @export
#' @param tabsetPanel_id string for tabsetPanel id prefix to use.
#' Default is "tabs". This is what controls the visibility of sidepanels.
#' @inheritParams shiny::moduleServer
#' @example man/examples/TabCondition_Example.R
#' @details
#' When deploying an app with multiple Hotgenes UIs,
#' don't 
#' 
#' @return string condition
TabCondition <- function(
    tabsetPanel_id = "tabs",
    id = NULL) {
  #out <- paste0("input.", tabsetPanel_id, ' =="', id, '"')
  out <- glue::glue("input.{tabsetPanel_id} ==\ '{id}\'")
 # paste0("input.", "tabsetPanel_id", ' =="', "id", '"')
  print(out)
  return(out)
}


#' shiny download handler as a module
#' @export
#' @inheritParams shiny::moduleServer
#' @inheritParams shiny::downloadButton
#' @param ... passed to [shiny::downloadButton()]
#' @md

UI_download_file <- function(id = NULL, 
                             label = "download",
                             ...) {
  ns <- shiny::NS(id)
  
  shiny::downloadButton(
    outputId = "DL_file" %>% ns(),
    label = label
  )
}
#' shiny plot download handler
#' @export
#' @inheritParams shiny::moduleServer
#' @inheritParams shiny::downloadHandler
#' @inheritParams Shiny_Hotgenes_Server
#' @param ... additional arguments for [shiny::downloadHandler()]
#' @rdname UI_download_file
#' @md

Server_download_file <- function(
    id = NULL,
    filename = NULL,
    content = NULL,
    input, output, session,
    ...) {
  
  
  shiny::moduleServer(id, {
    function(input, output, session) {
      output$DL_file <- shiny::downloadHandler(
        filename = filename,
        content = content,
        ...
      )
    }
  })
}


#' shiny plot download handler
#' @export
#' @inheritParams shiny::moduleServer
#' @inheritParams shiny::downloadButton
#' @rdname UI_interactive_plot

UI_download_plot <- function(id = NULL, label = "PDF") {
  ns <- shiny::NS(id)

  shiny::downloadButton(
    outputId = "DL_Plot" %>% ns(),
    label = label
  )
}
#' shiny plot download handler
#' @export
#' @inheritParams shiny::moduleServer
#' @inheritParams ggplot2::ggsave
#' @param filename string for file name without extension.
#' Default is"Figure"
#' @param extension string for file extension. Default is".png"
#' @param ... additional arguments for \code{\link[ggplot2]{ggsave}}
#' @rdname UI_interactive_plot

Server_download_plot <- function(
    id = NULL,
    reactive_plot = NULL,
    filename = "Figure",
    extension = ".pdf",
    units = "in",
    input, output, session,
    width = 10, height = 8, ...) {
  Final_name <- paste0(filename, extension)

  shiny::moduleServer(id, {
    function(input, output, session) {
      output$DL_Plot <- shiny::downloadHandler(
        filename = Final_name,
        content = function(file) {
         
          ggplot2::ggsave(
            filename = file,
            units = units,
            plot = reactive_plot(),
            width = width,
            height = height, ...
          )
        }
      )
    }
  })
}



#' interactive plot modules
#' @export
#' @inheritParams shiny::moduleServer
#' @inheritParams shiny::plotOutput
#' @param ... additional arguments for \code{\link[shiny]{plotOutput}}
#' @rdname UI_interactive_plot


UI_interactive_plot <- function(
    id,
    width = "100%", height = "400px",
    ...) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::plotOutput(
      outputId = "PlotOut" %>% ns(),
      width = width, height = height,
      click = "clickOut" %>% ns(),
      brush = shiny::brushOpts(
        id = "brushOut" %>% ns(),
        ...
      )
    ),
    shiny::fluidRow(
      shiny::column(
        width = 6,
        shiny::h4("Points near click"),
        shiny::verbatimTextOutput("click_info_plot" %>% ns())
      ),
      shiny::column(
        width = 6,
        shiny::h4("Brushed points"),
        shiny::verbatimTextOutput("brush_info_plot" %>% ns())
      )
    )
  )
}
#' interactive plot modules
#' @export
#' @inheritParams Shiny_Hotgenes_Server
#' @inheritParams shiny::moduleServer
#' @inheritParams shiny::plotOutput
#' @inheritParams shiny::renderPlot
#' @inheritParams shiny::brushedPoints
#' @param reactive_plot a reactive ggplot object.
#' @rdname UI_interactive_plot
#' @param refine_cols string indicating which columns
#' to show. If NULL (default), all cols will be shown.
Server_interactive_plot <- function(
    id,
    input, output, session,
    refine_cols = NULL,
    reactive_plot = NULL,
    xvar = reactive(NULL),
    yvar = reactive(NULL),
    execOnResize = FALSE) {
  shiny::moduleServer(id, {
    function(input, output, session) {
      output$PlotOut <- shiny::renderPlot(
        {
          reactive_plot()
        },
        execOnResize = execOnResize
      )

      df <- shiny::reactive({
        
        if(is.null(refine_cols)){
          df1 <- reactive_plot()$data
        } else {
          
          df1 <- reactive_plot()$data %>% 
            dplyr::select(dplyr::any_of(refine_cols))
        }
        return(df1)
      })


      ## points
      # click out
      output$click_info_plot <- shiny::renderPrint({
        # req(sel_plotData())
        # cells.located %>% print
        
        shiny::nearPoints(df = df(),
          coordinfo = input$clickOut,
          xvar = xvar(),
          yvar = yvar()
        )
      })


      # brush out
      output$brush_info_plot <- shiny::renderPrint({
        shiny::brushedPoints(df = df(),
          brush = input$brushOut,
          xvar = xvar(),
          yvar = yvar()
        )
      })
    }
  }) # end moduleServer
} # end function


# internal plot
#' @noRd
empty_plot <- function(
    message_ = paste("\n No data. Check inputs")) {
  
  ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::ggtitle(label = message_) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
    
}



#' @noMd
df_coldata_handler <- function(Hotgenes = NULL,
                               max_col_levels = Inf){
  
  if(is.null(max_col_levels)){
    stop("max_col_levels must be provided")
  }
  
  DF_coldata <- shiny::eventReactive(shiny::req(Hotgenes), {
   # t1 <- Sys.time()
    
    quali_distinct_raw <- coldata_(Hotgenes, mode = "quali") %>%
      dplyr::select_if(~ dplyr::n_distinct(.x,
                                           na.rm = TRUE
      ) > 1) %>%
      dplyr::select_if(~ dplyr::n_distinct(.x,
                                           na.rm = TRUE
      ) <= max_col_levels)
    
    sort_cols <- quali_distinct_raw %>%
      purrr::imap_int(~ .x %>% dplyr::n_distinct(na.rm = TRUE)) %>%
      sort() %>%
      names()
    
    #
    
    #
    
    
    Outlist <- list(
      quali = coldata_(Hotgenes, mode = "quali"),
      quanti = coldata_(Hotgenes, mode = "quanti"),
      all = coldata_(Hotgenes, mode = "all"),
      quali_distinct = quali_distinct_raw %>%
        dplyr::select(dplyr::any_of(sort_cols))
    )  %>% purrr::compact()
    
    Features_aliases <- Mapper_list(Hotgenes = Hotgenes) 
    
    Outlist$all_features <- list(
                 auxiliary_assays = auxiliary_assays_(Hotgenes,
              mode = "quanti") %>% names()) %>%
      append(Features_aliases) %>% 
      purrr::compact()
    
  #  t2 <- Sys.time()
    
 #   glue::glue("DF_coldata {signif(t2-t1)}") %>%
  #    print()
    
    return(Outlist)
  })
  
  return(DF_coldata)
}