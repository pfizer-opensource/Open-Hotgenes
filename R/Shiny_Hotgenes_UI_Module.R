# module as assembled tabPanel --------------------------------------------


#' [tabPanel_Hotgenes()] provides UI for loaded Hotgenes object
#' @export
#' @md
#' @inheritParams shiny::sidebarPanel
#' @inheritParams shiny::tabPanel
#' @rdname Shiny_Hotgenes_UI

tabPanel_Hotgenes <- function(
id = "hotgenes_load",
tabsetPanel_id = "tabset_Hotgenes",
title = "Hotgenes", 
value = title) {

# module start
ns <- shiny::NS(id)


outTab <-shiny::tabPanel(
title = title, 
value = value,

Hotgenes::Shiny_Hotgenes_UI(
id = id,
tabsetPanel_id = tabsetPanel_id)
) # tabPanel end

return(outTab)
}

#' Modular Hotgenes UI
#' @importFrom utils packageVersion
#' @importFrom shiny sidebarPanel
#' @export

Shiny_Hotgenes_UI <- function(
id = NULL,
tabsetPanel_id = "tabs",
width = 3) {
ns <- shiny::NS(id)

# sidePanel_list
sidePanel_params <- shiny::sidebarPanel(
  width = width,
  
  embed_Hotgenes_logo(),
  
  
  # BoxPlot_UI_inputs_support ------------------------------------------------
  
  
  BoxPlot_UI_inputs_support(
    id = "BoxPlot" %>% ns(),
    tabsetPanel_id = tabsetPanel_id
  ),
  
  
  # AuxAssays_sidepanel -----------------------------------------------------
  
  
  AuxAssays_sidepanel( id = "AuxAssays_A" %>% ns(),
                       tabsetPanel_id = tabsetPanel_id
  ),
  
  
  # DEstats_UI_inputs_support -----------------------------------------------
  
  
  DEstats_UI_inputs_support(
    id = "DEstats" %>% ns(),
    tabsetPanel_id = tabsetPanel_id
  ),
  
  # ExpsPlot_UI_inputs_support ----------------------------------------------
  
  
  ExpsPlot_UI_inputs_support(
    id = "ExpsPlot" %>% ns(),
    tabsetPanel_id = tabsetPanel_id
  ),
  
  
  # fgsea_UI_inputs_support -------------------------------------------------
  
  
  fgsea_UI_inputs_support(
    id = "GSEA" %>% ns(),
    tabsetPanel_id = tabsetPanel_id
  ),
  
  
  # VennDiagram_UI_inputs_support -------------------------------------------
  
  
  VennDiagram_UI_inputs_support(
    id = "VennDiag" %>% ns(),
    tabsetPanel_id = tabsetPanel_id
  ),
  
  
  # PCA_UI_inputs_support ---------------------------------------------------
  
  
  PCA_UI_inputs_support(
    id = "PCA" %>% ns(),
    tabsetPanel_id = tabsetPanel_id
  )
)

# mainPanel ---------------------------------------------------------------

mainPanel_params <- shiny::mainPanel(
  shiny::tabsetPanel(
    id = tabsetPanel_id,
    
    # BoxPlot_UI_main_support --------------------------------------------------
    
    
    BoxPlot_UI_main_support(
      id = "BoxPlot" %>% ns(),
      tabsetPanel_id = tabsetPanel_id
    ),
    
    # AuxAssays_mainpanel -----------------------------------------------------
    
    
    AuxAssays_mainpanel(
      id = "AuxAssays_A" %>% ns(),
      tabsetPanel_id = tabsetPanel_id
    ),
    
    # DEstats_UI_main_support -------------------------------------------------
    
    DEstats_UI_main_support(
      id = "DEstats" %>% ns(),
      tabsetPanel_id = tabsetPanel_id
    ),
    # fgsea_UI_main_support ---------------------------------------------------
    
    
    fgsea_UI_main_support(
      id = "GSEA" %>% ns(),
      tabsetPanel_id = tabsetPanel_id
    ),
    
    
    # VennDiagram_UI_main_support ---------------------------------------------
    
    
    VennDiagram_UI_main_support(
      id = "VennDiag" %>% ns(),
      tabsetPanel_id = tabsetPanel_id
    ),
    
    # PCA_UI_main_support -----------------------------------------------------
    
    
    PCA_UI_main_support(
      id = "PCA" %>% ns(),
      tabsetPanel_id = tabsetPanel_id
    ),
    
    # ExpsPlot_UI_main_support ------------------------------------------------
    
    
    ExpsPlot_UI_main_support(
      id = "ExpsPlot" %>% ns(),
      tabsetPanel_id = tabsetPanel_id
    ),
    
    
    # Metadata_UI_main_support ------------------------------------------------
    
    
    Metadata_UI_main_support(
      id = "Metadata" %>% ns(),
      tabsetPanel_id = tabsetPanel_id
    )
  )
)

#shiny::sidebarLayout()
UI_params <-shiny::sidebarLayout(
  
  sidebarPanel = sidePanel_params,
  
  mainPanel = mainPanel_params)


return(UI_params)
}


#' @importFrom shinybusy show_modal_spinner remove_modal_spinner
#' @importFrom shinyjs enable disable useShinyjs
#' @param input shiny server input parameter
#' @param output shiny server output parameter
#' @param session shiny server session parameter
#' @inheritParams BoxPlot
#' @inheritParams shiny::moduleServer
#' @inheritParams fgsea_
#' @inheritParams TabCondition
#' @inheritParams BoxPlot_Server_module
#' @param OntologyMethods output of OntologyMethods function.
#' @param Mapper_choices vector containing names of
#' Mapper_ columns that should be used to map expression data
#' to ontology gene identifiers. Default is names(Mapper_(Hotgenes)).
#' A named vector may be used, as well.
#' @param parallel.sz Number of threads of execution
#' to use when doing the calculations in parallel.
#' @param remove_modal_spinner logical, if TRUE remove_modal_spinner will
#' run at the end of initial server run.
#' @export
#' @example man/examples/Shiny_Hotgenes_Example.R
#' @rdname Shiny_Hotgenes_UI

Shiny_Hotgenes_Server <- function(
id = NULL,
input = NULL,
output = NULL,
session = NULL,
Hotgenes = NULL,
OntologyMethods = Hotgenes::OntologyMethods(),
Mapper_choices = names(Mapper_(Hotgenes)),
parallel.sz = 1L,
max_col_levels = Inf,
selected_fillby = "",
remove_modal_spinner = TRUE) {
shiny::moduleServer(
id = id,
function(input, output, session) {
# Turn off plotting device
if (!is.null(grDevices::dev.list())) {
grDevices::graphics.off()
}

# AuxAssays_server --------------------------------------------------------


shiny::observe({ 

Hotgenes_out <- AuxAssays_server(
id = "AuxAssays_A",
Hotgenes = Hotgenes,
session = session,
input = input,
output = output

) 


# BoxPlot_Server_module ----------------------------------------------------

BoxPlot_output <- BoxPlot_Server_module(
id = "BoxPlot",
input = input,
output = output,
session = session,
Hotgenes = Hotgenes_out,
max_col_levels = max_col_levels,
selected_fillby = selected_fillby
)


# DEstats_Server_module ---------------------------------------------------


DEstats_Server_module(
id = "DEstats",
input = input, output = output,
session = session,
Hotgenes = Hotgenes_out,
ExpressionSlots = BoxPlot_output$reactive_NormSlot,
SampleIDs = BoxPlot_output$SampleIDs_Sel
)

# fgsea_Server_module -----------------------------------------------------



fgsea_Server_module(
id = "GSEA",
input = input, output = output,
session = session,
Hotgenes = Hotgenes_out,
PCA_ranks = out,
ExpressionSlots = BoxPlot_output$reactive_NormSlot,
SampleIDs = BoxPlot_output$SampleIDs_Sel,
OntologyMethods = OntologyMethods,
Mapper_choices = Mapper_choices
)


# ExpsPlot_Server_module --------------------------------------------------


ExpsPlot_Server_module(
id = "ExpsPlot",
input = input, output = output,
session = session,
Hotgenes = Hotgenes_out,
ExpressionSlots = BoxPlot_output$reactive_NormSlot,
SampleIDs = BoxPlot_output$SampleIDs_Sel
)


# Metadata_Server_module --------------------------------------------------

Metadata_Server_module(
id = "Metadata",
input = input, output = output,
session = session,
Hotgenes = Hotgenes_out
)


# VennDiagram_Server_module -----------------------------------------------


VennDiagram_Server_module(
id = "VennDiag",
input = input, output = output,
session = session,
Hotgenes = Hotgenes_out,
ExpressionSlots = BoxPlot_output$reactive_NormSlot,
SampleIDs = BoxPlot_output$SampleIDs_Sel
)

# PCA_Server_module -------------------------------------------------------



out <-  PCA_Server_module(
id = "PCA",
input = input, output = output,
session = session,
Hotgenes = Hotgenes_out,
ExpressionSlots = BoxPlot_output$reactive_NormSlot,
SampleIDs = BoxPlot_output$SampleIDs_Sel
) 

serverOut <- shiny::reactive({
shiny::req(shiny::isTruthy(BoxPlot_output$reactive_NormSlot()))
BoxPlot_output$reactive_NormSlot()
})


# # verifies that everything is loaded
if (isTRUE(remove_modal_spinner)) {
shiny::observeEvent(shiny::req(serverOut() %in% ExpressionSlots_(Hotgenes_out)), {
print("remove_modal_spinner")
shinybusy::remove_modal_spinner(session = session)
})
}



})

}
) # moduleServer end

# adding return value
} # function end


# Generic shiny hotgenes methods ----------------------------------------------
#' R shiny app for visualizing DE analysis
#' @importFrom shiny renderText isTruthy updateSelectInput selectInput
#' updateSelectizeInput selectizeInput numericInput h5 bindEvent
#' observe reactive tabPanel fileInput conditionalPanel HTML
#' req actionButton tags textOutput tagList NS tagList checkboxInput
#' moduleServer verbatimTextOutput fluidRow column icon eventReactive
#' bindCache h3 h4 observeEvent radioButtons renderPlot tabsetPanel
#' updateRadioButtons
#' @importFrom RColorBrewer brewer.pal
#' @importFrom plyr ldply
#' @importFrom grDevices graphics.off dev.list
#' @importFrom openxlsx write.xlsx
#' @importFrom shinyWidgets pickerInput
#' @importFrom tidyr pivot_wider
#' @param Hotgenes R Hotgenes object
#' @rdname Shiny_Hotgenes_UI
#' @return a shiny object for visualization.
#' @export
#' @return a shiny UI and server for visualization.
Shiny_Hotgenes <- function(
Hotgenes = NULL,
OntologyMethods = Hotgenes::OntologyMethods(),
Mapper_choices = names(Mapper_(Hotgenes)),
theme = shinythemes::shinytheme("united"),
parallel.sz = 1L,
max_col_levels = Inf)
UseMethod("Shiny_Hotgenes", Hotgenes)

#' @rdname Shiny_Hotgenes_UI
#' @inheritParams  shiny::fluidPage
#' @importFrom shinythemes shinytheme
#' @export
Shiny_Hotgenes.default <- function(
Hotgenes = NULL,
OntologyMethods = Hotgenes::OntologyMethods(),
Mapper_choices = names(Mapper_(Hotgenes)),
theme = shinythemes::shinytheme("united"),
parallel.sz = 1L,
max_col_levels = Inf) {


# UI ------------------------------------------------------------------

ui <- shiny::fluidPage(
Shiny_Hotgenes_UI(
id = "Hotgenes_A", 
tabsetPanel_id = "tabs"),
theme = theme
)

server <- function(input, output, session) {


# this ensures that BoxPlot loads first   
shiny::observe({
shiny::updateTabsetPanel(
  session = session,
  inputId = "tabs",
  selected = "BoxPlot"
)
})

Shiny_Hotgenes_Server(
id = "Hotgenes_A",
input = input, output = output,
session = session,
Hotgenes = Hotgenes,
OntologyMethods = OntologyMethods,
Mapper_choices = Mapper_choices,
max_col_levels = max_col_levels
)
}

shiny::shinyApp(ui, server)
}



#' @rdname Shiny_Hotgenes_UI
#' @export
Shiny_Hotgenes.list <- function(
Hotgenes = NULL,
OntologyMethods = Hotgenes::OntologyMethods(),
Mapper_choices = names(Mapper_(Hotgenes)),
theme = shinythemes::shinytheme("united"),
parallel.sz = 1L,
max_col_levels = Inf) {

print("processing as list")

HotgenesList <- Hotgenes

# define UI
ui <- shiny::fluidPage(
theme =  theme,
shinyjs::useShinyjs(),

shiny::fluidRow(
shiny::column(
width = 2,
shiny::selectInput(
inputId = "StudyID",
label = "Choose Hotgenes Object", 
selectize = FALSE,
choices = ""),

offset = 0
),

shiny::column(
width = 1,
shiny::tags$div(

# style="padding-left: 20px;padding-top: 20px;align:left",
style="margin-left: 5px;padding-top: 20px",
align = "left",
# style = "padding-top: 20px;width:100%",
shiny::actionButton(
inputId = "loadObjbutton",
label = "Load"))
)

),

shiny::fixedRow(
Shiny_Hotgenes_UI(id = "Obj_A")
)

) # fluidPage close


# Set up server
server <- function(input, output, session) {


# this ensures that BoxPlot loads first   
shiny::observe({
shiny::updateTabsetPanel(
  session = session,
  inputId = "tabs",
  selected = "BoxPlot"
)
})

shiny::observe({

shiny::updateSelectInput(
session = session,
inputId = "StudyID",
selected = character(0),
choices = names(HotgenesList))

})

shiny::observe({

shinyjs::enable("loadObjbutton")

}) %>% 
shiny::bindEvent(shiny::req(input$StudyID))



##
Hotgenes_out <- shiny::eventReactive(
shiny::req(input$loadObjbutton), {
shiny::req(input$StudyID)

shinyjs::disable("loadObjbutton")

shinybusy::show_modal_spinner(
session = session,
spin = "cube-grid",
text = "loading")

print(input$StudyID)
ss <- HotgenesList[[input$StudyID]]



return(ss)

})

# loading object
shiny::observe({

shiny::req(input$loadObjbutton)

Hotgenes::Shiny_Hotgenes_Server(
input = input,
output = output,
session = session,
OntologyMethods = OntologyMethods,
max_col_levels = max_col_levels,
Hotgenes = Hotgenes_out() ,
parallel.sz = parallel.sz,
id = "Obj_A"
)



})



}

# run app
shiny::shinyApp(ui, server)

}