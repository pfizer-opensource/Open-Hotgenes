require(Hotgenes)

# load example data
dds_Hotgenes_dir <- system.file("extdata",
paste0("dds_Hotgenes", ".RDS"),
package = "Hotgenes",
mustWork = TRUE)

hotobj <- readRDS(dds_Hotgenes_dir)

# The example below works for all the Shiny_Tabs modules triplet functions
if (FALSE) {
  
options(device.ask.default = FALSE)

# Set UI
ui <- shiny::fluidPage(
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      
      # side panel module
      BoxPlot_UI_inputs_support() ),
    shiny::mainPanel( 
      shiny::tabsetPanel( 
        
        # main panel module
        BoxPlot_UI_main_support()


) ) )) 

# Set up server
server <- function(input, output, session) {
  
 BoxPlot_Server_module(
  Hotgenes = hotobj)


}

# run app

shiny::shinyApp(ui, server)
}
