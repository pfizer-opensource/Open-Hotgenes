#' Modules for Shiny_Tabs  
#' @name Shiny_Tabs
#' @description
#' All the tabs in [Shiny_Hotgenes()] are written as sets of functions.
#' 
#' Auxiliary Assays: 
#' [AuxAssays_sidepanel()], 
#' [AuxAssays_mainpanel()], and [AuxAssays_server()]
#' 
#' PCA tab:
#' [PCA_UI_inputs_support()], 
#' [PCA_UI_main_support()], and [PCA_Server_module()]
#' 
#' Metadata tab: 
#' [Metadata_UI_main_support()] and [Metadata_Server_module()]
#' 
#' GSEA tab: 
#' [fgsea_UI_inputs_support()], 
#' [fgsea_UI_main_support()], and [fgsea_Server_module()]
#' 
#' Exps Plot tab: 
#' [ExpsPlot_UI_inputs_support()], 
#' [ExpsPlot_UI_main_support()], and [ExpsPlot_Server_module()]
#' 
#' DEstats tab: 
#' [DEstats_UI_inputs_support()], 
#' [DEstats_UI_main_support()], and [DEstats_Server_module()]
#' 
#' Exps Data tab: 
#' [BoxPlot_UI_inputs_support()], 
#' [BoxPlot_UI_main_support()], and [BoxPlot_Server_module()]
#' 
#' VennDiagram tab: 
#' [VennDiagram_UI_inputs_support()], 
#' [VennDiagram_UI_main_support()], and [VennDiagram_Server_module()]
#' @md
NULL

#' Shiny module for VennDiagram
#' @export
#' @rdname Shiny_Tabs
#'
VennDiagram_UI_inputs_support <- function(
id = NULL,
tabsetPanel_id = "tabs") {
# module start
ns <- shiny::NS(id)

shiny::tagList(

# conditionalPanel start
shiny::conditionalPanel(
condition = TabCondition(tabsetPanel_id = tabsetPanel_id, id = id),
shinyWidgets::switchInput(inputId = "con_dir" %>% ns(),
label = "directionality",
value = FALSE),

# Input: contrasts
shiny::selectizeInput(
inputId = "Venn_Contrasts" %>% ns(),
label = "Select 2-4 contrasts:", multiple = TRUE,
options = list(maxItems = 4),
choices = ""
),
shiny::numericInput(
inputId = "padj_cut_DE_Venn" %>% ns(),
label = "padj cut off:", value = 0.1,
min = 0, max = 1, step = 0.01
),
shiny::numericInput(
inputId = "lfc_Venn" %>% ns(),
label = "absolute(log2fc) cut off:", value = 0,
min = 0, max = 20, step = 1
),
shiny::numericInput(
inputId = "name_size" %>% ns(),
label = "name size:", value = 6,
min = 0, max = 50, step = 1
),

# number of Features to show
shiny::numericInput(
inputId = "text_size" %>% ns(),
label = "text size:", value = 6,
min = 0, max = 50, step = 1
),
shiny::actionButton("execute_btn" %>% ns(),
"Execute",
class = "btn-primary"
)
) # condition end
) # tagList end
}

#' Shiny module for VennDiagram
#' @export
#' @inheritParams shiny::moduleServer
#' @inheritParams TabCondition
#' @inheritParams shiny::tabPanel
#' @rdname Shiny_Tabs
#' @example man/examples/Shiny_Tabs_Example.R
VennDiagram_UI_main_support <- function(
id = NULL,
tabsetPanel_id = "tabs",
title = "Venn Diagrams") {

ns <- shiny::NS(id)


shiny::tabPanel(
title = title,
value = id,
shiny::tabsetPanel(

# Venn Diagram
shiny::tabPanel(
title = "Venn Diagram",
value = "Venn_Report_tab" %>% ns(),
shiny::column(
width = 8,
UI_download_plot(id = "DL_VennDiagram" %>% ns()),
shiny::textOutput("Text_Help_Venn_Table" %>% ns()),
shiny::plotOutput(
outputId = "Venn_plot" %>% ns(),
width = "100%", height = "500px"
)
),
shiny::column(
width = 4,
shiny::selectInput(
inputId = "Cat_id" %>% ns(),
label = "Category selection:",
multiple = TRUE,
selectize = FALSE,

choices = "clust"
),
shiny::numericInput(
inputId = "Top_id" %>% ns(),
label = "Top selection:",
value = 20, min = 0, max = Inf, step = 1
)
),
shiny::column(
width = 8,
shiny::textOutput("Text_Venn_Table" %>% ns()),
DT::dataTableOutput("Venn_Table" %>% ns())
)
),

# VennD Heatmap
shiny::tabPanel(
title = "VennD Heatmap",
value = "Venn_Report_tab" %>% ns(),
shiny::fluidRow(DE_pheUI(id = "B" %>% ns()))
)
)
)
}


#' Shiny module for VennDiagram
#' @export
#' @inheritParams shiny::moduleServer
#' @inheritParams DE_pheServer
#' @inheritParams fgsea_
#' @inheritParams Shiny_Hotgenes_Server
#' @rdname Shiny_Tabs
VennDiagram_Server_module <- function(
id = NULL,
input = NULL,
output = NULL,
session = NULL,
Hotgenes = NULL,
ExpressionSlots = reactive(NULL),
SampleIDs = reactive(NULL)) {
shiny::moduleServer(
id,
function(input, output, session) {
# Venn Diagram tab updates ---------------------------------------------

shiny::observe({


if(isTRUE(input$con_dir)){
max <- 2
label <- "Select 2 contrasts:"
} else {
max <- 4
label <- "Select 2-4 contrasts:"
}



shiny::updateSelectizeInput(session, "Venn_Contrasts",
  choices = contrasts_(Hotgenes), 
  options = list(maxItems = max),
  label = label,
  server = FALSE
) 
})


shiny::observe({
shiny::req(length(input$Venn_Contrasts) > 1)
shiny::updateSelectInput(session, "Cat_id",
choices = VennCatChoices()
)
})

# Venn plot
Venn_p <- shiny::eventReactive(input$execute_btn,{
# verifies that more than contrasts available
shiny::req(length(input$Venn_Contrasts) > 1)

shiny::req(all(input$Venn_Contrasts %in% contrasts_(Hotgenes)))

shinybusy::show_modal_spinner(
session = session,
text = "Processing..."
)


if(isTRUE(input$con_dir)){
report <- "contrast_dir"

} else {
report <- "Features"
}

# Report = Features
V <- Hotgenes %>%
DE(
Report = report,
mapFeatures = FALSE,
contrasts = input$Venn_Contrasts,
padj_cut = input$padj_cut_DE_Venn,
.log2FoldChange = input$lfc_Venn
) %>%
Venn_Report(
set_name_size = input$name_size,
text_size = input$text_size
)

# Getting max absolute gene ranks for each gene
# This is to assign ranks for gene to simplify
# selection

All_V_DEGs <- V$Intsect %>%
unlist(use.names = FALSE) %>%
unique()

VennRanks <- rank_hotList(Hotgenes = Hotgenes,
hotList = All_V_DEGs,
contrasts = input$Venn_Contrasts)


#sum(grepl(pattern = ":", fixed = TRUE, x = .data$contrast))

V$df <- V$Intsect %>%
purrr::imap(~ dplyr::as_tibble(.x) %>%
purrr::set_names("Feature") %>%
dplyr::mutate("contrast" = .y) %>%

  
dplyr::mutate("Total" = length(.x)) %>%
dplyr::mutate("contrast" = gsub(":", "&", .data$contrast)) %>%
dplyr::mutate_at(c("contrast", "Total"), as.factor)) %>% 
purrr::list_rbind( )  %>%
dplyr::relocate(c("Total", "contrast"), .before = 1) %>%
# Adding ranks and setting descending order by contrast
dplyr::left_join(VennRanks, by = "Feature") %>%
dplyr::group_by(.data$contrast) %>%
dplyr::mutate(stat = signif(.data$stat, 3)) %>%
dplyr::arrange(dplyr::desc(.data$stat), .by_group = TRUE) %>%
dplyr::slice_head(n = input$Top_id) %>%
dplyr::rename("abs.stat" = "stat") %>%
dplyr::ungroup() %>%
# renaming contrast
dplyr::rename("Category" = .data$contrast)

# progress stop
shinybusy::remove_modal_spinner(session = session)

return(V)
})

# Prep for venn
VennCatChoices <- shiny::reactive({
shiny::req(length(input$Venn_Contrasts) > 1)

unique(Venn_p()$df$Category)
})


output$Venn_plot <- shiny::renderPlot({
print(Venn_p()$vennD)
})

# Venn_Table
output$Text_Help_Venn_Table <- shiny::renderText({
if (length(input$Venn_Contrasts) <= 1) {
"\n\nPlease select 2-4 contrasts"
}
})

output$Text_Venn_Table <- shiny::renderText({
if (length(input$Venn_Contrasts) > 1) {
"Any feature(s) selected from this table
will show up in the 'VennD Heatmap' tab."
}
})

Venn_df <- shiny::reactive({
Venn_p()$df %>%
# tibble::as_tibble() %>%
dplyr::filter(.data$Category %in% input$Cat_id)
})

# getting pre-selected
allRows <- shiny::reactive({
rownames(Venn_df()) %>% as.numeric()
})

output$Venn_Table <- DT::renderDataTable({
DT::datatable(
Venn_df(),
filter = "top",
rownames = FALSE,
selection = list(
mode = "multiple",
target = "row", selected = allRows()
),
extensions = "Buttons",
options = list(
dom = "Bt",
columnDefs = list(list(
className = "dt-left",
targets = "_all"
)),
pageLength = -1,
autoWidth = TRUE,
scrollY = 500,
scrollX = 300,
scrollCollapse = TRUE,
buttons = list(
list(extend = "copy", title = NULL),
"csv",
list(extend = "excel", title = NULL)
)
)
)
})


hotlist_Venn <- shiny::reactive({
shiny::req(input$Venn_Table_rows_selected)

Venn_df() %>%
# dplyr::slice(input$Venn_Table_rows_current) %>%
dplyr::slice(input$Venn_Table_rows_selected) %>%
dplyr::pull("Feature") %>%
unique()
})

Server_download_plot(
  id = "DL_VennDiagram",
  filename = "VennDiagram", bg = "white",
  height = 8, width = 10,
  reactive_plot = shiny::reactive(Venn_p()$vennD)
)


DE_pheServer(
id = "B",
Hotgenes = Hotgenes,
ExpressionSlots = ExpressionSlots,
hotList = hotlist_Venn,
contrasts = shiny::reactive(input$Venn_Contrasts),
padj_cut = shiny::reactive(input$padj_cut_DE_Venn),
.log2FoldChange = shiny::reactive(input$lfc_Venn),
SampleIDs = SampleIDs
) 
} # server end
)
}


#' Shiny module for PCA
#' @export
#' @rdname Shiny_Tabs
#'
PCA_UI_inputs_support <- function(
id = NULL,
tabsetPanel_id = "tabs") {
# module start
ns <- shiny::NS(id)


shiny::tagList(
# conditionalPanel start
shiny::conditionalPanel(
condition = TabCondition(tabsetPanel_id = tabsetPanel_id, id = id),

shiny::selectInput(
inputId = "PCA_contrasts" %>% ns(),
label = "Contrasts selection:",
multiple = TRUE,
choices = ""
),

# Quali
shiny::selectInput(
inputId = "Quali" %>% ns(),
label = "Qualitative variables:",
multiple = TRUE,
choices = ""
),

# Quanti
shiny::selectInput(
inputId = "Quanti" %>% ns(),
label = "Quantitative variables:",
multiple = TRUE,
choices = ""
),
shiny::numericInput(
inputId = "HCPC_max" %>% ns(),
label = "Cluster Max:",
value = 5, min = 0, max = 10, step = 1
),
shiny::numericInput(
inputId = "HCPC_min" %>% ns(),
label = "Cluster Min:",
value = 3, min = 0, max = 10, step = 1
),
shiny::numericInput(
inputId = "padj_cut_PCA" %>% ns(),
label = "padj cut off:", value = 0.1,
min = 0, max = 1, step = 0.01
),

# PCA_signif_
shiny::numericInput(
inputId = "PCA_signif_" %>% ns(),
label = "Significant digits:",
value = 3, min = 0, max = Inf, step = 1
),
shiny::actionButton("goButton2" %>% ns(),
"Execute",
class = "btn-primary"
)
) # condition end
) # tagList end
}

#' Shiny module for PCA
#' @export
#' @rdname Shiny_Tabs
#' @inheritParams shiny::moduleServer
#' @inheritParams TabCondition
#' @inheritParams shiny::tabPanel
PCA_UI_main_support <- function(
id = NULL,
tabsetPanel_id = "tabs",
title = "PCA") {
ns <- shiny::NS(id)


shiny::tabPanel(
title = title,
value = id,
shiny::tabsetPanel(
id = tabsetPanel_id ,

# PCA
shiny::tabPanel(
title = "PCA",
value = "FactoMiner_PCA", #%>% ns(),
shiny::column(
width = 8,

shiny::tags$div(
style="margin-bottom:10px;margin-top:10px;",
# shiny::column(width = 4,


shinyWidgets::materialSwitch(
inputId = "pca_to_gsea" %>% ns(),
label = "GSEA tab support",
inline = TRUE,
value = FALSE),
UI_download_plot(id = "DL_PCA" %>% ns())
#)
),


shiny::plotOutput(outputId = "PCA_plot" %>% ns()),
shiny::h4("Clustered Quali Vars: reports if a given qualitative variable is enriched or
absent from a cluster."),
shiny::h4("Clustered Quanti Vars: reports if a given quantitative variable is enriched or
absent from a cluster"),
shiny::h4("Clustered Feature: reports if a given feature is enriched or
absent from a cluster."),
shiny::h4("PCA Heatmap: any feature(s) selected from
'Clustered Feature' tab will show up here.")
),
shiny::column(
width = 4,
shiny::selectInput(
inputId = "habillage_id" %>% ns(),
label = "Color/ellipse by:",
multiple = FALSE,
choices = "clust"
),
shiny::selectInput(
inputId = "habillage_shape_id" %>% ns(),
label = "Shape by:",
multiple = FALSE,
choices = "clust"
),
shiny::sliderInput(
inputId = "point_size" %>% ns(),
label = "point size", value = 3,
min = 1, max = 10, step = 1
),
shiny::sliderInput(
inputId = "label_size" %>% ns(),
label = "label size:", value = 0,
min = 0, max = 10, step = 1
),
shiny::sliderInput(
inputId = "ellipse.level" %>% ns(),
label = "ellipse level:", value = 0.5,
min = 0, max = 1, step = 0.05
),
shiny::sliderInput(
inputId = "ellipse.alpha" %>% ns(),
label = "ellipse alpha:", value = 1,
min = 0, max = 1, step = 0.05
)
)
),

# Clustered Quali Vars
shiny::tabPanel(
title = "Clustered Quali Vars",
value = "FactoMiner_PCA" %>% ns(),
DT::dataTableOutput("Tab2_PCA_quali_sup" %>% ns()),
shiny::textOutput("Text_Tab2_PCA_quali_sup" %>% ns())
),

# Clustered Quanti Vars
shiny::tabPanel(
title = "Clustered Quanti Vars",
value = "FactoMiner_PCA" %>% ns(),
DT::dataTableOutput("Tab2_PCA_quanti_sup" %>% ns()),
shiny::textOutput("Text_Tab2_PCA_quanti_sup" %>% ns())
),

# Clustered Features
shiny::tabPanel(
title = "Clustered Features",
value = "FactoMiner_PCA" %>% ns(),
DT::dataTableOutput("Clustered_Features" %>% ns()),
shiny::textOutput("Text_Clustered_Features" %>% ns())
),



# PCA Heatmap
shiny::tabPanel(
title = "PCA Heatmap",
value = "FactoMiner_PCA" %>% ns(),
shiny::fluidRow(DE_pheUI(id = "C" %>% ns()))
)
)
)
}


#' Shiny module for PCA
#' @export
#' @inheritParams shiny::moduleServer
#' @inheritParams DE_pheServer
#' @inheritParams fgsea_
#' @inheritParams Shiny_Hotgenes_Server
#' @rdname Shiny_Tabs
PCA_Server_module <- function(
id = NULL,
input = NULL,
output = NULL,
session = NULL,
Hotgenes = NULL,
ExpressionSlots = reactive(NULL),
SampleIDs = reactive(NULL)) {
shiny::moduleServer(
id,
function(input, output, session) {
# onload reactive object --------------------------------------------------
  DF_coldata <- df_coldata_handler(Hotgenes = Hotgenes)

# updating vars -----------------------------------------------------------
shiny::observeEvent(Hotgenes, {
shiny::updateSelectInput(session, "PCA_contrasts",
choices = contrasts_(Hotgenes)
) # %>% observe


shiny::updateSelectInput(session, "Quanti",
choices = c(DF_coldata()$quanti %>% names(),
           auxiliary_assays_(Hotgenes,
                             mode = "quanti") %>% names())
) # %>% observe


shiny::updateSelectInput(session, "Quali",
choices = DF_coldata()$quali %>% names()
) # %>% observe
})


shiny::observe({
shiny::observeEvent(input$Quali, {
shiny::updateSelectInput(session, "habillage_id",
 choices = c("clust", input$Quali)
)
})


shiny::observeEvent(input$Quali, {
shiny::updateSelectInput(session, "habillage_shape_id",
 choices = c("clust", input$Quali)
)
})
})

# Server Tab PCA ----------------------------------------------------------

# resets the goButton2 button when ini() observed below
ini <- shiny::eventReactive(Hotgenes, {
input$goButton2
})

# must be observed
shiny:: observe({
paste0("Setting PCA ini() = ", ini()) # %>% print()
})

# FactoWrapper
output_temp <- shiny::eventReactive(input$goButton2, {
# works with resetting code above
# update the value, will be 1 upon change to Hotgenes
new_value <- ini() / input$goButton2

# print to console, just to check if it updated.
# print(new_value)
paste0("PCA runs when new_value is not 1. new_value = ", new_value) %>%
print()


shiny::req(new_value != 1)


# progress
shinybusy::show_modal_spinner(
session = session,
text = "Processing..."
)


PooledTraits <- c(input$Quanti, input$Quali)

PCAout <- FactoWrapper(
Hotgenes = Hotgenes,
max = input$HCPC_max,
min = input$HCPC_min,
padj_cut = input$padj_cut_PCA,
ExpressionSlots = ExpressionSlots(),
habillage_selection = input$habillage_id,
biplot = FALSE,
ellipse.level = input$ellipse.level,
ellipse.alpha = input$ellipse.alpha,
label_sel = c("all"),
contrasts = input$PCA_contrasts,
coldata_ids = PooledTraits,
aux_features = PooledTraits,

signif_ = input$PCA_signif_,
SampleIDs = SampleIDs()
)


# progress stop
shinybusy::remove_modal_spinner(session = session)


return(PCAout)
})

output$Text_Tab2_PCA_quali_sup <- shiny::renderText({


if (!is.null(output_temp()$TopGroups)) {
"Select a cluster here to see associated
features in the 'Clustered Features' tab."
} else if (is.null(output_temp()$TopGroups)) {
NULL
}
})

output$Tab2_PCA_quali_sup <- DT::renderDataTable({
shiny::req(nrow(output_temp()$TopGroups) > 0)
DT::datatable(output_temp()$TopGroups,
rownames = FALSE, filter = "top",
selection = list(
mode = "single",
target = "row"
),
extensions = "Buttons",
options = list(
dom = "Bt",
scrollY = 500,
scrollX = 300,
scrollCollapse = TRUE,
pageLength = -1,
buttons = list(
list(extend = "copy", title = NULL),
"csv",
list(extend = "excel", title = NULL)
)
)
)
})

# output_Quanti_sup_table
output$Text_Tab2_PCA_quanti_sup <- shiny::renderText({

if (!is.null(output_temp()$TopTibble_sup)) {
"Select a cluster here to see associated
features in the 'Clustered Features' tab."
} else if (is.null(output_temp()$TopTibble_sup)) {
NULL
}
})

# debug 
shiny::observe({
output_temp()$TopTibble_sup %>% print()
})

output$Tab2_PCA_quanti_sup <- DT::renderDataTable({


DT::datatable(output_temp()$TopTibble_sup,
rownames = FALSE, filter = "top",
selection = list(
mode = "single",
target = "row"
),
extensions = "Buttons",
options = list(
dom = "Bt",
scrollY = 500,
scrollX = 300,
scrollCollapse = TRUE,
pageLength = -1,
buttons = list(
list(extend = "copy", title = NULL),
"csv",
list(extend = "excel", title = NULL)
)
)
)
})

# reactive output from Categorical Table 1
# Influences Tab2 PCA Feature Table 2
Clust_Sel <- shiny::reactive({
# get/set Quali_Cluster
if (!is.null(input$Tab2_PCA_quali_sup_rows_selected)) {
Quali_Cluster <- output_temp()$TopGroups %>%
dplyr::slice(input$Tab2_PCA_quali_sup_rows_selected) %>%
dplyr::mutate_at("Cluster", as.character) %>%
dplyr::pull("Cluster")
} else if (is.null(input$Tab2_PCA_quali_sup_rows_selected)) {
Quali_Cluster <- NULL
}

# get/set Quanti_Cluster
if (!is.null(input$Tab2_PCA_quanti_sup_rows_selected)) {
Quanti_Cluster <- output_temp()$TopTibble_sup %>%
dplyr::slice(input$Tab2_PCA_quanti_sup_rows_selected) %>%
dplyr::mutate_at("Cluster", as.character) %>%
dplyr::pull("Cluster")
} else if (is.null(input$Tab2_PCA_quanti_sup_rows_selected)) {
Quanti_Cluster <- NULL
}

ClusterNumber <- c(Quali_Cluster, Quanti_Cluster)
ClusterNumber
})

Clustered_Features <- shiny::reactive({
if (!is.null(Clust_Sel())) {
output_temp()$TopTibble %>%
dplyr::mutate_at("Cluster", as.character) %>%
dplyr::filter(.data$Cluster %in% Clust_Sel())
} else if (is.null(Clust_Sel())) {
output_temp()$TopTibble
}
})



# Feature Table 2
output$Text_Clustered_Features <- shiny::renderText({
if (!is.null(Clustered_Features())) {
"Any feature(s) selected from this table
will show up in the 'PCA Heatmap' tab."
} else if (is.null(Clustered_Features())) {
NULL
}
})

output$Clustered_Features <- DT::renderDataTable(

server = TRUE, 
expr = Clustered_Features(),
rownames = FALSE, filter = "top",
selection = list(
mode = "multiple",
target = "row"
),
extensions = "Buttons",
options = list(
dom = "Bt",
scrollY = 500,
scrollX = 300,
scrollCollapse = TRUE,
pageLength = -1,
buttons = list(
list(extend = "copy", title = NULL),
"csv",
list(extend = "excel", title = NULL)
)
)
)



# output Tab 2 PCA plot -----
PCA_Plot <- shiny::reactive({
if (!"PCA" %in% class(output_temp()$res)) {
output_temp()$res_PPI_pa_1
} else if ("PCA" %in% class(output_temp()$res)) {
shiny::req(input$habillage_id %in% c("clust", output_temp()$quali_sup))

shiny::req(input$habillage_shape_id %in% c("clust", output_temp()$quali_sup))


Final_PCA_OUT <- factoExtra_DFs(
PCA_obj = output_temp(),
repel = TRUE,
label = c("all"),
col.quanti.sup = "red",
col.var = c("black"),
select.var = list(contrib = 10),
col.ind.sup = "black",
title = paste0("Expression: ", ExpressionSlots()),
ellipse.alpha = input$ellipse.alpha,
ellipse.level = input$ellipse.level,
habillage_id = input$habillage_id,
habillage_shape_id = input$habillage_shape_id,
point_size = input$point_size,
label_size = input$label_size
)

return(Final_PCA_OUT)
}
})


output$PCA_plot <- shiny::renderPlot({

final_plot <-  PCA_Plot()%>%
ForcePlot() %>%
grid::grid.draw()



return(final_plot)
})


# Server_download_plot
Server_download_plot(
  id = "DL_PCA",
  filename = "PCA", bg = "white",
  height = 10, width = 8,
  reactive_plot = shiny::reactive({
    PCA_Plot() +
      ggplot2::theme(plot.margin = ggplot2::unit(c(1, 1, 1, 1), "cm"))
  })
)


# PCA Features to Heatmap
DE_pheServer(
id = "C",
Hotgenes = Hotgenes,
ExpressionSlots = ExpressionSlots,
hotList = hotlist_PCA,
annotation_row = hotlist_annotation_row,
contrasts = shiny::reactive(input$PCA_contrasts),
padj_cut = shiny::reactive(input$padj_cut_PCA),
.log2FoldChange = shiny::reactive(0),
SampleIDs = SampleIDs
)


# hotlist_PCA
hotlist_PCA <- shiny::reactive({
anno <- shiny::req(hotlist_annotation_row())

anno %>%
rownames() 
})

# hotlist_PCA
hotlist_annotation_row <- shiny::reactive({
selected_rows <-  shiny::req(input$Clustered_Features_rows_selected)

Clustered_Features() %>%
tibble::as_tibble() %>%
dplyr::slice(selected_rows) %>%
dplyr::select(dplyr::any_of(c("clust"="Cluster", 
      "Feature"))) %>% 
tibble::column_to_rownames("Feature") %>% 
droplevels()
})


pca_out <- shiny::reactive({
if(shiny::isTruthy(input$pca_to_gsea)){

out_ranks <- shiny::req(output_temp()[["Ranks"]])
return(out_ranks)
} else {
NULL
}
})
return(pca_out)
} # server end
)
}


#' Shiny module for Metadata
#' @export
#' @rdname Shiny_Tabs
#' @inheritParams shiny::moduleServer
#' @inheritParams TabCondition
#' @inheritParams shiny::tabPanel
Metadata_UI_main_support <- function(
id = NULL,
tabsetPanel_id = "tabs",
title = "Metadata") {
ns <- shiny::NS(id)


shiny::tabPanel(
title = title,
value = id,
shiny::tabsetPanel(

# Metadata
  shiny::tabPanel(
title = "Sample metadata",
value = "coldata_tab" %>% ns(),
# column(width = 8,
DT::dataTableOutput("coldata_tab_Table" %>% ns())
# )
),

# Mapped names
shiny::tabPanel(
title = "Feature metadata",
value = "MappedFeatures" %>% ns(),
# column(width = 8,
DT::dataTableOutput("MappedFeatures_Table" %>% ns())
# )
)
)
)
}


#' Shiny module for Metadata
#' @export
#' @inheritParams shiny::moduleServer
#' @inheritParams DE_pheServer
#' @inheritParams fgsea_
#' @inheritParams Shiny_Hotgenes_Server
#' @rdname Shiny_Tabs
Metadata_Server_module <- function(
id = NULL,
input = NULL,
output = NULL,
session = NULL,
Hotgenes = NULL,
ExpressionSlots = reactive(NULL),
SampleIDs = reactive(NULL)) {
shiny::moduleServer(
id,
function(input, output, session) {
# onload reactive object --------------------------------------------------
  DF_coldata <- df_coldata_handler(Hotgenes = Hotgenes)



col_Meta <- shiny::reactive({
DT::datatable(DF_coldata()$all,
filter = "top",
extensions = "Buttons", rownames = TRUE,
options = list(
dom = "lBtp",
scrollY = 500,
scrollX = 300,
scrollCollapse = TRUE,
lengthMenu = list(c(10, 25, 100, -1), c("10", "25", "100", "All")),
pageLength = 10,
buttons = list(
list(extend = "copy", title = NULL),
"csv",
list(extend = "excel", title = NULL)
)
)
)
})

output$coldata_tab_Table <- DT::renderDataTable(
{
col_Meta()
},
server = TRUE
)




Mapper_Meta <- shiny::reactive({
DT::datatable(Mapper_(Hotgenes),
filter = "top",
extensions = "Buttons", rownames = FALSE,
options = list(
dom = "lBtp",
scrollY = 500,
scrollX = 300,
scrollCollapse = TRUE,
lengthMenu = list(c(10, 25, 100, -1), c("10", "25", "100", "All")),
pageLength = 10,
buttons = list(
list(extend = "copy", title = NULL),
"csv",
list(extend = "excel", title = NULL)
)
)
)
})

output$MappedFeatures_Table <- DT::renderDataTable(
{
Mapper_Meta()
},
server = TRUE
)
} # server end
)
}


#' Shiny module for DE_phe
#' @export
#' @rdname Shiny_Tabs
#' @inheritParams shiny::moduleServer
#' @importFrom rclipboard rclipboardSetup

# Module UI
DE_pheUI <- function(id) {
# `NS(id)` returns a namespace function, which was save as `ns` and will
# invoke later.
ns <- shiny::NS(id)

shiny::tagList(
  shiny::column(
width = 8,

UI_download_plot(id = "DL_heatmap" %>% ns()),
shiny::plotOutput(
outputId = ns("pheat_plot"),
width = "100%", height = "800px"
)
),
shiny::column(
width = 2,
shiny::h4("Data Selection"),
shiny::numericInput(
inputId = ns("TopLength"),
label = "Top Features:",
value = 50, min = 1,
max = Inf, step = 1
),

# rclipboard --------------------------------------------------------------

rclipboard::rclipboardSetup(),
shiny::uiOutput(outputId = ns("topFeatures")),

##

shiny::selectInput(
inputId = ns("annotations"),
label = "Design:", multiple = TRUE,
choices = ""
),
shiny::selectInput(
inputId = ns("arrangeby"),
label = "Sort by:", multiple = TRUE,
choices = ""
),


shiny::h4("Clustering"),
# Cut rows
shiny::sliderInput(
inputId = ns("col_cut"),
label = "# of col clusters:", value = 0,
min = 0, max = 10, step = 1
),

# Cut rows
shiny::sliderInput(
inputId = ns("row_cut"),
label = "# of row clusters:", value = 0,
min = 0, max = 10, step = 1
),
shiny::h4("Label by"),
# coldata_pal
shiny::selectInput(
inputId = ns("label_by"),
label = "Select alias:",
choices = "Feature"
)
),


# column(),
shiny::column(
width = 2,
shiny::h4("Cell Sizes"),
shiny::sliderInput(
inputId = ns("W"),
label = "Width:", value = 2,
min = 1, max = 30, step = 1
),
shiny::sliderInput(
inputId = ns("H"),
label = "Height:", value = 10,
min = 1, max = 30, step = 1
)
),
shiny::column(
width = 2,
shiny::h4("Color Scheme"),

# coldata_pal
shiny::selectInput(
inputId = ns("coldata_pal"),
label = "Design Palette:",
choices = ""
),
shiny::selectInput(
inputId = ns("col_pal"),
label = "Expression Palette:",
choices = "RdYlBu"
),

# scale range
shiny::sliderInput(
inputId = "rg" %>% ns(),
label = "Scale Range:", value = 2,
min = 1, max = 5, step = 1
),

shiny::h4("Labeling"),
# show_colnames
shiny::checkboxInput(
inputId = ns("show_colnames"),
label = "Columns",
value = FALSE
),

# show_rownames
shiny::checkboxInput(
inputId = ns("show_rownames"),
label = "Rows",
value = TRUE
)
)
)
}


#' @rdname Shiny_Tabs
#' @export
#' @md
#' @param Subtitle reactive object
#' containing a string, usually pathway name to append to
#' heatmap title. Default is NULL.
#' @inheritParams DEphe
#' @inheritParams shiny::moduleServer
#' @inheritParams Shiny_Hotgenes_Server
#' @inheritParams pheatmap::pheatmap
#' @details ExpressionSlots, hotList, contrasts, SampleIDs,
#' padj_cut, and .log2FoldChange must be reactive objects.
#'
#' For example:
#'
#' Correct: ExpressionSlots = NormalizedData_Sel
#'
#' Correct: ExpressionSlots = reactive(input$NormSlot)
#'
#' Incorrect: ExpressionSlots = NormalizedData_Sel()

DE_pheServer <- function(
id,
Hotgenes = NULL,
ExpressionSlots = shiny::reactive(NULL),
contrasts = shiny::reactive(NULL),
padj_cut = shiny::reactive(NULL),
hotList = shiny::reactive(NULL),
.log2FoldChange = shiny::reactive(NULL),
SampleIDs = shiny::reactive(NULL),
Subtitle = shiny::reactive(NULL),
annotation_row = shiny::reactive(NULL)) {
  shiny::moduleServer(
id,
function(input, output, session) {
# onload reactive object --------------------------------------------------

  DF_coldata <- df_coldata_handler(Hotgenes = Hotgenes)

# updates -----------------------------------------------------------------
shiny::observe({
  shiny::updateSelectInput(session, "annotations",
choices = DF_coldata()$quali %>% names()
)
})


shiny::observeEvent(input$annotations, {
  shiny::updateSelectInput(session, "arrangeby",
choices = input$annotations
)
})

# rclipbutton -------------------------------------------------------------

output$topFeatures <- shiny::renderUI({
rclipboard::rclipButton(
inputId = "clipbtn",
label = "Features",
clipText = stringr::str_c(pheat_p()$topHits, collapse = ", "),
icon = shiny::icon("clipboard")
)
})


# coldata_pal
Input_coldata_pal <- RColorBrewer::brewer.pal.info %>%
tibble::rownames_to_column("Name") %>%
dplyr::filter(.data$colorblind == TRUE) %>%
dplyr::filter(.data$category == "qual") %>%
dplyr::pull("Name")

shiny::updateSelectInput(session, "coldata_pal",
choices = Input_coldata_pal,
selected = "Paired"
)
# colors
pheaPal <- RColorBrewer::brewer.pal.info %>%
tibble::rownames_to_column("Name") %>%
dplyr::filter(.data$colorblind == TRUE) %>%
dplyr::filter(.data$category == "div") %>%
dplyr::pull("Name")


shiny::updateSelectInput(session, "col_pal",
choices = pheaPal, selected = "RdYlBu"
)

# label_by

shiny::updateSelectInput(session, "label_by",
choices = Mapper_(Hotgenes) %>% names(),
selected = "Feature"
)
# plots -------------------------------------------------------------------


pheat_p <- shiny::reactive({
# set up column clusters
if (input$col_cut == 0) {
cluster_cols <- FALSE
} else if (input$col_cut > 0) {
cluster_cols <- TRUE
}

# set up row clusters
if (input$row_cut == 0) {
cluster_rows <- FALSE
} else if (input$row_cut > 0) {
cluster_rows <- TRUE
}

# set scale
rg <- input$rg


chr_Subtitle <- as.character(Subtitle())

if (is.null(chr_Subtitle)) {
MainTitle <- paste0("Expression: ", ExpressionSlots())
} else if (!is.null(chr_Subtitle)) {
MainTitle <- paste0(chr_Subtitle, "\nExpression: ", ExpressionSlots())
}

if(!is.null(annotation_row())){
clust_col <- coldata_palettes(annotation_row())

annotation_colors = coldata_palettes(Hotgenes,
             brewer.pals = input$coldata_pal,
             SampleIDs = SampleIDs()
) %>% append(clust_col)
} else {
annotation_colors = coldata_palettes(Hotgenes,
             brewer.pals = input$coldata_pal,
             SampleIDs = SampleIDs()
)
}

Hotgenes %>%
DEphe(
SampleIDs = SampleIDs(),
Topn = input$TopLength,
breaks = seq(-rg, rg, length.out = 100),
hotList = hotList(),
ExpressionSlots = ExpressionSlots(),
contrasts = contrasts(),
padj_cut = padj_cut(),
.log2FoldChange = .log2FoldChange(),
annotations = input$annotations,
arrangeby = input$arrangeby,
annotation_colors = annotation_colors,
label_by = input$label_by,
scale = "row",
show_colnames = input$show_colnames,
show_rownames = input$show_rownames,
border_color = "white",
cutree_rows = input$row_cut,
cutree_cols = input$col_cut,
fontsize_row = input$H,
fontsize_col = input$W,
cellwidth = input$W,
cellheight = input$H,
treeheight_row = 10,
treeheight_col = 10,
main = MainTitle,
cluster_cols = cluster_cols,
cluster_rows = cluster_rows,
color = grDevices::colorRampPalette(
  base::rev(RColorBrewer::brewer.pal(
n = 9,  
name = input$col_pal 
)))(100),
annotation_row = annotation_row()
)
})


output$pheat_plot <- shiny::renderPlot(
print(pheat_p())
)

Server_download_plot(
  id = "DL_heatmap",
  filename = "heatmap", bg = "white",
  reactive_plot = shiny::reactive(pheat_p())
)


}
)
}


#' Shiny module for fgsea
#' @export
#' @rdname Shiny_Tabs
#' @inheritParams shiny::moduleServer
#' @param tabsetPanel_id string matching tabsetPanel id
#' for conditional panel setting. Default is "tabs",
#' but common ones include "tabs".

fgsea_UI_inputs_support <- function(id = "fgsea_tab",
  tabsetPanel_id = "dataset") {
ns <- shiny::NS(id)

shiny::tagList(
shiny::conditionalPanel(
condition = TabCondition(tabsetPanel_id = tabsetPanel_id, id = id),

# Input: contrasts
shiny::selectInput(
inputId = "fgsea_Contrasts" %>% ns(),
label = "Select contrasts:",
multiple = FALSE,
# options = list(maxItems = 4),
choices = ""
),

# MapperCol
shiny::selectInput("input_MapperCol" %>% ns(),
"Feature Identifier",
choices = '',
multiple = FALSE
),
shiny::selectInput("ontology_library" %>% ns(),
label = "Database:",
choices = ''
),

# ontology_sets
shiny::selectInput("ontology_sets" %>% ns(),
label = "Ontology Set:",
choices = '',
multiple = TRUE
),

# ontology_species
shiny::selectInput("ontology_species" %>% ns(),
label = "Species:",
choices = '',
multiple = FALSE
),


# ontology_gene_col
shiny::selectInput("ontology_gene_col" %>% ns(),
label = "Gene Identifier:",
choices = '',
multiple = FALSE
),


# sidepanel input row start -----------------------------------------------

shiny::h5(shiny::strong('Sizes:')),
shiny::fluidRow(
shiny::column(
width = 6,
offset = 0, 
shiny::tags$div(
align = "left",
shiny::numericInput(
inputId = "minSize_Onto" %>% ns(),
width = "400px",
label = "min:", value = 5,
min = 0, max = 100, step = 1
) )
),
shiny::column(
width = 6,
offset = 0,
shiny::tags$div( 

align = "left",
shiny::numericInput(
inputId = "maxSize_Onto" %>% ns(),
width = "400px",
label = "max:", value = 500,
min = 0, max = Inf, step = 1
) )
)
),

# sidepanel input row end -------------------------------------------------

shiny::actionButton("fgsea_Button" %>% ns(),
label = "Execute",
class = "btn-primary"
),
shiny::tags$hr(),
shiny::tags$a(
href = "https://bioconductor.org/packages/release/bioc/html/fgsea.html",
"This app uses fgsea",
target = "_blank"
),
shiny::tags$hr(),
UI_download_file(id = "fgsea_complete" %>% ns(),
                 label = "Unfiltered Results"),
shiny::tags$hr()
# ),
) # condition end
) # tagList end
}


#' Shiny module for GSEA
#' @export
#' @inheritParams shiny::moduleServer
#' @inheritParams fgsea_UI_inputs_support
#' @inheritParams shiny::tabPanel
#' @rdname Shiny_Tabs
fgsea_UI_main_support <- function(
id = "fgsea_tab",
tabsetPanel_id = "dataset",
title = "GSEA") {
ns <- shiny::NS(id)

shiny::tabPanel(
title = title, 
value = id,
shiny::tabsetPanel(
id = tabsetPanel_id ,

shiny::tabPanel(
title = "setup",
value = "setup" %>% ns(),

UI_download_plot(id = "DL_fgseaPlot" %>% ns()),

shiny::plotOutput("fgseaPlot" %>% ns()),
shiny::tabPanel(
"setup",
shiny::tags$a(
href = "https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm#_Interpreting_GSEA_Results",
"Read about interpreting GSEA Results",
target = "_blank"
),
shiny::tags$hr(),
shiny::h4("Refine results"),
shiny::fluidRow(
shiny::column(
width = 2,
shiny::numericInput(
inputId = "metabaseR_FDR_cutoff" %>% ns(),
label = "padj cut off:", value = 1,
min = 0, max = 1, step = 0.01
)
),
shiny::column(2, shiny::numericInput(
inputId = "metabaseR_top_num_paths" %>% ns(),
label = "Top n results:", value = 3,
min = 0, max = 50, step = 1
), offset = 1)
),
shiny::fluidRow(
shiny::column(
4,
shiny::sliderInput(
inputId = "MetaBasePlotFontsize" %>% ns(),
label = "Font size:", value = 20,
min = 0, max = 50, step = 1
)
),
shiny::column(4, shiny::sliderInput(
inputId = "MetaBasePlotLabelBreaks" %>% ns(),
label = "Line width:", value = 50,
min = 10, max = 200, step = 1
))
)
),
DT::dataTableOutput("fgseaTable" %>% ns())
),
shiny::tabPanel(
title = "Heatmap",
value = "phe_leadingFeatures",
shiny::fluidRow(DE_pheUI(id = "A" %>% ns()))
),
shiny::tabPanel(
title = "enrichment",
value = "enrichment",
UI_download_plot(id = "DL_enrichmentPlot" %>% ns()),

shiny::plotOutput("enrichmentPlot" %>% ns())
)
)
)
}


#' Shiny module for GSEA
#' @export
#' @inheritParams shiny::moduleServer
#' @inheritParams DE_pheServer
#' @inheritParams fgsea_
#' @inheritParams Shiny_Hotgenes_Server
#' @importFrom BiocParallel MulticoreParam
#' @param OntologyMethods output of OntologyMethods function.
#' @param Mapper_choices vector containing names of
#' Mapper_ columns that should be used to map expression data
#' to ontology gene identifiers. Default is names(Mapper_(Hotgenes)).
#' A named vector may be used, as well.
#' @param parallel.sz Number of threads of execution
#' to use when doing the calculations in parallel.
#' @param PCA_ranks reactive object goinging named list of ranks from
#' generated from PCA tab. 
#' @rdname Shiny_Tabs
fgsea_Server_module <- function(
id = "fgsea_tab",
input = NULL,
output = NULL,
session = NULL,
Hotgenes = NULL,
PCA_ranks = shiny::reactive(NULL),
ExpressionSlots = shiny::reactive(NULL),
SampleIDs = shiny::reactive(NULL),
OntologyMethods = OntologyMethods(),
Mapper_choices = names(Mapper_(Hotgenes)),
parallel.sz = 10) {
shiny::moduleServer(
id,
function(input, output, session) {

# Getting PCA ranks
PCA_ranks_reactive <- shiny::reactive({
rnks <- PCA_ranks()
return(rnks)
})


con_choices <- shiny::reactive({

choices_con <- contrasts_(Hotgenes)

if(length(choices_con) == 1){
choices<- list(choices_con) %>% 
purrr::set_names(~.x)
} else {
choices <- list(contrasts = choices_con) 
}



if(shiny::isTruthy(PCA_ranks_reactive())){
choices$from_pca_tab <- names(PCA_ranks_reactive())

}

return(choices)

})

# updating from OntologyMethods
db_sel <- shiny::reactive({
OntologyMethods[[input$ontology_library]]
}) %>% 
shiny::bindEvent(shiny::req(input$ontology_library))

shiny::observe({
# updating choices
shiny::updateSelectInput(
session = session,
inputId = "fgsea_Contrasts",
choices = con_choices() #%>% as.character()
)
})


shiny::observeEvent(Hotgenes, {
# input_MapperCol
shiny::updateSelectInput(
session = session, 
inputId = "input_MapperCol",
choices = Mapper_choices
)


new_choices <- OntologyMethods %>%
purrr::imap(function(x, y) {
tibble::as_tibble(glue::glue("{y} {x$Version}"))%>%
dplyr::mutate(name = y)
}) %>% 
purrr::list_rbind() %>%
tibble::deframe()

# ontology_library
shiny::updateSelectInput(
session = session, 
inputId = "ontology_library",
choices = new_choices,
selected = character(0)
)
})

shiny::observe({

# ontology_species
shiny::updateSelectInput(
session = session, inputId = "ontology_species",
choices = db_sel()$species_choices
)

# ontology_sets
shiny::updateSelectInput(
session = session, inputId = "ontology_sets",
choices = db_sel()$Choices
)

# ontology_gene_col
shiny::updateSelectInput(
session = session, inputId = "ontology_gene_col",
choices = db_sel()$gene_cols
)
})


# resets the fgsea_Button button when ini_fgsea() observed below

ini_fgsea <- shiny::eventReactive(Hotgenes, {
input$fgsea_Button
})

# must be observed
shiny::observe({
paste0(
"Setting GSEA ini_fgsea() = ",
ini_fgsea()
) 
})


##
Output_fgsea_reactive <- shiny::eventReactive(input$fgsea_Button, {
# works with resetting code above
# update the value
new_value <- ini_fgsea() / input$fgsea_Button


# print to console, just to check if it updated.
# print(new_value)

paste0("GSEA runs when new_value is not 1. new_value = ",
       new_value) %>% print()

shiny::req(new_value != 1)

# "" is when ontology_species is not available
ontology_species <-shiny::req(input$ontology_species)

shinybusy::show_modal_spinner(
session = session,
spin = "cube-grid",
text = "Processing..."
)

# converts uploaded data to InputRanks, which is required for analysis
list_Ranks <- Hotgenes %>%
DE(
Report = "Ranks",
Rank_name = input$input_MapperCol, # see above
padj_cut = 1
) %>% 
append(PCA_ranks_reactive())

list_Ranks <- list_Ranks[input$fgsea_Contrasts]

# getting pathways
pthyways <- OntologyFunctions(
Methods = OntologyMethods,
db = input$ontology_library,
species = ontology_species,
gene_col = input$ontology_gene_col,
set = input$ontology_sets
)


shinybusy::update_modal_spinner(
paste0(
"Checking ",
length(pthyways),
" pathways"
),
session = session
)


OutPut_fgsea <- fgsea_(
Ranks = list_Ranks,
pathways = pthyways,
BPPARAM = BiocParallel::MulticoreParam(workers = parallel.sz),
minSize = input$minSize_Onto,
maxSize = input$maxSize_Onto
)


shinybusy::remove_modal_spinner(session = session)

return(OutPut_fgsea)
})

# making reactive

summary_fgsea_reactive <- shiny::reactive({
# prevents error when switching between contrasts
shiny::req(input$fgsea_Contrasts %in% names(Output_fgsea_reactive()$Results))

fg_plot <- Output_fgsea_reactive() %>%
GSEA_Plots(
padj_cut = input$metabaseR_FDR_cutoff,
contrasts = input$fgsea_Contrasts,
width = input$MetaBasePlotLabelBreaks,
Topn = input$metabaseR_top_num_paths
)

fg_plot[[1]] +

ggplot2::theme_minimal(base_size = input$MetaBasePlotFontsize)
})

Server_download_plot(
id = "DL_fgseaPlot",
filename = "Summary_GSEA", bg = "white",
reactive_plot = shiny::reactive(summary_fgsea_reactive())
)


# Renders fgseaPlot
output$fgseaPlot <- shiny::renderPlot({
summary_fgsea_reactive()
})

# renders heatmap
DE_pheServer(
id = "A", Hotgenes = Hotgenes,
hotList = LeadingFeatures,
ExpressionSlots = ExpressionSlots,
SampleIDs = SampleIDs,
Subtitle = reactive_selectedPathway
)

LeadingFeatures <- shiny::reactive({
  shiny::req(input$fgseaTable_rows_selected)

LeadingEdgeGenes <- reactive_fgsea_results() %>%
dplyr::slice(input$fgseaTable_rows_selected) %>%
dplyr::pull("leadingEdge") %>%
unlist(use.names = FALSE)


LeadingFeatures <- Mapper_(Hotgenes) %>%
dplyr::filter(.data[[input$input_MapperCol]] %in% LeadingEdgeGenes) %>%
dplyr::pull("Feature")


return(LeadingFeatures)
})

# renders enrichment plot
reactive_selectedPathway <- shiny::reactive({
  shiny::req(input$fgseaTable_rows_selected)

reactive_fgsea_results() %>%
dplyr::slice(input$fgseaTable_rows_selected) %>%
dplyr::pull("pathway") %>%
unlist(use.names = FALSE)
})

# reactive plot
enrichmentPlot_p <- shiny::reactive({
shiny::req(reactive_selectedPathway())

# plotEnrichment_
plotEnrichment_(
Output_fgsea_reactive(), input$fgsea_Contrasts,
reactive_selectedPathway()
)
})

output$enrichmentPlot <- shiny::renderPlot({
enrichmentPlot_p()
})

# Download enrichmentPlot
Server_download_plot(
id = "DL_enrichmentPlot",
filename = "enrichmentPlot",
reactive_plot = shiny::reactive(enrichmentPlot_p())
)



# renders fgseaTable
reactive_fgsea_results <- shiny::reactive({
# prevents error when switching between contrasts
shiny::req(input$fgsea_Contrasts %in% names(Output_fgsea_reactive()$Results))

Output_fgsea_reactive() %>%
fgsea_Results(
padj_cut = input$metabaseR_FDR_cutoff,
Topn = input$metabaseR_top_num_paths,
contrasts = input$fgsea_Contrasts,
mode = "D"
) %>%
purrr::pluck(1) %>%
# set significant digits to simplify table
dplyr::mutate_at(c(
"pval", "padj", "log2err",
"ES", "NES"
), function(x) {
signif(x, digits = 6)
})
})

output$fgseaTable <- DT::renderDataTable({
reactive_fgsea_results() %>%
DT::datatable(
rownames = FALSE,
selection = list(
mode = "single",
target = "row"
),
extensions = "Buttons",
options = list(
dom = "Bt", buttons = c("copy", "csv", "excel"),
scrollY = 500,
scrollX = 300,
scrollCollapse = TRUE,
pageLength = -1
)
)
})

# fgsea_complete

Server_download_file(
  id = "fgsea_complete",
  filename = glue::glue("Unfiltered_fgsea.xlsx"),
  content = function(file){
    
    shinybusy::show_modal_spinner(
      session = session,
      text = "Preparing to download..."
    )
    
    # Complete unfiltered
    # Also removing name, which may be too long
    
   out_file <- fgsea_Results(Output_fgsea_reactive(),
                  padj_cut = 1
    ) %>%
      purrr::imap(function(xCon, yCon) {
        xCon %>%
          dplyr::mutate(Contrast = yCon,
                        .before = 1) 
      }) %>%
      purrr::set_names(NULL) %>%
      openxlsx::write.xlsx(file,
                           colNames = TRUE,
                           rowNames = FALSE,
                           firstCol = FALSE,
                           firstRow = TRUE
      )
    
   shinybusy::remove_modal_spinner(
     session = session
   )  
   
   return(out_file)
    
  }
) # end


}
)
}


#' Shiny module for ExpsPlot
#' @export
#' @importFrom grDevices pdfFonts
#' @rdname Shiny_Tabs
#'
ExpsPlot_UI_inputs_support <- function(
id = NULL,
tabsetPanel_id = "tabs") {
# module start
ns <- shiny::NS(id)

shiny::tagList(

# conditionalPanel start
shiny::conditionalPanel(
condition = TabCondition(tabsetPanel_id = tabsetPanel_id, id = id),

shiny::selectizeInput(
inputId = "yVar" %>% ns(),
label = "Feature:", multiple = TRUE,
options = list(delimiter = ",", create = TRUE),
choices = ""
),
shiny::selectInput(
inputId = "xVar" %>% ns(),
label = "xVar:", multiple = FALSE,
choices = ""
),

# fill
shiny::selectInput(
inputId = "fill" %>% ns(),
label = "Fill:", multiple = FALSE,
choices = ""
),

# color
shiny::selectInput(
inputId = "color" %>% ns(),
label = "Color:", multiple = FALSE,
choices = ""
),

# sidepanel input row start -----------------------------------------------

#shiny::h5(shiny::strong('Font:')),
shiny::fluidRow(
shiny::column(
width = 6,
offset = 0, 
shiny::tags$div(
align = "left",
shiny::selectInput(
inputId = "facet_wrap" %>% ns(),
label = "Facets:",
multiple = TRUE,
width = "700px",
choices = ""

) )
),
shiny::column(
width = 4,
offset = 1,
shiny::tags$div( 

align = "left",
shiny::numericInput(
inputId = "facet_ncol" %>% ns(),
width = "100px",
label = "# cols:", value = 4,
min = 1, max = 100, step = 1
) )
)
),

# sidepanel input row end -------------------------------------------------

shinyWidgets::pickerInput(
inputId = "Filter_Var_levels" %>% ns(),
label = "Reorder:",
multiple = TRUE,
options = list(`actions-box` = TRUE),
choices = ""
),
shiny::selectInput(
inputId = "Var_levels" %>% ns(),
label = "Set levels:",
choices = "",
selected = "",
multiple = TRUE
),

# sidepanel input row start -----------------------------------------------

shiny::h5(shiny::strong('Font:')),
shiny::fluidRow(
shiny::column(
width = 6,
offset = 0, 
shiny::tags$div(
align = "left",
shiny::selectInput(
inputId = "font" %>% ns(),
width = "700px",
label = "Name", 
multiple = FALSE,
choices = names(grDevices::pdfFonts()),
selected = "Helvetica"
) )
),
shiny::column(
width = 5,
offset = 1,
shiny::tags$div( 

align = "left",
shiny::numericInput(
inputId = "base_size" %>% ns(),
width = "100px",
label = "Size", 
value = 20,
min = 0, max = 100, step = 1
) )
)
)

# sidepanel input row end -------------------------------------------------

) # condition end
) # tagList end
}

#' Shiny module for ExpsPlot
#' @export
#' @rdname Shiny_Tabs
#' @inheritParams shiny::moduleServer
#' @inheritParams TabCondition
#' @inheritParams shiny::tabPanel
ExpsPlot_UI_main_support <- function(
id = NULL,
tabsetPanel_id = "tabs",
title = "Exps Plots") {
ns <- shiny::NS(id)


shiny::tabPanel(
title = title,
value = id,

# Row 1
shiny::fluidRow(
shiny::column(
width = 1,
shiny::tags$div(

style="margin-left: 0px;margin-right: 5px;padding-top: 20px",
align = "left",

shiny::actionButton(inputId = "makePlot" %>% ns(),
label = "Execute", 
class = "btn-primary",
icon = shiny::icon("paintbrush"))
)
),

# scales col
shiny::column(
width = 1,

shiny::tags$div(

style="margin-left: 20px;padding-top: 0px",
align = "left",
shiny::h5("Scales:"),
shiny::checkboxInput(
inputId = "free_x" %>% ns(),
label = "free_x",
value = FALSE )
) ),

shiny::column(
width = 1,
shiny::tags$div(

style="margin-left: 5px;padding-top: 25px",
align = "left",
shiny::checkboxInput(inputId = "free_y" %>% ns(),
label = "free_y",
value = FALSE )
) ),

# box plot cols
shiny::column(
width = 1,
shiny::tags$div(

style="margin-left: 5px;padding-top: 0px",
align = "left",
shiny::h5("Show:"),

shiny::checkboxInput(
inputId = "boxplot" %>% ns(),
label = "boxplot",
value = TRUE )
) ),

shiny::column(
width = 1,
shiny::tags$div(

style="margin-left: 5px;padding-top: 25px",
align = "left",
shiny::checkboxInput(
inputId = "pointplot" %>% ns(),
label = "points", 
value = FALSE )
) ),




# column to download button  
shiny::column(
width = 1,
shiny::tags$div(
style = "padding-top: 20px",
UI_download_plot(id = "DL_ScatterPlots" %>% ns(),
label = "PDF"))
)
), # fluid row ends


shiny::conditionalPanel( 
condition = glue::glue("input.yVar != '' "),
ns = ns,
style = "display: none;",
UI_interactive_plot(id = "Int_ScatterPlots" %>% ns(),
width = "100%", height = "700px") %>%
shinycssloaders::withSpinner()),
shiny::fluidRow(
shiny::column(
width = 6,
shiny::verbatimTextOutput("click_info" %>% ns())
),
shiny::column(
width = 6,
shiny::verbatimTextOutput("brush_info" %>% ns())
)
),
DT::dataTableOutput("Express_Dat_Table" %>% ns())
)
}


#' Shiny module for ExpsPlot
#' @export
#' @inheritParams shiny::moduleServer
#' @inheritParams DE_pheServer
#' @inheritParams fgsea_
#' @inheritParams Shiny_Hotgenes_Server
#' @rdname Shiny_Tabs
ExpsPlot_Server_module <- function(
id = NULL,
input = NULL,
output = NULL,
session = NULL,
Hotgenes = NULL,
ExpressionSlots = shiny::reactive(NULL),
SampleIDs = shiny::reactive(NULL)) {
shiny::moduleServer(
id,
function(input, output, session) {
# onload reactive object -------------------------------------------------

  DF_coldata <- df_coldata_handler(Hotgenes = Hotgenes)


# updating vars -----------------------------------------------------------
 shiny::observe({
   shiny::updateSelectInput(session, "facet_wrap",
                            choices = DF_coldata()$quali_distinct %>% names(),
                            selected = ""
   ) 
 })


shiny::observe({
  shiny::updateSelectizeInput(session, "yVar",
                              choices = DF_coldata()$all_features ,
                              selected = "",
                              server = TRUE
  )
})


shiny::observe({

shinyWidgets::updatePickerInput(
session = session,
inputId = "Filter_Var_levels" ,
choices = names(DF_coldata()$quali_distinct)
)

})

shiny::observe({

level_choices <- shiny::req(input$Filter_Var_levels)


shiny::updateSelectInput(session,
inputId = "Var_levels",
choices = DF_coldata()$quali_distinct[level_choices] %>%
 # getting levels
 purrr::imap(function(x, y) {
   level_raw <- x %>% levels()
   level_names <- paste(y, level_raw, sep = "_")
   
   purrr::set_names(level_names, level_raw)
 })
)
})



shiny::updateSelectInput(session, "xVar",
selected = "",
choices = DF_coldata()$quali_distinct %>% names()
) %>%
shiny::observe()

# fill
shiny::updateSelectInput(session, "fill",
selected = "",
choices = DF_coldata()$quali_distinct %>%
names()
) %>% shiny::observe()


# color
shiny::updateSelectInput(session, "color",
selected = "",
choices = DF_coldata()$quali_distinct %>%
names()
) %>% shiny::observe()


# lineVar
shiny::updateSelectInput(session, "lineVar",
selected = "",
choices = DF_coldata()$quali_distinct %>%
names()
) %>% 
shiny::observe()



# reactive plot -----------------------------------------------------------

# ggplot
reactive_scatter <- shiny::eventReactive(input$makePlot,{
# required objects
shiny::req(input$xVar)
shiny::req(input$yVar)

shiny::req(ExpressionSlots() %in% ExpressionSlots_(Hotgenes))

# reordering xVar
if (isTruthy(input$Var_levels)) {
print("reactive_input_vars")
level_choices <- shiny::req(input$Filter_Var_levels)

newMeta <- DF_coldata()$quali_distinct[level_choices] %>%
# getting levels
purrr::imap(function(x, y) {
# To ensure correct mapping between factor levels,
# a name must be assigned

level_raw <- x %>% levels()
level_names <- paste(y, level_raw, sep = "_")

# inverse format as the UI, since UI shows names and
# returns values

New_IDs <- purrr::set_names(level_raw, level_names)

# New output
New_IDs[input$Var_levels] %>% as.character()
}) %>%
purrr::compact()
} else {
newMeta <- NULL
}


scale_ref <- data.frame(
out = c("free","free_y", "free_x", "fixed"),
free_x = c(TRUE, FALSE, TRUE, FALSE),
free_y = c(TRUE, TRUE, FALSE, FALSE))

set_scales <- scale_ref %>% 
dplyr::filter(.data$free_x == input$free_x,
.data$free_y == input$free_y) %>% 
dplyr::pull("out")


# making plot
P <- Hotgenes %>%
ExpsPlot(
xVar = input$xVar,
named_levels = newMeta,
yVar = input$yVar,
fill = input$fill,
color = input$color,
boxplot = input$boxplot,
linevar = input$lineVar,
pointplot = input$pointplot,
basemean = FALSE,
ExpressionSlots = ExpressionSlots(),
SampleIDs = SampleIDs(),
facets = c(input$facet_wrap),
labeller = label_both,
ncol = input$facet_ncol,
scales = set_scales
) +


ggplot2::theme_classic(base_size = input$base_size) +

ggplot2::theme(
plot.title = ggplot2::element_text(hjust = 0.5),
text = ggplot2::element_text(family = input$font),
axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
)


return(P)
})


Server_interactive_plot(
id = "Int_ScatterPlots",
reactive_plot = shiny::reactive(reactive_scatter()),
xvar = shiny::reactive(input$xVar),
yvar = shiny::reactive("value"),
execOnResize = FALSE
)

Server_download_plot(
id = "DL_ScatterPlots",
filename = "Scatter",
reactive_plot = shiny::reactive(reactive_scatter())
)


# ggplot Express_Dat_Table
output$Express_Dat_Table <- DT::renderDataTable({
reactive_scatter() %>%
ExpsPlot_Table(
xVar = input$xVar,
facets = input$facet_wrap,
value_col = "value" 
) %>%
DT::datatable(

filter = "none",
rownames = FALSE,
extensions = "Buttons",
options = list(

dom = "Bfrtip",
pageLength = -1,
scrollY = 500,
scrollX = 300,
scrollCollapse = TRUE,
buttons = list(
list(extend = "copy", title = NULL),
"csv",
list(extend = "excel", title = NULL)
)
)
)
}) %>% shiny::bindEvent(input$makePlot)

}
)
}


#' Shiny module for DEstats
#' @export
#' @rdname Shiny_Tabs
#'
DEstats_UI_inputs_support <- function(
id = NULL,
tabsetPanel_id = "tabs") {
# module start
ns <- shiny::NS(id)

shiny::tagList(

# conditionalPanel start
shiny::conditionalPanel(
condition = TabCondition(tabsetPanel_id = tabsetPanel_id, id = id),

shinyWidgets::pickerInput(
inputId = "Filter_Contrasts" %>% ns(),
label = "Contrasts filter:",
multiple = TRUE,
options = list(`actions-box` = TRUE),
choices = ""
),

shiny::selectInput(
inputId = "DE_Contrasts" %>% ns(),
label = "Contrasts selection:",
multiple = FALSE,
choices = "",
selectize = FALSE,
selected = NULL
),
shiny::h5("'Empty' contrasts dropped"),

shiny::tags$hr(style = "border-color: black;"),
shiny::h5("Filter by stats"),
shiny::numericInput(
inputId = "padj_cut_DE_tables" %>% ns(),
label = "padj cut off:", value = 0.1,
min = 0, max = 1, step = 0.01
),
shiny::numericInput(
inputId = "lfc_Statistics" %>% ns(),
label = "absolute(log2fc) cut off:", value = 0,
min = 0, max = 20, step = 1
),

# signif_
shiny::numericInput(
inputId = "signif_" %>% ns(),
label = "Significant digits:",
value = 3, min = 0, max = Inf, step = 1
),
shiny::tags$hr(style = "border-color: black;"),
shiny::h5("Or pick from here"),

shiny::selectInput(
inputId = "mapper_col" %>% ns(),
label = "aliases selection:",
multiple = FALSE,
choices = "Feature",
selected = "Feature"
),


shiny::selectizeInput(
inputId = "DE_HotList" %>% ns(),
label = "Select Feature:",
options = list(delimiter = ",", create = TRUE),
multiple = TRUE, choices = ""
),


) # condition end
) # tagList end
}

#' Shiny module for DEstats
#' @export
#' @rdname Shiny_Tabs
#' @inheritParams shiny::moduleServer
#' @inheritParams TabCondition
#' @inheritParams shiny::tabPanel
#' @param v_plot_width numeric, width of interactive plot.
#' @param checkbox_width numeric, width of check box area on interactive plot.
DEstats_UI_main_support <- function(
id = NULL,
tabsetPanel_id = "tabs",
title = "DE Stats",
v_plot_width = 8,
checkbox_width = 2) {
ns <- shiny::NS(id)


shiny::tabPanel(
title = title,
value = id,
shiny::tabsetPanel(

# Statistics
shiny::tabPanel(
title = "Statistics",
value = "Statistics" %>% ns(),
UI_download_plot(id = "downloadDESummary_plot" %>% ns()),

shiny::conditionalPanel( 
condition = glue::glue("input.DE_Contrasts != '' "),
ns = ns,
style = "display: none;",
shiny::plotOutput(outputId = "DESummary_plot" %>% ns(),
width = "100%", height = "400px") %>%
shinycssloaders::withSpinner()),


DT::dataTableOutput("Tab3_Query" %>% ns()),
UI_download_file(id = "DE_complete" %>% ns(),
                 label = "Download Complete Unfiltered Results")
),


# Volcano Plots -----------------------------------------------------------

shiny::tabPanel(
title = "Volcano Plots",
value = "Statistics" %>% ns(),
shiny::fluidRow(
shiny::column(
width = v_plot_width,
UI_download_plot(id = "DL_Volcano" %>% ns()),

shiny::conditionalPanel(
condition = glue::glue("input.DE_Contrasts != '' "),
ns = ns,
style = "display: none;",
UI_interactive_plot( id = "v_plot" %>% ns() ) %>%
shinycssloaders::withSpinner()
)

),
shiny::column(
width = checkbox_width,

shiny::checkboxInput(
inputId = "Hide_labels" %>% ns(),
label = "Hide labels:",
value = TRUE
),



# base_size
shiny::numericInput(
inputId = "base_size" %>% ns(),
label = "Font size:",
value = 15,
min = 0, max = 100, step = 1
),
# pointSize
shiny::numericInput(
inputId = "pointSize" %>% ns(),
label = "Point size:",
value = 2,
min = 0, max = 100, step = 1
),
shiny::numericInput(
inputId = "label_point_size" %>% ns(),
label = "Point label size:",
value = 3,
min = 0, max = 100, step = 1
)
)
)
),

# Heatmap
shiny::tabPanel(
title = "Heatmap",
value = "Statistics" %>% ns(),
shiny::fluidRow(DE_pheUI(id = "A" %>% ns()))
)
)
)
}


#' Shiny module for DEstats
#' @export
#' @md
#' @inheritParams shiny::moduleServer
#' @inheritParams DE_pheServer
#' @inheritParams fgsea_
#' @inheritParams Shiny_Hotgenes_Server
#' @inheritParams Server_interactive_plot
#' @rdname Shiny_Tabs
DEstats_Server_module <- function(
id = NULL,
input = NULL,
output = NULL,
session = NULL,
Hotgenes = NULL,
ExpressionSlots = shiny::reactive(NULL),
SampleIDs = shiny::reactive(NULL)) {
shiny::moduleServer(
id,
function(input, output, session) {
# updating vars -----------------------------------------------------------

  shiny::observe({
# Filter_Contrasts
all_choices <- contrasts_(Hotgenes)



shinyWidgets::updatePickerInput(session = session, 
      inputId = "Filter_Contrasts",
      choices = all_choices,
      selected = all_choices) 

})

shiny::observe({

# now only contrasts that have any DEGs                        
contrast_choices <- input$Filter_Contrasts

shiny::updateSelectInput(session, "DE_Contrasts",
choices = contrast_choices,
selected = character(0)) 



})


shiny::observe({

shiny::req(input$mapper_col %in% names(Hotgenes@Mapper))
shiny::updateSelectizeInput(session, "DE_HotList",
  choices = Background_(Hotgenes,
                        Col = input$mapper_col),
  server = TRUE
) 
})


shiny::observe({

shiny::updateSelectInput(
inputId = "mapper_col" ,
choices = names(Hotgenes@Mapper),
selected = "Feature"
)
})


# DE_complete

Server_download_file(
  id = "DE_complete",
  filename = glue::glue("Unfiltered_DE.xlsx"),
  content = function(file){
    shinybusy::show_modal_spinner(
      session = session,
      text = "Preparing to download..."
    )
    
    # Complete unfiltered
    # Also removing name, which may be too long
  out_file  <- Output_DE_(Hotgenes, padj_cut = 1
    ) %>%
      plyr::dlply("contrast", identity) %>% 
      
      purrr::set_names(~base::abbreviate(.x)) %>%
      openxlsx::write.xlsx(file,
                           colNames = TRUE,
                           rowNames = FALSE,
                           firstCol = FALSE,
                           firstRow = TRUE
      )
    
  shinybusy::remove_modal_spinner(
    session = session
  )  
  
  return(out_file)
  }
)


DE_HotList_re <- shiny::reactive({

if(input$mapper_col == "Feature"){
out <- input$DE_HotList
} else {
out <- hotList_Feature(Hotgenes,
hotList = input$DE_HotList)
}

return(out)
})
# DESummary_plot
DEPlot_reactive <- shiny::reactive({


P <- Hotgenes %>%
DEPlot(
contrasts = input$Filter_Contrasts,
padj_cut = input$padj_cut_DE_tables,
.log2FoldChange = input$lfc_Statistics,
hotList = DE_HotList_re()
)

# set asp
deno <- P$data$contrast %>%
unique() %>%
length() / 4 %>% round()

if (deno < 1) {
deno <- 1
}

P <- P +
  ggplot2::theme(
aspect.ratio = 1 / deno,
legend.position = "top"
)

P <- ForcePlot(P)


return(P)
})

output$DESummary_plot <- shiny::renderPlot(
{
DEPlot_reactive() %>%
grid::grid.draw()
},
res = 72,
width = 800,
height = 400
)

# Download summary plot
Server_download_plot(
id = "downloadDESummary_plot",
filename = "SummaryPlot",
reactive_plot = shiny::reactive(DEPlot_reactive())
)


DE_coeff <- shiny::reactive({
  shiny::req(input$DE_Contrasts %in% contrasts_(Hotgenes))

# shinybusy::show_modal_spinner(session = session,
# text = "Loading...")

# DE Statistics
#t1 <- Sys.time()


DE_Output <- DE(Hotgenes,
contrasts = input$DE_Contrasts,
padj_cut = input$padj_cut_DE_tables,
Report = "Details",
hotList = DE_HotList_re(),
.log2FoldChange = input$lfc_Statistics,
mapFeatures = TRUE,
signif_ = input$signif_
)%>% 
purrr::chuck(1)

# progress stop
# shinybusy::remove_modal_spinner(session = session)

#t2 <- Sys.time()

# glue::glue("DE_coeff {signif(t2-t1)}") %>%
# print()
return(DE_Output)
})


output$Tab3_Query <- DT::renderDataTable(
{
DT::datatable(
DE_coeff(),
filter = "none",
extensions = "Buttons", rownames = FALSE,
options = list(
dom = "lpBt",
lengthMenu = list(
c(10, 25, 100, -1),
c("10", "25", "100", "All")
),
pageLength = 10,
scrollY = 500,
scrollX = 300,
scrollCollapse = TRUE,
# pageLength = -1,
# pageLength = 5,
# lengthMenu = c(5, 10, 15, 20),

buttons = list(
list(extend = "copy", title = NULL),
"csv",
list(extend = "excel", title = NULL)
)
)
)
},
server = TRUE
) 


# Server Tab Volcano Plot -------------------------------------------------

Vplot_dat <- shiny::reactive({
#       shinymode = TRUE,

DE_Contrasts <- shiny::req(input$DE_Contrasts)

print("Vplot_dat start")


#t1 <- Sys.time()

VPlot_Out <- VPlot(
shinymode = TRUE,
interactive = FALSE,
repel_labels = FALSE,
Hotgenes = Hotgenes,
contrasts = DE_Contrasts
)


VPlot_Out$refined_cols <- c(names(VPlot_Out$mapper_df),
  "log2FoldChange", "FDR"
)
#t2 <- Sys.time()


# glue::glue("Vplot_dat {signif(t2-t1)}") %>%
# print()


return(VPlot_Out)

})

Vplot_p <- shiny::reactive({
vList <- shiny::req(Vplot_dat())
vList$refined_cols <- NULL

input_List <- 
list( Hide_labels = input$Hide_labels,
padj_cut = input$padj_cut_DE_tables,
.log2FoldChange = input$lfc_Statistics,
#hotList = input$DE_HotList,
hotList = DE_HotList_re(),
base_size = input$base_size,
pointSize = input$pointSize,
#  max.overlaps = input$max.overlaps,
point_label_size = input$label_point_size,
col = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C"))

vList[names(vList) %in% names(input_List)] <-NULL

vList <- vList %>% append(input_List)

#t1 <- Sys.time()


VPlot_Out <- base::do.call(ggplot2_Volcano_prep,vList )


# 
# t2 <- Sys.time()
# 
# glue::glue("Vplot_p {signif(t2-t1)}") %>%
# print()


return(VPlot_Out)
}) 



shiny::observe({

Server_interactive_plot(
id = "v_plot",
refine_cols = Vplot_dat()$refined_cols ,
reactive_plot = shiny::reactive(Vplot_p()),
xvar = shiny::reactive("log2FoldChange"),
yvar = shiny::reactive("FDR"),
execOnResize = FALSE
) 
})


# Download VolcanoPlot
Server_download_plot(
id = "DL_Volcano",
filename = "VolcanoPlot",
reactive_plot = shiny::reactive(Vplot_p())
)




# Tab 5 pheatmap Plot ---------------------------------

DE_pheServer(
id = "A", Hotgenes = Hotgenes,
hotList = DE_HotList_re,
ExpressionSlots = ExpressionSlots,
contrasts = shiny::reactive({
  shiny::req(input$DE_Contrasts) 
}),
padj_cut = shiny::reactive(input$padj_cut_DE_tables),
.log2FoldChange = shiny::reactive(input$lfc_Statistics),
SampleIDs = SampleIDs
)
} # server end
)
}


#' Shiny module for BoxPlot
#' @export
#' @rdname Shiny_Tabs
#'
BoxPlot_UI_inputs_support <- function(id = NULL,
    tabsetPanel_id = "tabs") {
# module start
ns <- shiny::NS(id)

shiny::tagList(
# conditionalPanel start
  shiny::conditionalPanel(
condition = TabCondition(tabsetPanel_id = tabsetPanel_id, id = id),

embed_Hotgenes_logo(),


shiny::tags$hr(style = "border-color: black;"),
shiny::radioButtons(
inputId = "NormSlot" %>% ns(),
label = "Expression slot:",
inline = FALSE,
selected = "",
choices = ""
),


# SampleGroups
shiny::selectInput(
inputId = "SampleGroups" %>% ns(),
label = "Subset by:",
multiple = FALSE,
choices = ""
),


# Show_SampleIDs
shiny::checkboxInput(
inputId = "Show_SampleIDs" %>% ns(),
label = "Sample names",
value = FALSE
),
shinyWidgets::pickerInput(
inputId = "Exps_Samples" %>% ns(),
label = "Select levels:",
multiple = TRUE,
options = list(`actions-box` = TRUE),
choices = ""
),
shiny::h4("Values here applied to all tabs")
) # condition end
) # tagList end
}

#' Shiny module for BoxPlot
#' @export
#' @rdname Shiny_Tabs
#' @importFrom shinycssloaders withSpinner
#' @inheritParams shiny::moduleServer
#' @inheritParams TabCondition
#' @inheritParams shiny::tabPanel
BoxPlot_UI_main_support <- function(id = NULL,
  tabsetPanel_id = "tabs",
  title = "Exps Data") {
ns <- shiny::NS(id)


shiny::tabPanel(
title = title,
value = id,
UI_download_plot(id = "DL_ScatterPlots" %>% ns()),
shiny::conditionalPanel(
condition = glue::glue("input.NormSlot !== '' "),
ns = ns,
style = "display: none;",
shiny::plotOutput(outputId = "Tab1" %>% ns()) %>%
shinycssloaders::withSpinner()
),

UI_download_file(id = "Exps_complete" %>% ns(),
                 label = "Unfiltered Expression data")
)
}


#' Shiny module for BoxPlot
#' @export
#' @importFrom glue glue
#' @importFrom grid grid.draw
#' @param max_col_levels integer setting the maximum
#' number of levels a variable, from the coldata, can have for plots.
#' @param selected_fillby string for to be passed to BoxPlot as initial selection.
#' Default is NULL.
#' @inheritParams BoxPlot
#' @inheritParams shiny::moduleServer
#' @inheritParams DE_pheServer
#' @inheritParams fgsea_
#' @inheritParams Shiny_Hotgenes_Server
#' @rdname Shiny_Tabs
BoxPlot_Server_module <- function(id = NULL,
input = NULL,
output = NULL,
session = NULL,
Hotgenes = NULL,
max_col_levels = 12,
selected_fillby = "") {
  shiny::moduleServer(
id,
function(input, output, session) {
# # onload reactive object --------------------------------------------------
  
  DF_coldata <- df_coldata_handler(Hotgenes = Hotgenes,
                                   max_col_levels = max_col_levels)
  
  

# updating tab 1 on load --------------------------------------------------
reactive_NormSlot <- shiny::reactive({
  shiny::req(input$NormSlot %in% ExpressionSlots_(Hotgenes))

return(input$NormSlot)
})

# getting sampleids for everything else
SampleIDs_Sel <- shiny::reactive({
if (shiny::isTruthy(input$SampleGroups) &
    shiny::isTruthy(input$Exps_Samples)) {
fil_col() %>%
dplyr::filter(dplyr::if_any(
.cols = dplyr::any_of(input$SampleGroups),
.fns = ~ .x %in% input$Exps_Samples
)) %>%
rownames()
} else {
NULL
}
})



# update Exps_Samples for subsetting data
# requires input$SampleGroups
shiny::observeEvent(input$SampleGroups, {
  shiny::req(isTruthy(input$SampleGroups))

choices <- DF_coldata()$quali_distinct %>%
# input$SampleGroups is the column for sub setting
dplyr::select(dplyr::any_of(input$SampleGroups)) %>%
purrr::map(~ .x %>% levels())


shinyWidgets::updatePickerInput(session,
      "Exps_Samples",
      choices = choices,
      selected = ""
)
})



shiny::updateSelectInput(
session,
"SampleGroups",
choices = DF_coldata()$quali_distinct %>% names(),
selected = DF_coldata()$quali_distinct %>%
dplyr::select(dplyr::any_of(selected_fillby)) %>%
names()
) %>%
  shiny::observe()


# QC tab updates ----------------------------------------------------------
shiny::observeEvent(Hotgenes, {
  shiny::updateRadioButtons(session, "NormSlot",
choices = ExpressionSlots_(Hotgenes)
)
})


# col
fil_col <- shiny::reactive({
DF_coldata()$quali_distinct %>%
dplyr::select(dplyr::any_of(input$SampleGroups))
})

# Initial_QC_tab1 ---------------------------------------------------------


Initial_QC_tab1 <- shiny::reactive({
t1 <- Sys.time()

BoxPlot_out <- BoxPlot(
Hotgenes = Hotgenes,
ExpressionSlots = input$NormSlot,
SampleIDs = SampleIDs_Sel()
) +

ggplot2::theme(legend.position = "top")

BoxPlot_out_no_names <- BoxPlot_out +
  ggplot2::theme(
axis.title.x = ggplot2::element_blank(),
axis.text.x = ggplot2::element_blank(),
axis.ticks.x = ggplot2::element_blank()
)

Pout <- list(
BoxPlot_out = BoxPlot_out,
BoxPlot_out_no_names = BoxPlot_out_no_names
)
t2 <- Sys.time()

glue::glue("Initial_QC_tab1 {signif(t2-t1)}") %>%
print()

return(Pout)
}) %>%
  shiny::bindCache(
input$NormSlot,
Hotgenes,
SampleIDs_Sel()
) %>%
  shiny::bindEvent(
  shiny::req(shiny::isTruthy(input$NormSlot)),

shiny::req(input$NormSlot %in% ExpressionSlots_(Hotgenes)),
Hotgenes,
SampleIDs_Sel()
)



# This allows for the final PDF to match the screen
QC_tab1 <- shiny::reactive({


t1 <- Sys.time()

if(base::isFALSE(input$Show_SampleIDs)) {

Pout <- "BoxPlot_out_no_names"
} else {

Pout <- "BoxPlot_out"
}


OutPlot <- Pout

t2 <- Sys.time()

glue::glue("QC_tab1 {signif(t2-t1)}") %>%
print()

return(OutPlot)
}) %>%

shiny::bindEvent(
input$Show_SampleIDs
)


# Download Expression data plot
Server_download_plot(
id = "DL_ScatterPlots",
filename = "BoxPlot",

reactive_plot = shiny::reactive({
BoxPlot_out <- Initial_QC_tab1()[[QC_tab1()]]
if (shiny::isTruthy(input$SampleGroups)) {


c_breaks <- coldata_palettes(coldata_(Hotgenes),
       coldata_ids = input$SampleGroups
)[[input$SampleGroups]]

BoxPlot_out <- BoxPlot_out +
  ggplot2::geom_boxplot(
  ggplot2::aes(fill = .data[[input$SampleGroups]]),
outlier.shape = NA,
stat = "identity",
orientation = "x"
) +
  ggplot2::scale_fill_manual(values = c_breaks, name = input$SampleGroups)
}

return(BoxPlot_out)
})
)

# renderplot
output$Tab1 <- shiny::renderPlot(
{
t1 <- Sys.time()
BoxPlot_out <- Initial_QC_tab1()[[QC_tab1()]]

if (shiny::isTruthy(input$SampleGroups)) {
c_breaks <- coldata_palettes(coldata_(Hotgenes),
       coldata_ids = input$SampleGroups
)[[input$SampleGroups]]

BoxPlot_out <- BoxPlot_out +
  ggplot2::geom_boxplot(
  ggplot2::aes(fill = .data[[input$SampleGroups]]),
outlier.shape = NA,
stat = "identity",
orientation = "x"
) +
ggplot2::scale_fill_manual(values = c_breaks, name = input$SampleGroups)
}


final_p <- BoxPlot_out %>%
ForcePlot() %>%
grid::grid.draw()

t2 <- Sys.time()

glue::glue("renderPlot Tab1 {signif(t2-t1)}") %>%
print()

return(final_p)
},
res = 72,
width = 800,
height = 400,
execOnResize = FALSE
) %>% shiny::bindCache(
QC_tab1(),
Initial_QC_tab1(),
input$SampleGroups
)



Server_download_file(
  id = "Exps_complete",
  filename = glue::glue("Complete_ExpsDat.xlsx"),
  content = function(file){
    
    shinybusy::show_modal_spinner(
      session = session,
      text = "Preparing to download..."
    )
    
    
  out_file <- ExpressionSlots_(Hotgenes) %>%
      purrr::set_names(~.x) %>%
      purrr::imap(
        ~ DExps(
          Hotgenes,
          Query_set = TRUE,
          Q = TRUE,
          ExpressionSlots = .y
        ) %>%
          dplyr::relocate(
            dplyr::any_of(names(
              DF_coldata()$all
            )),
            .before = 1
          ) %>%
          tibble::rownames_to_column() %>%
          base::t()
      ) %>%
      openxlsx::write.xlsx(
        file,
        colNames = FALSE,
        rowNames = TRUE,
        firstCol = TRUE,
        firstRow = TRUE
      )
  
  shinybusy::remove_modal_spinner(
    session = session
  )  
  
  return(out_file)
  
  }
)



# final output ------------------------------------------------------------

Final_output <- list(
DF_coldata = DF_coldata,
SampleIDs_Sel = SampleIDs_Sel,
reactive_NormSlot = reactive_NormSlot
)

return(Final_output)
}
)
}




#
#' @noRd
# @importFrom grid grid.draw
# @importFrom ggplot2 ggplot_build ggplot_gtable
# @param P ggplot object
#'
ForcePlot <- function(P = NULL) {
TF <- ggplot2::ggplot_build(P)
TF <- ggplot2::ggplot_gtable(TF)
# TF <- grid::grid.draw(TF)


return(TF)
}


#' Shiny module for Auxiliary Assays
#' @export
#' @rdname Shiny_Tabs
AuxAssays_sidepanel <- function(
id = NULL,
tabsetPanel_id = "tabs") {

# module start
ns <- shiny::NS(id)


shiny::tagList(

# conditionalPanel start
shiny::conditionalPanel(
condition = TabCondition(tabsetPanel_id = tabsetPanel_id, id = id),

shiny::textOutput(ns("num_aux_features")),

shiny::tags$br(),


shiny::selectInput(inputId = "var_include" %>% ns(),
label = "Select metadata to include",
selected = "",
choices = "",
multiple = TRUE),

shiny::fileInput("raw_auxfile" %>% ns(), 
multiple = TRUE,
"Choose aux data file", 
accept = ".csv"
)


) # condition end
) # tagList end}
}

#' Shiny module for Auxiliary Assays
#' @export
#' @rdname Shiny_Tabs
#' @inheritParams shiny::moduleServer
#' @inheritParams TabCondition
#' @inheritParams shiny::tabPanel
AuxAssays_mainpanel <- function(
id = NULL,
tabsetPanel_id = "tabs",
title = "Auxiliary Assays") {
ns <- shiny::NS(id)


shiny::tabPanel(
title = title,
value = id,
shiny::tabsetPanel(

# Available Auxiliary_Assays
  shiny::tabPanel(
title = "Available Auxiliary Assays",
value = "auxassay_tab" %>% ns(),


DT::dataTableOutput("auxassay_tab_Table" %>% ns())

),

# template generation
shiny::tabPanel(
title = "Make template",
value = "MakeTemplate" %>% ns(),

shiny::h3("Use this tab to generate a template for upload."),

shiny::h3("Any new must be alighed with the SampleIDs column"),

# column(width = 8,
DT::dataTableOutput("MakeTemplate_Table" %>% ns())
),

# uploaded
shiny::tabPanel(title = "Uploaded Auxiliary Assays",
value = "raw_auxassay_tab" %>% ns(),

shiny::textOutput(ns("removed_features")),

shiny::tags$br(),

shiny::actionButton("Update_Button" %>% ns(),
 label = "Update",
 class = "btn-primary"
),
shiny::tags$hr(),

DT::dataTableOutput("raw_auxassay_tab_Table" %>% ns())
)

)
)
}


#' Shiny module for Auxiliary Assays
#' @export
#' @importFrom readr read_csv cols
#' @inheritParams shiny::moduleServer
#' @inheritParams DE_pheServer
#' @inheritParams fgsea_
#' @inheritParams Shiny_Hotgenes_Server
#' @rdname Shiny_Tabs
AuxAssays_server <- function(
id = NULL,
input = NULL,
output = NULL,
session = NULL,
Hotgenes = NULL,
ExpressionSlots = shiny::reactive(NULL),
SampleIDs = shiny::reactive(NULL)) {

  shiny::moduleServer(
id,
function(input, output, session) {


# check feature names -----------------------------------------------------

feature_drop_re <- shiny::reactive({

shiny::req(raw_aux_table_re())

parse_features(Hotgenes = Hotgenes, 
features = names(raw_aux_table_re()))
})      


# update
Hotgenes_re <- shiny::reactive({


if(shiny::isTruthy(input$Update_Button)){


print("updating")


auxiliary_assays_(Hotgenes) <- Final_auxdata_re() 


print("done")
}

return(Hotgenes)


})

output$removed_features <- shiny::renderText({
drp_these <- feature_drop_re() %>% 
unlist(use.names = FALSE) %>% 
stringr::str_c( collapse = ', ') 

return(glue::glue("Removed: {drp_these}") )
})







# number of available features --------------------------------------------
shiny::observeEvent( Hotgenes_re(),{
output$num_aux_features <- shiny::renderText({ 


Hotgenes <- Hotgenes_re()

glue::glue("Total auxiliary assays: {length(auxiliary_assays_features(Hotgenes))}") 
})


})



# raw uploaded data --------------------

raw_aux_table_re <- shiny::reactive({
  shiny::req(shiny::isTruthy(input$raw_auxfile))

uploadedData <- readr::read_csv(file = input$raw_auxfile$datapath,
      col_types = readr::cols())
return(uploadedData)
})

Final_auxdata_re <- shiny::reactive({

exlcudeThese <-feature_drop_re() %>%  
  unlist(use.names = FALSE)

raw_aux_table_re() %>% 
dplyr::select(-dplyr::any_of(exlcudeThese))
})

raw_auxdata_re <- shiny::reactive({

Final_auxdata_re()  %>% 

DT::datatable(
filter = "top",
#extensions = "Buttons", 
rownames = FALSE,
options = list(
dom = "ltp",
scrollY = 500,
scrollX = 300,
scrollCollapse = TRUE,
lengthMenu = list(c(10, 25, 100, -1), c("10", "25", "100", "All")),
pageLength = 10
)
)

})


# available data
DF_auxassays <- shiny::reactive({

#if(shiny::isTruthy(Hotgenes_re())){
Hotgenes <- Hotgenes_re()
#}

# get aux_assays
aux_data_wide <- Hotgenes %>% 
auxiliary_assays_() %>% 
tibble::rownames_to_column("SampleIDs")

if(shiny::isTruthy(input$var_include)){
# merging with phenodata
coldata_Out <- coldata_(Hotgenes, coldata_ids = input$var_include) %>% 
tibble::rownames_to_column("SampleIDs") %>% 
dplyr::left_join(aux_data_wide, by = "SampleIDs") 

return(coldata_Out)
} else {
return(aux_data_wide)
}

# var_include
})

# update inputs
shiny::observe({

shiny::updateSelectInput(session = session,
inputId = "var_include",
choices = coldata_names(Hotgenes_re())) 
})

col_AuxAssay <- shiny::reactive({
DT::datatable(DF_auxassays(),
filter = "top",
extensions = "Buttons", rownames = FALSE,
options = list(
dom = "lBtp",
scrollY = 500,
scrollX = 300,
scrollCollapse = TRUE,
lengthMenu = list(c(10, 25, 100, -1), c("10", "25", "100", "All")),
pageLength = 10,
buttons = list(
list(extend = "copy", title = NULL),
"csv"
)
)
)
})

output$auxassay_tab_Table <- DT::renderDataTable(
{
col_AuxAssay()
},
server = FALSE
)




AuxAssayUpdateTemplate <- shiny::reactive({
DT::datatable(coldata_(Hotgenes_re(), coldata_ids = input$var_include) %>% 
tibble::rownames_to_column("SampleIDs"),
filter = "top",
extensions = "Buttons", rownames = FALSE,
options = list(
dom = "lBtp",
scrollY = 500,
scrollX = 300,
scrollCollapse = TRUE,
lengthMenu = list(c(10, 25, 100, -1), c("10", "25", "100", "All")),
pageLength = 10,
buttons = list(
list(extend = "copy", title = NULL),
"csv"
)
)
)
})

output$MakeTemplate_Table <- DT::renderDataTable(
{
AuxAssayUpdateTemplate()
},
server = FALSE
)


output$raw_auxassay_tab_Table <- DT::renderDataTable(
{
raw_auxdata_re()
},
server = FALSE
)
return(Hotgenes_re())
} # server end
)

# return(outObj)
}
