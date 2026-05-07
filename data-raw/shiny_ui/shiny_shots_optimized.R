#library(Hotgenes)
library(shinytest2)
#library(magick)
#library(cli)

funs_path <- file.path(getwd(), "data-raw/shiny_ui/funs_helpers.R")
source(funs_path)

if(FALSE){
  rm(list = ls()); gc()
  devtools::document(); devtools::load_all()
  
  file.edit(funs_path)
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

all_DE_features <-DE(HotgenesObj, Topn = 1) |> unlist(use.names = FALSE) |> unique()
all_contrasts <- contrasts_(HotgenesObj)
all_expr_slots <- ExpressionSlots_(HotgenesObj)
all_coldata_cols <- coldata_names(HotgenesObj)
all_features <- Background_(HotgenesObj, Col = "Feature")
all_mapper_cols <- names(Mapper_(HotgenesObj))

default_contrast <- all_contrasts[1]
default_expr_slot <- all_expr_slots[1]
default_coldata_col <- all_coldata_cols[3]
default_feature <- all_DE_features[1]
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



# setup_output_directory -----------------------------------------------------------

#fig_dir <- file.path(tempdir(), "hotgenes_figures")
fig_dir <- file.path(getwd(), "vignettes", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cli::cli_h2("Output directory")
cli::cli_inform("Figures will be saved to: {.file {fig_dir}}")

# launch_app -----------------------------------------------------------

cli::cli_h2("Launching Shiny app")

# OntologyMethods <- Hotgenes::OntologyMethods()

OntologyMethods <- Hotgenes::OntologyMethods(
  Ontology_Function = list(
    msigdbr  = Hotgenes::msigdbr_wrapper),
  InputChoices = list(
    # msigdbr  = Hotgenes::msigdbr_wrapper_choices()$set
    
    msigdbr  =   c("C1", "CGP", "CP", "CP:BIOCARTA", "CP:KEGG_LEGACY", "CP:KEGG_MEDICUS", 
      "CP:PID", "CP:REACTOME", "CP:WIKIPATHWAYS", "MIR:MIRDB", "MIR:MIR_LEGACY", 
      "TFT:GTRD", "TFT:TFT_LEGACY", "3CA", "CGN", "CM", "GO:BP", "GO:CC", 
      "GO:MF", "HPO", "C6", "IMMUNESIGDB", "VAX", "C8", "H")
  ),
  gene_col_choices = list(
    msigdbr  = c("gene_symbol",  "ensembl_gene")
  ),
  species_choices = list(
    msigdbr  = c("human", "mouse", "rat", "dog")
    
  ),
  versions = list(
    msigdbr  = utils::packageVersion("msigdbr"))
)

app <- AppDriver$new(
  app = Shiny_Hotgenes(HotgenesObj, OntologyMethods = OntologyMethods),
  height = 1000,
  width = 1600
)

# Wait for first render at 1x
app$wait_for_value(output = "Hotgenes_A-BoxPlot-Tab1", timeout = 15000)

cli::cli_alert_success("App launched successfully")

# set_retina_quality -----------------------------------------------------------

cli::cli_h2("Configuring device pixel ratio (retina 2x)")

app$get_chromote_session()$Emulation$setDeviceMetricsOverride(
  width = 1600,
  height = 1000,
  deviceScaleFactor = 2,
  mobile = FALSE
)
app$get_js("$(window).trigger('resize')")

cli::cli_alert_success("Device pixel ratio set to 2x (retina quality)")

# capture_boxplot_tab -----------------------------------------------------------

cli::cli_h2("Capturing BoxPlot tab")

if(FALSE) {
  browseURL(fig_dir)
}


# Poll for 2x render
check_v <- app$get_js("$('#Hotgenes_A-BoxPlot-Tab1').find('img').prop('naturalWidth')")
while (is.null(check_v) || check_v <= 1195) {
  Sys.sleep(0.5)
  check_v <- app$get_js("$('#Hotgenes_A-BoxPlot-Tab1').find('img').prop('naturalWidth')")
  cli::cli_inform("naturalWidth: {check_v}")
}


###
take_screenshot(
  app = app,
  tab_id = "Hotgenes_A-BoxPlot",
  filename = "shiny-01-boxplot_raw.png",
  output_dir = fig_dir,
  wait_time = 1
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

# capture_destats_volcano -----------------------------------------------------------

cli::cli_h2("Capturing DEstats Volcano tab")

app$get_js("$('[role=\"tab\"]:contains(\"Volcano Plots\")').click()")

app$wait_for_idle()

take_screenshot(
  app = app,
  tab_id = "Hotgenes_A-DEstats",
  filename = "shiny-02b-volcano_raw.png",
  inputs = list("Hotgenes_A-DEstats-DE_Contrasts" = default_contrast),
  output_dir = fig_dir,
  wait_time = 1.5
)

annotate_screenshot(
  input_path = file.path(fig_dir, "shiny-02b-volcano_raw.png"),
  label = paste0("B2  Volcano Plot — contrast: ", default_contrast),
  output_path = file.path(fig_dir, "shiny-02b-volcano.png")
)

# capture_destats_heatmap -----------------------------------------------------------

cli::cli_h2("Capturing DEstats Heatmap tab")
#app$get_js("$('[role=\"tab\"]:contains(\"Heatmap\")').click()")
#app$get_js("$('#Hotgenes_A-DEstats').find('[role=\"tab\"]:contains(\"Heatmap\")').first().click()")
app$get_js("$('.tab-pane[data-value=\"Hotgenes_A-DEstats\"] [role=\"tab\"]:contains(\"Heatmap\")').first().click()")
app$wait_for_idle()


take_screenshot(
  app = app,
  tab_id = "Hotgenes_A-DEstats",
  filename = "shiny-02c-heatmap_raw.png",
  inputs = construct_inputs("DEstats",
                            DE_Contrasts = default_contrast,
                            #DE_HotList = NULL,
                            padj_cut_DE_tables = 0.1),
  output_dir = fig_dir,
  wait_time = 2
)

annotate_screenshot(
  input_path = file.path(fig_dir, "shiny-02c-heatmap_raw.png"),
  label = paste0("B3  DE Heatmap — contrast: ", default_contrast),
  output_path = file.path(fig_dir, "shiny-02c-heatmap.png")
)

# capture_expsplot_tab -----------------------------------------------------------

cli::cli_h2("Capturing ExpsPlot tab")


take_screenshot(
  app = app,
  tab_id = "Hotgenes_A-ExpsPlot",
  filename = "shiny-03-expsplot_raw.png",
  inputs = construct_inputs("ExpsPlot",
                             yVar = default_feature,
                            fill = "sh",
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

# Set ontology_library first and wait for observer to populate ontology_sets
app$set_inputs("Hotgenes_A-GSEA-ontology_library" = "msigdbr", wait_ = FALSE)
app$wait_for_idle()

take_screenshot(
  app = app,
  tab_id = "Hotgenes_A-GSEA",
  filename = "shiny-06-gsea_raw.png",
  inputs = construct_inputs("GSEA",
                             fgsea_Contrasts = default_contrast,
                            
                            ontology_sets = "H",
                           
                             input_MapperCol = all_mapper_cols[1]),
  output_dir = fig_dir,
  button_id = "Hotgenes_A-GSEA-fgsea_Button",
  wait_time = 10
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
  "shiny-02b-volcano.png",
  "shiny-02c-heatmap.png",
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
# unlink(fig_dir, recursive = TRUE)
n_files <- length(all_files)

cli::cli_h2("Generated files ({n_files} total)")

file_info <- lapply(all_files, function(f) {
  size_kb <- file.size(file.path(fig_dir, f)) / 1024
  paste0("  • {.file ", f, "} ({format(size_kb, digits = 0)} KB)")
})

cli::cli_bullets("{unlist(file_info)}")

cli::cli_h2("Output directory")
cli::cli_inform("Path: {.file {fig_dir}}")

# browseURL(fig_dir)

#app$get_logs()
