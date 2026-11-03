#' DeeDee
#'
#' `DeeDee` is a package that allows for the visual comparison of two or more
#'   DEA results.
#'
#' `DeeDee` contains functions for the preprocessing of raw DEA result data and
#'   multiple functions that create graphical output and give an impression of
#'   similarities and differences in the expression profiles. The functions are
#'   contained individually as well as included in an interactive Shiny Web App.
#'
#' @importFrom AnnotationDbi keytypes
#' @importFrom checkmate assert_choice assert_data_frame assert_list
#'   assert_logical assert_number assert_subset assertChoice assertClass
#'   assertDataFrame test_subset
#' @importFrom clusterProfiler enrichGO
#' @importFrom ComplexHeatmap draw Heatmap
#' @importFrom ComplexUpset intersection_size upset
#' @importFrom DESeq2 results DESeq
#' @importFrom dplyr bind_rows full_join inner_join
#' @importFrom enrichplot emapplot pairwise_termsim
#' @importFrom ggplot2 aes annotate coord_cartesian geom_line geom_point ggplot
#'   labs scale_fill_manual theme_light xlab ylab
#' @importFrom ggvenn ggvenn
#' @importFrom InteractiveComplexHeatmap makeInteractiveComplexHeatmap
#' @importFrom writexl write_xlsx
#' @importFrom readxl excel_sheets read_excel
#' @importFrom shiny actionButton brushedPoints brushOpts checkboxGroupInput
#'   checkboxInput column conditionalPanel downloadButton downloadHandler
#'   fileInput fluidRow includeMarkdown navbarPage need numericInput
#'   observeEvent plotOutput reactive reactiveValues renderPlot renderTable
#'   renderText renderUI req selectInput shinyApp tableOutput tabPanel
#'   textOutput uiOutput validate
#' @importFrom shinyBS bsCollapse bsCollapsePanel bsModal
#' @importFrom shinycssloaders withSpinner
#' @importFrom shinythemes shinytheme
#' @importFrom stats approx complete.cases
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom tools file_ext
#' @importFrom utils read.table
#' @importFrom viridis scale_color_viridis viridis
#'
#' @name DeeDee-pkg
#' @docType package
"_PACKAGE"

globalVariables(
  c("x.pval", "y.pval", "logFC1", "logFC2", "rowname")
)
