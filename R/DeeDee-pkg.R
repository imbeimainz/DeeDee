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
#' @importFrom checkmate assert_choice assert_data_frame assert_list assert_logical
#' assert_number assert_subset assertChoice assertClass assertDataFrame test_subset
#' @importFrom circlize colorRamp2
#' @importFrom clusterProfiler enrichGO
#' @importFrom ComplexHeatmap draw Heatmap rowAnnotation
#' @importFrom ComplexUpset intersection_size upset
#' @importFrom DESeq2 results DESeq vst
#' @importFrom dplyr bind_rows full_join inner_join
#' @importFrom edgeR topTags
#' @importFrom enrichplot emapplot pairwise_termsim
#' @importFrom ggplot2 aes annotate coord_cartesian geom_bar geom_line geom_point
#' geom_text ggplot labs position_dodge scale_fill_manual theme_bw theme_light
#' xlab ylab ggplot_build
#' @importFrom ggvenn ggvenn
#' @importFrom InteractiveComplexHeatmap makeInteractiveComplexHeatmap
#' InteractiveComplexHeatmapOutput
#' @importFrom limma topTable
#' @importFrom writexl write_xlsx
#' @importFrom readxl excel_sheets read_excel
#' @importFrom rlang .data
#' @importFrom rmarkdown render
#' @importFrom rvest html_node
#' @importFrom shiny actionButton brushedPoints brushOpts checkboxGroupInput
#' checkboxInput column conditionalPanel downloadButton downloadHandler fileInput
#' fluidPage fluidRow HTML includeHTML includeMarkdown isTruthy navbarPage need
#' numericInput observeEvent plotOutput reactive reactiveValues removeNotification
#' renderPlot renderPrint renderTable renderText renderUI req selectInput
#' shinyApp showNotification tableOutput tabPanel tagList tags textOutput
#' uiOutput validate verbatimTextOutput
#' @importFrom shinydashboard tabBox
#' @importFrom shinyBS bsCollapse bsCollapsePanel bsModal
#' @importFrom shinycssloaders withSpinner
#' @importFrom shinythemes shinytheme
#' @importFrom stats approx complete.cases na.omit
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importFrom S4Vectors metadata metadata<- DataFrame SimpleList
#' @importFrom SummarizedExperiment rowData rowData<- mcols assays
#' rowData rowData<- SummarizedExperiment
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom tools file_ext
#' @importFrom utils read.table data packageVersion browseURL
#' @importFrom viridis scale_color_viridis viridis
#' @importFrom xml2 read_html write_html
#' @importFrom methods show as callNextMethod is new validObject
#'
#' @name DeeDee-pkg
#' @docType package
NULL

globalVariables(
  c("x.pval", "y.pval", "logFC1", "logFC2", "rowname")
)
