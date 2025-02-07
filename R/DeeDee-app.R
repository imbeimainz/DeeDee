#' DeeDee App
#'
#' @description `deedee_app` opens the DeeDee Shiny web application, combining
#' the functionalities of all other DeeDee functions with a user-friendly
#' graphical user interface.
#'
#' @param deedee_obj An object of the class DeeDeeExperiment to be analyzed.
#' @param dde An object of the class DeeDeeExperiment to be analyzed.
#'
#' @return A shiny app
#' @export
#'
#' @examples
#' data("de_named_list", package = "DeeDee")
#' library("SummarizedExperiment")
#'
#' rd_macrophage <- DataFrame(
#'   gene_id = rownames(de_named_list$ifng_vs_naive))
#' rownames(rd_macrophage) <- rownames(de_named_list$ifng_vs_naive)
#' se_macrophage_noassays <- SummarizedExperiment(
#'   assays = SimpleList(),
#'   rowData = rd_macrophage
#' )
#' dde <- DeeDeeExperiment(
#'   se_macrophage_noassays,
#'   de_results = de_named_list
#' )
#' dde
#'
#' if (interactive()) {
#'   deedee_app(dde = dde)
#' }
#'
#' # one can provide the DeeDeeExperiment object also at runtime
#' if (interactive()) {
#'   deedee_app()
#' }
deedee_app <- function(deedee_obj = NULL,
                       dde = NULL) {


  # ui definition -----------------------------------------------------------
  deedee_ui <- shiny::navbarPage(
    title = "DeeDee",
    id = "tabs",
    theme = shinythemes::shinytheme("flatly"),

    # ui - data input ----------------------------------------------------------
    shiny::tabPanel(
      title = "Input",
      tags$head(
        tags$style(
          shiny::HTML(
            ".shiny-output-error-validation {
              font-size: 15px;
              color: forestgreen;
              text-align: center;
              }
            "
          )
        )
      ),

      uiOutput("ui_panel_datainput"),

      # shiny::downloadButton("vignette",
      #                "Download DeeDee Package vignette (.html)")
    ),

    # ui - scatter -------------------------------------------------------------
    shiny::tabPanel(
      title = "Scatterplot",

      uiOutput("ui_panel_scatter"),
    ),


    # ui - heatmap -------------------------------------------------------------
    shiny::tabPanel(
      title = "Heatmap",
      id = "heatmap",

      uiOutput("ui_panel_heatmap"),

    ),

    # ui - venn ----------------------------------------------------------------
    shiny::tabPanel(
      title = "Venn Diagram",

      uiOutput("ui_panel_venn"),
    ),


    # ui - upset ---------------------------------------------------------------
    shiny::tabPanel(
      title = "UpSet Plot",

      uiOutput("ui_panel_upset"),
    ),


    # ui - qq ------------------------------------------------------------------
    shiny::tabPanel(
      title = "Quantile-Quantile Plot",

      uiOutput("ui_panel_qq"),
    ),


    # ui - cat -----------------------------------------------------------------
    shiny::tabPanel(
      title = "Concordance At the Top Plot",

      uiOutput("ui_panel_cat"),
    ),

    # ui - summary -------------------------------------------------------------
    shiny::tabPanel(
      title = "Summary",

      uiOutput("ui_panel_summary")
    )
  )




  # server definition ----------------------------------------------------------
  # nocov start
  deedee_server <- function(input, output, session) {

    # initializing reactives ---------------------------------------------------
    reactive_values <- reactiveValues()

    if (!is.null(dde)) {
      message("DeeDeeExperiment object provided")

      reactive_values$dde <- dde
    } else {
      message("no DeeDeeExperiment data provided, upload at runtime expected")
      reactive_values$dde <- NULL
    }

    # server - data input ------------------------------------------------------
    output$ui_panel_datainput <- renderUI({
      tagList(
        shiny::fluidRow(
          shiny::column(
            width = 8,
            shiny::fileInput(
              inputId = "upload_de",
              label = "Upload your DEA results or DeeDee objects",
              multiple = TRUE,
              accept = c(".rds", ".txt", ".xlsx"),
              placeholder = "No files selected"
            ),

            shiny::fileInput(
              inputId = "upload_deedee",
              label = "Upload your DeeDeeExperiment object",
              # multiple = TRUE,
              multiple = FALSE,
              accept = c(".rds"),
              placeholder = "no DeeDeeExperiment provided"
            ),

            verbatimTextOutput("print_dde"),
            shiny::tableOutput("inp_infobox"),
            verbatimTextOutput("print_dde_touse")
          ),
          shiny::column(
            width = 4,
            shiny::selectInput(
              inputId = "in_organism",
              label = "Organism",
              choices = list(
                "Human" = "org.Hs.eg.db",
                "Mouse" = "org.Mm.eg.db",
                "Fly" = "org.Dm.eg.db",
                "Rat" = "org.Rn.eg.db"
              )
            ),
            shiny::uiOutput("ui_key_inp"),
            shiny::uiOutput("ui_datasets"),
            shiny::uiOutput("ui_datasets_deedee"),
            shiny::conditionalPanel(
              "output.inp_infobox",
              shiny::downloadButton(
                "btn_inp_download",
                "Download DeeDee object (.RDS)"
              )
            )
          )
        ),
        shinyBS::bsCollapse(
          shinyBS::bsCollapsePanel(
            title = "INFO - input",
            shiny::includeMarkdown(
              system.file("extdata", "input.md", package = "DeeDee")
            ),
            style = "primary"
          )
        )
      )
    })
    output$ui_key_inp <- shiny::renderUI({
      shiny::req(input$in_organism)
      anno <- input$in_organism
      require(anno, character.only = TRUE)
      shiny::selectInput(
        inputId = "in_key_type",
        label = "Key type of gene IDs",
        choices = keytypes(get(anno))
      )
    })

    mydata <- shiny::reactive({
      shiny::req(
        shiny::isTruthy(input$upload_de) || shiny::isTruthy(deedee_obj)
      )

      ext <- c()
      res <- list()

      if (!is.null(deedee_obj)) {
        if (is(deedee_obj, "DeeDeeObject")) {
          ext[1] <- "arg"
          # df <- data.frame(filename, type, contrast, genes)
          res[[1]] <- deedee_obj@DeeDeeList
        }
        else {
          stop("Your argument input is not of type 'DeeDeeObject'. Please check your input and re-open the application.")
        }
      }

      k <- 0

      # reading out input files
      if (length(input$upload_de[, 1] > 0)) {
        for (i in (length(ext) + 1):(length(input$upload_de[, 1]) + length(ext))) {
          k <- k + 1
          ext[i] <- tools::file_ext(input$upload_de[k, "datapath"])
          shiny::validate(
            shiny::need(
              ext[[i]] == "rds" ||
                ext[[i]] == "RDS" ||
                ext[[i]] == "xlsx" ||
                ext[[i]] == "txt",
              message = "Please upload only .RDS, .xlsx or .txt files"
            )
          )

          # .RDS input
          if (ext[[i]] == "rds" || ext[[i]] == "RDS") {
            res[[i]] <- readRDS(input$upload_de[[k, "datapath"]])
            res[[i]] <- rds_input(obj = res[[i]],
                                  nm = input$upload_de[[k, "name"]])

            # .xlsx input
          } else if (ext[[i]] == "xlsx") {
            sheets <- readxl::excel_sheets(input$upload_de[[k, "datapath"]])
            res[[i]] <- xlsx_input(obj = sheets,
                                   path= input$upload_de[[k, "datapath"]],
                                   nm = input$upload_de[k, "name"])

            # .txt input
          } else if (ext[[i]] == "txt") {
            temp <- utils::read.table(input$upload_de[[k, "datapath"]])
            res[[i]] <- list(temp)
            names(res[[i]]) <-
              unlist(strsplit(input$upload_de[k, "name"], split = ".", fixed = TRUE))[1]
          } else {
            return(NULL)
          }
        }
      }

      # merging input data structures
      dat <- list()
      for (i in seq_len(length(res))) {
        for (j in seq_len(length(res[[i]]))) {
          checkmate::assert_subset(names(res[[i]][[j]]), c("logFC", "pval"))
          dat[[names(res[[i]])[[j]]]] <- res[[i]][[j]]
        }
      }

      ###### to-replace ####### obj <- DeeDeeLegacy::DeeDeeObject(DeeDeeList = dat)
      obj <- dat #### TODO

      return(obj)
    })

    my_deedee <- reactive({
      # TODO: from the list to the DeeDeeExperiment

      # TODO: or directly from the uploads


      # TODO: or as of now, passed directly from the argument
      reactive_values$dde
    })

    output$ui_datasets <- shiny::renderUI({
      shiny::req(
        shiny::isTruthy(input$upload_de) || shiny::isTruthy(deedee_obj)
      )
      shiny::checkboxGroupInput(
        inputId = "select_datasets",
        label = "Select datasets to be used",
        choices = names(mydata()@DeeDeeList),
        selected = names(mydata()@DeeDeeList)
      )
    })

    output$ui_datasets_deedee <- shiny::renderUI({
      # shiny::req(
      #   shiny::isTruthy(input$upload_de) || shiny::isTruthy(deedee_obj)
      # )
      shiny::req(
        shiny::isTruthy(reactive_values$dde)
      )
      shiny::checkboxGroupInput(
        inputId = "select_datasets_deedee",
        label = "Select datasets to be used",
        choices = names(dea(my_deedee())),
        selected = names(dea(my_deedee()))
      )
    })

    mydata_use <- shiny::reactive({
      shiny::req(
        shiny::isTruthy(input$upload_de) || shiny::isTruthy(deedee_obj)
      )
      use <- input$select_datasets
      dat2 <- list()
      for (i in use) {
        dat2[i] <- mydata()@DeeDeeList[i]
      }
      return(dat2)
    })


    output$print_dde <- renderPrint({
      my_deedee()
    })

    my_deedee_use <- reactive({
      to_use <- input$select_datasets_deedee
      to_remove <- setdiff(names(dea(my_deedee())), to_use)

      dde_to_use <- remove_dea(my_deedee(), dea_name = to_remove)
      dde_to_use
    })

    output$print_dde_touse <- renderPrint({
      my_deedee_use()
    })

    output$btn_inp_download <- shiny::downloadHandler(
      filename = "DeeDee_object.RDS",
      content = function(file) {
        # shiny::req(
        #   shiny::isTruthy(input$upload_de) || shiny::isTruthy(deedee_obj)
        # )
        shiny::req(shiny::isTruthy(reactive_values$dde))
        ###### to-replace ####### dl <- DeeDeeLegacy::DeeDeeObject(DeeDeeList = mydata_use())
        ###### to-replace ####### saveRDS(dl, file)
      }
    )

    output$inp_infobox <- shiny::renderTable({
      shiny::validate(shiny::need(
        !is.null(mydata()@DeeDeeList),
        message = "Faulty input data provided."
      ))

      shiny::req(
        shiny::isTruthy(input$upload_de) || shiny::isTruthy(deedee_obj)
      )

      input_infobox(
        deedee_obj = deedee_obj,
        sets = input$upload_de,
        md = mydata()
      )
    })


    # server - scatter ---------------------------------------------------------
    output$ui_panel_scatter <- renderUI({
      tagList(
        shiny::fluidRow(
          shiny::column(
            width = 4,
            shiny::uiOutput("ui_scatter_choices1"),
            shiny::uiOutput("ui_scatter_choices2"),
            shiny::selectInput(
              inputId = "in_scatter_color_by",
              label = "Color by",
              choices = list(
                "1st p-value" = "pval1",
                "2nd p-value" = "pval2"
              ),
              selected = "pval1"
            ),
            shiny::numericInput(
              inputId = "in_scatter_pthresh",
              label = "P-value threshold",
              value = 0.05, min = 0.01, max = 1, step = 0.01
            ),
            shiny::actionButton(
              inputId = "btn_ora_button",
              label = "Over-representation analysis"
            )
          ),
          shiny::column(
            width = 8,
            shinycssloaders::withSpinner(
              shiny::plotOutput(
                outputId = "plot_deedee_scatter",
                # dblclick = "scatter_dblclick",
                brush = shiny::brushOpts(
                  id = "scatter_brush",
                  resetOnNew = FALSE
                )
              )
            )
          )
        ),
        shinyBS::bsCollapse(
          shinyBS::bsCollapsePanel(
            title = "INFO - scatter plot",
            shiny::includeMarkdown(
              system.file("extdata", "scatter.md", package = "DeeDee")
            ),
            style = "primary"
          )
        ),
        shinyBS::bsModal(
          id = "modalExample",
          title = "Gene Ontology over-representation analysis",
          trigger = "btn_ora_button",
          size = "large",
          shinycssloaders::withSpinner(
            shiny::plotOutput("scatter_ora")
          ),
          shiny::downloadButton(
            outputId = "btn_ora_download",
            label = "Download enrichment result object (.RDS)"
          )
        ),
        shiny::downloadButton(
          outputId = "btn_scatter_brush_download",
          label = "Download brushed genes (.xlsx)"
        ),
        shiny::tableOutput("scatter_brush_info")
      )
    })

    output$ui_scatter_choices1 <- shiny::renderUI({
      # shiny::req(
      #   shiny::isTruthy(input$upload_de) || shiny::isTruthy(deedee_obj)
      # )
      shiny::req(shiny::isTruthy(reactive_values$dde))
      shiny::selectInput(
        inputId = "in_scatter_select1",
        label = "1st data set",
        choices = names(dea(my_deedee_use()))
      )
    })

    output$ui_scatter_choices2 <- shiny::renderUI({
      # shiny::req(
      #   shiny::isTruthy(input$upload_de) || shiny::isTruthy(deedee_obj)
      # )
      shiny::req(shiny::isTruthy(reactive_values$dde))
      shiny::selectInput(
        inputId = "in_scatter_select2",
        label = "2nd data set",
        selected = names(dea(my_deedee_use()))[2],
        choices = names(dea(my_deedee_use()))
      )
    })

    # --- plot output ---
    ranges <- shiny::reactiveValues(
      x = NULL,
      y = NULL
    )

    output$plot_deedee_scatter <- shiny::renderPlot({
      shiny::validate(
        shiny::need(
          # length(mydata()@DeeDeeList) >= 2,
          length(dea(my_deedee())) >= 2,
          message = "Please upload at least two contrasts."
        )
      )
      shiny::validate(
        shiny::need(
          # length(mydata_use()) >= 2,
          length(dea(my_deedee_use())) >= 2,
          message = "Please select at least two contrasts."
        )
      )

      sel1 <- match(input$in_scatter_select1, names(dea(my_deedee_use())))
      sel2 <- match(input$in_scatter_select2, names(dea(my_deedee_use())))
      shiny::req(sel1)
      shiny::req(sel2)
      res <- deedee_scatter(
        my_deedee_use(),
        select1 = sel1,
        select2 = sel2,
        color_by = input$in_scatter_color_by,
        pthresh = input$in_scatter_pthresh
      )

      shiny::validate(
        shiny::need(
          !is.null(res),
          message = "No common genes in input datasets.")
      )
      res +
        ggplot2::coord_cartesian(
          xlim = ranges$x,
          ylim = ranges$y,
          expand = FALSE
        )
    })

    # --- brushing ---
    # scatter_download_button <- shiny::renderUI(
    #     conditionalPanel("output.scatter_brushed()",
    #     downloadButton("btn_scatter_brush_download",
    #                    "Download brushed genes (.txt)")))

    scatter_brushed <- shiny::reactive({
      shiny::req(input$scatter_brush)
      df <- data.frame(
        x = mydata_use()[[input$in_scatter_select1]],
        y = mydata_use()[[input$in_scatter_select2]]
      )

      df <- subset(df, x.pval < 0.05 & y.pval < 0.05)

      names(df) <- c(
        paste(input$in_scatter_select1, ".logFC", sep = ""),
        paste(input$in_scatter_select1, ".pval", sep = ""),
        paste(input$in_scatter_select2, ".logFC", sep = ""),
        paste(input$in_scatter_select2, ".pval", sep = "")
      )
      shiny::brushedPoints(
        df,
        input$scatter_brush,
        xvar = paste(input$in_scatter_select1, ".logFC", sep = ""),
        yvar = paste(input$in_scatter_select2, ".logFC", sep = "")
      )
    })

    output$scatter_brush_info <- shiny::renderTable(
      {
        shiny::req(scatter_brushed())
        scatter_brushed()
      },
      rownames = TRUE
    )

    output$btn_scatter_brush_download <- shiny::downloadHandler(
      filename = "scatter_brushed_genes.xlsx",
      content = function(file) {
        shiny::req(input$scatter_brush)
        bru <- tibble::rownames_to_column(scatter_brushed())
        first <- data.frame(bru[1], bru[2], bru[3])
        nm1 <- unlist(strsplit(names(first)[2],
                               split = ".",
                               fixed = TRUE
        ))[1]
        names(first) <- c("rowname", "logFC", "pval")
        second <- data.frame(bru[1], bru[4], bru[5])
        nm2 <- unlist(strsplit(names(second)[2],
                               split = ".",
                               fixed = TRUE
        ))[1]
        names(second) <- c("rowname", "logFC", "pval")
        l <- list(first, second)
        names(l) <- c(nm1, nm2)
        writexl::write_xlsx(l, file, col_names = TRUE)
      }
    )

    # --- enrich ---
    enrich <- shiny::reactive({
      # shiny::req(
      #   shiny::isTruthy(input$upload_de) || shiny::isTruthy(deedee_obj)
      # )
      shiny::req(shiny::isTruthy(reactive_values$dde))
      shiny::validate(
        shiny::need(
          scatter_brushed(),
          message = "No brushed genes."
        )
      )

      sel1 <- match(input$in_scatter_select1, names(mydata_use()))
      sel2 <- match(input$in_scatter_select2, names(mydata_use()))
      shiny::req(sel1)
      shiny::req(sel2)
      data <- list(mydata_use()[[sel1]], mydata_use()[[sel2]])

      res <- ora(
        geneList = scatter_brushed(),
        universe = data,
        orgDB = input$in_organism,
        key_type = input$in_key_type
      )
      shiny::validate(
        shiny::need(
          is(res, "enrichResult"),
          message = "Not working."
        )
      )

      return(res)
    })

    output$scatter_ora <- shiny::renderPlot({
      shiny::validate(
        shiny::need(
          !is.null(enrich()),
          message = "Something went wrong..."
        )
      )
      en <- enrich()
      en_df <- as.data.frame(en)
      shiny::validate(
        shiny::need(
          nrow(en_df) > 0,
          message = "No enriched terms found."
        )
      )

      options(ggrepel.max.overlaps = Inf)
      plt <- enrichplot::emapplot(enrichplot::pairwise_termsim(en))
      shiny::validate(
        shiny::need(
          !is.null(plt),
          message = "No enriched terms found."
        )
      )
      return(plt)
    })

    output$btn_ora_download <- shiny::downloadHandler(
      filename = "enrichment_results.RDS",
      content = function(file) {
        shiny::req(!is.null(enrich()))
        saveRDS(enrich(), file)
      }
    )

    # server - heatmap ---------------------------------------------------------
    output$ui_panel_heatmap <- renderUI({
      tagList(
        shiny::fluidRow(
          shiny::column(
            width = 4,
            shiny::numericInput(
              inputId = "in_heatmap_show_first",
              label = "Show first",
              value = 25,
              min = 1
            ),
            shiny::checkboxInput(
              inputId = "in_heatmap_show_gene_names",
              label = "Show gene names",
              value = FALSE
            ),
            shiny::checkboxInput(
              inputId = "in_heatmap_showNA",
              label = "Show NA",
              value = FALSE
            ),
            shiny::selectInput(
              inputId = "in_heatmap_dist",
              label = "Distance measure",
              choices = list(
                "Euclidean" = "euclidean",
                "Manhattan" = "manhattan",
                "Pearson" = "pearson",
                "Spearman" = "spearman"
              ),
              selected = "euclidean"
            ),
            shiny::selectInput(
              inputId = "in_heatmap_clust",
              label = "Clustering method",
              choices = list(
                "Single" = "single",
                "Complete" = "complete",
                "Average" = "average",
                "Centroid" = "centroid"
              ),
              selected = "average"
            ),
            shiny::numericInput(
              inputId = "in_heatmap_pthresh",
              label = "P-value threshold",
              value = 0.05, min = 0.01, max = 1, step = 0.01
            ),
            shiny::actionButton(
              inputId = "btn_heatmap_action",
              label = "Create heatmap")
          ),
          shiny::column(
            width = 8,
            shiny::textOutput("heatmap_errors"),
            shiny::conditionalPanel(
              "output.heatmap_errors == ''",
              shinycssloaders::withSpinner(
                InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput()
              )
            )
          )
        ),
        shinyBS::bsCollapse(
          shinyBS::bsCollapsePanel(
            title = "INFO - heatmap",
            shiny::includeMarkdown(
              system.file("extdata", "heatmap.md",package = "DeeDee")
            ),
            style = "primary"
          )
        )
      )
    })

    output$heatmap_errors <- shiny::renderText({
      shiny::validate(
        shiny::need(
          # length(mydata()@DeeDeeList) >= 2,
          length(dea(my_deedee())) >= 2,
          message = "Please upload at least two contrasts."
        )
      )
      shiny::validate(
        shiny::need(
          # length(mydata_use()) >= 2,
          length(dea(my_deedee_use())) >= 2,
          message = "Please select at least two contrasts."
        )
      )
      shiny::validate(
        shiny::need(
          !is.null(heatmap_output()),
          message = "No common genes in input datasets."
        )
      )
      ""
    })

    heatmap_output <- shiny::reactive({
      # shiny::req(
      #   shiny::isTruthy(input$upload_de) || shiny::isTruthy(deedee_obj)
      # )
      shiny::req(shiny::isTruthy(reactive_values$dde))
      shiny::req(input$in_heatmap_show_first)
      # shiny::req(mydata_use())
      res <- deedee_heatmap(
        my_deedee_use(),
        show_first = input$in_heatmap_show_first,
        show_gene_names = input$in_heatmap_show_gene_names,
        dist = input$in_heatmap_dist,
        clust = input$in_heatmap_clust,
        pthresh = input$in_heatmap_pthresh,
        show_na = input$in_heatmap_showNA
      )
      shiny::validate(
        shiny::need(!is.null(res), message = "No common genes in input datasets.")
      )

      res <- ComplexHeatmap::draw(res)

      return(res)
    })

    listen <- shiny::reactive({
      list(
        input$in_heatmap_show_first,
        input$in_heatmap_show_gene_names,
        input$in_heatmap_dist,
        input$in_heatmap_clust,
        input$in_heatmap_pthresh,
        input$in_heatmap_showNA,
        mydata_use()
      )
    })

    global <- reactiveValues(notify = FALSE)

    shiny::observeEvent(listen(), {
      global$notify <- TRUE
    })

    shiny::observeEvent(input$tabs, {
      if (input$tabs == "Heatmap") {
        if (global$notify == TRUE) {
          shiny::showNotification(
            ui = "Something changed. Click 'Create heatmap' to reload.",
            duration = NULL,
            type = "warning",
            id = "heatmap_warning"
          )
        }
      } else {
        shiny::removeNotification(id = "heatmap_warning")
      }
    })

    shiny::observeEvent(listen(), {
      req(input$tabs == "Heatmap")
      shiny::showNotification(
        ui = "Something changed. Click 'Create heatmap' to reload.",
        duration = NULL,
        type = "warning",
        id = "heatmap_warning"
      )
    })

    shiny::observeEvent(input$btn_heatmap_action, {
      shiny::req(
        # length(mydata_use()) >= 2
        length(dea(my_deedee_use())) >= 2,
      )
      global$notify <- FALSE
      InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
        input,
        output,
        session,
        heatmap_output()
      )
      shiny::removeNotification("heatmap_warning")
    })


    # server - venn ------------------------------------------------------------
    output$ui_panel_venn <- renderUI({
      tagList(
        shiny::fluidRow(
          shiny::column(
            width = 4,
            shiny::selectInput(
              inputId = "in_venn_mode",
              label = "Mode",
              choices = list(
                "Up" = "up",
                "Down" = "down",
                "Both" = "both"
              ),
              selected = "both"
            ),
            shiny::numericInput(
              inputId = "in_venn_pthresh",
              label = "P-value threshold",
              value = 0.05, min = 0.01, max = 1, step = 0.01
            )
          ),
          shiny::column(
            width = 8,
            shinycssloaders::withSpinner(
              shiny::plotOutput("plot_deedee_venn")
            )
          )
        ),
        shinyBS::bsCollapse(
          shinyBS::bsCollapsePanel(
            title = "INFO - venn",
            shiny::includeMarkdown(
              system.file("extdata", "venn.md", package = "DeeDee")
            ),
            style = "primary"
          )
        )
      )
    })

    output$plot_deedee_venn <- shiny::renderPlot({
      shiny::validate(
        shiny::need(
          # length(mydata()@DeeDeeList) >= 2,
          length(dea(my_deedee())) >= 2,
          message = "Please upload at least two contrasts."
        )
      )
      shiny::validate(
        shiny::need(
          # length(mydata_use()) >= 2,
          length(dea(my_deedee_use())) >= 2,
          message = "Please select at least two contrasts."
        )
      )

      # shiny::req(
      #   shiny::isTruthy(input$upload_de) || shiny::isTruthy(deedee_obj)
      # )
      shiny::req(shiny::isTruthy(reactive_values$dde))
      res <- deedee_venn(
        my_deedee_use(),
        mode = input$in_venn_mode,
        pthresh = input$in_venn_pthresh
      )

      shiny::validate(
        shiny::need(
          !is.null(res),
          message = "No genes in your datasets. Maybe your specified p-value threshold is too low?"
        )
      )

      res
    })


    # server - upset -----------------------------------------------------------
    output$ui_panel_upset <- renderUI({
      tagList(
        shiny::fluidRow(
          shiny::column(
            width = 4,
            shiny::selectInput(
              inputId = "in_upset_mode",
              label = "Mode",
              choices = list(
                "Up" = "up",
                "Down" = "down",
                "Both" = "both"
              ),
              selected = "both"
            ),
            shiny::conditionalPanel(
              condition = "input.in_upset_mode == 'both'",
              shiny::checkboxInput(
                inputId = "in_upset_colored",
                label = "Coloring",
                value = TRUE
              )
            ),
            shiny::numericInput(
              inputId = "in_upset_minset",
              label = "Minimum set size",
              value = 10, min = 0, step = 1
            ),
            shiny::numericInput(
              inputId = "in_upset_pthresh",
              label = "P-value threshold",
              value = 0.05, min = 0.01, max = 1, step = 0.01
            )
          ),
          shiny::column(
            width = 8,
            shinycssloaders::withSpinner(
              shiny::plotOutput("plot_deedee_upset")
            )
          )
        ),
        shinyBS::bsCollapse(
          shinyBS::bsCollapsePanel(
            title = "INFO - upset",
            shiny::includeMarkdown(
              system.file("extdata", "upset.md", package = "DeeDee")
            ),
            style = "primary"
          )
        )
      )
    })

    output$plot_deedee_upset <- shiny::renderPlot({
      shiny::validate(
        shiny::need(
          # length(mydata()@DeeDeeList) >= 2,
          length(dea(my_deedee())) >= 2,
          message = "Please upload at least two contrasts."
        )
      )
      shiny::validate(
        shiny::need(
          # length(mydata_use()) >= 2,
          length(dea(my_deedee_use())) >= 2,
          message = "Please select at least two contrasts."
        )
      )

      # shiny::req(
      #   shiny::isTruthy(input$upload_de) || shiny::isTruthy(deedee_obj)
      # )
      shiny::req(shiny::isTruthy(reactive_values$dde))
      if (input$in_upset_mode == "both" && input$in_upset_colored) {
        mode <- "both_colored"
      } else {
        mode <- input$in_upset_mode
      }

      res <- deedee_upset(
        my_deedee_use(),
        mode = mode,
        pthresh = input$in_upset_pthresh,
        min_setsize = input$in_upset_minset
      )

      shiny::validate(
        shiny::need(
          !is.null(res),
          message = "No genes in your datasets. Maybe your specified p-value threshold is too low?"
        )
      )
      res
    })


    # server - qq --------------------------------------------------------------
    output$ui_panel_qq <- renderUI({
      tagList(
        shiny::fluidRow(
          shiny::column(
            width = 4,
            shiny::checkboxInput(
              inputId = "in_qq_multiple",
              label = "Multiple",
              value = FALSE
            ),
            shiny::conditionalPanel(
              condition = "!input.in_qq_multiple",
              shiny::uiOutput("ui_qq_choices1"),
              shiny::uiOutput("ui_qq_choices2"),
              shiny::selectInput(
                inputId = "in_qq_color_by",
                label = "Color by",
                choices = list(
                  "1st p-value" = "pval1",
                  "2nd p-value" = "pval2"
                ),
                selected = "pval1"
              ),
              shiny::checkboxInput(
                inputId = "in_qq_line",
                "As line",
                value = FALSE)
            ),
            shiny::conditionalPanel(
              condition = "input.in_qq_multiple",
              shiny::uiOutput("ui_qq_ref"),
            ),
            shiny::numericInput(
              inputId = "in_qq_pthresh",
              label = "P-value threshold",
              value = 0.05, min = 0.01, max = 1, step = 0.01
            )
          ),
          shiny::column(
            width = 8,
            shiny::conditionalPanel(
              condition = "!input.in_qq_multiple",
              shinycssloaders::withSpinner(
                shiny::plotOutput(
                  outputId = "plot_deedee_qq",
                  brush = "qq_brush"
                )
              )
            ),
            shiny::conditionalPanel(
              condition = "input.in_qq_multiple",
              shinycssloaders::withSpinner(
                shiny::plotOutput("plot_deedee_qq_mult")
              )
            )
          )
        ),
        shinyBS::bsCollapse(
          shinyBS::bsCollapsePanel(
            title = "INFO - qq plot",
            shiny::includeMarkdown(
              system.file("extdata", "qq.md", package = "DeeDee")
            ),
            style = "primary"
          )
        ),
        shiny::conditionalPanel(
          condition = "!input.in_qq_multiple",
          shiny::downloadButton(
            outputId = "btn_qq_brush_download",
            label = "Download brushed genes (.xlsx)"
          ),
          shiny::tableOutput("qq_brush_info")
        )
      )
    })

    output$ui_qq_choices1 <- shiny::renderUI({
      # shiny::req(
      #   shiny::isTruthy(input$upload_de) || shiny::isTruthy(deedee_obj)
      # )
      shiny::req(shiny::isTruthy(reactive_values$dde))
      shiny::selectInput(
        inputId = "in_qq_select1",
        label = "1st data set",
        choices = names(dea(my_deedee_use()))
      )
    })

    output$ui_qq_choices2 <- shiny::renderUI({
      # shiny::req(
      #   shiny::isTruthy(input$upload_de) || shiny::isTruthy(deedee_obj)
      # )
      shiny::req(shiny::isTruthy(reactive_values$dde))
      shiny::selectInput(
        inputId = "in_qq_select2",
        label = "2nd data set",
        selected = names(dea(my_deedee_use()))[2],
        choices = names(dea(my_deedee_use()))
      )
    })

    output$ui_qq_ref <- shiny::renderUI({
      # shiny::req(
      #   shiny::isTruthy(input$upload_de) || shiny::isTruthy(deedee_obj)
      # )
      shiny::req(shiny::isTruthy(reactive_values$dde))
      shiny::selectInput(
        inputId = "in_qq_reference",
        label = "Reference",
        choices = names(dea(my_deedee_use()))
      )
    })

    output$plot_deedee_qq <- shiny::renderPlot({
      shiny::validate(
        shiny::need(
          # length(mydata()@DeeDeeList) >= 2,
          length(dea(my_deedee())) >= 2,
          message = "Please upload at least two contrasts."
        )
      )
      shiny::validate(
        shiny::need(
          # length(mydata_use()) >= 2,
          length(dea(my_deedee_use())) >= 2,
          message = "Please select at least two contrasts."
        )
      )

      # shiny::req(
      #   shiny::isTruthy(input$upload_de) || shiny::isTruthy(deedee_obj)
      # )
      shiny::req(shiny::isTruthy(reactive_values$dde))
      sel1 <- match(input$in_qq_select1, names(dea(my_deedee_use())))
      sel2 <- match(input$in_qq_select2, names(dea(my_deedee_use())))
      shiny::req(sel1)
      shiny::req(sel2)
      res <- deedee_qq(
        my_deedee_use(),
        select1 = sel1,
        select2 = sel2,
        color_by = input$in_qq_color_by,
        pthresh = input$in_qq_pthresh,
        as_line = input$in_qq_line
      )

      shiny::validate(
        shiny::need(
          !is.null(res),
          message = "No genes in your datasets. Maybe your specified p-value threshold is too low?"
        )
      )

      res
    })

    output$plot_deedee_qq_mult <- shiny::renderPlot({
      shiny::validate(
        shiny::need(
          # length(mydata()@DeeDeeList) >= 2,
          length(dea(my_deedee())) >= 2,
          message = "Please upload at least two contrasts."
        )
      )
      shiny::validate(
        shiny::need(
          # length(mydata_use()) >= 2,
          length(dea(my_deedee_use())) >= 2,
          message = "Please select at least two contrasts."
        )
      )

      # shiny::req(
      #   shiny::isTruthy(input$upload_de) || shiny::isTruthy(deedee_obj)
      # )
      shiny::req(shiny::isTruthy(reactive_values$dde))
      ref <- match(input$in_qq_reference, names(dea(my_deedee_use())))
      shiny::req(ref)
      res <- deedee_qqmult(
        my_deedee_use(),
        ref = ref,
        pthresh = input$in_qq_pthresh
      )

      shiny::validate(
        shiny::need(
          !is.null(res),
          message = "No genes in your datasets. Maybe your specified p-value threshold is too low?"
        )
      )

      res
    })

    qq_brushed <- shiny::reactive({
      shiny::req(input$qq_brush)
      x <- mydata()@DeeDeeList[[input$in_qq_select1]]$logFC
      y <- mydata()@DeeDeeList[[input$in_qq_select2]]$logFC
      pval1 <- mydata()@DeeDeeList[[input$in_qq_select2]]$pval
      pval2 <- mydata()@DeeDeeList[[input$in_qq_select2]]$pval
      names(x) <- row.names(mydata()@DeeDeeList[[input$in_qq_select2]])
      names(y) <- row.names(mydata()@DeeDeeList[[input$in_qq_select2]])

      sx_idx <- order(x)
      sy_idx <- order(y)

      sx <- x[sx_idx]
      sy <- y[sy_idx]
      pval1 <- pval1[sx_idx]
      pval2 <- pval2[sy_idx]

      lenx <- length(sx)
      leny <- length(sy)

      if (leny < lenx) {
        sx <- stats::approx(1L:lenx, sx, n = leny)$y
        pval1 <- stats::approx(1L:lenx, pval1, n = leny)$y
      }
      if (leny > lenx) {
        sy <- stats::approx(1L:leny, sy, n = lenx)$y
        pval2 <- stats::approx(1L:leny, pval2, n = lenx)$y
      }

      sx <- tibble::rownames_to_column(as.data.frame(sx))
      sx[3] <- pval1
      sy <- tibble::rownames_to_column(as.data.frame(sy))
      sy[3] <- pval2

      qq <- data.frame(
        x = sx,
        y = sy
      )

      names(qq) <- c(
        paste(input$in_qq_select1, "gene", sep = "."),
        paste(input$in_qq_select1, "logFC", sep = "."),
        paste(input$in_qq_select1, "pval", sep = "."),
        paste(input$in_qq_select2, "gene", sep = "."),
        paste(input$in_qq_select2, "logFC", sep = "."),
        paste(input$in_qq_select2, "pval", sep = ".")
      )

      temp1 <- paste(input$in_qq_select1, "logFC", sep = ".")
      temp2 <- paste(input$in_qq_select2, "logFC", sep = ".")

      shiny::brushedPoints(
        df = qq,
        brush = input$qq_brush,
        xvar = temp1,
        yvar = temp2
      )
    })

    output$qq_brush_info <- shiny::renderTable(
      {
        shiny::req(qq_brushed())
        qq_brushed()
      },
      rownames = FALSE
    )

    output$btn_qq_brush_download <- shiny::downloadHandler(
      filename = "qq_brushed_genes.xlsx",
      content = function(file) {
        shiny::req(input$qq_brush)
        bru <- qq_brushed()
        first <- data.frame(bru[1], bru[2], bru[3])
        nm1 <- unlist(strsplit(names(first)[2],
                               split = ".",
                               fixed = TRUE
        ))[1]
        names(first) <- c("rowname", "logFC", "pval")
        second <- data.frame(bru[4], bru[5], bru[6])
        nm2 <- unlist(strsplit(names(second)[2],
                               split = ".",
                               fixed = TRUE
        ))[1]
        names(second) <- c("rowname", "logFC", "pval")
        l <- list(first, second)
        names(l) <- c(nm1, nm2)
        writexl::write_xlsx(l, file, col_names = TRUE)
      }
    )


    # server - cat -------------------------------------------------------------
    output$ui_panel_cat <- renderUI({
      tagList(
        shiny::fluidRow(
          shiny::column(
            width = 4,
            shiny::selectInput(
              inputId = "in_cat_mode",
              label = "Mode",
              choices = list(
                "Up" = "up",
                "Down" = "down",
                "Both" = "both"
              ),
              selected = "up"
            ),
            shiny::numericInput(
              inputId = "in_cat_maxrank",
              label = "Max rank",
              value = 1000,
              min = 1
            ),
            shiny::uiOutput("ui_cat_choice"),
            shiny::numericInput(
              inputId = "in_cat_pthresh",
              label = "P-value threshold",
              value = 0.05, min = 0.01, max = 1, step = 0.01
            )
          ),
          shiny::column(
            width = 8,
            shinycssloaders::withSpinner(
              shiny::plotOutput("plot_deedee_cat")
            )
          )
        ),
        shinyBS::bsCollapse(
          shinyBS::bsCollapsePanel(
            title = "INFO - cat plot",
            shiny::includeMarkdown(
              system.file("extdata", "cat.md", package = "DeeDee")
            ),
            style = "primary"
          )
        )
      )
    })

    output$ui_cat_choice <- shiny::renderUI({
      # shiny::req(
      #   shiny::isTruthy(input$upload_de) || shiny::isTruthy(deedee_obj)
      # )
      shiny::req(shiny::isTruthy(reactive_values$dde))
      shiny::selectInput(
        inputId = "in_cat_ref",
        label = "Reference contrast",
        selected = names(dea(my_deedee_use()))[1],
        choices = names(dea(my_deedee_use()))
      )
    })

    output$plot_deedee_cat <- shiny::renderPlot({
      shiny::validate(
        shiny::need(
          # length(mydata()@DeeDeeList) >= 2,
          length(dea(my_deedee())) >= 2,
          message = "Please upload at least two contrasts."
        )
      )
      shiny::validate(
        shiny::need(
          # length(mydata_use()) >= 2,
          length(dea(my_deedee_use())) >= 2,
          message = "Please select at least two contrasts."
        )
      )

      # shiny::req(
      #   shiny::isTruthy(input$upload_de) || shiny::isTruthy(deedee_obj)
      # )
      shiny::req(shiny::isTruthy(reactive_values$dde))
      shiny::req(input$in_cat_maxrank)
      shiny::req(input$in_cat_ref)

      ref <- match(input$in_cat_ref, names(dea(my_deedee_use())))
      res <- deedee_cat(
        my_deedee_use(),
        ref = ref,
        maxrank = input$in_cat_maxrank,
        mode = input$in_cat_mode,
        pthresh = input$in_cat_pthresh
      )

      shiny::validate(
        shiny::need(
          !is.null(res),
          message = "No genes in your datasets. Maybe your specified p-value threshold is too low?"
        )
      )

      res
    })

    # server - summary ---------------------------------------------------------
    output$ui_panel_summary <- renderUI({
      tagList(
        shiny::numericInput(
          inputId = "sum_pthresh",
          label = "P-value threshold",
          value = 0.05, min = 0.01, max = 1, step = 0.01
        ),
        shiny::actionButton(
          inputId = "sum_button",
          label = "Create Summary"
        ),
        shinyBS::bsModal(
          id = "sum_modal",
          title = "DeeDee Summary",
          trigger = "sum_button",
          size = "large",
          shiny::fluidPage(
            shiny::downloadButton(
              outputId = "sum_download",
              label = "Download your DeeDee Summary (.html)"
            ),
            shinycssloaders::withSpinner(
              shiny::uiOutput("ui_show_html_summary")
            )
          )
        ),
        shiny::fluidRow(
          shinydashboard::tabBox(
            title = "",
            width = 12,
            id = "sum_params",
            shiny::tabPanel(
              "Scatterplot",
              shiny::uiOutput("ui_sum_scatter_choices1"),
              shiny::uiOutput("ui_sum_scatter_choices2"),
              shiny::selectInput(
                inputId = "sum_scatter_color_by",
                label = "Color by",
                choices = list(
                  "1st p-value" = "pval1",
                  "2nd p-value" = "pval2"
                ),
                selected = "pval1"
              )
            ),
            shiny::tabPanel(
              title = "Heatmap",
              shiny::numericInput(
                inputId = "sum_heatmap_show_first",
                label = "Show first",
                value = 25,
                min = 1
              ),
              shiny::checkboxInput(
                inputId = "sum_heatmap_show_gene_names",
                label = "Show gene names",
                value = FALSE
              ),
              shiny::checkboxInput(
                inputId = "sum_heatmap_show_na",
                label = "Show NA",
                value = FALSE
              ),
              shiny::selectInput(
                inputId = "sum_heatmap_dist",
                label = "Distance measure",
                choices = list(
                  "Euclidean" = "euclidean",
                  "Manhattan" = "manhattan",
                  "Pearson" = "pearson",
                  "Spearman" = "spearman"
                ),
                selected = "euclidean"
              ),
              shiny::selectInput(
                inputId = "sum_heatmap_clust",
                label = "Clustering method",
                choices = list(
                  "Single" = "single",
                  "Complete" = "complete",
                  "Average" = "average",
                  "Centroid" = "centroid"
                ),
                selected = "average"
              )
            ),
            shiny::tabPanel(
              title = "Venn Diagram",
              shiny::selectInput(
                inputId = "sum_venn_mode",
                label = "Mode",
                choices = list(
                  "Up" = "up",
                  "Down" = "down",
                  "Both" = "both"
                ),
                selected = "both"
              )
            ),
            shiny::tabPanel(
              title = "UpSet Plot",
              shiny::selectInput(
                inputId = "sum_upset_mode",
                label = "Mode",
                choices = list(
                  "Up" = "up",
                  "Down" = "down",
                  "Both" = "both"
                ),
                selected = "both"
              ),
              shiny::conditionalPanel(
                condition = "input.sum_upset_mode == 'both'",
                shiny::checkboxInput(
                  inputId = "sum_upset_colored",
                  label = "Coloring",
                  value = TRUE
                )
              ),
              shiny::numericInput(
                inputId = "sum_upset_min_setsize",
                label = "Minimum set size",
                value = 10, min = 0, step = 1
              )
            ),
            shiny::tabPanel(
              title = "Quantile-Quantile Plot",
              shiny::uiOutput("ui_sum_qq_ref"),
            ),
            shiny::tabPanel(
              title = "Concordance At the Top Plot",
              shiny::selectInput(
                inputId = "sum_cat_mode",
                label = "Mode",
                choices = list(
                  "Up" = "up",
                  "Down" = "down",
                  "Both" = "both"
                ),
                selected = "up"
              ),
              shiny::numericInput(
                inputId = "sum_cat_maxrank",
                label = "Max rank",
                value = 1000,
                min = 1
              ),
              shiny::uiOutput("ui_sum_cat_choice")
            )
          )
        ),
        shinyBS::bsCollapse(
          shinyBS::bsCollapsePanel(
            title = "INFO - summary",
            shiny::includeMarkdown(
              system.file("extdata", "summary.md", package = "DeeDee")
            ),
            style = "primary"
          )
        )
      )
    })

    output$ui_sum_scatter_choices1 <- shiny::renderUI({
      # shiny::req(
      #   shiny::isTruthy(input$upload_de) || shiny::isTruthy(deedee_obj)
      # )
      shiny::req(shiny::isTruthy(reactive_values$dde))
      shiny::selectInput(
        inputId = "sum_scatter_select1",
        label = "1st data set",
        choices = names(dea(my_deedee_use()))
      )
    })

    output$ui_sum_scatter_choices2 <- shiny::renderUI({
      # shiny::req(
      #   shiny::isTruthy(input$upload_de) || shiny::isTruthy(deedee_obj)
      # )
      shiny::req(shiny::isTruthy(reactive_values$dde))
      shiny::selectInput(
        inputId = "sum_scatter_select2",
        label = "2nd data set",
        selected = names(dea(my_deedee_use()))[2],
        choices = names(dea(my_deedee_use()))
      )
    })

    output$ui_sum_qq_ref <- shiny::renderUI({
      # shiny::req(
      #   shiny::isTruthy(input$upload_de) || shiny::isTruthy(deedee_obj)
      # )
      shiny::req(shiny::isTruthy(reactive_values$dde))
      shiny::selectInput(
        inputId = "sum_qqmult_ref",
        label = "Reference",
        selected = names(mydata_use())[1],
        choices = names(mydata_use())
      )
    })

    output$ui_sum_cat_choice <- shiny::renderUI({
      # shiny::req(
      #   shiny::isTruthy(input$upload_de) || shiny::isTruthy(deedee_obj)
      # )
      shiny::req(shiny::isTruthy(reactive_values$dde))
      shiny::selectInput(
        inputId = "sum_cat_ref",
        label = "Reference contrast",
        selected = names(mydata_use())[1],
        choices = names(mydata_use())
      )
    })

    summary <- reactive({
      shiny::validate(
        shiny::need(
          # length(mydata()@DeeDeeList) >= 2,
          length(dea(my_deedee())) >= 2,
          message = "Please upload at least two contrasts."
        )
      )
      shiny::validate(
        shiny::need(
          # length(mydata_use()) >= 2,
          length(dea(my_deedee_use())) >= 2,
          message = "Please select at least two contrasts."
        )
      )

      # TODO: update summary code from here
      shiny::req(mydata_use())

      outfile <- tempfile(fileext = ".html")

      sc_sel1 <- match(input$sum_scatter_select1, names(mydata_use()))
      if (is.null(sc_sel1)) {
        sc_sel1 <- 1
      }
      sc_sel2 <- match(input$sum_scatter_select2, names(mydata_use()))
      if (is.null(sc_sel2)) {
        sc_sel2 <- 2
      }

      qq_ref <- match(input$sum_qqmult_ref, names(mydata_use()))
      if (length(qq_ref) == 0) {
        qq_ref <- 1
      }

      cat_ref <- match(input$sum_cat_ref, names(mydata_use()))
      if (length(cat_ref) == 0) {
        cat_ref <- 1
      }

      if (input$sum_upset_mode == "both" && input$sum_upset_colored) {
        ups_mode <- "both_colored"
      } else {
        ups_mode <- input$sum_upset_mode
      }

      shiny::req(sc_sel1)
      shiny::req(sc_sel2)
      shiny::req(qq_ref)
      shiny::req(cat_ref)

      deedee_summary(
        # mydata_use(),
        dde = my_deedee_use(),
        output_path = outfile,
        overwrite = TRUE,
        pthresh = input$sum_pthresh,
        scatter_select1 = sc_sel1,
        scatter_select2 = sc_sel2,
        scatter_color_by = input$sum_scatter_color_by,
        heatmap_show_first = input$sum_heatmap_show_first,
        heatmap_show_gene_names = input$sum_heatmap_show_gene_names,
        heatmap_dist = input$sum_heatmap_dist,
        heatmap_clust = input$sum_heatmap_clust,
        heatmap_show_na = input$sum_heatmap_show_na,
        venn_mode = input$sum_venn_mode,
        upset_mode = ups_mode,
        upset_min_setsize = input$sum_upset_min_setsize,
        qqmult_ref = qq_ref,
        cat_ref = cat_ref,
        cat_maxrank = input$sum_cat_maxrank,
        cat_mode = input$sum_cat_mode,
        silent = TRUE,
        open_file = FALSE
      )

      return(outfile)
    })

    output$sum_download <- shiny::downloadHandler(
      filename = "DeeDee_Summary.html",
      content = function(file) {
        shiny::req(summary())

        file.copy(summary(), file)
      }
    )

    output$ui_show_html_summary <- shiny::renderUI({
      shiny::req(summary())

      out <- tempfile(fileext = ".html")

      xml2::write_html(rvest::html_node(xml2::read_html(summary()), "body"),
                       file = out
      )

      shiny::includeHTML(path = out)
    })



    # observers ----------------------------------------------------------------

    observeEvent(input$upload_deedee, {
      tryCatch(
        {
          showNotification(
            "Trying to upload a DeeDeeExperiment object as an RDS file...",
            id = "upload_try"
          )

          reactive_values$in_dde <- readRDS(input$upload_deedee$datapath)

          if (is(reactive_values$in_dde, "DeeDeeExperiment")) {
            reactive_values$dde <- reactive_values$in_dde
            showNotification(
              ui = "Upload complete, enjoy exploring your dataset with DeeDee! ",
              type = "message"
            )
          }
        },
        error = function(e) {
          showNotification(
            "Wrong/unexpected file format uploaded! Please check the content of your RDS file to be a DeeDeeExperiment.",
            type = "error"
          )
        }
      )
    })

  }
  # nocov end

  # launch app -----------------------------------------------------------------
  # shiny::shinyOptions(deedee_obj = deedee_obj)
  shiny::shinyApp(ui = deedee_ui, server = deedee_server)
}
