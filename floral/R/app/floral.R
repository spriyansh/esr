setwd("/home/vega/James/esr/floral/R/app")

# libs
suppressPackageStartupMessages(require(shiny))
suppressPackageStartupMessages(require(shinythemes))
suppressPackageStartupMessages(require(shinycssloaders))
suppressPackageStartupMessages(require(flowCore))
spinType <- 6

# Basic Structure

# Main UI
ui <- fluidPage(
  theme = shinytheme("darkly"),
  navbarPage(
    title = "Floral",
    tabPanel(
      title = "About",
      icon = icon("info")
    ),
    tabPanel(
      title = "Upload",
      icon = icon("upload"),

      # Download example
      downloadButton("downloadExampleData", "Download Example"),


      # Input fcs
      fileInput("fcsFile", "Select a FCS file",
        multiple = F,
        accept = c(".FCS", ".fcs", "fcs", "FCS")
      ),
      DT::dataTableOutput("summary"),
      textOutput("message")
    ),
    tabPanel(
      title = "Tabulate",
      icon = icon("table")
    ),
    tabPanel(
      title = "Visualize",
      icon = icon("chart-bar")
    ),
    tabPanel(
      title = "Contact",
      icon = icon("envelope")
    )
  )
)

# Main Server
server <- function(input, output, session) {
  fcsObjectReactive <- reactive({
    req(input$fcsFile$datapath)
    tryCatch(
      expr = {
        fcsObject <- read.FCS(input$fcsFile$datapath, transformation = FALSE)
        print(fcsObject)
        output$message <- renderText("File read sucsessfully")
        return(fcsObject)
      },
      error = function(e) {
        output$message <- renderText("Error: Unable to read the file")
      },
      warning = function(w) {
        output$message <- renderText("Warning: Caught an warning!")
      }
    )
  })


  output$summary <- DT::renderDataTable({
    if (!is.null(fcsObjectReactive())) {
      emptyDf <- data.frame()
      tryCatch(
        expr = {
          fcsObject <- fcsObjectReactive()
          return(summary(fcsObject)) # Set up in table
        }, error = function(e) {
          output$message <- renderText("Error: parse table")
          return(emptyDf)
        },
        warning = function(w) {
          output$message <- renderText(as.character(w))
          return(emptyDf)
        }
      )
    } else {
      output$message <- renderText("No file Uploaded yet! Please Upload file first")
    }
  })

  output$downloadExampleData <- downloadHandler(
    filename <- function() {
      paste("example", "fcs", sep = ".")
    },
    content <- function(file) {
      file.copy("../../testData/exampleInputFile.fcs", file)
    },
    contentType = ".fcs"
  )
}

shinyApp(ui, server)
