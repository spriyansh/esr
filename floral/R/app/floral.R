#setwd("/home/vega/James/esr/floral/R/app")

# libs
library(BiocManager)
options(repos = BiocManager::repositories())
suppressPackageStartupMessages(require(shiny))
suppressPackageStartupMessages(require(shinythemes))
suppressPackageStartupMessages(require(shinycssloaders))
suppressPackageStartupMessages(require(reshape2))
suppressPackageStartupMessages(require(flowCore))
suppressPackageStartupMessages(require(viridis))
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(ggridges))

spinType <- 6

# Custom Function
readSlots <- function(fcsObject) {
  # Check and Load
  for (i in slotNames(fcsObject)) {
    if (!is.null(slot(fcsObject, name = i))) {
      if (i == "exprs") {
        exprsDataSlot <- as.data.frame(slot(fcsObject, name = i))
      } else if (i == "parameters") {
        parametersSlot <- slot(fcsObject, name = i)
        parametersDataSlot <- slot(parametersSlot, name = "data")
        parametersVarMetaDataSlot <- slot(parametersSlot, name = "varMetadata")
      } else if (i == "description") {
        descriptionDataSlot <- as.data.frame(slot(fcsObject, name = i))
      }
    } else {
      output$message <- renderText(paste("Floral says:", i, " slot is empty, some visuals may not work", sep = ""))
    }
  }

  return(list(
    exprsData = exprsDataSlot, parameters = parametersSlot,
    parametersData = parametersDataSlot, descriptionData = descriptionDataSlot,
    parametersVarMetaData = parametersVarMetaDataSlot
  ))
}

# Main UI
ui <- fluidPage(titlePanel("Floral"),
  theme = shinytheme("united"),
  sidebarLayout(
    sidebarPanel(
      h4("About"), h5("Visualize Flow Cytometry Data"), br(),
      h4("Download a sample input file"),
      downloadButton("downloadExampleData", "Download Example"), br(),
      fileInput("fcsFile", h4("Download a sample input file"),
        multiple = F, accept = c(".FCS", ".fcs", "fcs", "FCS")
      ), selectInput("dynamicColumns", h4("Select a column from table"), choices = "", selected = ""),
      verbatimTextOutput("message"),
      br(), hr(), br(),
      actionButton(
        inputId = "support", label = "Contact Support",
        icon = icon("envelope"),
        onclick = "window.open('https://www.metapriyansh.com/', '_blank')"
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Tabulate",
          icon = icon("table"),
          h4(textOutput("sumTable")),
          br(), DT::dataTableOutput("summary")
        ),
        tabPanel("Visualize", icon = icon("chart-bar"),
                 
                 plotOutput("hist"))
      ),
    )
  )
)

# Main Server
server <- function(input, output, session) {
  output$sumTable <- renderText("Summary Table will appear once data file is uploaded")


  fcsObjectReactive <- reactive({
    req(input$fcsFile$datapath)
    tryCatch(
      expr = {
        fcsObject <- read.FCS(input$fcsFile$datapath, transformation = FALSE)
        output$message <- renderText(paste("Floral says", "File read sucsessfully", sep = ": "))
        return(fcsObject)
      },
      error = function(e) {
        output$message <- renderText(paste("Floral says", "Unable to read the file", sep = ": "))
      },
      warning = function(w) {
        output$message <- renderText(paste("Floral says", "Caught an warning!", sep = ": "))
      }
    )
  })

  output$summary <- DT::renderDataTable({
    if (!is.null(fcsObjectReactive())) {
      emptyDf <- data.frame()
      tryCatch(
        expr = {
          fcsObject <- fcsObjectReactive()
          summaryTable <- DT::datatable(summary(fcsObject),
            options = list(paging = F, searching = F, scrollX = TRUE)
          )
          output$sumTable <- renderText("Statistical Summary Table")
          return(summaryTable) # Set up in table
        }, error = function(e) {
          output$message <- renderText(paste("Floral says", "Unable to parse file", sep = ": "))
          return(emptyDf)
        },
        warning = function(w) {
          output$message <- renderText(as.character(w))
          return(emptyDf)
        }
      )
    } else {
      output$message <- renderText(paste("Floral says", "No file Uploaded yet! Please Upload file first", sep = ": "))
    }
  })


  observe({
    updateSelectInput(session,
      inputId = "dynamicColumns", label = NULL,
      if (!is.null(fcsObjectReactive())) {
        tryCatch(
          expr = {
            fcsObject <- fcsObjectReactive()
           # choices <- colnames(fcsObject)
            
            choices <- colnames(fcsObject)[! colnames(fcsObject) %in% c("Time")]
          },
          error = function(e) {
            output$message <- renderText(paste("Floral says", "Error in parsing file", sep = ": "))
          },
          warning = function(w) {
            output$message <- renderText(paste("Floral says", w, sep = ": "))
          }
        )
      }else {
        output$message <- renderText(paste("Floral says", "No file Uploaded yet! Please Upload file first", sep = ": "))
      }
    )
  })

  fcsListReactive <- reactive({
    if (!is.null(fcsObjectReactive())) {
    tryCatch(
      expr = {
        fcsList <- readSlots(fcsObjectReactive())
        output$message <- renderText(paste("Floral says", "Slots read sucsessfully", sep = ": "))
        return(fcsList)
      },
      error = function(e) {
        output$message <- renderText(paste("Floral says", "Unable to read the Slots", sep = ": "))
      },
      warning = function(w) {
        output$message <- renderText(paste("Floral says", "Caught an warning!", sep = ": "))
      }
    )}else{
        output$message <- renderText(paste("Floral says", "No file Uploaded yet! Please Upload file first", sep = ": "))
    }
  })
  
  
  output$hist <- renderPlot({
    exprsData <- fcsListReactive()$exprsData
    exprsData <- as.data.frame(melt(exprsData, id = "Time"))
    
    exprsDataSubset <- exprsData[exprsData$variable == input$dynamicColumns, ]
    
    K <- round(1 + 3.322 * log(nrow(exprsDataSubset)))
    R <- round(max(exprsDataSubset$value)) - round(min(exprsDataSubset$value))
    B <- round(R / K)
    
    histPlot <- ggplot(data = exprsDataSubset, aes(x = value)) +
      geom_histogram(binwidth = B, fill = "#69b3a2", color = "#e9ecef", alpha = 0.9) +
      ggtitle("FSC_A") +
      theme(
        plot.title = element_text(size = 15)
      )
    histPlot
  })
  
  output$downloadExampleData <- downloadHandler(
    filename <- function() {
      paste("example", "fcs", sep = ".")
    },
    content <- function(file) {
      file.copy("exampleInputFile.fcs", file)
      output$message <- renderText(paste("Floral says", "You just download the example file", sep = ": "))
    },
    contentType = ".fcs"
  )
}

shinyApp(ui, server)
