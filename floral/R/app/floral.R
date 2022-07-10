# setwd("/home/vega/James/esr/floral/R/app")
# styler::style_file("floral.R")

# libs
library(BiocManager)
options(repos = BiocManager::repositories())
suppressPackageStartupMessages(require(shiny))
suppressPackageStartupMessages(require(shinythemes))
suppressPackageStartupMessages(require(shinycssloaders))
suppressPackageStartupMessages(require(reshape2))
suppressPackageStartupMessages(require(flowCore))
suppressPackageStartupMessages(require(ggcyto))
suppressPackageStartupMessages(require(flowStats))
suppressPackageStartupMessages(require(viridis))
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(flowViz))

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
      output$message <- renderText(paste("Floral says:", i, " slot is empty, some visuals may not work(e1)", sep = ""))
    }
  }

  return(list(
    exprsData = exprsDataSlot, parameters = parametersSlot,
    parametersData = parametersDataSlot, descriptionData = descriptionDataSlot,
    parametersVarMetaData = parametersVarMetaDataSlot
  ))
}

# Nothing computed
nothingComputed <- function() {
  suppressPackageStartupMessages(require(ggplot2))

  df <- data.frame()
  p <- ggplot(df) +
    geom_point() +
    xlim(0, 10) +
    ylim(0, 10) +
    annotate("text", x = 3.9, y = 5.0, size = 30, col = "#f58a53", label = "(") +
    annotate("text", x = 5, y = 5.6, size = 10, col = "#f58a53", label = "o  o") +
    annotate("text", x = 6.1, y = 5.0, size = 30, col = "#f58a53", label = ")") +
    annotate("text", x = 5, y = 5.1, size = 10, col = "#f58a53", label = "|") +
    geom_segment(aes(x = 4.7, xend = 5.3, y = 4.4, yend = 4.4), size = 2, color = "#f58a53") +
    annotate("text", x = 5, y = 3, size = 6, col = "#f58a53", label = "Nothing to plot!") +
    theme_bw()
  return(p)
}

# Main UI
ui <- fluidPage(titlePanel("Floral"),
  theme = shinytheme("united"),
  sidebarLayout(
    sidebarPanel(
      h4("About"), h5("Visualize Flow Cytometry Data"), br(),
      downloadButton("downloadExampleData", "Download Example File"), br(), br(),
      fileInput("fcsFile", "Upload a fcs file",
        multiple = F, accept = c(".FCS", ".fcs", "fcs", "FCS")
      ), selectInput("dynamicColumns", "Select a column for histogram and density plot", choices = "", selected = ""),
      sliderInput("histBins", "Adjust Bins for Histogram",
        min = -100, max = 100, value = 0
      ),
      sliderInput("hexBin", "Adjust Bins for Flow Cell Chart",
        min = -100, max = 100, value = 0
      ),
      h4("Adjust Parameters for Flow plot"),
      selectInput("xCoord", h5("Select for x-axis"), choices = "", selected = ""),
      selectInput("yCoord", h5("Select for y-axis"), choices = "", selected = ""),
      sliderInput("scale", "Adjust scale for flow Cell Chart",
        min = 1, max = 10, value = 2
      ),
      hr(),
      verbatimTextOutput("message"),
      actionButton(
        inputId = "support", label = "Contact Support",
        icon = icon("envelope"),
        onclick = "window.open('https://www.metapriyansh.com/', '_blank')"
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Visuals",
          icon = icon("chart-bar"),
          h4(textOutput("sumTable")),
          br(), fluidRow(column(12, DT::dataTableOutput("summary"))),br(),br(),
          fluidRow(
            column(4, plotOutput("hist") %>% withSpinner(type = spinType, color = "#e65320")),
            column(4, plotOutput("dens") %>% withSpinner(type = spinType, color = "#e65320")),
            column(4, plotOutput("hex") %>% withSpinner(type = spinType, color = "#e65320"))
          )
        )
      )
    )
  )
)

# Main Server
server <- function(input, output, session) {
  output$sumTable <- renderText("Plots and tables will appear once correct data file is uploaded")


  fcsObjectReactive <- reactive({
    req(input$fcsFile$datapath)
    tryCatch(
      expr = {
        fcsObject <- read.FCS(input$fcsFile$datapath, transformation = FALSE)
        output$message <- renderText(paste("Floral says", "File read sucsessfully", sep = ": "))
        return(fcsObject)
      },
      error = function(e) {
        output$message <- renderText(paste("Floral says", "Unable to read the file (e2)", sep = ": "))
      },
      warning = function(w) {
        output$message <- renderText(paste("Floral says", "Caught an warning! (e3)", sep = ": "))
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
          output$message <- renderText(paste("Floral says", "Unable to parse file (e4)", sep = ": "))
          return(NULL)
        },
        warning = function(w) {
          output$message <- renderText(as.character(w))
          return(NULL)
        }
      )
    } else {
      output$message <- renderText(paste("Floral says", "No file Uploaded yet! Please Upload file first (e5)", sep = ": "))
      return(NULL)
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

            choices <- colnames(fcsObject)[!colnames(fcsObject) %in% c("Time")]
          },
          error = function(e) {
            output$message <- renderText(paste("Floral says", "Error in parsing file (e6)", sep = ": "))
            return(NULL)
          },
          warning = function(w) {
            output$message <- renderText(paste("Floral says", w, sep = ": "))
            return(NULL)
          }
        )
      } else {
        output$message <- renderText(paste("Floral says", "No file Uploaded yet! Please Upload file first (e7)", sep = ": "))
        return(NULL)
      }
    )
    updateSelectInput(session,
      inputId = "xCoord", label = NULL,
      if (!is.null(fcsObjectReactive())) {
        tryCatch(
          expr = {
            fcsObject <- fcsObjectReactive()
            choices <- colnames(fcsObject)[!colnames(fcsObject) %in% c("Time")]
          },
          error = function(e) {
            output$message <- renderText(paste("Floral says", "Error in parsing file (e6)", sep = ": "))
            return(NULL)
          },
          warning = function(w) {
            output$message <- renderText(paste("Floral says", w, sep = ": "))
            return(NULL)
          }
        )
      } else {
        output$message <- renderText(paste("Floral says", "No file Uploaded yet! Please Upload file first (e7)", sep = ": "))
        return(NULL)
      }
    )
    updateSelectInput(session,
      inputId = "yCoord", label = NULL,
      if (!is.null(fcsObjectReactive())) {
        tryCatch(
          expr = {
            fcsObject <- fcsObjectReactive()
            choices <- rev(colnames(fcsObject)[!colnames(fcsObject) %in% c("Time")])
          },
          error = function(e) {
            output$message <- renderText(paste("Floral says", "Error in parsing file (e6)", sep = ": "))
            return(NULL)
          },
          warning = function(w) {
            output$message <- renderText(paste("Floral says", w, sep = ": "))
            return(NULL)
          }
        )
      } else {
        output$message <- renderText(paste("Floral says", "No file Uploaded yet! Please Upload file first (e7)", sep = ": "))
        return(NULL)
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
          output$message <- renderText(paste("Floral says", "Incorrect file uploaded (e8)", sep = ": "))
          return(NULL)
        },
        warning = function(w) {
          output$message <- renderText(paste("Floral says", "Caught an warning! (e9)", sep = ": "))
          return(NULL)
        }
      )
    } else {
      output$message <- renderText(paste("Floral says", "No file Uploaded yet! Please Upload file first (e10)", sep = ": "))
      return(NULL)
    }
  })

  output$hist <- renderPlot({
    if (!is.null(fcsListReactive())) {
      tryCatch(
        expr = {
          exprsData <- fcsListReactive()$exprsData
          exprsData <- as.data.frame(melt(exprsData, id = "Time"))
          exprsDataSubset <- exprsData[exprsData$variable == input$dynamicColumns, ]
          K <- round(1 + 3.322 * log(nrow(exprsDataSubset)))
          R <- round(max(exprsDataSubset$value)) - round(min(exprsDataSubset$value))
          B <- round(R / K)

          B <- B + input$histBins

          histPlot <- ggplot(data = exprsDataSubset, aes(x = value)) +
            geom_histogram(binwidth = B, fill = "#69b3a2", color = "#e9ecef", alpha = 1) +
            ggtitle(input$dynamicColumns) +
            theme(
              plot.title = element_text(size = 15),
              panel.background = element_blank(),
              panel.grid.major = element_line(size = rel(0.3), linetype = 2, colour = "#e65320"),
              panel.grid.minor = element_line(size = rel(0.1), linetype = 1, colour = "#e65320")
            )
          histPlot
        }, error = function(e) {
          output$message <- renderText(paste("Floral says", "Please report error code as (e16)", sep = ": "))
          n <- nothingComputed()
          n
        }, warning = function(w) {
          output$message <- renderText(paste("Floral says", "Please report error code as (e17)", sep = ": "))
          n <- nothingComputed()
          n
        }
      )
    } else {
      n <- nothingComputed()
      n
    }
  })

  output$dens <- renderPlot({
    if (!is.null(fcsListReactive())) {
      tryCatch(
        expr = {
          exprsData <- fcsListReactive()$exprsData
          exprsData <- as.data.frame(melt(exprsData, id = "Time"))
          exprsDataSubset <- exprsData[exprsData$variable == input$dynamicColumns, ]
          densPlot <- ggplot(data = exprsDataSubset, aes(x = value)) +
            geom_density(fill = "#69b3a2", color = "#e9ecef", alpha = 1) +
            ggtitle(input$dynamicColumns) +
            theme(
              plot.title = element_text(size = 15),
              panel.background = element_blank(),
              panel.grid.major = element_line(size = rel(0.3), linetype = 2, colour = "#e65320"),
              panel.grid.minor = element_line(size = rel(0.1), linetype = 1, colour = "#e65320")
            )
          densPlot
        }, error = function(e) {
          output$message <- renderText(paste("Floral says", "Please report error code as (e14)", sep = ": "))
          n <- nothingComputed()
          n
        }, warning = function(w) {
          output$message <- renderText(paste("Floral says", "Please report error code as (e15)", sep = ": "))
          n <- nothingComputed()
          n
        }
      )
    } else {
      n <- nothingComputed()
      n
    }
  })


  output$hex <- renderPlot({
    if (!is.null(fcsObjectReactive())) {
      
      tryCatch(expr = {
      if (input$xCoord != input$yCoord) {
        fcsObject <- fcsObjectReactive()
        B <- 200 + input$hexBin
        flowChart <- ggcyto(fcsObject, aes_string(
          x = paste0("`", input$xCoord, "`"),
          y = paste0("`", input$yCoord, "`")
        )) +
          geom_hex(bins = B, alpha = 1, size = 2) +
          theme(
            plot.title = element_text(size = 15),
            strip.background = element_blank(),
            strip.text = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_line(size = rel(0.3), linetype = 2, colour = "#e65320"),
            panel.grid.minor = element_line(size = rel(0.1), linetype = 1, colour = "#e65320")
          )
        lg <- flowStats::lymphGate(fcsObject,
          channels = c(input$xCoord, input$yCoord),
          scale = input$scale
        )
        fres <- filter(fcsObject, lg)
        flowChart <- flowChart + geom_gate(lg) + geom_stats()
        output$message <- renderText(paste("Floral says: ", "Plotted for", input$xCoord, "vs", input$yCoord, sep = ""))

        return(flowChart)
      
        
      } else {
        output$message <- renderText(paste("Floral says", "Both x-coord and y-coord cannot be same! (e11)", sep = ": "))
        nothingComputed()
      }},
      error = function(e) {
        output$message <- renderText(paste("Floral says", "Please report error code as (e12)", sep = ": "))
        n <- nothingComputed()
        n
      }, warning = function(w) {
        output$message <- renderText(paste("Floral says", "Please report error code as (e13)", sep = ": "))
        n <- nothingComputed()
        n
      }
      )
    } else {
      nothingComputed()
    }
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
