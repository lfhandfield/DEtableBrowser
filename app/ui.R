library(shiny)
library(shinyjs)
library(DT)
library(Matrix)
library(rgl)
library(shinyRGL)
library(ggplot2)
library(stats)
library(shinydashboard)
library(shinyBS)

ui <- dashboardPage(dashboardHeader(disable = T),
  dashboardSidebar(
    shiny::tags$head(
      shiny::tags$style(shiny::HTML("
                    body {
                    background-image: url('https://www.sanger.ac.uk/sites/default/files/wellcomesangerinstitutelogo1.png');
                    background-size: 200px;
                    background-attachment: fixed;
                    background-repeat: no-repeat;
                    background-position: center;
                    }
  
                   .multicol { 
                     column-count: 1; 
                   } 
                   
                   .main-sidebar {color: #000000; background-color: #444444;}
                  
                   
                   #outertabBox {color: #000000; background-color: #888888;}  
                   #tabBox {color: #000000; background-color: #AAAAAA;}  

                     ")
        )
    ), box(width=24,id = "outertabBox",tabBox(
      width = 12,
      id = "tabBox",
      tabPanel("Select Table",fluidRow(
        selectInput(inputId = "dataset",
          label = "Celltype and Sample calling:",
          choices = c("Fine Celltypes / Multinomial" ,  "Broad Celltypes / Multinomial", "Scmap Celltypes / Multinomial" ,  "Fine Celltypes / Clustering" ,  "Broad Celltypes / Clustering", "Scmap Celltypes / Clustering"),
          selected = "gene"
        ),
        bsTooltip(
          "dataset",
          "Choose a combinasion of Method for cell-type identification and Method for sample calling from citeseq tags",
          "right",
          options = list(container = "body")
        ),
        selectInput(
          inputId = "resfield",
          label = "Choose a result type:",
          choices = c("gene" ,  "gene_preFDR", "consensus_gene" ,  "consensus_gene_preFDR", "go", "consensus_go"),
          selected = "gene"
        )
      )),tabPanel("Table options",fluidRow(
        selectInput(inputId = "obs",
          label = "Extented Annotation:",
          choices = c()
        ),
        bsTooltip("obs", "Displays an additionnal annotation below the table when a row is selected", "right", options = list(container = "body")),
        checkboxGroupInput("showCols", "Visible Columns:",c(), selected=c())
      )),tabPanel("Heatmap options",fluidRow(
        selectInput(inputId = "comtype",
          label = "Heatmap includes:",
          choices = c("All", "Paired Comparisons", "Pooled Comparisons")),
        bsTooltip("comtype", "Filters some collumn on heatmap, containing comparison solely based on 1vs1 sample comparison and/or pooled 2vs2 samples comparisons", "right", options = list(container = "body")),
        sliderInput("nbhistcols", label = "Nb colunm for histogram:", min = 3, max = 100, value = 15),
        bsTooltip("nbhistcols", "Maximum number of collumn displayed in heatmap", "right", options = list(container = "body"))
      )),tabPanel("Filters",fluidRow(
        selectInput(inputId = "filter",
          label = "Filter Field:",
          choices = c("Comparison" ,  "Celltype", "NAME"),
          selected = "Comparison"
        ),
        selectInput(inputId = "filterchoice", label = "Filter Value:", choices = NULL),
        actionButton(inputId = "fltaddbutton",label = "Add/Remove Filter"),
        bsTooltip("fltaddbutton",
          "Add Filter, or remove filter if it is already present or if rows from displayed filter are selected",
          "right",
          options = list(container = "body")
        )
  ))))
), dashboardBody(
      shiny::tags$h4(uiOutput("selectedQuery")),
      dataTableOutput("currentfilters"),
      dataTableOutput("results"),
      uiOutput("help"),
      uiOutput("help2"),
      plotOutput("map")
      )
)
