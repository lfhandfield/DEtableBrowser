
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
                  
                   .progress-text {
                            text-align: center;
                            font-size: 18px;
                            font-style: italic; 
                   }
                   .shiny-notification {position: fixed; top: 0% ;left: 50%;}

                   #outertabBox {color: #000000; background-color: #888888;}  
                   #tabBox {color: #000000; background-color: #AAAAAA;}  

                     ")
        )
    ), shinyjs::useShinyjs(), box(width=24,id = "outertabBox",tabBox(  # Table Selection Tab
      width = 12,
      id = "tabBox",
      tabPanel("Simple Browse",fluidRow(
        selectInput(inputId = "simplecondition",
          label = "Condition Investigated:",
          choices = c("APP V717I in Neurons", "APP V717I in Microglia", "PSEN1 M146I in Neurons", "PSEN1 M146I in Microglia", "PSEN1 Intron4 mutation in Neurons", "PSEN1 Intron4 mutation in Microglia", "LPS", "TREM2 knock-out Microglia")),
        selectInput(inputId = "simplecelltype",
          label = "Effect observed in:",
          choices = c("Microglia","Neurons", "Neuron Precursor","Neuron and Microglia","All")),
        selectInput(inputId = "simpledetype",
          label = "Effect:",
          choices = c("Differentially Expressed","Higher Expression in Disease", "Lower Expression in Disease")),
      #  bsTooltip("simplecondition", "Filters some collumn on heatmap, containing comparison solely based on 1vs1 sample comparison and/or pooled 2vs2 samples comparisons", "right", options = list(container = "body")),
        actionButton(inputId = "simplebutton",label = "Execute Query")
      )),tabPanel("Select Table",fluidRow(   
        selectInput(inputId = "dataset",
          label = "Celltype and Sample calling:",
          choices = c("Fine Celltypes / Multinomial" ,  "Broad Celltypes / Multinomial", "Scmap Celltypes / Multinomial" ,  "Fine Celltypes / Clustering" ,  "Broad Celltypes / Clustering", "Scmap Celltypes / Clustering"),
          selected = "Broad Celltypes / Multinomial"
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
        checkboxGroupInput("showCols", "Visible Columns:",c(), selected=c()),
        bsTooltip("obs", "Displays an additionnal annotation below the table when a row is selected", "right", options = list(container = "body"))
      )),tabPanel("Heatmap options",fluidRow(
        selectInput(inputId = "comtype",
          label = "Comparison Types:",
          choices = c("All", "Paired Comparisons", "Pooled Comparisons"), selected = "Paired Comparisons"),
        bsTooltip("comtype", "Filters some collumn on heatmap, containing comparison solely based on 1vs1 sample comparison and/or pooled 2vs2 samples comparisons", "right", options = list(container = "body")),
        selectInput(inputId = "samexcl",
          label = "Sample Exclude:",
          choices = c("non-Indel mutation", "Negatives","Include All"), selected = "Indel mutation"),
        bsTooltip("comtype", "Filters some collumn on heatmap, containing comparison solely based on 1vs1 sample comparison and/or pooled 2vs2 samples comparisons", "right", options = list(container = "body")),
        selectInput(inputId = "ctpexcl",
          label = "Celltype Exclude:",
          choices = c("All but Microglia", "All but Neurons", "All but Microglia or Neurons","Match Filters","Include All"), selected = "Match Filters"),
        bsTooltip("comtype", "Filters some collumn on heatmap, containing comparison solely based on 1vs1 sample comparison and/or pooled 2vs2 samples comparisons", "right", options = list(container = "body")),
        sliderInput("nbhistcols", label = "Nb colunm for histogram:", min = 3, max = 100, value = 15),
        bsTooltip("nbhistcols", "Maximum number of collumn displayed in heatmap", "right", options = list(container = "body"))
      )),tabPanel("Filters",fluidRow( # Filter Addition Tab
        selectInput(inputId = "filter",
          label = "Field:",
          choices = c("Comparison" ,  "Celltype", "Gene"),
          selected = "Comparison"
        ),
        selectInput(inputId = "filterchoice", label = "Criterion:", choices = NULL),
        numericInput(inputId = "filtervalue", label = "Value:", value=0),
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
