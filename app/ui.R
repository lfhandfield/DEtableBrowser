library(shiny)
library(shinyjs)
library(DT)
library(Matrix)
library(rgl)
library(shinyRGL)
library(ggplot2)
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

#                  .col-sm-24 {color: #0##00000; background-color: #444444;}  
#                   .box {color: #000000; background-color: #444444;}  
#                   #siderbar background
#                    '.skin-blue .main-sidebar { background-color: white;}',

#                  #siderbar text color
#                  '.skin-blue .main-sidebar .sidebar{
#                    color: #red;}'
                     
# mount-farm is required to access /lustre
# Define server logic to summarize and view selected dataset ----

ui <- dashboardPage(
 dashboardHeader(disable = T),
 dashboardSidebar(  shiny::tags$head(
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

                     "))), box(width=24,id = "outertabBox",tabBox(
        width = 12,
        id = "tabBox",
        tabPanel("Select Table",fluidRow(
     
      selectInput(inputId = "dataset",
            label = "Celltype and Sample calling:",
            choices = c("Fine Celltypes / Multinomial" ,  "Broad Celltypes / Multinomial", "Scmap Celltypes / Multinomial" ,  "Fine Celltypes / Clustering" ,  "Broad Celltypes / Clustering", "Scmap Celltypes / Clustering" ),
            selected = "gene"),
      bsTooltip("dataset", "Choose a combinasion of Method for cell-type identification and Method for sample calling from citeseq tags",
         "right", options = list(container = "body")), # ,  ,trigger="Click"
      
      selectInput(inputId = "resfield",
            label = "Choose a result type:",
            choices = c("gene" ,  "gene_preFDR", "consensus_gene" ,  "consensus_gene_preFDR", "go", "consensus_go"),
            selected = "gene")
      )),tabPanel("Table options",fluidRow(
              selectInput(inputId = "obs",
             label = "Extented Annotation:",
            choices = c()),
            bsTooltip("obs", "Displays an additionnal annotation below the table when a row is selected", "right", options = list(container = "body")),
            checkboxGroupInput("showCols", "Visible Columns:",c(), selected=c())
    #  selectizeInput( "filterchoice", "Choose Value:", choices = NULL, options = list(create = TRUE)),
      )),tabPanel("Heatmap options",fluidRow(
        selectInput(inputId = "comtype",
          label = "Heatmap includes:",
          choices = c("All", "Paired Comparisons", "Pooled Comparisons")),
        bsTooltip("comtype", "Filters some collumn on heatmap, containing comparison solely based on 1vs1 sample comparison and/or pooled 2vs2 samples comparisons", "right", options = list(container = "body")),
        sliderInput("nbhistcols", label = "Nb colunm for histogram:",
                  min = 3, max = 100, value = 15),
        bsTooltip("nbhistcols", "Maximum number of collumn displayed in heatmap",
         "right", options = list(container = "body"))
      )),tabPanel("Filters",fluidRow(
          selectInput(inputId = "filter",
                  label = "Filter Field:",
                  choices = c("Comparison" ,  "Celltype", "NAME"),
                  selected = "Comparison"),
            selectInput(inputId = "filterchoice",
            label = "Filter Value:",
            choices = NULL),
  
        actionButton(inputId = "fltaddbutton",label = "Add/Remove Filter"),
        bsTooltip("fltaddbutton", "Add Filter, or remove filter if it is already present or if rows from displayed filter are selected", "right", options = list(container = "body"))
    
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
library(stats)

# mount-farm is required to access /lustre
# Define server logic to summarize and view selected dataset ----
ui <- fluidPage(
  useShinyjs(),
#  extendShinyjs(text = "shinyjs.hideFilt = document.getElementById('currentfilters').style.display='block'", functions=c("hideFilt")),
#  extendShinyjs(text = "shinyjs.showFilt = document.getElementById('currentfilters').style.display='none'" , functions=c("showFilt")),
  shiny::tags$head(
    shiny::tags$style(shiny::HTML("
                    body {
                    background-image: url('https://www.sanger.ac.uk/sites/default/files/wellcomesangerinstitutelogo1.png');
                    background-size: 200px;
                    background-attachment: fixed;
                    background-repeat: no-repeat;
                    background-position: center;
                    }"))),
  
  # App title ----
  titlePanel("Browse DE genes"),

  # Sidebar layout with a input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      # Input: Selector for choosing dataset ----
      selectInput(inputId = "dataset",
            label = "Choose a dataset:",
            choices = c("MH" ,  "MHBR", "JaJn" ,  "MHS", "MHBRS", "JaJnS"),
            selected = "gene"),
      
      selectInput(inputId = "resfield",
            label = "Choose a result type:",
            choices = c("gene" ,  "gene_preFDR", "consensus_gene" ,  "consensus_gene_preFDR", "go", "consensus_go"),
            selected = "gene"),
      selectInput(inputId = "obs",
             label = "Extented Annotation:",
            choices = c()),

      selectInput(inputId = "filter",
                  label = "Filter Field:",
                  choices = c("Comparison" ,  "Celltype", "NAME"),
                  selected = "Comparison"),
    #  selectizeInput( "filterchoice", "Choose Value:", choices = NULL, options = list(create = TRUE)),
      selectInput(inputId = "filterchoice",
            label = "Filter Value to Add:",
            choices = NULL),

    
      # Input: Numeric entry for number of obs to view ----
      actionButton(inputId = "fltaddbutton",label = "Add/Remove Filter"),
      
      checkboxGroupInput("showCols", "Visible Columns:",c(), selected=c()),
      selectInput(inputId = "comtype",
        label = "Heatmap includes:",
        choices = c("All", "Paired Comparisons", "Pooled Comparisons")),
      sliderInput("nbhistcols", label = "Nb colunm for histogram:",
                  min = 3, max = 100, value = 15)
    ),
    
    

    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Verbatim text for data summary ----
      #verbatimTextOutput("summary"),

      # Output: HTML table with requested number of observations ----
      #tableOutput("view"),
      shiny::tags$h4(uiOutput("selectedQuery")),
      dataTableOutput("currentfilters"),
      dataTableOutput("results"),
      uiOutput("help"),
      uiOutput("help2"),
      plotOutput("map")
      
      #webGLOutput("myWebGL")
    )
  )
)
