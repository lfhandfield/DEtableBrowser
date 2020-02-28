
library(shiny)
library(shinyjs)
library(DT)
library(Matrix)
library(rgl)
library(shinyRGL)
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
      
      checkboxGroupInput("showCols", "Visible Columns:",c(), selected=c())
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
      plotOutput("map"),
      webGLOutput("myWebGL")
    )
  )
)
