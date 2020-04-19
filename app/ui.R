# 
# . overflow-x: scroll;

library(shiny)
library(shinyjs)
library(DT)
library(Matrix)
#library(rgl)
#library(shinyRGL)
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
                        tabPanel("Simple Browsing Tab",fluidRow(
                          selectInput(inputId = "simplecondition",
                                      label = "Condition Investigated:",
                                      choices = c("APP V717I in Neurons", "APP V717I in Microglia", "PSEN1 M146I in Neurons", "PSEN1 M146I in Microglia", "PSEN1 Intron4 mutation in Neurons", "PSEN1 Intron4 mutation in Microglia", "LPS", "TREM2 knock-out Microglia")),
                          bsTooltip("simplecondition", "Selects the characteristic that discriminates the condition of test samples organoids to control samples organoids.", "right", options = list(container = "body")),
                          selectInput(inputId = "simplecelltype",
                                      label = "Effect observed in:", # , "Neuron Precursor"
                                      choices = c("Microglia","Neurons","Neuron and Microglia","All")),
                          bsTooltip("simplecelltype", "Selects the cell-type within organoids in which gene are tested for differential expression, using both DEseq2 and Wilcoxon tests.", "right", options = list(container = "body")),
                          selectInput(inputId = "simpledetype",
                                      label = "Effect:",
                                      choices = c("Significant for DEseq2","Significant for Wilcox test","Higher Expression in Disease", "Lower Expression in Disease", "Upregulated pathways in Disease", "Downregulated pathways in Disease")),
                          bsTooltip("simpledetype", "Selects the filtering and ordering criterion for genes detected as differentilly expressed.", "right", options = list(container = "body")),
                          selectInput(inputId = "simpleextra",
                                      label = "Contextual Display",
                                      choices = c("Heatmap with other conditions", "Heatmap with other all celltypes", "Volcano Plot of DEseq", "Tsne Overlay of Wilcox test")),
                          bsTooltip("simpleextra", "Show Fold Change for more celltypes and/or conditions. Genes on the heatmap interactively match what the table reports.", "right", options = list(container = "body")),
                          actionButton(inputId = "simplebutton",label = "Execute Query"),
                          bsTooltip("simplebutton", "Execute the query matching fields from this tab. Other tabs can be used to refine the results displayed.", "right", options = list(container = "body"))
                        )),tabPanel("Tables",fluidRow(   
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
                            choices = c("genes (within batches)" , "genes (consensus)" , "pathways/annotations (within batches)", "pathways/annotations (consensus)"),
                            selected = "genes (consensus)" 
                          )
                        )),tabPanel("& Columns",fluidRow(
                          selectInput(inputId = "obs",
                                      label = "Extended Annotation:",
                                      choices = c()
                          ),
                          checkboxGroupInput("showCols", "Visible Columns:",c(), selected=c()),
                          bsTooltip("obs", "Displays an additionnal annotation below the table when a row is selected", "right", options = list(container = "body"))
                        )),tabPanel("Heatmap options",fluidRow(
                          selectInput(inputId = "comtype", label = "Comparison Types:",
                                      choices = c("All", "1-to-1 Comparisons", "Pooled Comparisons"), selected = "1-to-1 Comparisons"),
                          bsTooltip("comtype", "Filters columns on heatmap based on whether 1-to-1 sample comparison and/or pooled 2vs2 samples comparisons are performed", "right", options = list(container = "body")),
                          selectInput(inputId = "samexcl", label = "Sample Included:",
                                      choices = c("Match ConsensusGroup","Match Comparison","point-mutation conditions", "other disease conditions", "other neutral conditions","Include All"), selected = "point-mutation conditions"),
                          bsTooltip("samexcl", "Filters columns corresponding to condition tested. May match filters used for the table.", "right", options = list(container = "body")),
                          selectInput(inputId = "ctpexcl", label = "Celltype Included:",
                                      choices = c("Microglia", "Neurons", "Microglia and Neurons","Match Filters","All"), selected = "Match Filters"),
                          bsTooltip("ctpexcl", "Filters columns corresponding to celltype from each organoid in which fold change expression is reported on the heatmap", "right", options = list(container = "body")),
                          selectInput(inputId = "clusterheat", label = "rows & cols clustering:",
                                      choices = c("No ordering", "Cluster Genes", "Cluster Columns","Cluster Both"), selected = "No ordering"),
                          sliderInput("nbhistcols", label = "Nb column for histogram:", min = 3, max = 100, value = 30),
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
                                    "Add Filter, or remove filter if it is already present or if any rows from filter table are selected",
                                    "right",
                                    options = list(container = "body")
                          )
                        ), id= "heatmapopttab"),
                        selectInput(
                          inputId = "contextfield",
                          label = "Contextual display:",
                          choices = c("Heatmap" , "Volcano Plot" , "Tsne Overlay"),
                          selected = "genes (consensus)" 
                        ),downloadButton("downloadData", "Download Table")
                        ))
                    ), dashboardBody(
                      shiny::tags$h4(uiOutput("selectedQuery")),
                      DT::dataTableOutput("currentfilters"),
                      DT::dataTableOutput("results"),
                      uiOutput("help"),
                      uiOutput("help2"),
                      plotOutput("map")
                    )
)
