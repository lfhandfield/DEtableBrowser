#' Server handler for detablebrowser
#'
#' @importFrom DT renderDataTable datatable
server <- function(input, output, session) {
  last.query.state <- reactiveVal("genelist")
  debug.state <- reactiveVal(0)
  value <- reactiveVal("")
  data <- reactiveVal("")
  dataclean <- reactiveVal(0)
  curflt <- reactiveVal(data.frame(comp= c(), value=character())) 
  observe({
      if (dataclean() == 0){
       
      }else{
      updateSelectInput(session,"filter", choices = colnames(data()), selected = colnames(data())[1])

      updateSelectInput(session,"obs",choices = setdiff(colnames(data()), c("NAME", "DE", "Log2FC", "LogitAuroc", "Comparison", "Celltype", "Archtype", "TPMmean", "DEseq_Log10pval", "Wilcox_Log10pval", "DEseq_adj_Log10pval", "Wilcox_adj_Log10pval", "DESeq_basemean", "FAD_coverage", "Ctrl_coverage", "FAD_Log2FC_toEmpty", "Ctrl_Log2FC_toEmpty", "MeanLog2FC", "MeanLog2FC", "MeanLogitAuroc", "Nbgenes","ID", "Domain","Tail", "pvalue", "Test" )))
      helpstr <- "GOSlim"
      ext <- c("FullName", "GO", "GOslim", "Description", "Intersection")
      danames <- setdiff(colnames(data()), helpstr)
      updateCheckboxGroupInput(session,"showCols", choices = danames, selected = setdiff(danames, c("FullName", "GO", "GOslim", "Description")))
      }
    })
  observe({
    dataclean(0)
    progress <- Progress$new(session, min=0)
    on.exit(progress$close())
    progress$set(message = 'Loading Table...',
    detail = 'This may take a few seconds')
    data(readRDS(paste("/lustre/scratch117/cellgen/team218/lh20/results/table_",input$dataset,"_",input$resfield, ".rds", sep="")))
   #value(paste("/lustre/scratch117/cellgen/team218/lh20/results/table_NO_",input$dataset, ".rds", sep=""))

      updateSelectInput(session,"filter", choices = colnames(data()), selected = colnames(data())[1])

      updateSelectInput(session,"obs",choices = setdiff(colnames(data()), c("NAME", "DE", "Log2FC", "LogitAuroc", "Comparison", "Celltype", "Archtype", "TPMmean", "DEseq_Log10pval", "Wilcox_Log10pval", "DEseq_adj_Log10pval", "Wilcox_adj_Log10pval", "DESeq_basemean", "FAD_coverage", "Ctrl_coverage", "FAD_Log2FC_toEmpty", "Ctrl_Log2FC_toEmpty", "MeanLog2FC", "MeanLog2FC", "MeanLogitAuroc", "Nbgenes","ID", "Domain","Tail", "pvalue", "Test" )))
      helpstr <- "GOSlim"
      ext <- c("FullName", "GO", "GOslim", "Description", "Intersection")
      danames <- setdiff(colnames(data()), helpstr)
      updateCheckboxGroupInput(session,"showCols", choices = danames, selected = setdiff(danames, c("FullName", "GO", "GOslim", "Description")))

      dataclean(1)
    })
  

  observeEvent(
                input$results_rows_selected,{   
                    last.query.state("results")
                }
 )
              
  observe({
    # get all character or factor columns
    if (class(data()[[input$filter]]) == "factor"){
      updateSelectInput(session, "filterchoice", choices =unique(data()[[input$filter]]))
    }else{
      updateSelectInput(session, "filterchoice", choices =c("greater than", "less than", "equal to"))
    }
      #selected = NULL) # remove selection
  })
  
  #observe({
    # get all character or factor columns
  #  updateActionButton(session, "fltaddbutton", choices =unique(data$gene))
      #selected = NULL) # remove selection
  #})
    observe({
      if (nchar(input$filterchoice) == 0){
        shinyjs::disable("fltaddbutton")
      }else{
        shinyjs::enable("fltaddbutton")
      }
    
    })
    
  observeEvent(input$fltaddbutton, {
    if ((!is.null(input$currentfilters_rows_selected))&&(length(input$currentfilters_rows_selected) != 0)){
      tmp <- curflt()
      daflt <- (setdiff(1:nrow(tmp), input$currentfilters_rows_selected))
      value(daflt)
      if (length(daflt) == 0) curflt(data.frame(comp= c(), value=character()))
      else curflt(tmp[daflt,])
      if (nrow(curflt()) == 0) runjs("document.getElementById('currentfilters').style.display='none'")
    }else{
    #debugstate <- "button pressed"
    if (input$filter %in% rownames(curflt())){
        tmp <- curflt()
        tmp[input$filter, "value"] <- paste(sort(strsplit(paste(tmp[input$filter, "value"], as.character(input$filterchoice), sep=" ; "), "[[:space:]];[[:space:]]")[[1]]), collapse = " ; ")
        #if (nrow(curflt()) == 0) runjs("var today = new Date(); alert(today);")

      #  tmp[input$filter, "value"] <- paste(strsplit(as.character(tmp[input$filter, "value"]), "[[:space:]];[[:space:]]")[[1]], as.character(input$filterchoice), sep=" ; ", collapse = '')
        #value(class(tmp[input$filter, "value"]))
        #value( paste(strsplit(as.character(tmp[input$filter, "value"]), "[[:space:]];[[:space:]]")[[1]], "test", sep=" ; "))
        #value( paste(strsplit(as.character(tmp[input$filter, "value"]), "[[:space:]];[[:space:]]")[[1]], "test", sep=" ; "))
        curflt(tmp)
    }else{
      that <- rbind(curflt(),data.frame(row.names = c(input$filter), comp=ifelse(class(data()[[input$filter]]) == "factor", "Is among", input$filterchoice), value= as.character(input$filterchoice) )) #
        that$value <- as.character(that$value)
        runjs("document.getElementById('currentfilters').style.display='block'")
        curflt(that)
    }
    }
  })
  
  # Return the requested resfield ----
#  datasetInput <- reactive({
#    switch(input$resfield,
#           "rock" = rock,
#           "pressure" = pressure,
#           "cars" = cars,
#           "10" = runif(10))
#  })

  # Generate a summary of the resfield ----
 # output$summary <- renderPrint({
#    dataset <- datasetInput()
#    summary(dataset)
#  })

  # Show the first "n" observations ----
#  output$view <- renderTable({
#    head(datasetInput(), n = input$obs)
#  })

 # output$map <- renderPlot({
#    plot(1:3,c(dim(data$gene),input$obs))
#  })
  
 # recommended.queries <- reactive({
                #selected.genes <- gene.list()
                #selected.datasets <- input$datasetCheckbox
                #if (length(selected.genes) > 1 && length(grep(TRUE, (caseCorrect(object, selected.genes) %in% object@index$genes()))) != 0 && length(selected.datasets) != 0){ 
                #    available.queries <-  as.data.table(markerGenes(object, selected.genes, selected.datasets))
                #}else{
                #    available.queries <- data.table()
                #}
                #if(length(available.queries) != 0){
                #    available.queries$tfidf <- signif(available.queries$tfidf, digits = 5) 
                #}
 #   available.queries <- data.frame(row.names = c("a","b", "c"))
 #     available.queries$day <- c("Monday", "Tuesday", "Friday")
  #    available.queries$value <- c(5:7)
  #return(available.queries)})
              
  output$currentfilters <- renderDataTable({
      if (is.null(curflt())) datatable(data.frame())
      else if (nrow(curflt()) == 0) datatable(data.frame())
      else datatable(curflt(), options = list(pageLength = nrow(curflt()), dom = 'tip'), caption = 'Current Filters Applied:')
    })
  
  output$results <- renderDataTable({
                if (dataclean() == 0){
                  data.frame(row.names=c("Nothing"))
                }else{
                # if(all(startsWith(object@index$genes(), "chr") == T)) "Peaks" else
           
                #input$selectInput$choices <- unique(data$gene$Celltype)
                fltrow <- rep(T, nrow(data()))
                if (nrow(curflt())>0) {
                  for(i in 1:nrow(curflt())){
                    if (curflt()$comp[i] == "Is among"){
                      fltrow <- fltrow &  data()[[rownames(curflt())[i]]] %in% strsplit(curflt()$value[i] , "[[:space:]];[[:space:]]")[[1]]
                    }else if (curflt()$comp[i] == "greater than"){
                      fltrow <- fltrow &  (data()[[rownames(curflt())[i]]] > curflt()$value[i])
                    }else if (curflt()$comp[i] == "less than"){
                      fltrow <- fltrow &  (data()[[rownames(curflt())[i]]] <  curflt()$value[i])
                    }else{
                      fltrow <- fltrow &  (data()[[rownames(curflt())[i]]] == curflt()$value[i])
                    }
                  }
                }
                #fltrow
                  #value(colnames(data))
#                  value(match(input$showCols, colnames(data)))

                datatable(data()[fltrow,input$showCols], selection = 'single',
                          #options = list(columnDefs = list(list(width = '70px', targets = c(2, 3, 4)), list(width = '10px', targets = c(0))), pageLength = 5, autoWidth = TRUE, dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel')),
                         extensions = 'Scroller', colnames = input$showCols,
                        rownames = F)
                }
})
  
  #observeEvent(output$results, {
  #  if 
    
  #input$results_rows_selected
  #})
    
output$help <- renderText({
      #input$results_rows_selected
  
      ifelse(length(input$results_rows_selected) == 0, "" ,as.character(data()[input$results_rows_selected, input$obs]))
        #ifelse(last.query.state() == "genelist", 'not right', 'tight')
    })
output$help2 <- renderText({
    value()
    #if (debug.state == "button pressed") "button pressed"
    #else "something"
  })

output$myWebGL <- renderWebGL({
    points3d(1:10, 1:10, 1:10)
    axes3d()
  })

}

