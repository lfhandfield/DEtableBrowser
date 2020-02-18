
server <- function(input, output, session) {
  last.query.state <- reactiveVal("genelist")
  debug.state <- reactiveVal(0)
  value <- reactiveVal("")
  curflt <- reactiveVal(data.frame(comp= c(), value=character())) 
  observe({

      updateSelectInput(session,"filter", choices = colnames(data[[input$resfield]]), selected = colnames(data[[input$resfield]])[1])

      updateSelectInput(session,"obs",choices = setdiff(colnames(data[[input$resfield]]), c("NAME", "DE", "Log2FC", "LogitAuroc", "Comparison", "Celltype", "Archtype", "TPMmean", "DEseq_Log10pval", "Wilcox_Log10pval", "DEseq_adj_Log10pval", "Wilcox_adj_Log10pval", "DESeq_basemean", "FAD_coverage", "Ctrl_coverage", "FAD_Log2FC_toEmpty", "Ctrl_Log2FC_toEmpty", "MeanLog2FC", "MeanLog2FC", "MeanLogitAuroc", "Nbgenes","ID", "Domain","Tail", "pvalue", "Test" )))
      helpstr <- "GOSlim"
      ext <- c("FullName", "GO", "GOslim", "Description", "Intersection")
      danames <- setdiff(colnames(data[[input$resfield]]), helpstr)
      updateCheckboxGroupInput(session,"showCols", choices = danames, selected = setdiff(danames, c("FullName", "GO", "GOslim", "Description")))
  })
  
  #observe({
    #data <- readRDS(paste("table_NO_",input$dataset, ".rds", sep=""))
   # })
  

  observeEvent(
                input$results_rows_selected,{   
                    last.query.state("results")
                }
 )
              
  observe({
    # get all character or factor columns
    if (class(data[[input$resfield]][[input$filter]]) == "factor"){
      updateSelectInput(session, "filterchoice", choices =unique(data[[input$resfield]][[input$filter]]))
    }else{
      updateSelectInput(session, "filterchoice", choices =c("greater than", "less than", "equal to"))
    }
      #selected = NULL) # remove selection
  })
  
  #observe({
    # get all character or factor columns
  #  updateActionButton(session, "fltaddbutton", choices =unique(data$gene[[input$resfield]]))
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
    if (length(input$currentfilters_rows_selected) != 0){
      tmp <- curflt()
      daflt <- (setdiff(1:nrow(tmp), input$currentfilters_rows_selected))
      value(daflt)
      if (length(daflt) == 0) curflt(data.frame(comp= c(), value=character()))
      else curflt(tmp[daflt,])
    }else{
    #debugstate <- "button pressed"
    if (input$filter %in% rownames(curflt())){
        tmp <- curflt()
        tmp[input$filter, "value"] <- paste(sort(strsplit(paste(tmp[input$filter, "value"], as.character(input$filterchoice), sep=" ; "), "[[:space:]];[[:space:]]")[[1]]), collapse = " ; ")
        
      #  tmp[input$filter, "value"] <- paste(strsplit(as.character(tmp[input$filter, "value"]), "[[:space:]];[[:space:]]")[[1]], as.character(input$filterchoice), sep=" ; ", collapse = '')
        #value(class(tmp[input$filter, "value"]))
        #value( paste(strsplit(as.character(tmp[input$filter, "value"]), "[[:space:]];[[:space:]]")[[1]], "test", sep=" ; "))
        #value( paste(strsplit(as.character(tmp[input$filter, "value"]), "[[:space:]];[[:space:]]")[[1]], "test", sep=" ; "))
        curflt(tmp)
    }else{
      that <- rbind(curflt(),data.frame(row.names = c(input$filter), comp=ifelse(class(data[[input$resfield]][[input$filter]]) == "factor", "Is among", input$filterchoice), value= as.character(input$filterchoice) )) #
        that$value <- as.character(that$value)
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
  
  recommended.queries <- reactive({
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
    available.queries <- data.frame(row.names = c("a","b", "c"))
      available.queries$day <- c("Monday", "Tuesday", "Friday")
      available.queries$value <- c(5:7)
  return(available.queries)})
              
  output$currentfilters <- renderDataTable({
     curflt()
    })
  
  output$results <- renderDataTable({
                
                if(is.null(recommended.queries())) return(data.table(row.names=c("Nothing")))
                # if(all(startsWith(object@index$genes(), "chr") == T)) "Peaks" else
                selected.colnames = input$showCols            
                #input$selectInput$choices <- unique(data$gene$Celltype)
                
                fltrow <- rep(T, nrow(data[[input$resfield]]))
                if (nrow(curflt())>0) {
                  for(i in 1:nrow(curflt())){
                    if (curflt()$comp[i] == "Is among"){
                      fltrow <- fltrow &  data[[input$resfield]][[rownames(curflt())[i]]] %in% strsplit(curflt()$value[i] , "[[:space:]];[[:space:]]")[[1]]
                    }else if (curflt()$comp[i] == "greater than"){
                      fltrow <- fltrow &  (data[[input$resfield]][[rownames(curflt())[i]]] > curflt()$value[i])
                    }else if (curflt()$comp[i] == "less than"){
                      fltrow <- fltrow &  (data[[input$resfield]][[rownames(curflt())[i]]] <  curflt()$value[i])
                    }else{
                      fltrow <- fltrow &  (data[[input$resfield]][[rownames(curflt())[i]]] == curflt()$value[i])
                    }
                  }
                }
                
                datatable(data[[input$resfield]][fltrow ,selected.colnames], selection = 'single',
                          #options = list(columnDefs = list(list(width = '70px', targets = c(2, 3, 4)), list(width = '10px', targets = c(0))), pageLength = 5, autoWidth = TRUE, dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel')),
                          extensions = 'Scroller', colnames = selected.colnames,
                          rownames = F)
})
  
  #observeEvent(output$results, {
  #  if 
    
  #input$results_rows_selected
  #})
    
output$help <- renderText({
      #input$results_rows_selected
  
      ifelse(length(input$results_rows_selected) == 0, "" ,as.character(data[[input$resfield]][input$results_rows_selected, input$obs]))
        #ifelse(last.query.state() == "genelist", 'not right', 'tight')
    })
  output$help2 <- renderText({
      value()
      #if (debug.state == "button pressed") "button pressed"
      #else "something"
    })
}
