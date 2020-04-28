
source("render.R")
options(DT.fillContainer = FALSE)
options(DT.autoHideNavigation = FALSE)
options(shiny.trace=TRUE)
#' Server handler for detablebrowser
#'
#' @importFrom DT renderDataTable datatable
server <- function(input, output, session) {
  last.query.state <- reactiveVal("genelist")
  debug.state <- reactiveVal(0); dataclean <- reactiveVal(c())
  value <- reactiveVal("");      data <- reactiveVal(""); overlay <- reactiveVal("")
  mat <- reactiveVal("");        simplesort <- reactiveVal("")
  plotgenes <- reactiveVal(c(""));
  shownrows <- reactiveVal(c(""));
  sortcol <- reactiveVal(c());
  curflt <- reactiveVal(data.frame(criterion= c(), value=character())) 
  filtrow <- reactiveVal(c(T))

  
  observe({
      if (length(dataclean()) == 0){
       
      }else{
      updateSelectInput(session,"filter", choices = colnames(data()), selected = colnames(data())[1])

      updateSelectInput(session,"obs",choices = setdiff(colnames(data()), c("Gene", "DE", "Log2FC", "LogitAuroc", "Comparison", "Celltype", "Archtype", "TPMmean", "DEseq_Log10pval", "Wilcox_Log10pval", "DEseq_adj_Log10pval", "Wilcox_adj_Log10pval", "DESeq_basemean", "FAD_coverage", "Ctrl_coverage", "FAD_Log2FC_toEmpty", "Ctrl_Log2FC_toEmpty", "MeanLog2FC", "MeanLog2FC", "MeanLogitAuroc", "Nbgenes","ID", "Domain","Tail", "pvalue", "Test" )))
      ext <- c("FullName", "GO", "GOslim", "Description", "Intersection")
      danames <- colnames(data())
      
      
      # set tablecol filter, with default values
      defaultselect <- setdiff(danames, c("GO", "GOslim", "GOSLIM", "FAD_coverage", "Ctrl_coverage", "Description", "DESCRIPTION","Fullname", "FULLNAME", "Intersection", "sample_testIds", "sample_ctrlIds", "ALIAS", "Alias", "biotype", "DE_concat", "Wilcox_adj_Log10pval"))
      if (grepl("genes", input$resfield)) {
        defaultselect <- setdiff(defaultselect, c("DEseq_Log10pval", "Wilcox_Log10pval"))
        if (grepl("consensus",input$resfield)) defaultselect <- setdiff(defaultselect, c("DE"))
      }else defaultselect <- setdiff(defaultselect, c("DEseq_adj_Log10pval", "Wilcox_adj_Log10pval"))
      logjs(defaultselect);
      logjs(setdiff(danames,dataclean()));
      defaultselect <- c(intersect(defaultselect, setdiff(danames,dataclean())), intersect(input$showCols,danames))
      logjs(defaultselect);
      updateCheckboxGroupInput(session,"showCols", choices = danames, selected = defaultselect)
      
      }
      
    })

  observe({
    if (input$tabContext == "Heatmap"){
      shinyjs::showElement("comtype");
      shinyjs::showElement("samexcl");
      shinyjs::showElement("ctpexcl");
      shinyjs::showElement("nbhistcols");
      shinyjs::showElement("clusterheat");
    }else{
      shinyjs::hideElement("comtype");
      shinyjs::hideElement("samexcl");
      shinyjs::hideElement("ctpexcl");
      shinyjs::hideElement("nbhistcols");
      shinyjs::hideElement("clusterheat");
    }
    
    if (input$tabContext == "Tsne Overlay"){
      if (class(overlay()) == "character"){
        shinyjs::disable("resfield"); shinyjs::disable("dataset"); shinyjs::disable("simplebutton");  shinyjs::disable("downloadData")
        progress <- Progress$new(session, min=0)
        on.exit(progress$close())
        progress$set(message = 'Loading Overlay data...',
                     detail = 'This may take a few seconds')    
        dastr <- switch(input$dataset, "FINE", "Fine Celltypes / Clustering"  = "FINES" , "Broad Celltypes / Clustering" = "FINES", "Scmap Celltypes / Clustering" = "FINES")
        logjs(paste("opening /lustre/scratch117/cellgen/team218/lh20/SnakeFolderEv4/shinydata/overlay_NO_",dastr, ".rds", sep=""))
        overlay(readRDS(paste("/lustre/scratch117/cellgen/team218/lh20/SnakeFolderEv4/shinydata/overlay_NO_",dastr, ".rds", sep="")))
        shinyjs::enable("resfield"); shinyjs::enable("dataset") ; shinyjs::enable("simplebutton"); shinyjs::enable("downloadData")
      }
      shinyjs::showElement("atlas");
    }else{
      shinyjs::hideElement("atlas");
    }
    })
  
  observe({ # Load the current table data
    
    shinyjs::disable("resfield"); shinyjs::disable("dataset"); shinyjs::disable("simplebutton");  shinyjs::disable("downloadData")
    progress <- Progress$new(session, min=0)
    on.exit(progress$close())
    progress$set(message = 'Loading Table...',
    detail = 'This may take a few seconds')
    dastr <- switch(input$dataset, "MHBR", "Fine Celltypes / Multinomial"  = "MH" , "Broad Celltypes / Multinomial" = "MHBR", "Scmap Celltypes / Multinomial" = "JaJn", "Fine Celltypes / Clustering"  = "MHS" , "Broad Celltypes / Clustering" = "MHBRS", "Scmap Celltypes / Clustering" = "JaJnS")
    
    restr <- switch(input$resfield, "gene_table",  "genes (consensus)"="gene_consensus", "pathways/annotations (within batches)"= "annot_table", "pathways/annotations (consensus)"= "annot_consensus" )
    if (class(data()) == "character") oldcols <- c()
    else oldcols <- colnames(data())
    
    data(readRDS(paste("/lustre/scratch117/cellgen/team218/lh20/SnakeFolderEv4/shinydata/NO_",dastr,"_",restr, ".rds", sep="")))
    
   #value(paste("/lustre/scratch117/cellgen/team218/lh20/results/table_NO_",input$dataset, ".rds", sep=""))
# 
      updateSelectInput(session,"filter", choices = colnames(data()), selected = colnames(data())[1])

      updateSelectInput(session,"obs",choices = setdiff(colnames(data()), c("Gene", "DE", "Log2FC", "LogitAuroc", "Comparison", "Celltype", "Archtype", "TPMmean", "DEseq_Log10pval", "Wilcox_Log10pval", "DEseq_adj_Log10pval", "Wilcox_adj_Log10pval", "DESeq_basemean", "FAD_coverage", "Ctrl_coverage", "FAD_Log2FC_toEmpty", "Ctrl_Log2FC_toEmpty", "MeanLog2FC", "MeanLog2FC", "MeanLogitAuroc", "Nbgenes","ID", "Domain","Tail", "pvalue", "Test" )))
      dataclean(oldcols)
      
      
      #tmp <- curflt()
      #cnt <- is.na(match(rownames(tmp), colnames(data()))) 
      #if (sum(cnt) != 0){
      #  if (sum(cnt) != nrow(tmp)) curflt(tmp[!cnt,])
      #  else curflt(data.frame(criterion= factor(c(),levels= c("is among","greater than", "less than", "equal to", "norm greater than", "norm less than")), value=character()))
      #}
      
      
      
      shinyjs::enable("resfield"); shinyjs::enable("dataset") ; shinyjs::enable("simplebutton"); shinyjs::enable("downloadData")
      dataclean(1)
    })
  
  
  observe({ # Load the current matrix for histogram plot
        progress <- Progress$new(session, min=0)
        shinyjs::disable("resfield"); shinyjs::disable("dataset");shinyjs::disable("simplebutton"); shinyjs::disable("downloadData")
        on.exit(progress$close())
        progress$set(message = 'Loading Table...',
        detail = 'This may take a few seconds')
            dastr <- switch(input$dataset, "MHBR",  "Fine Celltypes / Multinomial"  = "MH" , "Broad Celltypes / Multinomial" = "MHBR", "Scmap Celltypes / Multinomial" = "JaJn", "Fine Celltypes / Clustering"  = "MHS" , "Broad Celltypes / Clustering" = "MHBRS", "Scmap Celltypes / Clustering" = "JaJnS")
        mat(readRDS(paste("/lustre/scratch117/cellgen/team218/lh20/SnakeFolderEv4/shinydata/NO_",dastr,"_matrix.rds", sep="")))
        shinyjs::enable("resfield");shinyjs::enable("dataset");shinyjs::enable("simplebutton"); shinyjs::enable("downloadData")
    }) # Load the current matrix for histogram plot

  observe({ #### Update gene list for histogram
    if (sum(is.na(match(c("start", "length"), names(input$results_state)))) == 0){ # attribute might be missing at times
    maxo = input$results_state$start + input$results_state$length
    dalist <- which(filtrow())
    if (length(dalist) != 0){
      value("")
      if (length(input$results_state$order) != 0) {
        if ((input$showCols[as.numeric(input$results_state$order[[1]][1]) +1]) %in% colnames(data())) {
        #value(paste(data()[dalist,input$showCols[as.numeric(input$results_state$order[[1]][1]) +1]][1:10]))
        toord <- data()[dalist,input$showCols[as.numeric(input$results_state$order[[1]][1]) +1]]
        #value(ifelse(input$results_state$order[[1]][2]== "decr", T,F))
        if (class(toord) == "factor") dalist <- dalist[ order(as.character(toord, decreasing=ifelse(as.character(input$results_state$order[[1]][2])== "asc", F,T)))] 
        else dalist <- dalist[order(toord, decreasing=ifelse(as.character(input$results_state$order[[1]][2])== "asc", F,T))]
        }
      } # 
      #value(length(which(filtrow())[(input$results_state$start+1):maxo]))
      shownrows(dalist[(input$results_state$start+1): maxo])
      if (grepl("genes", input$resfield)) {
        plotgenes(unique(as.character(data()[dalist[(input$results_state$start+1): maxo], "Gene"])))
      }else{
        toplot <-strsplit(as.character(data()[ ifelse(length(input$results_rows_selected) == 0, dalist[input$results_state$start+1], which(filtrow())[input$results_rows_selected]), "Intersection"]), "," )[[1]]
        plotgenes(toplot)
}}}})  # Update gene list for histogram
  
  
  observe({ # Update filter value sets
    # get all character or factor columns
    if ((class(data()) != "data.frame")||(!input$filter %in% colnames(data()))){
    }else if (class(data()[[input$filter]]) == "factor"){
      shinyjs::disable("filtervalue")
      updateSelectInput(session, "filterchoice", choices =unique(data()[[input$filter]]))
    }else{
      shinyjs::enable("filtervalue")
      updateSelectInput(session, "filterchoice", choices =c("greater than", "less than", "equal to", "norm greater than", "norm less than"), selected = "greater than")
    }
  })
  
  observeEvent(input$simplebutton, {  # Simple Query is being executed
    print("query simple");
    if (input$dataset != "Broad Celltypes / Multinomial") updateSelectInput(session, "dataset", label = "Celltype and Sample calling:",
          choices = c("Fine Celltypes / Multinomial" ,  "Broad Celltypes / Multinomial", "Scmap Celltypes / Multinomial" ,  "Fine Celltypes / Clustering" ,  "Broad Celltypes / Clustering", "Scmap Celltypes / Clustering"),selected = "Broad Celltypes / Multinomial")

    tmp <- ifelse(input$simpledetype %in% c("Upregulated pathways in Disease", "Downregulated pathways in Disease"), "pathways/annotations (consensus)", "genes (consensus)")
    if (input$resfield != tmp){
      updateSelectInput(session, "resfield", label = "Choose a result type:",
      choices = c("genes (within batches)" , "genes (consensus)" , "pathways/annotations (within batches)", "pathways/annotations (consensus)"),selected = tmp)
    }
    
    tlvl <-c("is among","greater than", "less than", "equal to", "norm greater than", "norm less than")
    compset <- switch(input$simplecondition, "APP V717I in Neurons" = "V717IHtNeuro", "APP V717I in Microglia"= "V717IHtMicro", "PSEN1 M146I in Neurons" = "M146IHtNeuro", "PSEN1 M146I in Microglia" = "M146IHtMicro", "PSEN1 Intron4 mutation in Neurons"= "Intr4HtNeuro", "PSEN1 Intron4 mutation in Microglia"= "Intr4HtMicro", "LPS"="LPS", "TREM2 knock-out Microglia" = "TREM2KO")
    ctset <- switch(input$simplecelltype, "Microglia" = "Microglia", "Neurons" = "Neuron_cortical", "Neuron Precursor" = "IPC", "Neuron and Microglia"= "Neuron_cortical ; Microglia", "All" ="")


    
    tmp <- data.frame(row.names =c("ConsensusGroup"), criterion=factor(c("is among"), levels= tlvl), value = as.character(compset))
    if (nchar(ctset) != 0) tmp <- rbind(tmp,data.frame(row.names = c("Celltype") , criterion=factor(c("is among"), levels= tlvl), value=as.character(ctset)))
    if (input$simpledetype == "Higher Expression in Disease"){
      simplesort("pfc")
    tmp <- rbind(tmp,data.frame(row.names = c("Log2FC") , criterion=factor(c("greater than"), levels= tlvl), value=as.character(0)))
    }else if (input$simpledetype == "Lower Expression in Disease"){
      simplesort("nfc")
      tmp <- rbind(tmp,data.frame(row.names = c("Log2FC") , criterion=factor(c("less than"), levels= tlvl), value=as.character(0)))
    }else if (input$simpledetype == "Significant for DEseq2"){
      simplesort("signifD")
      tmp <- rbind(tmp,data.frame(row.names = c("DEseq_adj_Log10pval") , criterion=factor(c("less than"), levels= tlvl), value=as.character(-1.3)))
    }else if (input$simpledetype == "Significant for Wilcox test"){
      simplesort("signifW")
      tmp <- rbind(tmp,data.frame(row.names = c("Wilcox_adj_Log10pval") , criterion=factor(c("less than"), levels= tlvl), value=as.character(-1.3)))
    }else if (input$simpledetype == "Upregulated pathways in Disease"){
      simplesort("signifP")
      tmp <- rbind(tmp,data.frame(row.names = c("Domain") , criterion=factor(c("is among"), levels= tlvl), value="keg ; rea"))
      tmp <- rbind(tmp,data.frame(row.names = c("pvalue") , criterion=factor(c("less than"), levels= tlvl), value="0.001"))
    }else{
      simplesort("signifN")
      tmp <- rbind(tmp,data.frame(row.names = c("Domain") , criterion=factor(c("is among"), levels= tlvl), value="keg ; rea"))
      tmp <- rbind(tmp,data.frame(row.names = c("pvalue") , criterion=factor(c("less than"), levels= tlvl), value="0.001"))
    }
    tmp$value <- as.character(tmp$value)
    curflt(tmp)
    
    tmp <- c("Microglia", "Neurons", "Microglia and Neurons","Match Filters","All")
    tmp2 <- c("Match ConsensusGroup","Match Comparison","point-mutation conditions", "other disease conditions", "other neutral conditions","Include All")
#    tmp3 <- c("Heatmap" , "Volcano Plot" , "Tsne Overlay")

#    if (input$simpleextra == "Heatmap with other all celltypes"){
#      updateSelectInput(session, "samexcl",choices =  tmp2,selected = "Match ConsensusGroup")
#      updateSelectInput(session, "ctpexcl",choices =  tmp,selected = "All")
#      updateSelectInput(session, "contextfield",choices =  tmp3,selected ="Heatmap")
#    }else if (input$simpleextra == "Heatmap with other conditions"){
#      updateSelectInput(session, "samexcl",choices =  tmp2,selected = "Include All")
#      updateSelectInput(session, "ctpexcl",choices =  tmp,selected = "Match Filters")
#      updateSelectInput(session, "contextfield",choices =  tmp3,selected ="Heatmap")
#    }else if (input$simpleextra == "Tsne Overlay of Wilcox test"){
#      updateSelectInput(session, "samexcl",choices =  tmp2,selected = tmp2[3])
#      updateSelectInput(session, "ctpexcl",choices =  tmp,selected = "Match Filters")
#      updateSelectInput(session, "contextfield",choices =  tmp3,selected ="Tsne Overlay")
#    }else { # "Volcano Plot of DEseq"
#      updateSelectInput(session, "samexcl",choices =  tmp2,selected = "Include All")
#      updateSelectInput(session, "ctpexcl",choices =  tmp,selected = "All")
#      updateSelectInput(session, "contextfield",choices =  tmp3,selected ="Volcano Plot")
#    }

    
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("DEtable_", gsub(":","-",gsub(" ", "_", Sys.time())) , ".csv", sep = "")
    },
    content = function(file) {
      write.csv(data()[filtrow(),input$showCols], file, row.names = FALSE)
    }
  )
  
  observeEvent(input$fltaddbutton, {  # Filter being Added or Removed
    simplesort("")
    if ((!is.null(input$currentfilters_rows_selected))&&(length(input$currentfilters_rows_selected) != 0)){
      tmp <- curflt()
      daflt <- (setdiff(1:nrow(tmp), input$currentfilters_rows_selected))
      if (length(daflt) == 0) curflt(data.frame(criterion= factor(c(),levels= c("is among","greater than", "less than", "equal to", "norm greater than", "norm less than")), value=character()))
      else curflt(tmp[daflt,])
      #if (nrow(curflt()) == 0) runjs("document.getElementById('currentfilters').style.display='none'")
    }else{
    #debugstate <- "button pressed"
      if (input$filter %in% rownames(curflt())){
          tmp <- curflt()
          if ( input$filterchoice %in% c("greater than", "less than", "equal to", "norm greater than", "norm less than")){
            tmp[input$filter, "criterion"] = as.character(input$filterchoice)
            tmp[input$filter, "value"] = as.character(input$filtervalue)
            curflt(tmp)
          }else{
            curcur <- strsplit(tmp[input$filter, "value"], "[[:space:]];[[:space:]]")[[1]]
            if (is.na(match(as.character(input$filterchoice), curcur))) tmp[input$filter, "value"] <- paste(sort(c(curcur, as.character(input$filterchoice))), collapse = " ; ")
            else tmp[input$filter, "value"] <- paste(sort(setdiff(curcur, as.character(input$filterchoice))), collapse = " ; ")
            
            if (nchar(tmp[input$filter, "value"]) == 0){
              if (nrow(tmp) == 1) curflt(data.frame(criterion= factor(c(),levels= c("is among","greater than", "less than", "equal to", "norm greater than", "norm less than")), value=character()))
              else curflt(tmp[setdiff(1:nrow(tmp), match(input$filter, rownames(tmp))),])
            }else curflt(tmp)
          }
      }else{
        if ( input$filterchoice %in% c("greater than", "less than", "equal to", "norm greater than", "norm less than")){
          that <- rbind(curflt(),data.frame(row.names = c(input$filter), criterion=factor(c(as.character(input$filterchoice)),levels= c("is among","greater than", "less than", "equal to", "norm greater than", "norm less than")), value=as.character(input$filtervalue)))
          that$value <- as.character(that$value)
          #runjs("document.getElementById('currentfilters').style.display='block'")
          curflt(that)
        }else{
          that <- rbind(curflt(),data.frame(row.names = c(input$filter), criterion=factor(c("is among"),levels= c("is among","greater than", "less than", "equal to", "norm greater than", "norm less than")), value= as.character(input$filterchoice) )) #
          #runjs("document.getElementById('currentfilters').style.display='block'")
          that$value <- as.character(that$value)
          curflt(that)
        }
      }
    }
  })
  
  
              
  output$currentfilters <- renderDataTable({ #update filter table
      if (is.null(curflt())) datatable(data.frame())
      else if (nrow(curflt()) == 0) datatable(data.frame())
      else {
      datatable(curflt(), options = list(pageLength = nrow(curflt()), dom = 't'), caption = 'Filters currently applied on rows of the table:', colnames = c("", ""))
        
        }
    })
  
  output$results <- renderDataTable({ #update result table
                if (dataclean() == 0){
                  DT::datatable(data.frame(row.names=c("Nothing")))
                }else{
                # if(all(startsWith(object@index$genes(), "chr") == T)) "Peaks" else
           
                #input$selectInput$choices <- unique(data$gene$Celltype)
                fltrow <- rep(T, nrow(data()))
                if (nrow(curflt())>0) {
                  for(i in 1:nrow(curflt())){
                    if (rownames(curflt())[i] %in% colnames(data())){
                      if (curflt()$criterion[i] == "is among"){
                        fltrow <- fltrow & data()[[rownames(curflt())[i]]] %in% strsplit(curflt()$value[i] , "[[:space:]];[[:space:]]")[[1]]
                      }else{
                        if (curflt()$criterion[i] == "greater than") compres <- (as.numeric(data()[[rownames(curflt())[i]]]) > as.numeric(as.character(curflt()$value[i])))
                        else if (curflt()$criterion[i] == "less than") compres <- (as.numeric(data()[[rownames(curflt())[i]]]) < as.numeric(as.character(curflt()$value[i])))
                        else if (curflt()$criterion[i] == "norm greater than") compres <- (abs(as.numeric(data()[[rownames(curflt())[i]]])) > as.numeric(as.character(curflt()$value[i])))
                        else if (curflt()$criterion[i] == "norm less than") compres <- (abs(as.numeric(data()[[rownames(curflt())[i]]])) < as.numeric(as.character(curflt()$value[i])))
                        else compres <- (as.numeric(data()[[rownames(curflt())[i]]]) == as.numeric(as.character(curflt()$value[i])))
                        compres[is.na(compres)] <- F
                        fltrow <- fltrow & compres
                      }
                    }
                  }
                }
                #fltrow
                  #value(colnames(data))
#                  value(match(input$showCols, colnames(data)))
                
                filtrow(fltrow)
                if ((sum(fltrow) == 0)||(sum(is.na(match(input$showCols,colnames(data())))) != 0)){ DT::datatable(data.frame(row.names=c("Nothing")))
                  #value(setdiff(input$showCols, colnames(data())))
                }else{
                  lengthlist = c(5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100, 120)
                  if (sum(fltrow) < 120) lengthlist = c(lengthlist[lengthlist < sum(fltrow)], sum(fltrow))
                  
                  defsort <- switch(simplesort(), list(NA, NA),"pfc" = list(match("Log2FC", input$showCols), "desc"), "nfc" = list(match("Log2FC", input$showCols), "asc"), "signifD"= list(match("DEseq_adj_Log10pval", input$showCols), "asc"), "signifW"= list(match("Wilcox_adj_Log10pval", input$showCols), "asc"), "signifP"= list(match("MeanLog2FC", input$showCols), "desc"), "signifN"= list(match("MeanLog2FC", input$showCols), "asc"))
                  #defsort <- list(NA, NA)
                  #toroundlist <- c()
                  #for(elem in input$showCols) if    (class(data()[[input$filter]]) == "factor")
                  optstr <- list(scrollX = TRUE, dom = 'lpt', stateSave=T, lengthMenu = lengthlist)
                  if (!is.na(defsort[1])) optstr <- c(optstr, list(order = list(defsort)))
                  #value(intersect(input$showCols, c("Log2FC", "MeanLog2FC", "LogitAuroc","TPMmean","DEseq_adj_Log10pval")))
                  DT::datatable(data()[fltrow,input$showCols], selection = 'single',
                  extensions = 'Scroller', colnames = input$showCols, options = optstr, rownames = F)
                  # %>% DT::formatRound(columns=intersect(input$showCols, c("Log2FC", "MeanLog2FC", "LogitAuroc","TPMmean","DEseq_adj_Log10pval")), digits=3)

                  }
                }
})
  
  
  observe({ #draw Heatmap / Vocano / overlay
    if (length(input$results_rows_selected) == 0) {
      comps = ""
      gsize <- c(1,1)
    }else if (grepl("consensus", input$resfield)){
      comps <- mat()$cons[[as.character(data()[which(filtrow())[input$results_rows_selected], "ConsensusGroup"])]]
      gsize <- c(trunc((2+length(comps)) / 3) , ifelse(length(comps) <3, length(comps) ,3))
      if (length(comps) == 2) gsize <- c(2,2)
    }else{
      comps <- as.character(data()[which(filtrow())[input$results_rows_selected], "Comparison"])
      gsize <- c(1,1)
    } 
  
  output$map <- renderPlot({
      if (input$tabContext == "Heatmap"){

        if (length(plotgenes()) > 1){
              curcolnames <- colnames(mat()$deseq$logpval)
              
              if (input$comtype == "All") colfilt <- rep(T, length(curcolnames))
              else if (input$comtype == "Pooled Comparisons") colfilt <- mat()$ispool[mat()$coltotest]
              else colfilt <- !mat()$ispool[mat()$coltotest]
              
              if (input$ctpexcl == "Microglia") tmp <- grepl("[Mm]icroglia", levels(mat()$celltypes))
              else if (input$ctpexcl == "Neurons") tmp <- grepl("[Nn]euron", levels(mat()$celltypes))
              else if (input$ctpexcl == "Microglia and Neurons") tmp <- grepl("[Nn]euron", levels(mat()$celltypes)) | grepl("[Mm]icroglia", levels(mat()$celltypes))
              else if (input$ctpexcl == "Match Filters"){
                tmp <- match("Celltype", rownames(curflt()))
                if (is.na(tmp)) tmp <- rep(T, length(levels(mat()$celltypes)))
                else tmp <- !is.na(match(levels(mat()$celltypes) , strsplit(curflt()$value[tmp] , "[[:space:]];[[:space:]]")[[1]]))
              }else tmp <- rep(T, length(levels(mat()$celltypes)))
              #value(tmp)
             colfilt <- colfilt & tmp[mat()$coltoct]
              tmp <- rep(T, length(mat()$comparisons))
              if (input$samexcl == "Match ConsensusGroup"){
                tmp <- match("ConsensusGroup", rownames(curflt()))
                if (is.na(tmp)) tmp <- rep(T, length(mat()$comparisons))
                else tmp <- !is.na(match(mat()$archt, strsplit(curflt()$value[tmp] , "[[:space:]];[[:space:]]")[[1]]))
              }else if (input$samexcl == "Match Comparison"){
                tmp <- match("Comparison", rownames(curflt()))
                if (is.na(tmp)) tmp <- rep(T, length(mat()$comparisons))
                else tmp <- !is.na(match(mat()$comparisons , strsplit(curflt()$value[tmp] , "[[:space:]];[[:space:]]")[[1]]))
              }else if (input$samexcl == "point-mutation conditions"){  tmp <- mat()$archt
                tmp <- !(grepl("LPS", tmp) | grepl("H9_vs_KOLF2", tmp) | grepl("TREM2KO",tmp) | grepl("TrueNegative",tmp))
              }else if (input$samexcl == "other disease conditions") {  tmp <- mat()$archt
                tmp <- grepl("LPS", tmp) | grepl("TREM2KO",tmp)
              }else if (input$samexcl =="Include All") tmp <- rep(T, length(mat()$comparisons))
              else { tmp <- mat()$archt
                tmp <-  grepl("Replicate_for_WT",tmp) | grepl("Replicate_for_Mutation",tmp) | grepl("H9Micro", tmp)
              }
              colfilt <- colfilt & tmp[mat()$coltotest]
              
              names(colfilt) <- NULL
              
              if (sum(colfilt) > input$nbhistcols) {
                colselect <- order(colSums(mat()$deseq$logpval[plotgenes(),colfilt,drop=F] < -1.3), decreasing=T)
                colselect <- colselect[1:input$nbhistcols]
                colselect <- sort(match(curcolnames[colfilt], curcolnames)[colselect])
              }else colselect <- which(colfilt)
  
              
              c1mat <- matrix("#AAAAAA", nrow= length(plotgenes()), ncol = length(colselect))
              c2mat <- c1mat
              
              override <- list(top = rep("", length(colselect)), bot = rep("", length(colselect)), axenames=c("Comparisons", "Cell Types"))
              whiteCMP <- rgb(col2rgb(mat()$color_CMP)/1020 + 0.75)
              whiteCT <- rgb(col2rgb(mat()$color_CT)/1020 + 0.75)
              colcolors <- rep("#AAAAAA", length(colselect)) 
              for(j in 1:length(colselect)) {
                override$top[j] <- levels(mat()$celltypes)[mat()$coltoct[colselect[j]]]
                override$bot[j] <- mat()$comp_titles[mat()$coltotest[colselect[j]]]
                c1mat[,j] <- rep(whiteCMP[mat()$coltotest[colselect[j]]], nrow(c1mat))
                c2mat[,j] <- rep(whiteCT[mat()$coltoct[colselect[j]]], nrow(c1mat)) 
                colcolors[j] <- mat()$color_CT[mat()$coltoct[colselect[j]]]
              }
              return(plotDataGrid(list(data = mat()$deseq$log2FC[rev(plotgenes()),colselect,drop=F], w=mat()$deseq$logpval[rev(plotgenes()),colselect,drop=F], c1 = c1mat, c2 = c2mat), colcolors = colcolors, do.cluster = c(input$clusterheat %in% c("Cluster Genes","Cluster Both"),input$clusterheat %in% c("Cluster Columns","Cluster Both")), transform=list(w="log10pval"), override.colnames = override, plot.attribs =list(xlabel = "Cell-type x Comparison", ylabel= "Genes")))
        }else if ((length(plotgenes()) == 0)||(! plotgenes() %in% rownames(mat()$deseq$logpval))) ggplot()
        else {#comtype
            rnam = levels(mat()$celltypes)
            cnam = names(mat()$ispool)[colfilt]
  
            dmat <- matrix(0, nrow= length(rnam), ncol = length(cnam), dimnames = list(rnam,cnam))
            wmat <- dmat
  
            c1mat <- matrix("#AAAAAA", nrow= length(rnam), ncol = length(cnam), dimnames = list(rnam,cnam))
            c2mat <- matrix("#888888", nrow= length(rnam), ncol = length(cnam), dimnames = list(rnam,cnam))
  
  
            for(j in 1:length(cnam)) {
              for(i in 1:length(rnam)) {
              k <- match(paste(rnam[i], cnam[j],sep="_"), colnames(mat()$deseq$log2FC))
              if (!is.na(k)) dmat[i,j] <- mat()$deseq$log2FC[plotgenes(),k]
              k <- match(paste(rnam[i], cnam[j],sep="_"), colnames(mat()$deseq$logpval))
              if (!is.na(k)) wmat[i,j] <- mat()$deseq$logpval[plotgenes(),k]
              }
              c1mat[,j] <- rep(mat()$color_CMP[j], nrow(c1mat))
              c2mat[,j] <- rep(mat()$color_CMP[j], nrow(c1mat)) 
            }
  
  
            return(plotDataGrid(list(data = dmat , w=wmat, c1 = c1mat, c2 = c2mat), transform=list(w="log10pval")))
        }
  }else if (length(input$results_rows_selected) == 0){
    return(ggplot() + ggtitle("Select a row in the table above to view a contextual display here"))
  } else {
    currow <- which(filtrow())[input$results_rows_selected]
    if (input$tabContext == "Volcano Plot"){
      gglist <- list();
      
      dact <- as.character(data()[currow, "Celltype"])
      colsel <- paste(dact, comps, sep = "_")
      if (length(colsel) == 0) stop("nolength!")
      
      danames <- rownames(mat()$deseq$logpval)
      danames[! (danames %in% plotgenes()) ] <- ""

      daalpha <- rep(1, nrow(mat()$deseq$log2FC))
      daalpha[!(danames %in% plotgenes()) ] <- daalpha[!(danames %in% plotgenes()) ] * 0.125
      for(i in 1:length(colsel)){
        dacolor <- rep("#000000", nrow(mat()$dese$log2FC))
        dacolor[(mat()$deseq$logpval[, colsel[i]] < -1.30103) & (mat()$deseq$log2FC[, colsel[i]] < 0) ] <- "#0088FF"
        dacolor[(mat()$deseq$logpval[, colsel[i]] < -1.30103) & (mat()$deseq$log2FC[, colsel[i]] > 0) ] <- "#FF6600"
        dacolor[(mat()$deseq$logpval[, colsel[i]] < -1.30103) & (mat()$deseq$log2FC[, colsel[i]] < 0) & (danames %in% plotgenes()) ] <- "#0000FF"
        dacolor[(mat()$deseq$logpval[, colsel[i]] < -1.30103) & (mat()$deseq$log2FC[, colsel[i]] > 0) & (danames %in% plotgenes()) ] <- "#DD0000"
        gglist <- c(gglist, list(plotLabels(mat()$deseq$log2FC[, colsel[i]], -mat()$deseq$logpval[, colsel[i]], danames, color = dacolor, alpha = daalpha, filter = (mat()$deseq$logpval[, colsel[i]] < -1.0), point.size = 3, plot.attribs = list(xlabel = "Log2FC", ylabel= "-log10 Pvalue", title = mat()$comp_titles[match(comps[i], mat()$comparisons)]))))
      }
      #value(paste(colsel , length(gglist)))
      return(grid_arrange_shared_legend(gglist, do.share.legend=F, nrow=gsize[1], ncol=gsize[2],position = "right", main.title = paste("Deseq DE genes in ", dact, sep="")) )
    }else{
      value(colnames(data()))
      if ("Gene" %in% colnames(data())){
        dagene <- data()[currow, "Gene"]
        dacol <- match(dagene, colnames(overlay()$dropout))
        subr <- overlay()$dropout@p[c(dacol,dacol+1)]
      return(makeOverlay(overlay(), overlay()$dematrices[[match(dagene, names(overlay()$dematrices))]], overlay()$dropout@i[(subr[1]+1):(subr[2])], dagene, comps, gridsize = gsize, titles = mat()$comp_titles[match(comps, mat()$comparisons)]))
      }else{
        return(ggplot() + ggtitle("Tsne Overlay not available for pathways/gene sets"))
        #return(makeOverlay(overlay(), overlay()$dematrices[[match( plotgenes() , names(overlay()$dematrices))]], c(), data()[currow, "Term"], comps, gridsize = gsize, titles = mat()$comp_titles[match(comps, mat()$comparisons)]))
      }
    }
  }

  
  
  }, height = ifelse(input$tabContext == "Heatmap" , 300 + ifelse(length(plotgenes()) == 1, length(unique(data()[["Celltype"]])) , length(plotgenes()))* 24,gsize[2] * 500)
 # ,width = ifelse(gsize[1] == 3, "100%", ifelse(gsize[1] == 2 , "75%","50%"))
  )})

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

  
  #observeEvent(output$results, {
  #  if 
    
  #input$results_rows_selected
  #})

output$help <- renderText({
      #input$results_rows_selected
      ifelse(length(input$results_rows_selected) == 0, "" ,as.character(data()[which(filtrow())[input$results_rows_selected], input$obs]))
        #ifelse(last.query.state() == "genelist", 'not right', 'tight')
    })
output$help2 <- renderText({
    value()
    #if (debug.state == "button pressed") "button pressed"
    #else "something"
  })

observe({ #draw tsne overlay
  output$atlas <- renderPlot({
    if (is.null(names(overlay()))) return(ggplot())
    else{
      return(makeTsne(overlay()$coords, mat()$celltype, mat()$color_CT))
    }
})})
    
    
    
#output$myWebGL <- renderWebGL({
#    points3d(1:10, 1:10, 1:10)
#    axes3d()
#  })

}

# Create Shiny app ----
# shinyApp(ui = ui, server = server)




