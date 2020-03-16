source("render.R")
options(DT.fillContainer = FALSE)
options(DT.autoHideNavigation = FALSE)
#' Server handler for detablebrowser
#'
#' @importFrom DT renderDataTable datatable
server <- function(input, output, session) {
  last.query.state <- reactiveVal("genelist")
  debug.state <- reactiveVal(0); dataclean <- reactiveVal(0)
  value <- reactiveVal("");      data <- reactiveVal("")
  mat <- reactiveVal("");        simplesort <- reactiveVal("")
  plotgenes <- reactiveVal(c(""))
  
  curflt <- reactiveVal(data.frame(criterion= c(), value=character())) 
  filtrow <- reactiveVal(c(T))
  
  
  observe({
      if (dataclean() == 0){
       
      }else{
      updateSelectInput(session,"filter", choices = colnames(data()), selected = colnames(data())[1])

      updateSelectInput(session,"obs",choices = setdiff(colnames(data()), c("Gene", "DE", "Log2FC", "LogitAuroc", "Comparison", "Celltype", "Archtype", "TPMmean", "DEseq_Log10pval", "Wilcox_Log10pval", "DEseq_adj_Log10pval", "Wilcox_adj_Log10pval", "DESeq_basemean", "FAD_coverage", "Ctrl_coverage", "FAD_Log2FC_toEmpty", "Ctrl_Log2FC_toEmpty", "MeanLog2FC", "MeanLog2FC", "MeanLogitAuroc", "Nbgenes","ID", "Domain","Tail", "pvalue", "Test" )))
      helpstr <- "GOSlim"
      ext <- c("FullName", "GO", "GOslim", "Description", "Intersection")
      danames <- setdiff(colnames(data()), helpstr)
      # set tablecol filter, with default values
      if (input$resfield == "gene") tofilt <- c("DEseq_Log10pval", "Wilcox_Log10pval")
      else tofilt <- c("DEseq_adj_Log10pval", "Wilcox_adj_Log10pval")
      updateCheckboxGroupInput(session,"showCols", choices = danames, selected = setdiff(danames, c(c("FullName", "GO", "GOslim", "Description", "Archtype", "FAD_Log2FC_toEmpty", "Ctrl_Log2FC_toEmpty","Alias", "LogitAuroc", "sample_ctrlIds", "sample_testIds"), tofilt)))
      }
      
    })
  
  
  observe({ # Load the current table data
    dataclean(0)
    shinyjs::disable("resfield"); shinyjs::disable("dataset"); shinyjs::disable("simplebutton")
    progress <- Progress$new(session, min=0)
    on.exit(progress$close())
    progress$set(message = 'Loading Table...',
    detail = 'This may take a few seconds')
    dastr <- switch(input$dataset, "MHBR", "Fine Celltypes / Multinomial"  = "MH" , "Broad Celltypes / Multinomial" = "MHBR", "Scmap Celltypes / Multinomial" = "JaJn", "Fine Celltypes / Clustering"  = "MHS" , "Broad Celltypes / Clustering" = "MHBRS", "Scmap Celltypes / Clustering" = "JaJnS")
    
    
    data(readRDS(paste("/lustre/scratch117/cellgen/team218/lh20/results/table_",dastr,"_",input$resfield, ".rds", sep="")))
   #value(paste("/lustre/scratch117/cellgen/team218/lh20/results/table_NO_",input$dataset, ".rds", sep=""))

      updateSelectInput(session,"filter", choices = colnames(data()), selected = colnames(data())[1])

      updateSelectInput(session,"obs",choices = setdiff(colnames(data()), c("Gene", "DE", "Log2FC", "LogitAuroc", "Comparison", "Celltype", "Archtype", "TPMmean", "DEseq_Log10pval", "Wilcox_Log10pval", "DEseq_adj_Log10pval", "Wilcox_adj_Log10pval", "DESeq_basemean", "FAD_coverage", "Ctrl_coverage", "FAD_Log2FC_toEmpty", "Ctrl_Log2FC_toEmpty", "MeanLog2FC", "MeanLog2FC", "MeanLogitAuroc", "Nbgenes","ID", "Domain","Tail", "pvalue", "Test" )))
      helpstr <- "GOSlim"
      ext <- c("FullName", "GO", "GOslim", "Description", "Intersection")
      danames <- setdiff(colnames(data()), helpstr)
      updateCheckboxGroupInput(session,"showCols", choices = danames, selected = setdiff(danames, c("FullName", "GO", "GOslim", "Description")))
      shinyjs::enable("resfield"); shinyjs::enable("dataset") ; shinyjs::enable("simplebutton")
      dataclean(1)
    })
  
  
  observe({ # Load the current matrix for histogram plot
        progress <- Progress$new(session, min=0)
        shinyjs::disable("resfield"); shinyjs::disable("dataset");shinyjs::disable("simplebutton")
        on.exit(progress$close())
        progress$set(message = 'Loading Table...',
        detail = 'This may take a few seconds')
            dastr <- switch(input$dataset, "MHBR",  "Fine Celltypes / Multinomial"  = "MH" , "Broad Celltypes / Multinomial" = "MHBR", "Scmap Celltypes / Multinomial" = "JaJn", "Fine Celltypes / Clustering"  = "MHS" , "Broad Celltypes / Clustering" = "MHBRS", "Scmap Celltypes / Clustering" = "JaJnS")
        mat(readRDS(paste("/lustre/scratch117/cellgen/team218/lh20/results/matrix_",dastr,".rds", sep="")))
        shinyjs::enable("resfield");shinyjs::enable("dataset");shinyjs::enable("simplebutton")
    }) # Load the current matrix for histogram plot

  observe({ #### Update gene list for histogram
    if (sum(is.na(match(c("start", "length"), names(input$results_state)))) == 0){ # attribute might be missing at times
    maxo = input$results_state$start + input$results_state$length
    dalist <- which(filtrow())
    if (length(dalist) == 0) value("no row selected")
    else{
      if (length(input$results_state$order) != 0) {
        #value(paste(data()[dalist,input$showCols[as.numeric(input$results_state$order[[1]][1]) +1]][1:10]))
        toord <- data()[dalist,input$showCols[as.numeric(input$results_state$order[[1]][1]) +1]]
        #value(ifelse(input$results_state$order[[1]][2]== "decr", T,F))
        #value(as.character(input$results_state$order[[1]][2]))
        if (class(toord) == "factor") dalist <- dalist[ order(as.character(toord, decreasing=ifelse(as.character(input$results_state$order[[1]][2])== "asc", F,T)))] 
        else dalist <- dalist[order(toord, decreasing=ifelse(as.character(input$results_state$order[[1]][2])== "asc", F,T))] 
      } # 
      #value(length(which(filtrow())[(input$results_state$start+1):maxo]))
      
      if (is.na(match(input$resfield ,c("go", "consensus_go")))) plotgenes(unique(as.character(data()[dalist[(input$results_state$start+1): maxo], "Gene"])))
      else{
        
        toplot <-strsplit(as.character(data()[dalist[input$results_state$start+ 1 + ifelse(length(input$results_rows_selected) == 0, 0, input$results_rows_selected)], "Intersection"]), "," )[[1]]
        value(length(toplot))
        if (length(toplot) > 30) toplot <- toplot[1:30]
        plotgenes(toplot)
}}}})  # Update gene list for histogram
  
  
  observe({ # Update filter value sets
    # get all character or factor columns
    if (class(data()[[input$filter]]) == "factor"){
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
    }else{
      simplesort("signifW")
      tmp <- rbind(tmp,data.frame(row.names = c("Wilcox_adj_Log10pval") , criterion=factor(c("less than"), levels= tlvl), value=as.character(-1.3)))
    }
    tmp$value <- as.character(tmp$value)
    curflt(tmp)
    
    tmp <- c("Microglia", "Neurons", "Microglia and Neurons","Match Filters","All")
    tmp2 <- c("Match ConsensusGroup","Match Comparison","point-mutation conditions", "other disease conditions", "other neutral conditions","Include All")
    
    if (input$simpleheat == "All celltypes"){
      updateSelectInput(session, "samexcl",choices =  tmp2,selected = "Match ConsensusGroup")
      updateSelectInput(session, "ctpexcl",choices =  tmp,selected = "All")
    }else if (input$simpleheat == "All point mutations"){
      updateSelectInput(session, "samexcl",choices =  tmp2,selected = tmp2[3])
      updateSelectInput(session, "ctpexcl",choices =  tmp,selected = "Match Filters")
    }else if (input$simpleheat == "All conditions"){
      updateSelectInput(session, "samexcl",choices =  tmp2,selected = "Include All")
      updateSelectInput(session, "ctpexcl",choices =  tmp,selected = "Match Filters")
    }else{
      updateSelectInput(session, "samexcl",choices =  tmp2,selected = "Include All")
      updateSelectInput(session, "ctpexcl",choices =  tmp,selected = "All")
    }

    
  })
  
  observeEvent(input$fltaddbutton, {  # Filter being Added or Removed
    simplesort("")
    if ((!is.null(input$currentfilters_rows_selected))&&(length(input$currentfilters_rows_selected) != 0)){
      tmp <- curflt()
      daflt <- (setdiff(1:nrow(tmp), input$currentfilters_rows_selected))
      value(daflt)
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
  
  observe({ #draw Heatmap
  output$map <- renderPlot({
      if (length(plotgenes()) > 1){
            curcolnames <- colnames(mat()$deseq$logpval)
            
            if (input$comtype == "All") colfilt <- rep(T, length(curcolnames))
            else if (input$comtype == "Pooled Comparisons") colfilt <- mat()$ispool[mat()$coltotest]
            else colfilt <- !mat()$ispool[mat()$coltotest]
            
            if (input$ctpexcl == "Microglia") tmp <- grepl("[Mm]icroglia", mat()$celltype)
            else if (input$ctpexcl == "Neurons") tmp <- grepl("[Nn]euron", mat()$celltype)
            else if (input$ctpexcl == "Microglia and Neurons") tmp <- grepl("[Nn]euron", mat()$celltype) | grepl("[Mm]icroglia", mat()$celltype)
            else if (input$ctpexcl == "Match Filters"){
              tmp <- match("Celltype", rownames(curflt()))
              if (is.na(tmp)) tmp <- rep(T, length(mat()$celltype))
              else tmp <- !is.na(match(mat()$celltype , strsplit(curflt()$value[tmp] , "[[:space:]];[[:space:]]")[[1]]))
            }else tmp <- rep(T, length(mat()$celltype))

           colfilt <- colfilt & tmp[mat()$coltoct]
            tmp <- rep(T, length(mat()$comparisons))
            if (input$samexcl == "Match ConsensusGroup"){
              tmp <- match("ConsensusGroup", rownames(curflt()))
              if (is.na(tmp)) tmp <- rep(T, length(mat()$comparisons))
              else tmp <- !is.na(match(mat()$archt, mat()$ConsensusGroup))
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
              tmp <- grepl("H9_vs_KOLF2", tmp) | grepl("TrueNegative",tmp)
            }
            colfilt <- colfilt & tmp[mat()$coltotest]
            
            names(colfilt) <- NULL

            colselect <- order(colSums(mat()$deseq$logpval[plotgenes(),colfilt,drop=F] < -1.3), decreasing=T)
            if (length(colselect) > input$nbhistcols) colselect <- colselect[1:input$nbhistcols]
            colselect <- match(curcolnames[colfilt], curcolnames)[colselect]
            c1mat <- matrix("#AAAAAA", nrow= length(plotgenes()), ncol = length(colselect))
            c2mat <- c1mat
            for(j in 1:length(colselect)) {
              curname <- curcolnames
              c1mat[,j] <- rep(mat()$colA[mat()$coltotest[colselect[j]]], nrow(c1mat))
              c2mat[,j] <- rep(mat()$colB[mat()$coltotest[colselect[j]]], nrow(c1mat)) 
            }
            plotDataGrid(list(data = mat()$deseq$log2FC[plotgenes(),colselect,drop=F], w=mat()$deseq$logpval[plotgenes(),colselect,drop=F], c1 = c1mat, c2 = c2mat), transform=list(w="log10pval"),plot.attribs =list(xlabel = "Cell-type x Comparison", ylabel= "Genes"))
          
      }else if ((length(plotgenes()) == 0)||(! plotgenes() %in% rownames(mat()$deseq$logpval))) ggplot()
      else {#comtype
          rnam = mat()$celltypes
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
            c1mat[,j] <- rep(mat()$colA[j], nrow(c1mat))
            c2mat[,j] <- rep(mat()$colB[j], nrow(c1mat)) 
          }

          plotDataGrid(list(data = dmat , w=wmat, c1 = c1mat, c2 = c2mat), transform=list(w="log10pval"))
      }
  }, height = 300 + ifelse(length(plotgenes()) == 1, length(unique(data()[["Celltype"]])) , length(plotgenes()))* 24)
  })
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
                #fltrow
                  #value(colnames(data))
#                  value(match(input$showCols, colnames(data)))
                
                filtrow(fltrow)
                if (sum(fltrow) == 0) DT::datatable(data.frame(row.names=c("Nothing")))
                else{
                  lengthlist = c(5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100, 120)
                  if (sum(fltrow) < 120) lengthlist = c(lengthlist[lengthlist < sum(fltrow)], sum(fltrow))
                  
                  defsort <- switch(simplesort(), list(NA, NA),"pfc" = list(match("Log2FC", input$showCols), "desc"), "nfc" = list(match("Log2FC", input$showCols), "asc"), "signifD"= list(match("DEseq_adj_Log10pval", input$showCols), "asc"), "signifW"= list(match("Wilcox_adj_Log10pval", input$showCols), "asc"))
                  #defsort <- list(NA, NA)
                  if (is.na(defsort[1])) {
                          DT::datatable(data()[fltrow,input$showCols], selection = 'single',
                            extensions = 'Scroller', colnames = input$showCols,
                           options = list(dom = 'lpt', stateSave=T, lengthMenu = lengthlist),
                          rownames = F)
                  }else {
                    defsort[1] <- as.numeric(defsort[1]) - 1
                  DT::datatable(data()[fltrow,input$showCols], selection = 'single',
                            extensions = 'Scroller', colnames = input$showCols,
                           options = list(dom = 'lpt', stateSave=T, lengthMenu = lengthlist,order = list(defsort)),
                          rownames = F)                  
                  }
                  # ,
                  
                  # 

                  }
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

#output$myWebGL <- renderWebGL({
#    points3d(1:10, 1:10, 1:10)
#    axes3d()
#  })

}
