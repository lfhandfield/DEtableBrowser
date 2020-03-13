#' Server handler for detablebrowser
#'
#' @importFrom DT renderDataTable datatable
server <- function(input, output, session) {
  last.query.state <- reactiveVal("genelist")
  debug.state <- reactiveVal(0)
  value <- reactiveVal("")
  data <- reactiveVal("")
  mat <- reactiveVal("")
  simplesort <- reactiveVal("")
  plotgenes <- reactiveVal(c(""))
  dataclean <- reactiveVal(0)
  curflt <- reactiveVal(data.frame(criterion= c(), value=character())) 
  filtrow <- reactiveVal(c(T))
  
changeStyle <- function(p, plot.attribs, classprefix=""){
	library(ggplot2)
	if ("flags" %in% names(plot.attribs)) flags.plot <- plot.attribs[["flags"]]
	else flags.plot <- c()
	if ("title" %in% names(plot.attribs)) p <- p + ggtitle(plot.attribs[["title"]])
	if ("xlabel" %in% names(plot.attribs)) p <- p + xlab(plot.attribs[["xlabel"]])
	else if (("no.xlabel" %in% flags.plot)||(paste(classprefix,"no.xlabel",sep=".") %in% flags.plot))  p <- p + xlab(NULL)
	if ("ylabel" %in% names(plot.attribs)) p <- p + ylab(plot.attribs[["ylabel"]])
	else if (("no.ylabel" %in% flags.plot)||(paste(classprefix,"no.ylabel",sep=".") %in% flags.plot))  p <- p + ylab(NULL)
	
	if ("ylabel" %in% names(plot.attribs)) p <- p + ylab(plot.attribs[["ylabel"]])
	else if (("no.ylabel" %in% flags.plot)||(paste(classprefix,"no.ylabel",sep=".") %in% flags.plot))  p <- p + ylab(NULL)

	if ("no.legend" %in% flags.plot) p <- p + theme(legend.position="none")
	else{
		if ("legnames.size" %in% names(plot.attribs))  p <- p + theme(legend.text=element_text(size=plot.attribs[["legnames.size"]]))
	}

	# building xticks style
	themearg <- list()
	if (("no.xnamedticks" %in% flags.plot)||(paste(classprefix,"no.xnamedticks",sep=".") %in% flags.plot )) {
		p <- p + theme(axis.ticks.x=element_blank()) 
		p <- p + theme(axis.text.x=element_blank())
	}else{
		if ("no.xticks" %in% flags.plot) p <- p + theme(axis.ticks.x=element_blank())
		if ("no.xnames" %in% flags.plot) p <- p + theme(axis.text.x=element_blank())
		else if (!(("rot.xnames" %in% flags.plot)||(paste(classprefix,"rot.xnames",sep=".") %in% flags.plot))) themearg <- c(themearg,list(angle = 90, hjust = 1))
		if ("xnames.size" %in% names(plot.attribs)) themearg <- c(themearg,list(size= plot.attribs[["xnames.size"]]))
		if ("xnames.color" %in% names(plot.attribs)) themearg <- c(themearg,list(colour= plot.attribs[["xnames.color"]]))
		if ("xnames.bold" %in% flags.plot) themearg <- c(themearg,list(face="bold"))
	}
	if (length(themearg) > 0) p <- p + theme(axis.text.x= do.call(element_text, themearg))
	themearg <- list()
	if (("no.ynamedticks" %in% flags.plot)||(paste(classprefix,"no.ynamedticks",sep=".") %in% flags.plot )) {
		p <- p + theme(axis.ticks.y=element_blank()) 
		p <- p + theme(axis.text.y=element_blank())
	}else{
		if ("no.yticks" %in% flags.plot) p <- p + theme(axis.ticks.y=element_blank())
		if ("no.ynames" %in% flags.plot) p <- p + theme(axis.text.y=element_blank())
		else if ((("rot.ynames" %in% flags.plot)||(paste(classprefix,"rot.ynames",sep=".") %in% flags.plot))) themearg <- c(themearg,list(angle = 90, hjust = 1))
		if ("ynames.size" %in% names(plot.attribs)) themearg <- c(themearg,list(size= plot.attribs[["ynames.size"]]))
		if ("ynames.color" %in% names(plot.attribs)) themearg <- c(themearg,list(colour= plot.attribs[["ynames.color"]]))
		if ("ynames.bold" %in% flags.plot) themearg <- c(themearg,list(face="bold"))
	}
	if (length(themearg) > 0) p <- p + theme(axis.text.y= do.call(element_text, themearg))

	if ("scale.xrange" %in% names(plot.attribs)) p <- p + scale_x_continuous(limits = plot.attribs$scale.xrange)
	if ("scale.yrange" %in% names(plot.attribs)) p <- p + scale_y_continuous(limits = plot.attribs$scale.yrange)	

	if ("nb.col.legend" %in% names(plot.attribs)) p <- p + guides(color=guide_legend(ncol = plot.attribs$nb.col.legend),fill=guide_legend(ncol = plot.attribs$nb.col.legend))

	return(p)
}
plotDataGrid <- function(data, wdata= c(), xdata = c(), ydata =c(), transform=c(), plot.attribs=c(),do.zero.center=T, bgcolor = "#BBBBBB", do.cluster = c(T,T) ){
	library(ggplot2)
	if (class(data) != "list") {
		data <- list(data=data)
		if (is.null(wdata)) data$w <-matrix(1, dim(data$data)[1], dim(data$data)[2])
		else data$w <- dataw
		if (is.null(xdata)) data$x <-matrix(0.5, dim(data$data)[1], dim(data$data)[2])
		else data$x <- datax
		if (is.null(ydata)) data$y <-matrix(0.5, dim(data$data)[1], dim(data$data)[2])
		else data$y <- datay
	}else{
		if (!"data" %in% names(data)) stop("input list has 'data' as mandatory field")
		for( i in names(data)){
			if (!i %in% c("data", "x", "y", "w", "c1", "c2")) print(paste("The nknown field",i," is ignored, valid fields are \"data\", \"x\", \"y\", \"w\", \"c1\", \"c2\" only"))
		}
		
		if (!"w" %in% names(data)) data$w <- matrix(1, dim(data$data)[1], dim(data$data)[2])
		if (!"x" %in% names(data)) data$x <- matrix(0.5, dim(data$data)[1], dim(data$data)[2])
		if (!"y" %in% names(data)) data$y <- matrix(0.5, dim(data$data)[1], dim(data$data)[2])
	}
 	dd <- dim(data$data)
	cliprect <- c(0,0,dd[2],dd[1])

	if ((do.cluster[1])&&(nrow(data$data) > 1)){
	   dtable <- as.matrix(data$data); dtable[is.na(dtable)] <- 0 ; dtable[is.infinite(dtable)] <- 0
	   dareorder <- hclust(dist(data$data), method="complete")$order
	   for(i in names(data)) data[[i]] <- data[[i]][dareorder,,drop=F]
	}
	if ((do.cluster[2])&&(ncol(data$data) > 1)){
	   dtable <- as.matrix(data$data); dtable[is.na(dtable)] <- 0 ; dtable[is.infinite(dtable)] <- 0
	   dareorder <- hclust(dist(t(data$data)), method="complete")$order
	   for(i in names(data)) data[[i]] <- data[[i]][,dareorder,drop=F]
	}

	if (is.null(plot.attribs)) plot.attribs <- list(flags=c())
	if ("flags" %in% names(plot.attribs)) pflags <- plot.attribs[["flags"]]
	else pflags <- c()
	
	trformval <-c()
 	aurange <- range(as.vector(data$data),na.rm=T);
        if (do.zero.center){
                if (abs(aurange[1]) > abs(aurange[2])) aurange[2] <- -aurange[1]  
                else aurange[1] <- -aurange[2]
        }
	for(i in names(transform)){
	if (!i %in% c("data", "x", "y", "w")) stop("allowed transform fiels are \"data\", \"x\", \"y\", \"w\" only")
		if (transform[i] == "pval"){
			data[[i]] <- -0.8 + 0.09 / (0.05 + data[[i]])
			data[[i]][data[[i]] < 0.1] <- 0.1
		}else if (transform[i] == "log10pval"){
			data[[i]] <- 0.05 / (0.05 + exp(data[[i]]/log(10)))
		}else if (transform[i] == "lerfp1"){
			tmp <- as.matrix(log(data[[i]] + 1))
			m <- mean(tmp,na.rm=T)
			v <- var(as.vector(tmp),na.rm=T) ^-0.5
			if (is.na(v)){
				data[[i]] <- 0.5
			}else{
				data[[i]]  <- pnorm((tmp - m)*v)
			}
			if (i == "data") aurange <- c(0,1)
		}else if (transform[i] == "colwise.erf"){
			trformval <- matrix(0, 5, dd[2])
			for(j in 1:dd[2]){
				tmp <- data[[i]][,j]
				m <- mean(tmp,na.rm=T)
				v <- var(tmp,na.rm=T) ^-0.5
				if (is.na(v)){
					data[[i]][,j] <- 0.5
					trformval[,j] <- rep(m,5)
				}else{
					data[[i]][,j]  <- pnorm((tmp - m)*v)
					trformval[4,j] <- m - 0.5244005 / v
					trformval[5,j] <- m - 1.281552 / v
					trformval[3,j] <- m
					trformval[2,j] <- m + 0.5244005 / v
					trformval[1,j] <- m + 1.281552 / v
				}
			}
			if (i == "data") aurange <- c(0,1)
		}else if (transform[i] == "rowwise.erf"){
			trformval <- matrix(0, 5, dd[1])
			for(j in 1:dd[2]){
				tmp <- data[[i]][j,]
				m <- mean(tmp,na.rm=T)
				v <- var(tmp,na.rm=T) ^-0.5
				if (is.na(v)){
					data[[i]][j,] <- 0.5
					trformval[,j] <- rep(m,5)
				}else{
					data[[i]][j,]  <- pnorm((tmp - m)*v)
					trformval[4,j] <- m - 0.5244005 / v
					trformval[5,j] <- m - 1.281552 / v
					trformval[3,j] <- m
					trformval[2,j] <- m + 0.5244005 / v
					trformval[1,j] <- m + 1.281552 / v
				}
			}
			if (i == "data") aurange <- c(0,1)
		}else if (transform[i] == "threshold"){
			flt <- (data[[i]] <= 0.05)
flt[is.na(flt)] <- FALSE
			data[[i]][flt] <- 1;
			data[[i]][!flt] <- 0;
		}else if (transform[i] == "threshold2"){
			flt <- (data[[i]] <= 0.05)
			flt[is.na(flt)] <- FALSE
			data[[i]][flt] <- 0;
			data[[i]][!flt] <- 1;
		}

	}

	fgdata <- data.frame(row.names = 1:(dd[1] * dd[2] * 8))
	bgdata <- data.frame(row.names = 1:(dd[1] * dd[2] * 8))
  for(j in 1:dd[2]){
    for(i in 1:dd[1]){
	offset <- (i-1+ (j-1) * dd[1])*8
	fgdata$Log2FC[(offset+1):(offset+8)] <- rep(data$data[i,j],8)
	bgdata$Log2FC[(offset+1):(offset+4)] <- rep(0,4)
	bgdata$Log2FC[(offset+5):(offset+8)] <- rep(1,4)

	
	fgdata$I[(offset+1):(offset+8)] <- rep(offset/8,8)
	bgdata$I[(offset+1):(offset+4)] <- rep(offset/4,4)
	bgdata$I[(offset+5):(offset+8)] <- rep(1+offset/4,4)
	radata <- c(1 - sqrt(data$w[i,j]), 1 - sqrt(data$w[i,j]/2)) /2
	fgdata$X[offset+1] <- j - 1.0 + radata[1];		fgdata$Y[offset+1] <- i - 0.5;
 	fgdata$X[offset+2] <- j - 1.0 + radata[2];		fgdata$Y[offset+2] <- i - radata[2];
	fgdata$X[offset+3] <- j - 0.5;			fgdata$Y[offset+3] <- i - radata[1];
	fgdata$X[offset+4] <- j - radata[2];	fgdata$Y[offset+4] <- i - radata[2];
	fgdata$X[offset+5] <- j - radata[1];	fgdata$Y[offset+5] <- i - 0.5;
	fgdata$X[offset+6] <- j - radata[2]; 	fgdata$Y[offset+6] <- i - 1.0 + radata[2];
	fgdata$X[offset+7] <- j - 0.5;			fgdata$Y[offset+7] <- i -1.0 + radata[1];
	fgdata$X[offset+8] <- j - 1.0 + radata[2];		fgdata$Y[offset+8] <- i - 1.0 + radata[2];
	
	bgdata$X[offset+1] <- j -1.0;			bgdata$Y[offset+1] <- i - 1.0;
 	bgdata$X[offset+2] <- j -1.0 + data$x[i,j];		bgdata$Y[offset+2] <- i - 1.0;
	bgdata$X[offset+3] <- j -1.0 + data$x[i,j];		bgdata$Y[offset+3] <- i - 1.0 + data$y[i,j];
	bgdata$X[offset+4] <- j -1.0;			bgdata$Y[offset+4] <- i -1.0 + data$y[i,j];
	bgdata$X[offset+5] <- j -1.0 + data$x[i,j];		bgdata$Y[offset+5] <- i -1.0 + data$y[i,j];
	bgdata$X[offset+6] <- j;		 	bgdata$Y[offset+6] <- i -1.0 + data$y[i,j];
	bgdata$X[offset+7] <- j;			bgdata$Y[offset+7] <- i;
	bgdata$X[offset+8] <- j -1.0 + data$x[i,j];		bgdata$Y[offset+8] <- i;
    }
  }
 
	dabgcol <- rep(c(bgcolor), dd[1] * dd[2]*8)
  if ("c1" %in% names(data)){
    for(j in 1:dd[2]){
      for(i in 1:dd[1]){
        offset <- (i-1+ (j-1) * dd[1])*8
        dabgcol[offset+1] = data$c1[i,j]
        dabgcol[offset+2] = data$c1[i,j]
        dabgcol[offset+3] = data$c1[i,j]
        dabgcol[offset+4] = data$c1[i,j]
      }
    }
  }
  if ("c2" %in% names(data)){
    for(j in 1:dd[2]){
      for(i in 1:dd[1]){
        offset <- (i-1+ (j-1) * dd[1])*8
        dabgcol[offset+5] = data$c2[i,j]
        dabgcol[offset+6] = data$c2[i,j]
        dabgcol[offset+7] = data$c2[i,j]
        dabgcol[offset+8] = data$c2[i,j]
      }
    }
  }

	p <- ggplot(data = fgdata,mapping=aes(fill= Log2FC, group = I, y=Y, x=X) )
	p <- p + theme(axis.text.x=element_text(angle=90,vjust=0.5))
	p <- p + scale_x_discrete(limits= (1:dd[2])-0.5, labels= colnames(data$data) )# + xlab(NULL)
	p <- p + scale_y_discrete(limits= (1:dd[1])-0.5, labels= rownames(data$data) )# + ylab(NULL) 
	p <- p + geom_polygon(data=bgdata, mapping=aes(group = I, y=Y, x=X), fill = dabgcol)
	if (!is.null(trformval)){
		newbot = cliprect[2]  - (cliprect[4] / 10)
		cbdata <- data.frame(row.names = 1:164)
		for(i in 1:41){
			cbdata$X[i *4 -3] = 0; cbdata$X[i *4 -2] = 0; cbdata$X[i *4 -1] = cliprect[3]; cbdata$X[i *4] = cliprect[3];
			cbdata$Y[i *4 -3] = (newbot * i)/41; cbdata$Y[i *4 -2] = (newbot * (i-1))/41; cbdata$Y[i *4 -1] = (newbot * (i-1))/41; cbdata$Y[i *4] = (newbot * i)/41;
			cbdata$I[i *4 -3] = i; cbdata$I[i *4 -2] = i; cbdata$I[i *4 -1] = i; cbdata$I[i *4] = i;
			cbdata$Log2FC[i *4 -3] = (41-i)/40; cbdata$Log2FC[i *4 -2] = (41-i)/40; cbdata$Log2FC[i *4 -1] = (41-i)/40; cbdata$Log2FC[i *4] = (41-i)/40;
		}
		p <- p + geom_polygon(data=cbdata, mapping=aes(group = I, y=Y, x=X, fill=C))
		cliprect[2] = newbot
		cbdata <- data.frame(row.names = 1:(5*dd[2]))
		cbdata$Log2FC = 1:(5*dd[2])
		for(i in 1:dd[2]){
			cbdata$X[i *5 -4] = i-0.5; cbdata$X[i *5 -3] = i-0.5; cbdata$X[i *5 -2] = i-0.5; cbdata$X[i *5 -1] = i-0.5; cbdata$X[i *5] = i-0.5;
			cbdata$Y[i *5 -4] = newbot * 0.1; cbdata$Y[i *5 -3] = newbot * 0.3; cbdata$Y[i *5 -2] = newbot * 0.5; cbdata$Y[i *5 -1] = newbot * 0.7; cbdata$Y[i *5] = newbot * 0.9;
			cbdata$I[i *5 -4] = myFlFormat(trformval[1,i]); cbdata$I[i *5 -3] = myFlFormat(trformval[2,i]); cbdata$I[i *5 -2] = myFlFormat(trformval[3,i]); cbdata$I[i *5 -1] = myFlFormat(trformval[4,i]); cbdata$I[i *5] = myFlFormat(trformval[5,i]);
		}
		p <- p + geom_text(data=cbdata, mapping=aes(label = I, y=Y, x=X), color=rep(c("#FFFFFF"), dd[2]*5),size=rep(2,dd[2]*5))
	}


  p <- p + geom_polygon()
	p <- p + theme(axis.text=element_text(color="#000000",face="bold",size=10))
	p <- p + scale_fill_gradientn(colours=c("#00FFFF", "#00B0FF","#0079FF","#0000E8", "#000074","#000000","#4B0000","#960000","#E10000","#FF8000","#FFD600"),limits=aurange)
	p <- p + scale_size(limits=c(0,1))
	p <- p + coord_cartesian(xlim=c(cliprect[1],cliprect[3]),ylim=c(cliprect[2],cliprect[4]),expand=F)
	return(changeStyle(p,plot.attribs = plot.attribs))
}

  
  
  observe({
      if (dataclean() == 0){
       
      }else{
      updateSelectInput(session,"filter", choices = colnames(data()), selected = colnames(data())[1])

      updateSelectInput(session,"obs",choices = setdiff(colnames(data()), c("Gene", "DE", "Log2FC", "LogitAuroc", "Comparison", "Celltype", "Archtype", "TPMmean", "DEseq_Log10pval", "Wilcox_Log10pval", "DEseq_adj_Log10pval", "Wilcox_adj_Log10pval", "DESeq_basemean", "FAD_coverage", "Ctrl_coverage", "FAD_Log2FC_toEmpty", "Ctrl_Log2FC_toEmpty", "MeanLog2FC", "MeanLog2FC", "MeanLogitAuroc", "Nbgenes","ID", "Domain","Tail", "pvalue", "Test" )))
      helpstr <- "GOSlim"
      ext <- c("FullName", "GO", "GOslim", "Description", "Intersection")
      danames <- setdiff(colnames(data()), helpstr)
      # set tablecol filter, with default values
      updateCheckboxGroupInput(session,"showCols", choices = danames, selected = setdiff(danames, c(c("FullName", "GO", "GOslim", "Description", "Archtype", "FAD_Log2FC_toEmpty", "Ctrl_Log2FC_toEmpty","Alias", "LogitAuroc", "sample_ctrlIds", "sample_testIds"), ifelse(input$resfield == "gene", c("DEseq_Log10pval", "Wilcox_Log10pval"), c("DEseq_adj_Log10pval", "Wilcox_adj_Log10pval")))))
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
  
  # Load the current matrix for histogram plot
  observe({
        progress <- Progress$new(session, min=0)
        shinyjs::disable("resfield"); shinyjs::disable("dataset");shinyjs::disable("simplebutton")
        on.exit(progress$close())
        progress$set(message = 'Loading Table...',
        detail = 'This may take a few seconds')
            dastr <- switch(input$dataset, "MHBR",  "Fine Celltypes / Multinomial"  = "MH" , "Broad Celltypes / Multinomial" = "MHBR", "Scmap Celltypes / Multinomial" = "JaJn", "Fine Celltypes / Clustering"  = "MHS" , "Broad Celltypes / Clustering" = "MHBRS", "Scmap Celltypes / Clustering" = "JaJnS")
        mat(readRDS(paste("/lustre/scratch117/cellgen/team218/lh20/results/matrix_",dastr,".rds", sep="")))
        shinyjs::enable("resfield");shinyjs::enable("dataset");shinyjs::enable("simplebutton")
    })

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
       # plotgenes(strsplit(as.character(data()[dalist[input$results_state$start+ 1 + ifelse(length(input$results_rows_selected) == 0, 0, input$results_rows_selected)], "Intersection"]), "," )[1])
        #plotgenes(as.character(data()[dalist[(input$results_state$start+1): maxo], "Gene"]))
        #if length(input$results_rows_selected) == 0 
      }
    }
    }
    }) # filtrow()
  # 

  
  # Update filter value sets
  observe({
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
      }else{
        simplesort("signif")
      }
    tmp$value <- as.character(tmp$value)
    curflt(tmp)
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
  
  observe({
  output$map <- renderPlot({
    print(plotgenes())
      if (length(plotgenes()) > 1){
            curcolnames <- colnames(mat()$deseq$logpval)
            if (input$comtype == "All") colfilt <- rep(T, length(curcolnames))
            else if (input$comtype == "Pooled Comparisons") colfilt <- mat()$ispool[mat()$coltotest]
            else colfilt <- !mat()$ispool[mat()$coltotest]
            print(match(plotgenes(),rownames(mat()$deseq$logpval)))
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
                if (sum(fltrow) == 0) data.frame(row.names=c("Nothing"))
                else{
                  lengthlist = c(5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100, 120)
                  if (sum(fltrow) < 120) lengthlist = c(lengthlist[lengthlist < sum(fltrow)], sum(fltrow))
                  
                  defsort <- switch(simplesort(), c(NA, NA),"pfc" = c(match("Log2FC", input$showCols), "desc"), "nfc" = c(match("Log2FC", input$showCols), "asc"), "signif"= c(match("DEseq_adj_Log10pval", input$showCols), "asc"))

		  if (is.na(defsort[1])) defsort <- c()
                  else defsort[1] <- as.numeric(defsort[1]) + 1
 
                  datatable(data()[fltrow,input$showCols], selection = 'single',
                            #options = list(columnDefs = list(list(width = '70px', targets = c(2, 3, 4)), list(width = '10px', targets = c(0))), pageLength = 5, autoWidth = TRUE, dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel')),
                            extensions = 'Scroller', colnames = input$showCols,
                           options = list(dom = 'lpt', stateSave=T, lengthMenu = lengthlist,order = list(defsort)),
                          rownames = F)
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
