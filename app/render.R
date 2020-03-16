
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

  