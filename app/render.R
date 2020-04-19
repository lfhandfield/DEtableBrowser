
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
getGGlegend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  return(tmp$grobs[[leg]])}

grid_arrange_shared_legend <- function(plots, ncol = length(plots), nrow = 1, position = c("bottom", "right"), do.share.legend=T,do.newpage=T, main.title=c()) {
  library(ggplot2)
  library(gridExtra)
  library(grid)
  position <- match.arg(position)
  cur_legend <- getGGlegend(plots[[1]])
  if (do.share.legend){
    lheight <- sum(cur_legend$height)
    lwidth <- sum(cur_legend$width)
    gl <- lapply(plots, function(x) x + theme(legend.position="none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    combined <- switch(position,
                       "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                              cur_legend, top = main.title,
                                              ncol = 1,
                                              heights = unit.c(unit(1, "npc") - lheight, lheight)),
                       "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                             cur_legend, top = main.title, 
                                             ncol = 2,
                                             widths = unit.c(unit(1, "npc") - lwidth, lwidth))
    )
  }else{
    lheight <- sum(cur_legend$height)
    lwidth <- sum(cur_legend$width)
    gl <- lapply(plots, function(x) x )
    gl <- c(gl, ncol = ncol, nrow = nrow)
    combined <- arrangeGrob(gl, top= main.title)
  }
  if (do.newpage) grid.newpage()
  grid.draw(combined)
  # return gtable invisibly
  invisible(combined)
}

makeOverlay <- function(overdata, gene, compset, titles, gridsize){
  library(ggplot2)
  logjs(compset)
  aurange <- as.vector(overdata$dematrices[[gene]][,compset])
  logjs(dim(overdata$dematrices[[gene]][,compset]))
  logjs(aurange)
  frange <- range(aurange[!is.infinite(aurange)],na.rm=T)
  logjs(frange)
  aurange <- range(aurange,na.rm=T)
  if (is.infinite(aurange[1])) aurange[1] <- ifelse((frange[1] > 0), -1, frange[1]-1) 
  if (is.infinite(aurange[2])) aurange[2] <- ifelse((frange[2] < 0),  1, frange[2]+1) 
  logjs(aurange)
  
  if (abs(aurange[1]) > abs(aurange[2])) {
    daccrange = 1:(21+ floor(-20 *aurange[2] /aurange[1]))
    aurange[2] <- -aurange[1]
  }else if (aurange[1] != aurange[2]){
    daccrange = (21- floor(-20 *aurange[1] /aurange[2])):41
    aurange[1] <- -aurange[2]
  }else if (aurange[1] != 0){
    if (aurange[1] > 0){daccrange = 1:21; aurange[1] <- -aurange[2]}
    else{daccrange = 21:41; aurange[2] <- -aurange[1]}
  }else{daccrange = 1:41; aurange <- c(-1,1)}
  
  daccrange <- colorRampPalette(c("#00FFFF","#00B0FF","#0079FF","#0000E8","#000074","#000000","#4B0000","#960000","#E10000","#FF8000","#FFD600"))(41)[daccrange]
  
  gglist <- list()
  for(flist in 1:length(compset)){
    flt <- overdata$comptosmpls[, compset[flist]] != 0
  gdata <- data.frame(row.names = rownames(overdata$coords)[flt])
  gdata$X <- overdata$coords[flt,1]; gdata$Y <- overdata$coords[flt,2]
  frange <- overdata$dematrices[[gene]][,compset[flist]]
  logjs(paste(flist, compset[flist]))
  logjs(frange)
  frange[frange < aurange[1]] <- aurange[1]; frange[frange > aurange[2]] <- aurange[2]
#  tmp <- frange[overdata$partition@.Data]
  #tmp[(!overdata$dropout[, gene]) ] <- NA
  
  #    sampleset <- which(overdata[,compset[[flist]]])
  #    bg <- !(overdata$sample %in% sampleset)
  #    tmp[bg] <- NA
  degr <- sapply(overdata$dematrices[[gene]][,compset[flist]],function(x){return(ifelse(x==0,0.125,1))})
  gdata$C <- overdata$partition@.Data[flt] # tmp
  gdata$A <- degr[gdata$C]

  logjs(degr)
  logjs(table(gdata$A))
  p <- ggplot(gdata, aes(x=X,y=Y,color=C, alpha=A)) + geom_point();
#  p <- p + scale_color_gradientn(name=transform,colours=daccrange, na.value= "#BBBBBB")
   p <- p + scale_alpha_continuous(position=NULL,guide="none", na.value=0.25, range = c(0, 1))
   gglist <- c(gglist,list(changeStyle(p, list(title=titles[flist]))))
  }
return(grid_arrange_shared_legend(gglist, nrow =gridsize[1],ncol =gridsize[2], position = "right",main.title = paste("Cells supporting",gene,"as DE by Wilcox test")))}


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
plotDataGrid <- function(data, wdata= c(), xdata = c(), ydata =c(), transform=c(), colcolors=c() , plot.attribs=c(),do.zero.center=T, bgcolor = "#BBBBBB", do.cluster = c(T,T) ){
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
      if (!(i %in% c("data", "x", "y", "w", "c1", "c2"))) print(paste("The nknown field",i," is ignored, valid fields are \"data\", \"x\", \"y\", \"w\", \"c1\", \"c2\" only"))
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
  shapeno <- 0
  for(j in 1:dd[2]){
    for(i in 1:dd[1]){
      offset <- (i-1+ (j-1) * dd[1])*8
      fgdata$Log2FC[(offset+1):(offset+8)] <- rep(data$data[i,j],8)
      bgdata$Log2FC[(offset+1):(offset+4)] <- rep(0,4)
      bgdata$Log2FC[(offset+5):(offset+8)] <- rep(1,4)
      
      
      
      
      bgdata$I[(offset+1):(offset+4)] <- rep(offset/4,4)
      bgdata$I[(offset+5):(offset+8)] <- rep(1+offset/4,4)
      bgdata$X[offset+1] <- j -1.0;			bgdata$Y[offset+1] <- i - 1.0;
      bgdata$X[offset+2] <- j -1.0 + data$x[i,j];		bgdata$Y[offset+2] <- i - 1.0;
      bgdata$X[offset+3] <- j -1.0 + data$x[i,j];		bgdata$Y[offset+3] <- i - 1.0 + data$y[i,j];
      bgdata$X[offset+4] <- j -1.0;			bgdata$Y[offset+4] <- i -1.0 + data$y[i,j];
      bgdata$X[offset+5] <- j -1.0 + data$x[i,j];		bgdata$Y[offset+5] <- i -1.0 + data$y[i,j];
      bgdata$X[offset+6] <- j;		 	bgdata$Y[offset+6] <- i -1.0 + data$y[i,j];
      bgdata$X[offset+7] <- j;			bgdata$Y[offset+7] <- i;
      bgdata$X[offset+8] <- j -1.0 + data$x[i,j];		bgdata$Y[offset+8] <- i;
      
      if (data$w[i,j] < 0.5){
        fgdata$I[(offset+1):(offset+4)] <- shapeno
        fgdata$I[(offset+5):(offset+8)] <- shapeno +1
        shapeno <- shapeno + 2
        fgdata$X[offset+1] <- j - 0.23;	fgdata$Y[offset+1] <- i - 0.27;
        fgdata$X[offset+2] <- j - 0.27;	fgdata$Y[offset+2] <- i - 0.23;
        fgdata$X[offset+3] <- j - 0.77;	fgdata$Y[offset+3] <- i - 0.73;
        fgdata$X[offset+4] <- j - 0.73;	fgdata$Y[offset+4] <- i - 0.77;
        fgdata$X[offset+5] <- j - 0.23;	fgdata$Y[offset+5] <- i - 0.73;
        fgdata$X[offset+6] <- j - 0.27;	fgdata$Y[offset+6] <- i - 0.77;
        fgdata$X[offset+7] <- j - 0.77;	fgdata$Y[offset+7] <- i - 0.27;
        fgdata$X[offset+8] <- j - 0.73;	fgdata$Y[offset+8] <- i - 0.23;
      }else{
        fgdata$I[(offset+1):(offset+8)] <- shapeno
        shapeno <- shapeno + 1
        
        radata <- c(1 - sqrt(data$w[i,j]), 1 - sqrt(data$w[i,j]/2)) /2
        fgdata$X[offset+1] <- j - 1.0 + radata[1];		fgdata$Y[offset+1] <- i - 0.5;
        fgdata$X[offset+2] <- j - 1.0 + radata[2];		fgdata$Y[offset+2] <- i - radata[2];
        fgdata$X[offset+3] <- j - 0.5;			fgdata$Y[offset+3] <- i - radata[1];
        fgdata$X[offset+4] <- j - radata[2];	fgdata$Y[offset+4] <- i - radata[2];
        fgdata$X[offset+5] <- j - radata[1];	fgdata$Y[offset+5] <- i - 0.5;
        fgdata$X[offset+6] <- j - radata[2]; 	fgdata$Y[offset+6] <- i - 1.0 + radata[2];
        fgdata$X[offset+7] <- j - 0.5;			fgdata$Y[offset+7] <- i -1.0 + radata[1];
        fgdata$X[offset+8] <- j - 1.0 + radata[2];		fgdata$Y[offset+8] <- i - 1.0 + radata[2];
      }
    }
  }
  if (!is.null(colcolors)){
    agdata <- data.frame(row.names = 1:(dd[2] * 8))
    daagcol <- rep("#FFFFFF", dd[2] * 8)
    for(j in 1:dd[2]){
      offset <- (j-1) *8
      daagcol[(offset+1) : (offset+8)] <- colcolors[j] 
      agdata$I[(offset+1):(offset+4)] <- rep(offset/4,4)
      agdata$I[(offset+5):(offset+8)] <- rep(1+offset/4,4)
      agdata$X[offset+1] <- j -1.0;	agdata$Y[offset+1] <- -2.0;
      agdata$X[offset+2] <- j;		  agdata$Y[offset+2] <- -2.0;
      agdata$X[offset+3] <- j;      agdata$Y[offset+3] <- -1.0;
      agdata$X[offset+4] <- j -1.0; agdata$Y[offset+4] <- -1.0;
      agdata$X[offset+5] <- j -1.0;	agdata$Y[offset+5] <- dd[1];
      agdata$X[offset+6] <- j;		 	agdata$Y[offset+6] <- dd[1];
      agdata$X[offset+7] <- j;			agdata$Y[offset+7] <- dd[1]+1;
      agdata$X[offset+8] <- j -1.0;	agdata$Y[offset+8] <- dd[1]+1;
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
  #if (!is.null(colcolors)) p <- p + geom_polygon(data=bgdata, mapping=aes(group = I, y=Y, x=X), fill = daagcol)
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


plotLabels <- function(xdata, ydata, names = c(), color = c(), alpha = c(), filter = c(), xlim=c(), ylim=c(),plot.attribs=c(),label.size=4,point.size=4){
  library(ggplot2)
  if (is.null(plot.attribs)) plot.attribs <- list()
  if (!is.null(xlim)) plot.attribs$xlim = xlim
  if (!is.null(ylim)) plot.attribs$ylim = ylim
  if (is.null(filter)) filter <- rep(TRUE,length(xdata))
  
  if ("label.size" %in% names(plot.attribs)) label.size = plot.attribs$label.size
  if ("flags" %in% names(plot.attribs)) flags.plot <- plot.attribs[["flags"]]
  else flags.plot <- c()
  
  if (is.null(names)) names <- names(xdata);
  if (is.null(color)) color <- rep("#000000", length(xdata))
  if (is.null(alpha)) alpha <- rep(1, length(xdata))
  
  subflt <- filter & (names != "")
  hehe <- cbind(xdata[subflt],ydata[subflt],names[subflt], alpha[subflt])
  rownames(hehe) <- NULL
  colnames(hehe) <- c("X","Y","Label", "Alpha")
  hehe <- data.frame(hehe)
  hehe[,1] <- as.numeric(as.character(hehe[,1]))
  hehe[,2] <- as.numeric(as.character(hehe[,2]))
  hehe[,4] <- as.numeric(as.character(hehe[,4]))
  
  p <- ggplot(hehe, aes(x=X,y=Y))
  if ("xlim" %in% names(plot.attribs)) p <- p + scale_x_continuous(limits = plot.attribs$xlim)
  if ("ylim" %in% names(plot.attribs)) p <- p + scale_y_continuous(limits = plot.attribs$ylim)
  
  #colscale <- unique(color)
  #names(colscale) <- colscale
  #       p <- p + scale_color_manual(name="color",values=colscale)
  
  p <- p + geom_text(aes(label=Label,alpha=Alpha),size=label.size, color = color[subflt])
  
  subflt <- filter & (names == "")
  hehe <- cbind(xdata[subflt],ydata[subflt],names[subflt], alpha[subflt])
  rownames(hehe) <- NULL
  colnames(hehe) <- c("X","Y","Label", "Alpha")
  hehe <- data.frame(hehe)
  hehe[,1] <- as.numeric(as.character(hehe[,1]))
  hehe[,2] <- as.numeric(as.character(hehe[,2]))
  hehe[,4] <- as.numeric(as.character(hehe[,4]))
  p <- p + geom_point(mapping= aes(label=Label,alpha=Alpha), data=hehe,size=point.size, color = color[subflt])
return(changeStyle(p,plot.attribs))}

