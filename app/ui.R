
makeOverlay <- function(overdata, gene, compset){
  library(ggplot2)
  gglist <- list()
  aurange <- range(as.vector(overdata$dematrices[[gene]][,compset]),na.rm=T)
  if (abs(aurange[1]) > abs(aurange[2])) {
    daccrange = 1:(21+ floor(-20 *aurange[2] /aurange[1]))
    aurange[2] <- -aurange[1]
  }else{
    daccrange = (21- floor(-20 *aurange[1] /aurange[2])):41
    aurange[1] <- -aurange[2]
  }
  value(aurange)
  
daccrange <- colorRampPalette(c("#00FFFF", "#00B0FF","#0079FF","#0000E8", "#000074","#000000","#4B0000","#960000","#E10000","#FF8000","#FFD600"))(41)[daccrange]

  
  for(flist in 1:length(compset)){
    gdata <- data.frame(row.names = rownames(overdata$coords))
    gdata$X <- overdata$coords[,1]; gdata$Y <- overdata$coords[,2]
    tmp <- overdata$partition@.Data
    tmp <- overdata$dematrices[[gene]][tmp,compset[flist]] 
    tmp[(!overdata$dropout[, gene]) ] <- NA

    sampleset <- which(overdata[,compset[[flist]]])
    bg <- !(overdata$sample %in% sampleset)
    tmp[bg] <- NA
    gdata$C <- tmp
  
    p <- ggplot(gdata, aes(x=X,y=Y,fill=C, alpha=C)) + geom_point();
    p <- p + scale_color_gradientn(name=transform,colours=colpal, na.value= "#BBBBBB")
    p <- p + scale_alpha_continuous(position=NULL,guide="none", na.value=0.25, range = c(1, 1))
   
    gglist <- c(gglist,changeStyle(p, list(title=gene))) 
  }
return(gglist[[1]])}
