#QUEUE=yesterday
#NBCORES=2
#MEMGIGA=40


# Arguments 
# [1] this script name (arguments are shifted by one if used by Rscript...)
# [2] input selector string
# [3] output RDS file path
args <- commandArgs(TRUE)
argsall <- commandArgs(FALSE)
match <- grep("--file=", argsall)
if (length(match) > 0){
	path <- dirname(normalizePath(sub("--file=","",argsall[match])))
}else stop("run this using Rscript only")
rm(argsall)
rm(match)

s = as.numeric(args[1])
source("/nfs/users/nfs_l/lh20/Rcode.R")
library(Seurat)
library(methods)

print("reading")
print(args[1]);runargs <- readRDS(args[1])
print(args[2]);pair <- readRDS(args[2])
print(args[3]);wpair <- readRDS(args[3])
print(args[4]);base <- readRDS(args[4])

summary <- matrix(0,length(runargs$cel$celltypes) * length(names(runargs$cel$Comparisons)),12)
summary2 <- matrix(0,length(runargs$cel$celltypes) * length(names(runargs$cel$Comparisons)),12)

rownames(summary) <- gsub(" ", "_",paste(rep(runargs$cel$celltypes, length(names(runargs$cel$Comparisons)) ),rep(names(runargs$cel$Comparisons), each=length(runargs$cel$celltypes))))
colnames(summary) <- c("Q+","Q-","W+","W-", "Q+W+", "Q-W-", "Q+W-","Q-W+","Q+W_","Q-W_","Q_W+","Q_W-")
rownames(summary2) <- rownames(summary)
colnames(summary2) <- c("Q+","Q-","W+","W-", "Q+W+", "Q-W-", "Q+W-","Q-W+","Q+W_","Q-W_","Q_W+","Q_W-")

geneannot <- read.csv("/lustre/scratch117/cellgen/team218/lh20/geneannotation.tsv",sep="\t")
genecolnames = c("Comparison", "Archtype", "Celltype", "NAME", "DE", "Log2FC", "LogitAuroc", "DEseq_Log10pval", "Wilcox_Log10pval", "DEseq_adj_Log10pval", "Wilcox_adj_Log10pval","DESeq_basemean", "TPMmean", "FAD_coverage", "Ctrl_coverage", "FAD_Log2FC_toEmpty", "Ctrl_Log2FC_toEmpty", "Alias", "Description", "FullName", "GOslim", "GO")

genetable <- data.frame(row.names=genecolnames)
genetable2 <- data.frame(row.names=genecolnames)
genetable3 <- data.frame(row.names=c(genecolnames[c(1:15,18:22)], "DE_concat"))
genetable4 <- data.frame(row.names=c(genecolnames[c(1:15,18:22)], "DE_concat"))

log2FC_FAD_to_Empty <- data.frame(row.names = rownames(wpair$rawDE$wilcox.log10pval))
log2FC_CTRL_to_Empty <- data.frame(row.names = rownames(wpair$rawDE$wilcox.log10pval))


gotable <- data.frame(row.names=c("Domain", "ID", "Term", "Comparison", "Archtype","Celltype", "Tail", "pvalue", "Test", "Nbgenes", "MeanLog2FC", "MeanLogitAuroc", "DESeq_basemean", "TPMmean", "FAD_coverage", "Ctrl_coverage", "Intersection"))
gotable2 <- data.frame(row.names=c("Domain", "ID", "Term", "Comparison", "Archtype","Celltype", "Tail", "pvalue", "Test", "Nbgenes", "MeanLog2FC", "MeanLogitAuroc", "DESeq_basemean", "TPMmean", "FAD_coverage", "Ctrl_coverage", "Intersection"))

   for(ss in names(runargs$cel$Comparisons)){
	print(paste("Start to Process", ss))
	colnames(base[[paste(ss,"POS",sep='_')]]$deseq.log2FC) <- gsub(" ", "_",colnames(base[[paste(ss,"POS",sep='_')]]$deseq.log2FC))
	colnames(base[[paste(ss,"NEG",sep='_')]]$deseq.log2FC) <- gsub(" ", "_",colnames(base[[paste(ss,"NEG",sep='_')]]$deseq.log2FC))
	for(ct in runargs$cel$celltypes){
	  ct <- gsub(" ", "_", ct)
	  dn <- gsub(" ", "_",paste(ct,ss))

	  archt <- "Unkwown"
	  for(i in names(runargs$cel$archtype)){
		if (ct %in% gsub(" ", "_",runargs$cel$archtype[[i]])) archt <- i
	  }


	  print(paste("Processing", ct))
	  tmp <- paste("TOP", ct,sep="_")
	  if (tmp %in% names(pair$signifDE[[ss]])) {
		showlistQ <- pair$signifDE[[ss]][tmp][[1]][,"Corrected_deseq.log10pvalue",drop=F]
		print(paste(nrow(showlistQ),"+Q"))
		showlistQ <- rbind(showlistQ,pair$signifDE[[ss]][paste("BOT", ct,sep="_")][[1]][,"Corrected_deseq.log10pvalue",drop=F])
	  }else {showlistQ <- matrix(0,0,1); print("SoSad"); colnames(showlistQ) <- c("Corrected_deseq.log10pvalue")}
	  if (tmp %in% names(wpair$signifDE[[ss]])) {
		showlistW <- wpair$signifDE[[ss]][tmp][[1]][,"Corrected_wilcox.log10pval",drop=F]
		print(paste(nrow(showlistW),"+W"))
		showlistW <- rbind(showlistW,wpair$signifDE[[ss]][paste("BOT", ct,sep="_")][[1]][,"Corrected_wilcox.log10pval",drop=F])
	  }else {showlistW <- matrix(0,0,1); print("SoDas"); colnames(showlistW) <- c("Corrected_wilcox.log10pval")}
	  print("and now")
	 
	  tmp <- paste("TGO", ct,sep="_")
 	  gblock <- matrix(0,0,8); colnames(gblock) = c("pvalue","Domain","ID","Term","Enrich","intersection", "Tail", "Test")
	  if (tmp %in% names(pair$signifDE[[ss]])){
		  if (length(pair$signifDE[[ss]][[paste("TGO", ct,sep="_")]]) > 0){
	                  tmp <- pair$signifDE[[ss]][[tmp]]$list
			  gblock <- cbind(tmp, "+", "DESeq") ; colnames(gblock) <- c("pvalue","Domain","ID","Term","Enrich","intersection", "Tail", "Test")
		  }
		  
  		if (length(pair$signifDE[[ss]][[paste("BGO", ct,sep="_")]]) > 0){
		  tmp <- pair$signifDE[[ss]][[paste("BGO", ct,sep="_")]]$list
		  tmp <- cbind(tmp, "-", "DESeq");  
		  colnames(tmp) <- colnames(gblock)
		  gblock <- rbind(gblock, tmp);
		}
	 } 
	  print("and again")
	  if (paste("TGO", ct,sep="_") %in% names(wpair$signifDE[[ss]])){
	  	  if (length(wpair$signifDE[[ss]][[paste("TGO", ct,sep="_")]]) > 0){
                    tmp <- wpair$signifDE[[ss]][[paste("TGO", ct,sep="_")]]$list
		    tmp <- cbind(tmp, "+", "Wilcox");  
		    print(dim(tmp))
		    colnames(tmp) <- colnames(gblock)
		    gblock <- rbind(gblock, tmp);	  
		  }
		  
		  if (length(wpair$signifDE[[ss]][[paste("BGO", ct,sep="_")]]) > 0){
   		    tmp <- wpair$signifDE[[ss]][[paste("BGO", ct,sep="_")]]$list
		    tmp <- cbind(tmp, "-", "Wilcox"); 
		    print(dim(tmp));
	            colnames(tmp) <- colnames(gblock)
		    gblock <- rbind(gblock, tmp);
		  }
	  }

	  print("wellD") 

	  allnames <- unique(c(rownames(showlistQ),rownames(showlistW)))


	  if (!(dn %in% colnames(pair$rawDE$deseq.log10pvalue))){
		pair$rawDE$deseq.log10pvalue[,dn] <- rep(NA,nrow(pair$rawDE$deseq.log10pvalue))
		pair$rawDE$deseq.basemean[,dn] <- rep(NA,nrow(pair$rawDE$deseq.log10pvalue))
		pair$rawDE$deseq.log2FC[,dn] <- rep(NA,nrow(pair$rawDE$deseq.log10pvalue))
	  }

	  if (!(dn %in% colnames(wpair$rawDE$wilcox.log10pval))){
		wpair$rawDE$wilcox.log10pval[,dn] <- rep(NA,nrow(wpair$rawDE$wilcox.log10pval))
		wpair$rawDE$meanTPM[,dn] <- rep(NA,nrow(wpair$rawDE$wilcox.log10pval))
		wpair$rawDE$wilcox.logitAuroc[,dn] <- rep(NA,nrow(wpair$rawDE$wilcox.log10pval))
		wpair$rawDE$dropoutPosClass[,dn] <- rep(NA,nrow(wpair$rawDE$wilcox.log10pval))
		wpair$rawDE$dropoutNegClass[,dn] <- rep(NA,nrow(wpair$rawDE$wilcox.log10pval)) 
	  }

	  qmap <- pair$rawDE$deseq.log10pvalue[,dn] < -1.3; qmap[is.na(qmap)] <- F
	  wmap <- wpair$rawDE$wilcox.log10pval[,dn] < -1.3; wmap[is.na(wmap)] <- F
	  
	  print("wellE") 

	  allnames2 <- unique(c(rownames(pair$rawDE$deseq.log2FC)[qmap] ,rownames(wpair$rawDE$wilcox.log10pval)[wmap]))


	  ordterm <- rep(0, length(allnames2))
	  qmap <- match(allnames2 , rownames(pair$rawDE$deseq.log2FC))
	  wmap <- match(allnames2 , rownames(wpair$rawDE$wilcox.log10pval))
	  ordterm <- pair$rawDE$deseq.log2FC[allnames2,dn]
	  ordterm[is.na(ordterm)] <- 0
	  ordterm <- order(ordterm,decreasing=T)
#	  ordterm <- pair$signifDE[[ss]]$deseq.log2FC[qmap,ct] + wpair$signifDE[[ss]]$wilcox.logitAuroc[wmap,ct] * 3
#	  flt <- (allnames2 %in% rownames(showlistQ))&(allnames2 %in% rownames(showlistW))
#	  ordterm[flt] <- ordterm[flt] + sign(ordterm[flt]) * 10000
#	  ordterm <- order(ordterm,decreasing=T)
	  allnames2 <- allnames2[ordterm]
	  qmap <- qmap[ordterm]
	  wmap <- wmap[ordterm]

	  gmap <- match(allnames2 , geneannot$NAME)

	  print("wellE")
	  spmap <- match(rownames(wpair$rawDE$wilcox.log10pval), rownames(base[[paste(ss,"POS",sep='_')]]$deseq.log2FC))
	  snmap <- match(rownames(wpair$rawDE$wilcox.log10pval), rownames(base[[paste(ss,"NEG",sep='_')]]$deseq.log2FC))
  	  print("welldE") 
  	  print(sum(!is.na(spmap)))
          print(sum(!is.na(snmap)))
          print(length(allnames2))
 	  if (length(allnames2) > 0){ 
		print(paste(ss,"POS",sep='_'))
		print(dn)
		print(colnames(log2FC_FAD_to_Empty))
		log2FC_FAD_to_Empty[[dn]] <- NA
		log2FC_CTRL_to_Empty[[dn]] <- NA
		print(colnames(log2FC_FAD_to_Empty))

		print(length(base[[paste(ss,"POS",sep='_')]]$deseq.log2FC[spmap[!is.na(spmap)], ct]))
		print(length(base[[paste(ss,"NEG",sep='_')]]$deseq.log2FC[snmap[!is.na(snmap)], ct]))
		if (length(base[[paste(ss,"POS",sep='_')]]$deseq.log2FC[spmap[!is.na(spmap)], ct]) != 0) {
			log2FC_FAD_to_Empty[!is.na(spmap), dn ] <- base[[paste(ss,"POS",sep='_')]]$deseq.log2FC[spmap[!is.na(spmap)], ct]
			spmap <- match(allnames2 , rownames(base[[paste(ss,"POS",sep='_')]]$deseq.log2FC))
		}else spmap <- rep(NA, length(allnames2))
		
		if (length(base[[paste(ss,"NEG",sep='_')]]$deseq.log2FC[snmap[!is.na(snmap)], ct]) != 0) {
			log2FC_CTRL_to_Empty[!is.na(snmap), dn ] <- base[[paste(ss,"NEG",sep='_')]]$deseq.log2FC[snmap[!is.na(snmap)], ct]
			snmap <- match(allnames2 , rownames(base[[paste(ss,"NEG",sep='_')]]$deseq.log2FC))
		}else snmap <- rep(NA, length(allnames2))


	  

	  block <- data.frame(matrix(0,length(allnames2),length(genecolnames)))
	  colnames(block) <-  genecolnames
	 
	  print("welddlE") 
	  block$NAME <- allnames2;
	  block$Log2FC <- pair$rawDE$deseq.log2FC[qmap,dn]
	  block$LogitAuroc <- wpair$rawDE$wilcox.logitAuroc[wmap,dn]
	  block$DEseq_Log10pval <- pair$rawDE$deseq.log10pvalue[qmap,dn]
	  block$Wilcox_Log10pval <- wpair$rawDE$wilcox.log10pval[wmap,dn]

	  block$DESeq_basemean <- pair$rawDE$deseq.basemean[qmap,dn]
	  block$TPMmean <- wpair$rawDE$meanTPM[wmap,dn]
	  block$FAD_coverage <- wpair$rawDE$dropoutPosClass[wmap,dn]
	  block$Ctrl_coverage <- wpair$rawDE$dropoutNegClass[wmap,dn]

	  block$Alias[!is.na(gmap)] <- geneannot[gmap[!is.na(gmap)], "ALIAS"]
	  block$Description[!is.na(gmap)] <- geneannot[gmap[!is.na(gmap)], "DESCRIPTION"]
  	  block$FullName[!is.na(gmap)] <- geneannot[gmap[!is.na(gmap)], "FULLNAME"]
	  block$GOslim[!is.na(gmap)] <- geneannot[gmap[!is.na(gmap)], "GOSLIM"]
	  block$GO[!is.na(gmap)] <- geneannot[gmap[!is.na(gmap)], "GO"]

	  print("wellE") 
	  if (sum(!is.na(spmap)) != 0) block$FAD_Log2FC_toEmpty[!is.na(spmap)] <- base[[paste(ss,"POS",sep='_')]]$deseq.log2FC[spmap[!is.na(spmap)], ct]
	  if (sum(!is.na(snmap)) != 0) block$Ctrl_Log2FC_toEmpty[!is.na(snmap)] <- base[[paste(ss,"NEG",sep='_')]]$deseq.log2FC[snmap[!is.na(snmap)], ct]



	  block$DEseq_adj_Log10pval <- NA;  block$Wilcox_adj_Log10pval <- NA
	  map <- match(block$NAME , rownames(showlistQ))
	  block$DEseq_adj_Log10pval[!is.na(map)] <- showlistQ[map[!is.na(map)],1]
	  map <- match(block$NAME , rownames(showlistW))
	  block$Wilcox_adj_Log10pval[!is.na(map)] <- showlistW[map[!is.na(map)],1]
 	  block$DEseq_adj_Log10pval[block$DEseq_adj_Log10pval > -1.3] <- NA
	  block$Wilcox_adj_Log10pval[block$Wilcox_adj_Log10pval > -1.3] <- NA

	  block$DE <- ""
	  fakeDE <- rep("", length(block$DE))
	  flt <- (block$DEseq_adj_Log10pval < -1.3)&(block$Log2FC < 0); flt[is.na(flt)] <- F; block$DE[flt] <- "Q-"
	  flt <- (block$DEseq_adj_Log10pval < -1.3)&(block$Log2FC > 0); flt[is.na(flt)] <- F; block$DE[flt] <- "Q+"
	  flt <- (block$Wilcox_adj_Log10pval < -1.3)&(block$LogitAuroc < 0); flt[is.na(flt)] <- F; block$DE[flt] <- paste(block$DE[flt],"W-",sep="")
	  flt <- (block$Wilcox_adj_Log10pval < -1.3)&(block$LogitAuroc > 0); flt[is.na(flt)] <- F; block$DE[flt] <- paste(block$DE[flt],"W+",sep="") 
	  
	  flt <- (block$DEseq_Log10pval < -1.3)&(block$Log2FC < 0); flt[is.na(flt)] <- F; fakeDE[flt] <- "Q-"
	  flt <- (block$DEseq_Log10pval < -1.3)&(block$Log2FC > 0); flt[is.na(flt)] <- F; fakeDE[flt] <- "Q+"
	  flt <- (block$Wilcox_Log10pval < -1.3)&(block$LogitAuroc < 0); flt[is.na(flt)] <- F; fakeDE[flt] <- paste(fakeDE[flt],"W-",sep="")
	  flt <- (block$Wilcox_Log10pval < -1.3)&(block$LogitAuroc > 0); flt[is.na(flt)] <- F; fakeDE[flt] <- paste(fakeDE[flt],"W+",sep="")  

 
	  block$Celltype <- ct
	  block$Comparison <- ss
	  block$Archtype <- archt

	  showlistQ <- showlistQ[showlistQ[,1] <= -1.3 ,,drop=F]
	  showlistW <- showlistW[showlistW[,1] <= -1.3 ,,drop=F]
	  allnames <- unique(c(rownames(showlistQ),rownames(showlistW)))
	  flt <- block$NAME %in% allnames
	  print(ct)
	  concat <- paste(ct,ss, sep= '_')

	  print(sum(flt))
	  tmp2 <- table(fakeDE)
	   if (sum(flt) > 0){
	  	  tmp <- table(block$DE[flt])
		print(tmp)
	  summary[concat,"Q+"] = length(grep("Q\\+",block$DE[flt]))

	  summary[concat,"Q-"] = length(grep("Q\\-",block$DE[flt]))
	  summary[concat,"W+"] = length(grep("W\\+",block$DE[flt]))
	  summary[concat,"W-"] = length(grep("W\\-",block$DE[flt]))
	  summary[concat,"Q+"] = length(grep("Q\\+",block$DE[flt]))
	  summary[concat,"Q-"] = length(grep("Q\\-",block$DE[flt]))
	  summary[concat,"W+"] = length(grep("W\\+",block$DE[flt]))
	  summary[concat,"W-"] = length(grep("W\\-",block$DE[flt]))
  
	  summary[concat,"Q+W_"] = tmp["Q+"]
	  summary[concat,"Q-W_"] = tmp["Q-"]
	  summary[concat,"Q_W+"] = tmp["W+"]
	  summary[concat,"Q_W-"] = tmp["W-"]
	  }
		print(tmp2)	
	  summary2[concat,"Q+"] = length(grep("Q\\+",fakeDE))
	  summary2[concat,"Q-"] = length(grep("Q\\-",fakeDE))
	  summary2[concat,"W+"] = length(grep("W\\+",fakeDE))
	  summary2[concat,"W-"] = length(grep("W\\-",fakeDE))
	  summary2[concat,"Q+"] = length(grep("Q\\+",fakeDE))
	  summary2[concat,"Q-"] = length(grep("Q\\-",fakeDE))
	  summary2[concat,"W+"] = length(grep("W\\+",fakeDE))
	  summary2[concat,"W-"] = length(grep("W\\-",fakeDE))
  
	  summary2[concat,"Q+W_"] = tmp2["Q+"]
	  summary2[concat,"Q-W_"] = tmp2["Q-"]
	  summary2[concat,"Q_W+"] = tmp2["W+"]
	  summary2[concat,"Q_W-"] = tmp2["W-"]
	  for(j in c("Q+W+","Q+W-","Q-W+","Q-W-")) {
	        if (sum(flt) > 0) summary[concat,j] = tmp[j]
		summary2[concat,j] = tmp2[j]
		}


	
	  genetable <- rbind(genetable, block[block$NAME %in% allnames,,drop=F])
	  genetable2 <- rbind(genetable2, block[!(block$NAME %in% allnames),,drop=F])

	  }
	  if (nrow(gblock) > 0){
	  block2 <- data.frame(matrix(0,nrow(gblock),17))
	  colnames(block2) <-  c("Domain", "ID", "Term", "Comparison", "Archtype","Celltype", "Tail", "pvalue", "Test", "Nbgenes","MeanLog2FC", "MeanLogitAuroc", "DESeq_basemean", "TPMmean", "FAD_coverage", "Ctrl_coverage", "Intersection")


	  block2$Domain <- gblock$Domain
	  block2$ID <- gblock$ID
	  block2$Term <- gblock$Term
	  block2$Intersection <- gblock$intersection
	  block2$Test <- gblock$Test
	  block2$Tail <- gblock$Tail
	  block2$Comparison <- ss
	  block2$Celltype <- ct
	  block2$pvalue <- gblock$pvalue
	  for(j in 1:nrow(gblock)){
		curlist <- parselist(gblock[j,"intersection"])
		block2[j,c("Nbgenes")] <- length(curlist)
	  	block2[j,c("MeanLog2FC")] <- mean(pair$rawDE$deseq.log2FC[curlist,dn],na.rm=T)
		block2[j,c("MeanLogitAuroc")] <- mean(wpair$rawDE$wilcox.logitAuroc[curlist,dn],na.rm=T)
	  	block2[j,c("FAD_coverage")] <- mean(wpair$rawDE$dropoutPosClass[curlist,dn],na.rm=T)
		block2[j,c("Ctrl_coverage")] <- mean(wpair$rawDE$dropoutNegClass[curlist,dn],na.rm=T)
	  	block2[j,c("DESeq_basemean")] <- mean(pair$rawDE$deseq.basemean[curlist,dn],na.rm=T)
		block2[j,c("TPMmean")] <- mean(wpair$rawDE$meanTPM[curlist,dn],na.rm=T)
	 }
	
	block2$Archtype <- archt;

 
	  gotable<- rbind(gotable,block2) 
	   if (ncol(gotable) != 17) {
		print("Super Exit")
		print(colnames(gotable))

		print(ncol(gotable))
		print(colnames(block2))
		break
		}
	  }
	
	}
}


 for(ss in names(runargs$cel$Consensus)){
	print(paste("Start to Process", ss))
	for(ct in runargs$cel$celltypes){
	  ct <- gsub(" ", "_", ct)
	  dn <- gsub(" ", "_",paste(ct,runargs$cel$Consensus[[ss]]))
	  wn <- dn[!is.na(match(dn, colnames(wpair$rawDE$wilcox.logitAuroc)))] 
	  dn <- dn[!is.na(match(dn, colnames(pair$rawDE$deseq.log2FC)))]
	  print(wn)
	  print(dn)
	  archt <- "Unkwown"
	  for(i in names(runargs$cel$archtype)){
		if (ct %in% gsub(" ", "_",runargs$cel$archtype[[i]])) archt <- i
	  }


	  print(paste("Processing", ct))
	  tmp <- paste("TOP", ct,sep="_")
	  if (tmp %in% names(pair$consensusDE[[ss]])) {
		showlistQ <- pair$consensusDE[[ss]][tmp][[1]][,"corr_log10pval",drop=F]
		print(paste(nrow(showlistQ),"+Q"))
		showlistQ <- rbind(showlistQ,pair$consensusDE[[ss]][paste("BOT", ct,sep="_")][[1]][,"corr_log10pval",drop=F])
	  }else {showlistQ <- matrix(0,0,1); print("SoSad"); colnames(showlistQ) <- c("Corrected_deseq.log10pvalue")}
	  if (tmp %in% names(wpair$consensusDE[[ss]])) {
		showlistW <- wpair$consensusDE[[ss]][tmp][[1]][,"corr_log10pval",drop=F]
		print(paste(nrow(showlistW),"+W"))
		showlistW <- rbind(showlistW,wpair$consensusDE[[ss]][paste("BOT", ct,sep="_")][[1]][,"corr_log10pval",drop=F])
	  }else {showlistW <- matrix(0,0,1); print("SoDas"); colnames(showlistW) <- c("Corrected_wilcox.log10pval")}
	  print("and now")
	 
	  tmp <- paste("TGO", ct,sep="_")
	  gblock <- matrix(0,0,8); colnames(gblock) = c("pvalue","Domain","ID","Term","Enrich","intersection", "Tail", "Test")
	  if (tmp %in% names(pair$consensusDE[[ss]])){
		  if (length(pair$consensusDE[[ss]][paste("TGO", ct,sep="_")][[1]]) > 0){
		  	tmp <- pair$consensusDE[[ss]][tmp][[1]]$list
		  	gblock <- cbind(tmp, "+", "DESeq") ; colnames(gblock) <- c("pvalue","Domain","ID","Term","Enrich","intersection", "Tail", "Test")
		  }
		  if (length(pair$consensusDE[[ss]][paste("BGO", ct,sep="_")][[1]]) > 0){
			  tmp <- pair$consensusDE[[ss]][paste("BGO", ct,sep="_")][[1]]$list
			  tmp <- cbind(tmp, "-", "DESeq");  
			  colnames(tmp) <- colnames(gblock)
			  gblock <- rbind(gblock, tmp);
		  }
	  }
	  print("and again")
	  if (paste("TGO", ct,sep="_") %in% names(wpair$consensusDE[[ss]])){
	  	  if (length(wpair$consensusDE[[ss]][paste("TGO", ct,sep="_")][[1]]) > 0){
                    tmp <- wpair$consensusDE[[ss]][paste("TGO", ct,sep="_")][[1]]$list
		    tmp <- cbind(tmp, "+", "Wilcox");  
		    print(dim(tmp))
		    colnames(tmp) <- colnames(gblock)
		    gblock <- rbind(gblock, tmp);	  
		  }
		  
		  if (length(wpair$consensusDE[[ss]][paste("BGO", ct,sep="_")][[1]]) > 0){
   		    tmp <- wpair$consensusDE[[ss]][paste("BGO", ct,sep="_")][[1]]$list
		    tmp <- cbind(tmp, "-", "Wilcox"); 
		    print(dim(tmp));
	            colnames(tmp) <- colnames(gblock)
		    gblock <- rbind(gblock, tmp);
		  }
	  }

	  allnames <- unique(c(rownames(showlistQ),rownames(showlistW)))
	  ordterm <- rep(0, length(allnames))
	  qmap <- match(allnames , rownames(pair$rawDE$deseq.log2FC))
	  wmap <- match(allnames , rownames(wpair$rawDE$wilcox.log10pval))

	  ordterm <- pair$consensusDE[[ss]]$deseq.log2FC[qmap,ct] + wpair$consensusDE[[ss]]$wilcox.logitAuroc[wmap,ct] * 3
	  flt <- (allnames %in% rownames(showlistQ))&(allnames %in% rownames(showlistW))
	  ordterm[flt] <- ordterm[flt] + sign(ordterm[flt]) * 10000
	  ordterm <- order(ordterm,decreasing=T)
	  allnames <- allnames[ordterm]
	  qmap <- qmap[ordterm]
	  wmap <- wmap[ordterm]


	  gmap <- match(allnames , geneannot$NAME)
	  print(allnames[is.na(gmap)])

	  if (length(allnames) > 0){
	  block <- data.frame(matrix(0,length(allnames),ncol(genetable3)))
	  colnames(block) <-  colnames(genetable3)
	 
	  print("pop block") 

	  block$NAME <- allnames;
	  block$Log2FC <- Matrix::rowMeans(pair$rawDE$deseq.log2FC[qmap,dn,drop=F],na.rm=T)
	  block$LogitAuroc <- Matrix::rowMeans(wpair$rawDE$wilcox.logitAuroc[wmap,wn,drop=F],na.rm=T)
	  block$DEseq_Log10pval <- Matrix::rowMeans(pair$rawDE$deseq.log10pvalue[qmap,dn,drop=F],na.rm=T)
	  block$Wilcox_Log10pval <- Matrix::rowMeans(wpair$rawDE$wilcox.log10pval[wmap,wn,drop=F],na.rm=T)

	  block$DESeq_basemean <- Matrix::rowMeans(pair$rawDE$deseq.basemean[qmap,dn,drop=F],na.rm=T)
	  block$TPMmean <- Matrix::rowMeans(wpair$rawDE$meanTPM[wmap,wn,drop=F],na.rm=T)
	  block$FAD_coverage <- Matrix::rowMeans(wpair$rawDE$dropoutPosClass[wmap,wn,drop=F],na.rm=T)
	  block$Ctrl_coverage <- Matrix::rowMeans(wpair$rawDE$dropoutNegClass[wmap,wn,drop=F],na.rm=T)

	  block$Alias <- geneannot[gmap, "ALIAS"]
	  block$Description <- geneannot[gmap, "DESCRIPTION"]
  	  block$FullName <- geneannot[gmap, "FULLNAME"]
	  block$GOslim <- geneannot[gmap, "GOSLIM"]
	  block$GO <- geneannot[gmap, "GO"]

 
	  block$DEseq_adj_Log10pval <- NA;  block$Wilcox_adj_Log10pval <- NA
	  map <- match(block$NAME , rownames(showlistQ))
	  block$DEseq_adj_Log10pval[!is.na(map)] <- showlistQ[map[!is.na(map)],1]
	  map <- match(block$NAME , rownames(showlistW))
	  block$Wilcox_adj_Log10pval[!is.na(map)] <- showlistW[map[!is.na(map)],1]
 	  block$DEseq_adj_Log10pval[block$DEseq_adj_Log10pval > -1.3] <- NA
	  block$Wilcox_adj_Log10pval[block$Wilcox_adj_Log10pval > -1.3] <- NA

	  block$DE <- ""
	  fakeDE <- rep("", length(block$DE))
	  flt <- (block$DEseq_adj_Log10pval < -1.3)&(block$Log2FC < 0); flt[is.na(flt)] <- F; block$DE[flt] <- "Q-"
	  flt <- (block$DEseq_adj_Log10pval < -1.3)&(block$Log2FC > 0); flt[is.na(flt)] <- F; block$DE[flt] <- "Q+"
	  flt <- (block$Wilcox_adj_Log10pval < -1.3)&(block$LogitAuroc < 0); flt[is.na(flt)] <- F; block$DE[flt] <- paste(block$DE[flt],"W-",sep="")
	  flt <- (block$Wilcox_adj_Log10pval < -1.3)&(block$LogitAuroc > 0); flt[is.na(flt)] <- F; block$DE[flt] <- paste(block$DE[flt],"W+",sep="") 
	  
	  flt <- (block$DEseq_Log10pval < -1.3)&(block$Log2FC < 0); flt[is.na(flt)] <- F; fakeDE[flt] <- "Q-"
	  flt <- (block$DEseq_Log10pval < -1.3)&(block$Log2FC > 0); flt[is.na(flt)] <- F; fakeDE[flt] <- "Q+"
	  flt <- (block$Wilcox_Log10pval < -1.3)&(block$LogitAuroc < 0); flt[is.na(flt)] <- F; fakeDE[flt] <- paste(fakeDE[flt],"W-",sep="")
	  flt <- (block$Wilcox_Log10pval < -1.3)&(block$LogitAuroc > 0); flt[is.na(flt)] <- F; fakeDE[flt] <- paste(fakeDE[flt],"W+",sep="")  


	  block$DE_concat <- ""
	  for( subscr in runargs$cel$Consensus[[ss]]){
		map <- match(block$NAME, rownames((pair$signifDE[[subscr]][[paste("ORD",ct,sep="_")]])))
		map[!is.na(map)] <- sapply(pair$signifDE[[subscr]][[paste("ORD",ct,sep="_")]][map[!is.na(map)], "Corrected_deseq.log10pvalue"], function(x){if (x > -1.3) return(NA); return(T);})

		dadaname <- paste(ct, subscr, sep= '_')
		if (dadaname %in% names(pair$rawDE$deseq.log2FC)){
		map2 <- match(block$NAME, rownames(pair$rawDE$deseq.log10pvalue))
		print(length(map2))
		print(sum(!is.na(map2)))
		qres <- mapply(function(x,y,z){if (x > -1.3) return("__"); if (y < 0) { if (!is.na(z)) return("Q-"); return("_-")} else if (!is.na(z)) return("Q+"); return("_+")}, 
			pair$rawDE$deseq.log10pvalue[map2, dadaname], pair$rawDE$deseq.log2FC[map2,dadaname], map)
		}else{
		qres <- rep("na", nrow(block))
		}

		map <- match(block$NAME, rownames((wpair$signifDE[[subscr]][[paste("ORD",ct,sep="_")]])))
		map[!is.na(map)] <- sapply(wpair$signifDE[[subscr]][[paste("ORD",ct,sep="_")]][map[!is.na(map)], "Corrected_wilcox.log10pval"], function(x){if (x > -1.3) return(NA); return(T);})
		
		if (dadaname %in% names(pair$rawDE$deseq.log2FC)){
		map2 <- match(block$NAME, rownames(wpair$rawDE$wilcox.log10pval))
		print(length(map2))
		print(sum(!is.na(map2)))
		wres <- mapply(function(x,y,z){if (x > -1.3) return("__"); if (y < 0) { if (!is.na(z)) return("W-"); return("_-")} else if (!is.na(z)) return("W+"); return("_+")}, 
			wpair$rawDE$wilcox.log10pval[map2,dadaname], wpair$rawDE$wilcox.logitAuroc[map2,dadaname], map)
		}else{
		wres <- rep("na", nrow(block))
		}
		print(length(qres))
		print(length(wres))
		print(length(block$DE_concat))
		block$DE_concat <- mapply(function(x,y,z){return(paste(x,y,z,sep=""))}, block$DE_concat, qres, wres)
	  }

 
	  block$Celltype <- ct
	  block$Comparison <- ss
	  block$Archtype <- archt
	  showlistQ <- showlistQ[showlistQ[,1] <= -1.3 ,,drop=F]
	  showlistW <- showlistW[showlistW[,1] <= -1.3 ,,drop=F]
	  allnames <- unique(c(rownames(showlistQ),rownames(showlistW)))
	  flt <- block$NAME %in% allnames

	  genetable3 <- rbind(genetable3, block[block$NAME %in% allnames,,drop=F])
	  genetable4 <- rbind(genetable4, block[!(block$NAME %in% allnames),,drop=F])




	}

	  if (nrow(gblock) > 0){
	  block2 <- data.frame(matrix(0,nrow(gblock),17))
	  colnames(block2) <-  c("Domain", "ID", "Term", "Comparison", "Archtype","Celltype", "Tail", "pvalue", "Test", "Nbgenes","MeanLog2FC", "MeanLogitAuroc", "DESeq_basemean", "TPMmean", "FAD_coverage", "Ctrl_coverage", "Intersection")


	  block2$Domain <- gblock$Domain
	  block2$ID <- gblock$ID
	  block2$Term <- gblock$Term
	  block2$Intersection <- gblock$intersection
	  block2$Test <- gblock$Test
	  block2$Tail <- gblock$Tail
	  block2$Comparison <- ss
	  block2$Celltype <- ct
	  block2$pvalue <- gblock$pvalue
	  for(j in 1:nrow(gblock)){
		curlist <- parselist(gblock[j,"intersection"])
		block2[j,c("Nbgenes")] <- length(curlist)
		block2[j,c("MeanLog2FC")] <- mean(as.matrix(pair$rawDE$deseq.log2FC[curlist,dn,drop=F]),na.rm=T)
		block2[j,c("MeanLogitAuroc")] <- mean(as.matrix(wpair$rawDE$wilcox.logitAuroc[curlist,wn,drop=F]),na.rm=T)
	  	block2[j,c("FAD_coverage")] <- mean(as.matrix(wpair$rawDE$dropoutPosClass[curlist,wn,drop=F]),na.rm=T)
		block2[j,c("Ctrl_coverage")] <- mean(as.matrix(wpair$rawDE$dropoutNegClass[curlist,wn,drop=F]),na.rm=T)
	  	block2[j,c("DESeq_basemean")] <- mean(as.matrix(pair$rawDE$deseq.basemean[curlist,dn,drop=F]),na.rm=T)
		block2[j,c("TPMmean")] <- mean(as.matrix(wpair$rawDE$meanTPM[curlist,wn,drop=F]),na.rm=T)
	 }
	
	block2$Archtype <- archt;

 
	  gotable2<- rbind(gotable2,block2) 
	   if (ncol(gotable2) != 17) {
		print("Super Exit")
		print(colnames(gotable))

		print(ncol(gotable))
		print(colnames(block2))
		break
		}
	  }
	}
}

summary[is.na(summary)] <- 0
summary2[is.na(summary2)] <- 0

print("Saving output!")
saveRDS(list(gene=genetable,gene_preFDR=genetable2,go=gotable,summary= summary,summary_preFDR= summary2,consensus.gene=genetable3, consensus.gene_preFDR=genetable4,consensus.go=gotable2,DEseq= pair$rawDE, Wilcox= wpair$rawDE,log2FC_FAD_to_Empty=log2FC_FAD_to_Empty,log2FC_CTRL_to_Empty=log2FC_CTRL_to_Empty),args[5])

