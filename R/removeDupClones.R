#!/usr/bin/env Rscript
if (commandArgs()[1] != "RStudio") {
  
  ARGS <- c(
    "mutfile","character","file path of mutation file",
    "readfile","character","file path of summary file",
    "refseqfile","character","reference sequence file"
    
  )
  
  OPTS <- c(
    "tstart","numeric",0,"Start of reference to include",
    "tend","numeric",0,"End of reference to include",
    "j_thresh","numeric",0.9,"percent similarity threshold for clones to be repeats",
    "cores","numeric",4,""
  )
  
  source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
  }
  
  source_local("Rsub.R")
  source_local("SHMHelper.R")
  
  parseArgs("removeDupClones.R", ARGS, OPTS)
  
} else {
  mutfile <- "/Volumes//AltLab/SHM//Alt071-20140417/Results-Test/JKH058_Alt071//JKH058_Alt071_muts.txt"
  readfile <- "/Volumes//AltLab/SHM//Alt071-20140417/Results-Test/JKH058_Alt071//JKH058_Alt071_reads.txt"
  refseqfile <- "/Volumes/AltLab/SHM/Alt071-20140417/Reference/VB18_productive_reference.fas"
  tstart <- 0
  tend <- 0
  j_thresh <- .9
  cores <- 4
  
  source("~/SHMPipeline/R/Rsub.R")
  source("~/SHMPipeline/R/SHMHelper.R")
}

suppressPackageStartupMessages(library(plyr,quietly=TRUE))
suppressPackageStartupMessages(library(reshape2,quietly=TRUE))
suppressPackageStartupMessages(library(ggplot2,quietly=TRUE))
suppressPackageStartupMessages(library(Biostrings, quietly=TRUE))
suppressPackageStartupMessages(library(grid, quietly=TRUE))

refseq <- readDNAStringSet(refseqfile)

if (tstart == 0) {
  tstart <- 1
}
if (tend == 0) {
  tend <- nchar(refseq)
}

muts <- read.delim(mutfile,header=F,as.is=T,col.names=c("Expt","Read","Pos","Type","From","To","Size","End","Ins"))
reads <- read.delim(readfile,header=F,as.is=T,col.names=c("Expt","Read","Bp","Coords","Dup"))
muts$Type <- factor(muts$Type,levels=c("sub","del","ins"))



reads$Index <- sapply(1:nrow(reads),function(i,rs,ms) {
  read <- rs[i,]
  readMuts <- getMutsFromRead(read,ms)
  readMutsTable <- table(readMuts$Type)
  return(readMutsTable["sub"] + 2*readMutsTable["del"] + 2*readMutsTable["ins"])
}, reads, muts)

reads <- reads[rev(order(reads$Index)),]
reads$Dup <- ""

mutmat <- createMutationMatrix(reads,muts,refseq,tstart,tend)
insmat <- createInsertionMatrix(reads,muts,refseq,tstart,tend)



jaccard <- laply(mclapply(1:nrow(reads[reads$Index > 0,]),function(i,reads,mutmat,insmat) {
  jaccard <- rep(0,nrow(reads))
  if (i<2) return(jaccard)
  for (j in 1:(i-1)) {
    if (reads$Index[i] == 0 || reads$Index[j] == 0) {
      jaccard[j] <- 0
    } else {
      comparemuts <- compareMutations(rbind(mutmat[i,],mutmat[j,],insmat[i,],insmat[j,]))
      intersection <- comparemuts[1,] + comparemuts[2,] - comparemuts[3,]
      union <- comparemuts[3,]
      if (sum(intersection)==0) {
        jaccard[j] <- 0
      } else {
        jaccard[j] <- sum(union)/sum(intersection)
      }
    }
  }
  return(jaccard)
},reads[reads$Index > 0,], mutmat, insmat, mc.cores=cores),identity)

colnames(jaccard) <- reads$Read[reads$Index > 0]
rownames(jaccard) <- reads$Read[reads$Index > 0]

reads$Dup[reads$Index > 0] <- unlist(mclapply(1:nrow(reads[reads$Index > 0,]),function(i,reads,jaccard) {
  j_max <- max(jaccard[i,])
  if (j_max > j_thresh) {
    return(reads$Read[which.max(jaccard[i,])])
  } else {
    return("")
  } 
},reads[reads$Index > 0,],jaccard,mc.cores=cores))

for (i in 1:nrow(reads)) {
  if (reads$Dup[i] != "") {
    while (reads$Dup[reads$Read==reads$Dup[i]] != "") {
      reads$Dup[i] <- reads$Dup[reads$Read==reads$Dup[i]]
    }
  }
}


reads$Index <- 1:nrow(reads)

sorted_reads <- reads[0,]
clusters <- data.frame(id=numeric(),x=numeric(),y=numeric())



for (i in hclust(as.dist(1-jaccard),method="single")$order) {
  if (reads$Dup[i] != "") next
  dup_reads <- reads[reads$Dup == reads$Read[i],]
  
  sorted_reads <- rbind(sorted_reads,reads[i,],dup_reads)
  c2 <- nrow(sorted_reads) + 0.5
  c1 <- c2 - (nrow(dup_reads) + 1)
  if (nrow(dup_reads) > 0) clusters <- rbind(clusters,data.frame(id=rep(i,3),x=c(c1,c1,c2),y=c(c1,c2,c2)))
}



j_reflected <- t(jaccard)
j_reflected[lower.tri(j_reflected)] <- jaccard[lower.tri(jaccard)]

sorted_jaccard <- j_reflected[sorted_reads$Index,sorted_reads$Index]
colnames(sorted_jaccard) <- 1:nrow(sorted_jaccard)
rownames(sorted_jaccard) <- 1:nrow(sorted_jaccard)

sorted_reads <- rbind(sorted_reads,reads[! reads$Read %in% sorted_reads$Read,])
sorted_reads$Index <- NULL

write.table(sorted_reads,readfile,quote=F,sep="\t",na="",row.names=F,col.names=F)


pdf(sub(".txt","_similarity.pdf",readfile))

  jaccard.m <- melt(sorted_jaccard)
  jaccard.m$value <- ifelse(jaccard.m$Var2 >= jaccard.m$Var1,NA,jaccard.m$value)

  hm <- ggplot(jaccard.m, aes(y=Var1,x=Var2)) + scale_y_reverse() + geom_tile(aes(fill=value)) + scale_fill_gradient(low="blue",high="red") + guides(fill=F)
  if (nrow(clusters) > 0) hm <- hm + geom_polygon(aes(x=x,y=y,group=id),data=clusters,color="black",alpha=0,size=0.25)
  print(hm)
  
  breaks <- seq(0,1,length.out=20)
  rects <- data.frame(starts=breaks[1:(length(breaks)-1)])
  rects$ends <- rects$starts + diff(breaks)
  
  j_maxes <- data.frame(jmax=apply(jaccard,1,max))

 
  p_hist <- ggplot() + geom_rect(aes(xmin=starts,xmax=ends,ymin=-Inf,ymax=Inf,fill=starts),rects) + scale_fill_gradient(low="blue",high="red") + guides(fill=F) + geom_histogram(aes(x=jmax),j_maxes,breaks=breaks) 
  print(p_hist,vp=viewport(1,1,.45,.45,just=c("right","top")))
  
  
dev.off()


