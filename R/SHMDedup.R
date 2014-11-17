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
  
  parseArgs("SHMDedup.R", ARGS, OPTS)
  
} else {
  mutfile <- "/Volumes//AltLab/SHM/Alt102/results-robin//JKH109_Alt102/JKH109_Alt102_muts.txt"
  readfile <- "/Volumes//AltLab/SHM/Alt102/results-robin//JKH109_Alt102/JKH109_Alt102_reads.txt"
  refseqfile <- "/Volumes/AltLab/SHM/Alt102/ref//VB18_productive_reference.fas"
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

muts <- read.delim(mutfile,header=F,colClasses=mutsColClasses(),as.is=T,col.names=c("Expt","Read","Pos","Type","From","To","Size","End","Ins"))

if (nrow(muts) < 1) quit()

reads <- read.delim(readfile,header=F,as.is=T,col.names=c("Expt","Read","Bp","Coords","Dup"))
muts$Type <- factor(muts$Type,levels=c("sub","del","ins"))

reads <- reads[reads$Bp > 0,]

# The Index is simply a measure of how mutated the read is.
# Each sub adn insertion counts once, each del counts 2x.
# We will use it to sort on


indextable <- table(muts$Read,muts$Type)
reads$Index <- 0
reads$Index[match(rownames(indextable),reads$Read)] <- indextable[,"sub"] + 2*indextable[,"del"] + indextable[,"ins"]

reads <- reads[rev(order(reads$Index)),]
reads$Dup <- ""

cleanreads <- reads[reads$Index == 0,]
reads <- reads[reads$Index > 0,]

if (nrow(reads) < 2) quit()

mutmat <- createMutationMatrix(reads,muts,refseq,tstart,tend)

insmat <- createInsertionMatrix(reads,muts,refseq,tstart,tend)

t1 <- proc.time()
jaccard <- laply(mclapply(1:nrow(reads),function(i,nreads,scores,mutmat,insmat) {
  jaccard <- rep(0,nreads)
  if (i<2) return(jaccard)
  for (j in 1:(i-1)) {
    # Only consider reads that have mutation numbers within 50% of each other
    if (scores[i]/scores[j] > 2 || scores[j]/scores[i] > 2) {
      jaccard[j] <- 0
    } else {
      mut.comp <- rowSums(mapply(function(m1,m2,i1,i2) {
        mut.int <- 0
        mut.union <- 0
        ins.int <- 0
        ins.union <- 0
        if (any(grepl("(^$|^-$)",c(m1,m2)))) return(c(0,0))
        mut.union <- mut.union + sum(grepl("[ACGT<>]",c(m1,m2)))
        ins.union <- ins.union + sum(grepl("[ACGT]",c(i1,i2)))
        if (mut.union > 0 && m1==m2) {
          mut.union <- mut.union - 1
          mut.int <- 1
        }
        if (ins.union > 0 && i1==i2) {
          ins.union <- ins.union - 1
          ins.int <- 1
        }
        return(c(mut.union + ins.union,mut.int + ins.int))
      },mutmat[i,],mutmat[j,],insmat[i,],insmat[j,]))
      
      jaccard[j] <- ifelse(mut.comp[1]==0,0,mut.comp[2]/mut.comp[1])
    }
  }
  
  return(jaccard)
},nrow(reads), reads$Index, mutmat, insmat, mc.cores=cores),identity)
t2 <- proc.time() - t1

# 
# t1 <- proc.time()
# jaccard <- laply(mclapply(1:nrow(reads),function(i,reads,mutmat,insmat) {
#   jaccard <- rep(0,nrow(reads))
#   if (i<2) return(jaccard)
#   for (j in 1:(i-1)) {
#     if (reads$Index[i] == 0 || reads$Index[j] == 0) {
#       jaccard[j] <- 0
#     } else {
#       comparemuts <- compareMutations(rbind(mutmat[i,],mutmat[j,],insmat[i,],insmat[j,]))
#       intersection <- comparemuts[1,] + comparemuts[2,] - comparemuts[3,]
#       union <- comparemuts[3,]
#       if (sum(intersection)==0) {
#         jaccard[j] <- 0
#       } else {
#         jaccard[j] <- sum(union)/sum(intersection)
#       }
#     }
#   }
#   return(jaccard)
# },reads, mutmat, insmat, mc.cores=cores,mc.preschedule=F),identity)
# t2 <- proc.time() - t1

if (nrow(reads) <= 1) quit()

colnames(jaccard) <- reads$Read
rownames(jaccard) <- reads$Read

reads$Dup <- unlist(mclapply(1:nrow(reads),function(i,reads,jaccard) {
  j_max <- max(jaccard[i,])
  if (j_max > j_thresh) {
    return(reads$Read[which.max(jaccard[i,])])
  } else {
    return("")
  } 
},reads,jaccard,mc.cores=cores))

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

sorted_reads <- rbind(sorted_reads,cleanreads)
sorted_reads$Index <- NULL

write.table(sorted_reads,readfile,quote=F,sep="\t",na="",row.names=F,col.names=F)


pdf(sub(".txt","_similarity.pdf",readfile))

  jaccard.m <- melt(sorted_jaccard)
  jaccard.m$value <- ifelse(jaccard.m$Var2 >= jaccard.m$Var1,NA,jaccard.m$value)

  hm <- ggplot(jaccard.m, aes(y=Var1,x=Var2)) + scale_y_reverse() + geom_tile(aes(fill=value)) + scale_fill_gradient(low="blue",high="red",limits=c(0,1)) + guides(fill=F)
  if (nrow(clusters) > 0) hm <- hm + geom_polygon(aes(x=x,y=y,group=id),data=clusters,color="black",alpha=0,size=0.25)
  print(hm)
  
  breaks <- seq(0,1,length.out=20)
  rects <- data.frame(starts=breaks[1:(length(breaks)-1)])
  rects$ends <- rects$starts + diff(breaks)
  
  j_maxes <- data.frame(jmax=apply(jaccard,1,max))

 
  p_hist <- ggplot() + geom_rect(aes(xmin=starts,xmax=ends,ymin=-Inf,ymax=Inf,fill=starts),rects) + scale_fill_gradient(low="blue",high="red") + guides(fill=F) + geom_histogram(aes(x=jmax),j_maxes,breaks=breaks) 
  print(p_hist,vp=viewport(1,1,.45,.45,just=c("right","top")))
  
  
dev.off()


