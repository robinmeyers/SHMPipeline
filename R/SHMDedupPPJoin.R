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
  mutfile <- "/Volumes/AltLab/SHM/Alt053/results-new/JKH001_Alt053/JKH001_Alt053_muts.txt"
  readfile <- "/Volumes/AltLab/SHM/Alt053/results-new/JKH001_Alt053/JKH001_Alt053_reads.txt"
  refseqfile <- "/Volumes/AltLab/SHM/Alt053/ref//VRCPG04_UCA_VDJ_reference_sequence.fas"
  tstart <- 0
  tend <- 0
  j_thresh <- .8
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
muts$Type <- factor(muts$Type,levels=c("sub","del","ins"))

reads <- read.delim(readfile,header=F,as.is=T,col.names=c("Expt","Read","Bp","Coords","Dup"))

reads <- reads[reads$Bp > 0,]

rownames(reads) <- reads$Read
  
# First create tokens

tokens <- sapply(reads$Read,function(x) NULL)

for (m in 1:nrow(muts)) {
  if (muts[m,"Type"] == "sub") {
    if (muts[m,"Pos"] >= tstart && muts[m,"Pos"] <= tend) {
      tokens[[muts[m,"Read"]]] <- c(tokens[[muts[m,"Read"]]],paste(muts[m,"Pos"],muts[m,"To"],sep=""))
    }
  } else if (muts[m,"Type"] == "del") {
    if (muts[m,"Pos"] >= tstart && muts[m,"End"] <= tend) {
      tokens[[muts[m,"Read"]]] <- c(tokens[[muts[m,"Read"]]],paste(muts[m,"Pos"],"<",sep=""),paste(muts[m,"End"],">",sep=""))
    }
  } else {
    if (muts[m,"Pos"] >= tstart && muts[m,"Pos"] <= tend) {
      tokens[[muts[m,"Read"]]] <- c(tokens[[muts[m,"Read"]]],paste(muts[m,"Pos"],"i",muts[m,"Ins"],sep=""))
    }
  }
}

reads$Tokens <- unlist(lapply(tokens,length))

tokens.count <- table(unlist(tokens))
tokens.count <- sort(tokens.count)

reads.clean <- reads[reads$Tokens == 0,]
reads <- reads[reads$Tokens > 0,]

reads <- reads[order(reads$Tokens),]

tokens <- tokens[reads$Read]
tokens <- lapply(tokens, function(t) { sort(factor(t,levels=names(tokens.count),ordered=T)) } )

dup.pairs <- sapply(reads$Read,function(x) NULL)
inv.index <- sapply(names(tokens.count),function(x) NULL)

for (x in 1:length(tokens)) {
  overlap.map <- rep(0,length(tokens))
  names(overlap.map) <- names(tokens)
  x.rec.id <- names(tokens[x])
  x.rec <- tokens[[x]]
  prefix.length <- length(x.rec) - ceiling(j_thresh*length(x.rec)) + 1;
  
  for (i in 1:prefix.length) {
    token <- x.rec[i]
    y.recs <- inv.index[[token]]
    if (length(y.recs) > 0) {
      for (y in 1:length(y.recs)) {
        y.rec.id <- names(y.recs)[y]
        j <- y.recs[y]
        y.rec <- tokens[[y.rec.id]]
        if (length(y.rec) < j_thresh*length(x.rec)) next  
        alpha <- ceiling(j_thresh/(1+j_thresh)*(length(x.rec)+length(y.rec)))
        ubound <- 1 + min(length(x.rec)-i,length(y.rec)-j)
        if (overlap.map[y.rec.id] + ubound >= alpha) {
          overlap.map[y.rec.id] <- overlap.map[y.rec.id] + 1
        } else {
          overlap.map[y.rec.id] <- 0
        }
      }
    }
    inv.index[[token]][[x.rec.id]] <- i
  }
  
  for (y in 1:length(overlap.map)) {
    if (overlap.map[y] > 0) {
      p.x <- prefix.length
      y.rec.id <- names(overlap.map[y])
      y.rec <- tokens[[y.rec.id]]
      p.y <- length(y.rec) - ceiling(j_thresh*length(y.rec)) + 1;
      w.x <- x.rec[p.x]
      w.y <- y.rec[p.y]
      overlap <- overlap.map[y]
      alpha <- ceiling(j_thresh/(1+j_thresh)*(length(x.rec)+length(y.rec)))
      if (w.x < w.y) {
        ubound <- overlap.map[y] + length(x.rec) - p.x
        if (ubound >= alpha) {
          overlap <- overlap + length(intersect(tail(x.rec,n=length(x.rec)-p.x),tail(y.rec,n=length(y.rec)-overlap.map[y])))
        }
      } else {
        ubound <- overlap.map[y] + length(y.rec) - p.y
        if (ubound >= alpha) {
          overlap <- overlap + length(intersect(tail(x.rec,n=length(x.rec)-overlap.map[y]),tail(y.rec,n=length(y.rec)-p.y)))
        }
      }
      if (overlap >= alpha) {
        dup.pairs[[y.rec.id]] <- c(x.rec.id,dup.pairs[[y.rec.id]])
      }
    }
  }
    
}


reads <- rbind(reads,reads.clean)
reads$Dup <- ""

reads[names(dup.pairs),]$Dup <- unlist(lapply(dup.pairs,function(matches) {if (length(matches) > 0) return(matches[1]) else return("")}))

reads$Tokens <- NULL

write.table(reads,readfile,quote=F,sep="\t",na="",row.names=F,col.names=F)

