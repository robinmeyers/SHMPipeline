#!/usr/bin/Rscript

if (commandArgs()[1] != "RStudio") {
  
  ARGS <- c(
    "metafile","character","file path of meta file",
    "results","character","results directory"
  )
  
  OPTS <- c(
    "grouping","character","genotype,allele,tissue,pna","meta file variables to group by"
  )
  
  source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
  }
  
  source_local("Rsub.R")
  source_local("SHMHelper.R")
  
  parseArgs("SHMAggregate.R", ARGS, OPTS)
  
  
} else {
  metafile <- "/Volumes//AltLab/SHM//Alt102/meta/JKH109_meta.txt"
  results <- "/Volumes//AltLab/SHM//Alt102/testresults/"
  grouping <- "genotype,allele,tissue,pna"
  
  
  source("~/SHMPipeline/R/Rsub.R")
  source("~/SHMPipeline/R/SHMHelper.R")
  
}


meta <- read.delim(metafile,header=T,as.is=T)

grouping <- unlist(strsplit(grouping,","))
groups <- data.frame(unique(meta[,grouping]))
rownames(groups) <- 1:nrow(groups)

totalreads <- data.frame()
totalmuts <- data.frame()
totalclones <- data.frame()
totalbasex <- data.frame()

groupdir <-file.path(results,"groups")
dir.create(groupdir)

for (i in 1:nrow(groups)) {
  expts <- merge(meta,groups[i,])
  groupreads <- data.frame()
  groupclones <- data.frame()
  groupmuts <- data.frame()
  groupbasex <- data.frame()
  
  groupreadfile <- file.path(groupdir,paste(paste(groups[i,],collapse="-"),"_reads.txt",sep=""))
  groupclonefile <- file.path(groupdir,paste(paste(groups[i,],collapse="-"),"_clones.txt",sep=""))
  groupmutfile <- file.path(groupdir,paste(paste(groups[i,],collapse="-"),"_muts.txt",sep=""))
  groupbasexfile <- file.path(groupdir,paste(paste(groups[i,],collapse="-"),"_basex.txt",sep=""))
  
  for (j in 1:nrow(expts)) {

    readfile <- file.path(results,expts$experiment[j],paste(expts$experiment[j],"_reads.txt",sep=""))    
    clonefile <- file.path(results,expts$experiment[j],paste(expts$experiment[j],"_clones.txt",sep=""))
    mutfile <- file.path(results,expts$experiment[j],paste(expts$experiment[j],"_muts.txt",sep=""))
    basexfile <- file.path(results,expts$experiment[j],paste(expts$experiment[j],"_basex.txt",sep=""))
    
    if (!(file.exists(readfile) && file.exists(clonefile) && file.exists(mutfile) && file.exists(basexfile))) {
      next
    }
    
    exptreads <- read.delim(readfile,header=F,as.is=T,col.names=c("Expt","Read","Bp","Coords","Dup"))
    exptclones <- read.delim(clonefile,header=T,as.is=T)
    exptmuts <- read.delim(mutfile,header=F,as.is=T,col.names=c("Expt","Read","Pos","Type","From","To","Size","End","Ins"))
    exptbasex <- read.delim(basexfile,header=T,as.is=T)

    if (nrow(groupreads) < 1) {
      groupreads <- exptreads
    } else {
      groupreads <-rbind(groupreads,exptreads)
    }

    if (nrow(groupclones) < 1) {
      groupclones <- exptclones
    } else {
      groupclones <-rbind(groupclones,exptclones)
    }

    if (nrow(groupmuts) < 1) {
      groupmuts <- exptmuts
    } else {
      groupmuts <- rbind(groupmuts,exptmuts)
    }

    if (nrow(groupbasex) < 1) {
      groupbasex <- exptbasex
    } else {
      groupbasex <-rbind(groupbasex,exptbasex)
    }
  }
  
  if (nrow(totalreads) < 1) {
    totalreads <- groupreads
  } else {
    totalreads <-rbind(totalreads,groupreads)
  }

  if (nrow(totalclones) < 1) {
    totalclones <- groupclones
  } else {
    totalclones <-rbind(totalclones,groupclones)
  }

  if (nrow(totalbasex) < 1) {
    totalbasex <- groupbasex
  } else {
    totalbasex <-rbind(totalbasex,groupbasex)
  }
  
  write.table(groupreads,groupreadfile,sep="\t",quote=F,row.names=F,col.names=F,na="")
  write.table(groupclones,groupclonefile,sep="\t",quote=F,row.names=F,col.names=T,na="")
  write.table(groupmuts,groupmutfile,sep="\t",quote=F,row.names=F,col.names=F,na="")
  write.table(groupbasex,groupbasexfile,sep="\t",quote=F,row.names=F,col.names=T,na="")

  
}

meta$Reads <- 0
readnum <- table(totalreads$Expt)
meta$Reads[match(names(readnum),meta$experiment)] <- readnum

groups <- aggregate(as.formula(paste("Reads ~",paste(grouping,collapse=" + "))),meta,sum)

totalclones$Clones <- 1

exptstats <- aggregate(cbind(Clones,Bp,Subs,Dels,DelBp,Ins,InsBp) ~ Expt,totalclones,sum)
exptstats <- merge(meta,exptstats,by=1)

write.table(exptstats,file.path(results,"Expts.txt"),sep="\t",quote=F,row.names=F,col.names=T,na="")

groupstats <- aggregate(as.formula(paste("cbind(Reads,Clones,Bp,Subs,Dels,DelBp,Ins,InsBp) ~",paste(grouping,collapse=" + "))),exptstats,sum)

write.table(groupstats,file.path(results,"Groups.txt"),sep="\t",quote=F,row.names=F,col.names=T,na="")

exptbasex <- aggregate( as.formula( paste("cbind(",paste(colnames(totalbasex)[(ncol(totalbasex)-11):ncol(totalbasex)],collapse=","),") ~ Expt",sep="")), totalbasex, sum)

write.table(exptbasex,file.path(results,"Expt_BaseX.txt"),sep="\t",quote=F,row.names=F,col.names=T,na="")

exptbasex <- merge(meta,exptbasex,by=1)

groupbasex <- aggregate(as.formula( paste("cbind(",paste(colnames(exptbasex)[(ncol(exptbasex)-11):ncol(exptbasex)],collapse=","),") ~ ",paste(grouping,collapse=" + "),sep="")),exptbasex,sum)

write.table(groupbasex,file.path(results,"Group_BaseX.txt"),sep="\t",quote=F,row.names=F,col.names=T,na="")




