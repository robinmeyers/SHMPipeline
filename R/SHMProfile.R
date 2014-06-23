
if (commandArgs()[1] != "RStudio") {
  
  ARGS <- c(
    "mutfile","character","file path of mutation file",
    "readfile","character","file path of summary file",
    "refseqfile","character","reference sequence file",
    "outstub","character","file to write profile to"
    
  )
  
  OPTS <- c(
    "tstart","numeric",0,"Start of reference to include",
    "tend","numeric",0,"End of reference to include",
    "filtstart","numeric",0,"Start of reference to consider in filtering clones",
    "filtend","numeric",0,"End of reference to consider in filtering clones",
    "minsubs","numeric",0,"minimum number of substitutions for a clone to be included",
    "maxsubs","numeric",0,"maximum number of substitutions for a clone to be included",
    "mindels","numeric",0,"minimum number of deletions for a clone to be included",
    "maxdels","numeric",0,"maximum number of deletions for a clone to be included",
    "minins","numeric",0,"minimum number of insertions for a clone to be included",
    "maxins","numeric",0,"maximum number of insertions for a clone to be included",
    "rmdups","logical",TRUE,"remove clones marked as duplicates"
  )
  
  source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
  }
  
  source_local("Rsub.R")
  source_local("SHMHelper.R")
  
  parseArgs("SHMProfile.R", ARGS, OPTS)
  
} else {
  mutfile <- "/Volumes//AltLab/SHM//Alt071-20140417/Results-Test/JKH057_Alt071/JKH057_Alt071_muts.txt"
  readfile <- "/Volumes//AltLab/SHM//Alt071-20140417/Results-Test/JKH057_Alt071/JKH057_Alt071_reads.txt"
  refseqfile <- "/Volumes//AltLab/SHM//Alt071-20140417/Reference//VB18_productive_reference.fas"
  outstub <- "/Volumes//AltLab/SHM//Alt071-20140417/Results-Test/JKH057_Alt071/JKH057_Alt071"
  tstart <- 141
  tend <- 500
  filtstart <- 0
  filtend <- 0
  minsubs <- 0
  maxsubs <- 0
  mindels <- 0
  maxdels <- 0
  rmdups <- T
  
  source("~/SHMPipeline/R/Rsub.R")
  source("~/SHMPipeline/R/SHMHelper.R")
  
}

suppressPackageStartupMessages(library(plyr, quietly=TRUE))
suppressPackageStartupMessages(library(Biostrings, quietly=TRUE))

refseq <- readDNAStringSet(refseqfile)

reads <- read.delim(readfile,header=F,as.is=T,col.names=c("Expt","Read","Bp","Coords","Dup"),)
reads$Dup <- ifelse(is.na(reads$Dup),"",reads$Dup)


muts <- read.delim(mutfile,header=F,as.is=T,col.names=c("Expt","Read","Pos","Type","From","To","Size","End","Ins"))
# muts$From <- factor(muts$From,levels=c("A","C","G","T"))
# muts$To <- factor(muts$To,levels=c("A","C","G","T"))

if (!all(grepl("([[:digit:]]+-[[:digit:]]+)(,[[:digit:]]+-[[:digit:]]+)*",reads$Coords))) {
  stop ("Read coordinates not in correct format")
}


#### Exit if mutations exist that don't have a read in read data.frame
# if (any(is.na(merge(reads,muts,by=1:2,all.y=T)$Bp))) {
#   stop ("Mutations from clones not included in clonefile")
# }

if (anyDuplicated(reads[,c("Expt","Read")])) {
  stop ("Duplicate read IDs in clonefile")
}

if (tstart == 0) {
  tstart <- 1
}
if (tend == 0) {
  tend <- nchar(refseq)
}

if (filtstart == 0) {
  filtstart <- tstart
}
if (filtend == 0) {
  filtend <- tend
}


reads$filtsubs <- sapply(1:nrow(reads),function(i,rs,ms) {
    readMuts <- getMutsFromRead(rs[i,],muts)
    return(sum(readMuts$Type == "sub" & readMuts$Pos >= filtstart & readMuts$Pos <= filtend))
  },reads,muts)

reads$filtdels <- sapply(1:nrow(reads),function(i,rs,ms) {
  readMuts <- getMutsFromRead(rs[i,],muts)
  return(sum(readMuts$Type == "del" & readMuts$Pos >= filtstart & (readMuts$Pos + readMuts$Size - 1) <= filtend))
},reads,muts)

reads$filtins <- sapply(1:nrow(reads),function(i,rs,ms) {
  readMuts <- getMutsFromRead(rs[i,],muts)
  return(sum(readMuts$Type == "ins" & readMuts$Pos >= filtstart & readMuts$Pos <= filtend))
},reads,muts)

reads$filter <- 0

if (rmdups) {
  reads$filter <- ifelse(reads$Dup != "",1,reads$filter)
}

if (minsubs > 0) {
  reads$filter <- ifelse(reads$filtsubs < minsubs,1,reads$filter)
}

if (maxsubs > 0) {
  reads$filter <- ifelse(reads$filtsubs > maxsubs,1,reads$filter)
}

if (mindels > 0) {
  reads$filter <- ifelse(reads$filtdels < mindels,1,reads$filter)
}

if (maxdels > 0) {
  reads$filter <- ifelse(reads$filtdels > maxdels,1,reads$filter)
}

if (minins> 0) {
  reads$filter <- ifelse(reads$filtins < minins,1,reads$filter)
}

if (maxins > 0) {
  reads$filter <- ifelse(reads$filtins > maxins,1,reads$filter)
}

clones <- reads[reads$filter == 0,]
muts <- getMutsFromReads(clones,muts)
# 
# subs <- muts[muts$Type == "sub" & muts$Pos >= tstart & muts$Pos <= tend,]
# dels <- muts[muts$Type == "del" & muts$Pos >= tstart & muts$Pos + muts$Size - 1 <= tend,]
# ins <- muts[muts$Type == "ins" & muts$Pos >= tstart & muts$Pos <= tend,]


mut.mat <- createMutationMatrix(clones,muts,refseq,tstart,tend)
ins.mat <- createInsertionMatrix(clones,muts,refseq,tstart,tend)

profile <- calculateProfile(mut.mat,refseq)
delprofile <- calculateDeletionProfile(mut.mat,refseq)

clones$Subs <- apply(mut.mat,1,function(x) {sum(grepl("[ACGT]",x))})

clones[,c("Dels","DelBp","LargeDel")] <- t(apply(mut.mat,1,function(x) {
  delstarts <- which(grepl("<",x))
  delends <- which(grepl(">",x))
  delsizes <- delends - delstarts + 1
  if (length(delsizes) > 0) {
    return(c(length(delsizes),sum(delsizes),max(delsizes)))
  } else {
    return(c(0,0,0))
  }}))

clones[,c("Ins","InsBp","LargeIns")] <- t(apply(ins.mat,1,function(x) {
  inssizes <- nchar(x[grepl("[ACGTN]",x)])
  if (length(inssizes) > 0) {
    return(c(length(inssizes),sum(inssizes),max(inssizes)))
  } else {
    return(c(0,0,0))
  }}))



bases <- c("A","C","G","T")
bcols <- c()
for (i in bases) {
  for (j in bases) {
    if (i != j) bcols <- c(bcols,paste(i,">",j,sep=""))
  }
}

basex <- t(apply(mut.mat,1,function(x,ref) {
  basetable <- matrix(0,ncol=4,nrow=4,dimnames=list(bases,bases))
  pos <- which(grepl("[ACGT]",x))
  for (i in pos) {
    basetable[ref[i],x[i]] <- basetable[ref[i],x[i]] + 1
  }
  return(as.vector(t(basetable))[c(-1,-6,-11,-16)])
},unlist(strsplit(as.character(refseq),""))))



colnames(basex) <- bcols
basex <- cbind(clones[,1:2],basex)

write.table(profile,paste(outstub,"_profile.txt",sep=""),quote=F,sep="\t",na="",row.names=F,col.names=T)
write.table(clones,paste(outstub,"_clones.txt",sep=""),quote=F,sep="\t",na="",row.names=F,col.names=T)
write.table(basex,paste(outstub,"_basex.txt",sep=""),quote=F,sep="\t",na="",row.names=F,col.names=T)
write.table(mut.mat,paste(outstub,"_mutmatrix.txt",sep=""),quote=F,sep="\t",row.names=F,col.names=T)

