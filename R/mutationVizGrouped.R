#!/usr/bin/Rscript

ARGS <- c(
  "metafile","character","file path of meta file",
  "mutfile","character","file path of mutation file",
  "clonefile","character","file pat of clonefile",
  "refdir","character","reference sequence dir",
  "outdir","character","directory to send plots to"
)

OPTS <- c(
  "grouping","character","genotype,allele,tissue,pna","meta file variables to group by",
  "plotrows","numeric",4,"Rows on plot",
  "minsubs","numeric",0,"minimum number of substitutions for a clone to be included"
)

source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local("Rsub.R")
source_local("SHMHelper.R")

parseArgs("mutationVizGrouped.R", ARGS, OPTS)


library(Biostrings)
library(RColorBrewer)
basecolors <- getBasecolors()
bases <- getBases()
ascii <- getAscii()



meta <- read.delim(metafile,header=T,as.is=T)

muts <- read.delim(mutfile,header=F,as.is=T,col.names=c("Expt","Clone","Pos","Type","From","To","Size","End","Ins"))

clones <- read.delim(clonefile,header=T,as.is=T)
clones <- clones[order(clones$Bp),]


grouping <- unlist(strsplit(grouping,","))
groups <- data.frame(unique(meta[,grouping]))
rownames(groups) <- 1:nrow(groups)

meta$group <- match(lapply(1:nrow(meta),function(x){paste(meta[x,grouping],collapse=" ")}),lapply(1:nrow(groups),function(x){paste(groups[x,],collapse=" ")}))


clones$group <- meta$group[match(clones$Expt,meta$experiment)]


muts$group <- meta$group[match(muts$Expt,meta$experiment)]


for (group in 1:nrow(groups)) {
  
  refseqfile <- meta$reference[match(group,meta$group)[1]]
  refseq <- readDNAStringSet(paste(refdir,refseqfile,sep="/"))
  
  tstart <- meta$start[match(group,meta$group)[1]]
  tend <- meta$end[match(group,meta$group)[1]]
  # tstart <- 1
  # tend <- nchar(refseq)
  
  ref <- data.frame(Pos=1:nchar(as.character(refseq)),Base=strsplit(as.character(refseq),""))
  colnames(ref) <- c("Pos","Base")
  ref$color <- match(ref$Base,bases)
  ref$pch <- match(ref$Base,bases)
  
  group_muts <- muts[muts$group == group,]
  group_clones <- clones[clones$group == group,]
  group_subs <- group_muts[group_muts$Type == "sub",]
  group_dels <- group_muts[group_muts$Type == "del",]
  
  
  blocks <- data.frame(Clone=integer(),Start=integer(),End=integer())
  
  cloneIDs <- group_clones$Clone[group_clones$Bp > 0 & group_clones$Subs >= minsubs]
  
  if (length(cloneIDs) > 0) {
    for (i in 1:length(cloneIDs)) {
      coordlist <- unlist(strsplit(group_clones[group_clones$Clone == cloneIDs[i],"Coords"],","))
      for (j in 1:length(coordlist)) {
        coords <- as.integer(unlist(strsplit(coordlist[j],"-")))
        if (j == 1 && coords[1] > 1) {
          blocks[nrow(blocks)+1,] <- c(i,1,coords[1]-1)
        }
        if (j > 1) {
          coordslast <- as.integer(unlist(strsplit(coordlist[j-1],"-")))
          blocks[nrow(blocks)+1,] <- c(i,coordslast[2]+1,coords[1]-1)
        }
        if (j == length(coordlist) && coords[2] < nrow(ref)) {
          blocks[nrow(blocks)+1,] <- c(i,coords[2]+1,nrow(ref))
        }
      }
    }
  }  
  

  
  
  pdf(paste(outdir,"/",paste(groups[group,],collapse="_"),"_mutViz.pdf",sep=""),height=8,width=11)
  tictactoePlot(group_subs,group_dels,blocks,ref,tstart,tend,plotrows,cloneIDs)
  dev.off()
  
  pdf(paste(outdir,"/",paste(groups[group,],collapse="_"),"_subDens.pdf",sep=""),height=8,width=11)
  ref <- connectfourSubPlot(group_subs, blocks, ref, tstart, tend, plotrows,cloneIDs)
  write.table(ref,paste(outdir,"/",paste(groups[group,],collapse="_"),"_bases.txt",sep=""),quote=F,sep="\t",na="",row.names=F,col.names=T)
  dev.off()
  
  pdf(paste(outdir,"/",paste(groups[group,],collapse="_"),"_delDens.pdf",sep=""),height=8,width=11)
  connectfourDelPlot(group_dels, blocks, ref, tstart, tend, plotrows,cloneIDs)
  dev.off()
  
#   if (nrow(group_dels) > 0) { 
#     pdf(paste(outdir,"/",paste(groups[group,],collapse="_"),"_delHist.pdf",sep=""),height=8,width=11)
#     hist(group_dels$Size,breaks=20,freq=T,xlab="Deletion Size",ylab="Frequency")
#     dev.off()
#   }
#   
#   pdf(paste(outdir,"/",paste(groups[group,],collapse="_"),"_cloneSubHist.pdf",sep=""),height=8,width=11)
#   hist(group_clones$Subs,breaks=50,freq=T,xlab="Number of Substitutions",ylab="Frequency")
#   dev.off()
}

