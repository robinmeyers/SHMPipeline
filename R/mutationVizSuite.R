ARGS <- c(
  "mutfile","character","file path of mutation file",
  "clonefile","character","file path of summary file",
  "refseqfile","character","reference sequence file",
  "outstub","character","file to plot to"
  
)

OPTS <- c(
  "tstart","numeric",0,"Start of reference to view",
  "tend","numeric",0,"End of reference to view",
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

parseArgs("mutationViz.R", ARGS, OPTS)

library(Biostrings)
library(RColorBrewer)
basecolors <- getBasecolors()
bases <- getBases()
ascii <- getAscii()

refseq <- readDNAStringSet(refseqfile)


clones <- read.delim(clonefile,header=T,as.is=T)
clones <- clones[order(-clones$Bp),]

muts <- read.delim(mutfile,header=F,as.is=T,col.names=c("Expt","Clone","Pos","Type","From","To","Size","End","Ins"))

subs <- muts[ muts$Type == "sub",]

dels <- muts [ muts$Type == "del",]

cloneIDs <- clones$ID[clones$Bp > 0 & clones$Subs > minsubs]

if (tstart == 0) {
  tstart <- 1
}
if (tend == 0) {
  tend <- nchar(refseq)
}



ref <- data.frame(Pos=1:nchar(as.character(refseq)),Base=strsplit(as.character(refseq),""))
colnames(ref) <- c("Pos","Base")

ref$color <- match(ref$Base,bases)
ref$pch <- match(ref$Base,bases)



blocks <- data.frame(Clone=integer(),Start=integer(),End=integer())

if (length(cloneIDs) > 0) {
  for (i in 1:length(cloneIDs)) {
    coordlist <- unlist(strsplit(clones[clones$ID == cloneIDs[i],"Coords"],","))
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

pdf(paste(outstub,"_mutViz.pdf",sep=""),height=8,width=11)
tictactoePlot(subs,dels,blocks,ref,tstart,tend,plotrows,cloneIDs)
dev.off()

pdf(paste(outstub,"_subDens.pdf",sep=""),height=8,width=11)
ref <- connectfourSubPlot(subs, blocks, ref, tstart, tend, plotrows,cloneIDs)
write.table(ref,paste(outstub,"_bases.txt",sep=""),quote=F,sep="\t",na="",row.names=F,col.names=T)
dev.off()

pdf(paste(outstub,"_delDens.pdf",sep=""),height=8,width=11)
connectfourDelPlot(dels, blocks, ref, tstart, tend, plotrows,cloneIDs)
dev.off()

