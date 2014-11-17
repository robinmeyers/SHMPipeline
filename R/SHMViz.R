if (commandArgs()[1] != "RStudio") {
  
  ARGS <- c(
    "mutfile","character","file path of mutation file",
    "clonefile","character","file path of summary file",
    "refseqfile","character","reference sequence file",
    "output","character","file to write viz to"
    
  )
  
  OPTS <- c(
    "tstart","numeric",0,"Start of reference to include",
    "tend","numeric",0,"End of reference to include",
    "plotrows","numeric",4,"Rows on plot",
    "blankclones","logical",F,"",
    "showsubs","logical",T,"",
    "showdels","logical",T,"",
    "showins","logical",F,"",
    "figureheight","numeric",8,"height in inches",
    "showsequence","logical",T,"display sequence on plots",
    "regex1","character","AGCT","",
    "regex2","character","[AG]G[CT][AT]", ""
  )
  
  source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
  }
  
  source_local("Rsub.R")
  source_local("SHMHelper.R")
  
  parseArgs("SHMViz.R", ARGS, OPTS)
  
} else {
  mutfile <- "/Volumes//AltLab/SHM//Alt071-20140417/Results-New//JKH057_Alt071/JKH057_Alt071_muts.txt"
  clonefile <- "/Volumes//AltLab/SHM//Alt071-20140417/Results-New/JKH057_Alt071/JKH057_Alt071_clones.txt"
  refseqfile <- "/Volumes//AltLab/SHM//Alt071-20140417/Reference//VB18_productive_reference.fas"
  output <- "/Volumes//AltLab/SHM//Alt071-20140417/Results-New/viz//JKH057_Alt071_viz.pdf"
  tstart <- 141
  tend <- 500
  plotrows <- 3
  blankclones <- F
  showsubs <- T
  showdels <- T
  showins <- F
  figureheight <- 8
  showsequence <- T
  regex1 <- "AGCT"
  regex2 <- "[AG]G[CT][AT]"
  
  
  source("~/SHMPipeline/R/Rsub.R")
  source("~/SHMPipeline/R/SHMHelper.R")
  
}


suppressPackageStartupMessages(library(plyr, quietly=TRUE))
suppressPackageStartupMessages(library(RColorBrewer, quietly=TRUE))
suppressPackageStartupMessages(library(Biostrings, quietly=TRUE))
bases <- getBases()
basecolors <- getBasecolors()
ascii <- getAscii()

refseq <- readDNAStringSet(refseqfile)

ref <- data.frame(Pos=1:nchar(as.character(refseq)),Base=strsplit(as.character(refseq),""))
colnames(ref) <- c("Pos","Base")

ref$color <- match(ref$Base,bases)
ref$pch <- match(ref$Base,bases)

clones <- read.delim(clonefile,header=T,as.is=T)
muts <- read.delim(mutfile,header=F,colClasses=mutsColClasses(),as.is=T,col.names=c("Expt","Read","Pos","Type","From","To","Size","End","Ins"))

muts <- getMutsFromReads(clones,muts)

if (!all(grepl("([[:digit:]]+-[[:digit:]]+)(,[[:digit:]]+-[[:digit:]]+)*",clones$Coords))) {
  stop ("Clone coordinates not in correct format")
}

if (anyDuplicated(clones[,c("Expt","Read")])) {
  stop ("Duplicate read IDs in clonefile")
}

if (tstart == 0) {
  tstart <- 1
}
if (tend == 0) {
  tend <- nchar(refseq)
}

muttypes <- c()
if (showsubs) {
  muttypes <- c(muttypes,"sub")
}
if (showdels) {
  muttypes <- c(muttypes,"del")
}
if (showins) {
  muttypes <- c(muttypes,"ins")
}

muts <- muts[muts$Type %in% muttypes,]


if (! blankclones) {
  clones <- getReadsFromMuts(clones,muts)
}

blocks <- invertCoords(clones$Coords,refseq)


pdf(output,height=figureheight,width=11)

par(mai=c(0.2,0.75,0.2,0.75),omi=c(0.5,0,0,0))
  
layout(as.matrix(1:plotrows,ncol=1,nrow=plotrows))

ymax <- nrow(clones) + 0.5

rowwidth <- ceiling((tend-tstart+1)/plotrows)
  
tstarts <- tstart+rowwidth*0:(plotrows-1)
tends <- tstarts+rowwidth-1


refseq_rc <- reverseComplement(refseq)

regex1_plot <- data.frame(start=numeric(),end=numeric())
regex2_plot <- data.frame(start=numeric(),end=numeric())

if (regex1 != "") {
  
  regex1_match <- gregexpr(regex1,as.character(refseq))[[1]]
  regex1_match_rc <- gregexpr(regex1,as.character(refseq_rc))[[1]]
 
  if (regex1_match[1] > 0) {
    
    for (i in 1:length(regex1_match)) {

      start <- regex1_match[i]
      end <- start + attr(regex1_match,"match.length")[i] - 1
      regex1_plot[nrow(regex1_plot)+1,] <- c(start,end)

    }
  }
  
  if (regex1_match_rc[1] > 0) {
    for (i in 1:length(regex1_match_rc)) {

      end <- nchar(as.character(refseq)) - regex1_match_rc[i] + 1
      start <- end - attr(regex1_match,"match.length")[i] + 1
      regex1_plot[nrow(regex1_plot)+1,] <- c(start,end)

    }
  }
}

if (regex2 != "") {
  
  regex2_match <- gregexpr(regex2,as.character(refseq))[[1]]
  regex2_match_rc <- gregexpr(regex2,as.character(refseq_rc))[[1]]
 
  if (regex2_match[1] > 0) {
    
    for (i in 1:length(regex2_match)) {

      start <- regex2_match[i]
      end <- start + attr(regex2_match,"match.length")[i] - 1
      regex2_plot[nrow(regex2_plot)+1,] <- c(start,end)

    }
  }
  
  if (regex2_match_rc[1] > 0) {
    for (i in 1:length(regex2_match_rc)) {

      end <- nchar(as.character(refseq)) - regex2_match_rc[i] + 1
      start <- end - attr(regex2_match,"match.length")[i] + 1
      regex2_plot[nrow(regex2_plot)+1,] <- c(start,end)

    }
  }
}




for (i in 1:plotrows) {
  
  plot(c(),c(),ylab="",xlab="",xaxt="n",xlim=c(max(1,tstarts[i]-2),min(nchar(refseq),tends[i]+2)),ylim=c(0,ymax),xaxs="i",bty="o")
  axis(1,lwd=0,lwd.ticks=1)
  
  rect(xleft=regex2_plot$start-0.5,ybottom=-1,xright=regex2_plot$end+0.5,ytop=ymax,col=rgb(254,217,142,max=255),border=F)
  rect(xleft=regex1_plot$start-0.5,ybottom=-1,xright=regex1_plot$end+0.5,ytop=ymax,col=rgb(254,153,41,max=255),border=F)
    
  grid(ny=0,col=grey(0.5),lty=3)
  points(1:nrow(ref),rep(0,nrow(ref)),col=basecolors[ref$color],pch=ascii[ref$pch],cex=0.6)
  if (nrow(blocks) > 0) {
    rect(xleft=blocks$start-0.5,ybottom=blocks$i-0.5,xright=blocks$end+0.5,ytop=blocks$i+0.5,col=grey(0.1,0.25),border=F)
  }

  muts$y <- match(paste(muts$Expt,muts$Read),paste(clones$Expt,clones$Read))
  subs <- muts[muts$Type == "sub",]
  if (nrow(subs) > 0) {
    subs$color <- match(subs$To,bases)
    subs$pch <- match(subs$To,bases)
    points(subs$Pos,subs$y,col=basecolors[subs$color],pch=ascii[subs$pch],cex=0.5)
  }

  dels <- muts[muts$Type == "del",]
  if (nrow(dels) > 0) {
    segments(x0=dels$Pos-0.5,y0=dels$y,x1=dels$End+0.5,col=basecolors[7],lwd=1)
  }

  ins <- muts[muts$Type == "ins",]
  if (nrow(ins) > 0) {
    text(x=ins$Pos,y=ins$y,label=ins$Ins,cex=0.5)
  }

}
dev.off()
