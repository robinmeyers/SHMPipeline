ARGS <- c(
  "mutfile", "character", "file path of mutation file",
  "outfile","character","file to plot to",
  "refseq","character","reference sequence"
)

OPTS <- c(
  "tstart","numeric",0,"Start of reference to view",
  "tend","numeric",0,"End of reference to view",
  "plotrows","numeric",3,"Rows on plot"
)

source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local("Rsub.R")

parseArgs("mutationViz.R", ARGS, OPTS)

library(ShortRead)
library(RColorBrewer)
basecolors <- brewer.pal(5,"Set1")

refseq <- readDNAStringSet(refseqfile)

muts <- read.delim(mutfile,header=F,as.is=T)

colnames(muts)[1:3] <- c("ID","Pos","Type")

muts <- muts [ muts$Type != "", ]

subs <- muts[ muts$Type == "sub", ]
colnames(subs)[4:5] <- c("From","To")
subs$Pos <- as.integer(subs$Pos)

dels <- muts [ muts$Type == "del", ]
colnames(dels)[4:5] <- c("Size","End")
dels$Pos <- as.integer(dels$Pos)
dels$Size <- as.integer(dels$Size)
dels$End <- as.integer(dels$End)

clones <- unique(muts$ID[muts$Type != "ins"])

if (tstart == 0) {
  tstart <- min(c(subs$Pos,dels$Pos))
}
if (tend == 0) {
  tend <- max(c(subs$Pos,dels$End))
}

subs$y <- match(subs$ID,clones)
bases <- c("A","C","G","T","N")
ascii <- c(65,67,71,84)
subs$color <- match(subs$To,bases)
subs$pch <- match(subs$To,bases)

dels$y <- match(dels$ID,clones)

ref <- data.frame(Pos=1:nchar(as.character(refseq)),Base=strsplit(as.character(refseq),""))
colnames(ref) <- c("Pos","Base")

ref$color <- match(ref$Base,bases)
ref$pch <- match(ref$Base,bases)

agct <- unlist(gregexpr("AGCT",as.character(refseq)))


rowwidth <- ceiling((tend-tstart+1)/plotrows)

tstarts <- tstart+rowwidth*0:(plotrows-1)
tends <- tstarts+rowwidth-1




pdf(outfile,height=8,width=11)
par(mai=c(0.2,0.75,0.2,0.75),omi=c(0.5,0,0,0))

layout(as.matrix(1:plotrows,ncol=1,nrow=plotrows))

for (i in 1:plotrows) {

  plot(c(),c(),ylab="Clone",xlab="",yaxt="n",xlim=c(tstarts[i]-2,tends[i]+2),ylim=c(0,length(clones)),xaxs="i")
  axis(2,at=0:length(clones),labels=c("ref",clones),las=1)
  axis(4,at=0:length(clones),labels=c("ref",clones),las=1)
  rect(xleft=agct-0.5,ybottom=-0.5,xright=agct+3.5,ytop=length(clones)+0.5,col="LemonChiffon",border=F)
  grid(ny=length(clones)+1)
  points(1:nrow(ref),rep(0,nrow(ref)),col=basecolors[ref$color],pch=ascii[ref$pch],cex=0.6)
  points(subs$Pos,subs$y,col=basecolors[subs$color],pch=ascii[subs$pch],cex=0.8)
  segments(dels$Pos-0.5,y0=dels$y,x1=dels$End+0.5)

}

dev.off()
