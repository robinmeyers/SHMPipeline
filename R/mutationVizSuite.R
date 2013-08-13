ARGS <- c(
  "mutfile","character","file path of mutation file",
  "clonefile","character","file path of summary file",
  "refseqfile","character","reference sequence file",
  "outstub","character","file to plot to"
  
)

OPTS <- c(
  "tstart","numeric",0,"Start of reference to view",
  "tend","numeric",0,"End of reference to view",
  "plotrows","numeric",4,"Rows on plot"
)

source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local("Rsub.R")

parseArgs("mutationViz.R", ARGS, OPTS)

library(Biostrings)
library(RColorBrewer)
basecolors <- brewer.pal(5,"Set1")

refseq <- readDNAStringSet(refseqfile)


clones <- read.delim(clonefile,header=T,as.is=T)
clones <- clones[order(-clones$Bp),]

muts <- read.delim(mutfile,header=F,as.is=T,col.names=c("Expt","ID","Pos","Type","From","To","Size","End","Ins"))

subs <- muts[ muts$Type == "sub", c("ID","Pos","Type","From","To")]

dels <- muts [ muts$Type == "del", c("ID","Pos","Type","Size","End")]

cloneIDs <- clones$ID[clones$ID %in% muts$ID[muts$Type != "ins"] & clones$Bp > 0]

if (tstart == 0) {
  tstart <- 1
}
if (tend == 0) {
  tend <- nchar(refseq)
}

subs$y <- match(subs$ID,cloneIDs)
bases <- c("A","C","G","T","N")
ascii <- c(65,67,71,84,78)
subs$color <- match(subs$To,bases)
subs$pch <- match(subs$To,bases)

dels$y <- match(dels$ID,cloneIDs)

ref <- data.frame(Pos=1:nchar(as.character(refseq)),Base=strsplit(as.character(refseq),""))
colnames(ref) <- c("Pos","Base")

ref$color <- match(ref$Base,bases)
ref$pch <- match(ref$Base,bases)

agct <- unlist(gregexpr("AGCT",as.character(refseq)))
rgyw <- unique(c(unlist(gregexpr("[AG]G[CT][AT]",as.character(refseq))),unlist(gregexpr("[AT][AG]C[CT]",as.character(refseq)))))


rowwidth <- ceiling((tend-tstart+1)/plotrows)

tstarts <- tstart+rowwidth*0:(plotrows-1)
tends <- tstarts+rowwidth-1

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
par(mai=c(0.2,0.75,0.2,0.75),omi=c(0.5,0,0,0))

layout(as.matrix(1:plotrows,ncol=1,nrow=plotrows))
ymax <- max(5.5,length(cloneIDs))

for (i in 1:plotrows) {

  plot(c(),c(),ylab="",xlab="",xaxt="n",yaxt="n",xlim=c(max(1,tstarts[i]-2),min(nchar(refseq),tends[i]+2)),ylim=c(0,ymax),xaxs="i",bty="n")
  axis(1,lwd=0,lwd.ticks=1)
  axis(2,at=0:length(cloneIDs),labels=c("ref",cloneIDs),las=1,lwd=0)#,lwd.ticks=1)
  axis(4,at=0:length(cloneIDs),labels=c("ref",cloneIDs),las=1,lwd=0)#,lwd.ticks=1)
  rect(xleft=rgyw-0.5,ybottom=-0.5,xright=rgyw+3.5,ytop=ymax+0.5,col="LemonChiffon",border=F)
  rect(xleft=agct-0.5,ybottom=-0.5,xright=agct+3.5,ytop=ymax+0.5,col="LightGoldenrod",border=F)
  grid(ny=0,col=grey(0.5),lty=3)
  abline(h=(0:ymax)+0.5,col=grey(0.5),lty=3)
  points(1:nrow(ref),rep(0,nrow(ref)),col=basecolors[ref$color],pch=ascii[ref$pch],cex=0.6)
  if (length(cloneIDs) > 0) {
    rect(xleft=blocks$Start-0.5,ybottom=blocks$Clone-0.5,xright=blocks$End+0.5,ytop=blocks$Clone+0.5,col=grey(0.1,0.25),border=F)
    points(subs$Pos,subs$y,col=basecolors[subs$color],pch=ascii[subs$pch],cex=0.8)
    segments(dels$Pos-0.5,y0=dels$y,x1=dels$End+0.5,col=basecolors[5],lwd=2)
  } else {
    text(x=(tstarts[i]+tends[i])/2,y=ymax/2-0.5,"No Mutations To Display")
  }
  

}

dev.off()


pdf(paste(outstub,"_mutDens.pdf",sep=""),height=8,width=11)
par(mai=c(0.2,0.75,0.2,0.75),omi=c(0.5,0,0,0))

layout(as.matrix(1:plotrows,ncol=1,nrow=plotrows))

ref$dens <- 0
if (nrow(subs) > 0) {
  ref$dens <- hist(subs$Pos,c(0,1:nrow(ref)),plot=F)$counts
}
dens <- data.frame(x=c(ref$Pos-0.25,ref$Pos+0.25),y=c(ref$dens,ref$dens))
dens <- dens[order(dens$x),]
ymax <- max(5,max(dens$y))

for (i in 1:plotrows) {
  plot(c(),c(),ylab="",yaxt="n",xaxt="n",xlab="",xlim=c(max(1,tstarts[i]-2),min(nchar(refseq),tends[i]+2)),ylim=c(-0.5,ymax),xaxs="i",bty="n")
  axis(1,lwd=0,lwd.ticks=1)
  axis(2,at=c(-0.5,0:ymax),labels=c("ref",0:ymax),las=1,lwd=0)#,lwd.ticks=1)
  axis(4,at=c(-0.5,0:ymax),labels=c("ref",0:ymax),las=1,lwd=0)#,lwd.ticks=1)
  rect(xleft=rgyw-0.5,ybottom=-0.5,xright=rgyw+3.5,ytop=ymax+0.5,col="LemonChiffon",border=F)
  rect(xleft=agct-0.5,ybottom=-0.5,xright=agct+3.5,ytop=ymax+0.5,col="LightGoldenrod",border=F)
  grid(col=grey(0.5))
  points(1:nrow(ref),rep(-0.5,nrow(ref)),col=basecolors[ref$color],pch=ascii[ref$pch],cex=0.6)
  lines(dens$x,dens$y)

  
}

dev.off()
