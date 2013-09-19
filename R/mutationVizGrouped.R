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
  "plotrows","numeric",4,"Rows on plot"
)

source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local("Rsub.R")

parseArgs("mutationVizGrouped.R", ARGS, OPTS)


library(Biostrings)
library(RColorBrewer)
basecolors <- brewer.pal(7,"Set1")
bases <- c("A","C","G","T","N")
ascii <- c(65,67,71,84,78)


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
  
  group_muts <- muts[muts$group == group,]
  group_clones <- clones[clones$group == group,]
  
  cloneIDs <- group_clones$Clone[group_clones$Bp > 0]
  
  group_subs <- group_muts[group_muts$Type == "sub",]
  group_subs$y <- match(group_subs$Clone,cloneIDs)
  group_subs$color <- match(group_subs$To,bases)
  group_subs$pch <- match(group_subs$To,bases)
  
  group_dels <- group_muts[group_muts$Type == "del",]
  group_dels$y <- match(group_dels$Clone,cloneIDs)
  
  ref <- data.frame(Pos=1:nchar(as.character(refseq)),Base=strsplit(as.character(refseq),""))
  colnames(ref) <- c("Pos","Base")
  ref$color <- match(ref$Base,bases)
  ref$pch <- match(ref$Base,bases)
  
  agct <- unlist(gregexpr("AGCT",as.character(refseq)))
  rgyw <- unique(c(unlist(gregexpr("[AG]G[CT][AT]",as.character(refseq))),unlist(gregexpr("[AT][AG]C[CT]",as.character(refseq)))))
  dgyw <- unique(c(unlist(gregexpr("[AGT]G[CT][AT]",as.character(refseq))),unlist(gregexpr("[AT][AG]C[ACT]",as.character(refseq)))))
  ggg <- unique(c(unlist(gregexpr("GGG",as.character(refseq))),unlist(gregexpr("CCC",as.character(refseq)))))
  taa <- unique(c(unlist(gregexpr("TAA",as.character(refseq))),unlist(gregexpr("TTA",as.character(refseq)))))
  
  rowwidth <- ceiling((tend-tstart+1)/plotrows)
  
  tstarts <- tstart+rowwidth*0:(plotrows-1)
  tends <- tstarts+rowwidth-1
  
  blocks <- data.frame(Clone=integer(),Start=integer(),End=integer())
  
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
  par(mai=c(0.2,0.75,0.2,0.75),omi=c(0.5,0,0,0))
  
  layout(as.matrix(1:plotrows,ncol=1,nrow=plotrows))
  ymax <- max(5.5,length(cloneIDs))
  
  for (i in 1:plotrows) {
    
    plot(c(),c(),ylab="",xlab="",xaxt="n",xlim=c(max(1,tstarts[i]-2),min(nchar(refseq),tends[i]+2)),ylim=c(0,ymax),xaxs="i",bty="o")
    axis(1,lwd=0,lwd.ticks=1)
    #axis(2,at=0:length(cloneIDs),labels=c("ref",cloneIDs),las=1,lwd=0)#,lwd.ticks=1)
    #axis(4,at=0:length(cloneIDs),labels=c("ref",cloneIDs),las=1,lwd=0)#,lwd.ticks=1)
    rect(xleft=rgyw-0.5,ybottom=-1,xright=rgyw+3.5,ytop=ymax,col=rgb(254,217,142,max=255),border=F)
    rect(xleft=agct-0.5,ybottom=-1,xright=agct+3.5,ytop=ymax,col=rgb(254,153,41,max=255),border=F)
    rect(xleft=ggg-0.5,ybottom=-1,xright=ggg+2.5,ytop=0.15*ymax,col=rgb(204, 235, 197,max=255),border=F)
    rect(xleft=taa-0.5,ybottom=-1,xright=taa+2.5,ytop=0.15*ymax,col=rgb(179, 205, 227,max=255),border=F)
    
    grid(ny=0,col=grey(0.5),lty=3)
    points(1:nrow(ref),rep(0,nrow(ref)),col=basecolors[ref$color],pch=ascii[ref$pch],cex=0.6)
    if (length(cloneIDs) > 0) {
      rect(xleft=blocks$Start-0.5,ybottom=blocks$Clone-0.5,xright=blocks$End+0.5,ytop=blocks$Clone+0.5,col=grey(0.1,0.25),border=F)
      points(group_subs$Pos,group_subs$y,col=basecolors[group_subs$color],pch=ascii[group_subs$pch],cex=0.5)
      segments(x0=group_dels$Pos-0.5,y0=group_dels$y,x1=group_dels$End+0.5,col=basecolors[7],lwd=1)
    } else {
      text(x=(tstarts[i]+tends[i])/2,y=ymax/2-0.5,"No Mutations To Display")
    }
  }
  dev.off()
  
  pdf(paste(outdir,"/",paste(groups[group,],collapse="_"),"_mutDens.pdf",sep=""),height=8,width=11)
  par(mai=c(0.2,0.75,0.2,0.75),omi=c(0.5,0,0,0))
  
  layout(as.matrix(1:plotrows,ncol=1,nrow=plotrows))
  
  ref$dens_denom <- 0
  ref$dens_numer <- 0
  ref$dens <- 0
  
  if (nrow(group_subs) > 0) {
    ref$dens_numer <- hist(group_subs$Pos,0:nrow(ref),plot=F)$counts
  }
  ref$dens_denom <- length(cloneIDs) - unlist(lapply(1:nrow(ref),function(x) {
    return(nrow(blocks[blocks$Start <= ref$Pos[x] & blocks$End >= ref$Pos[x],]))
  }))
  
  ref$dens <- ifelse(ref$dens_denom==0,0,ref$dens_numer/ref$dens_denom)
  
  dens <- data.frame(x=c(ref$Pos-0.25,ref$Pos+0.25),y=c(ref$dens,ref$dens))
  dens <- dens[order(dens$x),]
  ymax <- max(0.01,dens$y)
  
  refy <- -ymax/20
  
  for (i in 1:plotrows) {
    plot(c(),c(),ylab="",xaxt="n",xlab="",xlim=c(max(1,tstarts[i]-2),min(nchar(refseq),tends[i]+2)),ylim=c(refy,ymax),xaxs="i",bty="o")
    axis(1,lwd=0,lwd.ticks=1)
    #axis(2,at=c(-0.5,0:ymax),labels=c("ref",0:ymax),las=1,lwd=0)#,lwd.ticks=1)
    #axis(4,at=c(-0.5,0:ymax),labels=c("ref",0:ymax),las=1,lwd=0)#,lwd.ticks=1)
# 	rect(xleft=dgyw-0.5,ybottom=refy,xright=dgyw+3.5,ytop=ymax,col=rgb(255,247,188,max=255),border=F)
    rect(xleft=rgyw-0.5,ybottom=2*refy,xright=rgyw+3.5,ytop=ymax,col=rgb(254,217,142,max=255),border=F)
    rect(xleft=agct-0.5,ybottom=2*refy,xright=agct+3.5,ytop=ymax,col=rgb(254,153,41,max=255),border=F)
    rect(xleft=ggg-0.5,ybottom=2*refy,xright=ggg+2.5,ytop=0.15*ymax,col=rgb(204, 235, 197,max=255),border=F)
    rect(xleft=taa-0.5,ybottom=2*refy,xright=taa+2.5,ytop=0.15*ymax,col=rgb(179, 205, 227,max=255),border=F)
    
    
    grid(col=grey(0.5))
    points(1:nrow(ref),rep(refy,nrow(ref)),col=basecolors[ref$color],pch=ascii[ref$pch],cex=0.6)
    lines(dens$x,dens$y)
  }
  dev.off()
  
  if (nrow(group_dels) > 0) { 
    pdf(paste(outdir,"/",paste(groups[group,],collapse="_"),"_delHist.pdf",sep=""),height=8,width=11)
    hist(group_dels$Size,breaks=20,freq=T,xlab="Deletion Size",ylab="Frequency")
    dev.off()
  }
  
  pdf(paste(outdir,"/",paste(groups[group,],collapse="_"),"_cloneSubHist.pdf",sep=""),height=8,width=11)
  hist(group_clones$Subs,breaks=50,freq=T,xlab="Number of Substitutions",ylab="Frequency")
  dev.off()
}

