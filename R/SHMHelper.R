getBases <- function () {
  bases <- c("A","C","G","T","N")
}
getAscii <- function () {
  ascii <- c(65,67,71,84,78)
}
getBasecolors <- function () {
  basecolors <- brewer.pal(7,"Set1")
}

getMutsFromRead <- function(read,muts) {
  return(muts[muts$Expt == read$Expt[1] & muts$Read == read$Read[1],])
}

getReadFromMut <- function(mut,reads) {
  return(reads[reads$Expt == mut$Expt[1] & reads$Read == mut$Read[1],])
}

getMutsFromReads <- function(reads,muts) {
  return(merge(reads,muts,by=1:2)[,colnames(muts)])
}


convertCoords <- function(coordlist) {
  coordlist <- strsplit(coordlist,",")
  coords <- ldply(lapply(1:length(coordlist),function(i,x) { 
                                              df <- ldply(strsplit(x[[i]],"-"))
                                              df <- cbind(i,df)
                                              return(df)},coordlist))
  colnames(coords) <- c("i","start","end")
  coords$start <- as.integer(coords$start)
  coords$end <- as.integer(coords$end)
  return(coords)
}

invertCoords <- function(coordlist,refseq) {
  coords <- convertCoords(coordlist)
  incoords <- coords[c(),]
  for (i in unique(coords$i)) {
    readcoords <- coords[coords$i==i,]
    
    if ( readcoords$start[1] > 1 ) incoords <- rbind(incoords,c(i,1,readcoords$start[1]-1))
    
    if (nrow(readcoords) > 1) {
      for(j in 2:nrow(readcoords)) {
        incoords <- rbind(incoords,c(i,readcoords$end[j-1]+1,readcoords$start[j]-1))
      }
    }
    
    if (nchar(as.character(refseq)) > readcoords$end[nrow(readcoords)]) incoords <- rbind(incoords,c(i,readcoords$end[nrow(readcoords)]+1,nchar(as.character(refseq))))
  }
  colnames(incoords) <- c("i","start","end")
  return(incoords)
}


createMutationMatrix <- function(reads,muts,refseq,tstart,tend) {
  mutmat <- matrix(".",ncol=nchar(as.character(refseq)),nrow=nrow(reads))
  incoords <- invertCoords(reads$Coords,refseq)
  
  for (i in 1:nrow(incoords)) {
    mutmat[incoords[i,"i"],incoords[i,"start"]:incoords[i,"end"]] <- ""
  }
  
  
  subs <- muts[muts$Type=="sub",]
  dels <- muts[muts$Type=="del",]
  
  if (nrow(subs) > 0) {
    for (i in 1:nrow(subs)) {
      sub <- subs[i,]
      read <- getReadFromMut(sub,reads)
      row <- which(reads$Expt == read$Expt[1] & reads$Read == read$Read[1])
      mutmat[row,sub$Pos] <- sub$To
    }
  }
  
  if (nrow(dels) > 0) {
    for (i in 1:nrow(dels)) {
      del <- dels[i,]
      read <- getReadFromMut(del,reads)
      row <- which(reads$Expt == read$Expt[1] & reads$Read == read$Read[1])
      if (del$Size > 1) {
        mutmat[row,del$Pos] <- "<"
        if (del$Size > 2) mutmat[row,(del$Pos+1):(del$Pos+del$Size-2)] <- "-"
        mutmat[row,del$Pos+del$Size-1] <- ">"
      } else {
        mutmat[row,del$Pos] <- "<>"
      }
    }
  }
  
  if (tstart > 1) {
    mutmat[,1:(tstart-1)] <- ""
  }
  
  if (tend < ncol(mutmat)) {
    mutmat[,(tend+1):ncol(mutmat)] <- ""
  }
  
  return(mutmat)
}

createInsertionMatrix <- function(reads,muts,refseq,tstart,tend) {
  insmat <- matrix(".",ncol=nchar(as.character(refseq)),nrow=nrow(reads))
  incoords <- invertCoords(reads$Coords,refseq)
  
  for (i in 1:nrow(incoords)) {
    insmat[incoords[i,"i"],incoords[i,"start"]:incoords[i,"end"]] <- ""
  }
  
  inss <- muts[muts$Type=="ins",]
  
  if (nrow(inss) > 0) {
    for (i in 1:nrow(inss)) {
      ins <- inss[i,]
      read <- getReadFromMut(ins,reads)
      row <- which(reads$Expt == read$Expt[1] & reads$Read == read$Read[1])
      insmat[row,ins$Pos] <- ins$Ins
    }
  }
  
  
  if (tstart > 1) {
    insmat[,1:(tstart-1)] <- ""
  }
  
  if (tend < ncol(insmat)) {
    insmat[,(tend+1):ncol(insmat)] <- ""
  }
  
  return(insmat)
}


compareMutations <- function(mutpair) {
  apply(mutpair,2,function(x) {
    mutchars <- c("A","C","G","T")
    delchars <- c("<",">","<>")
    result <- rep(0,3)
    
    if (all(x %in% c("","."))) return(result)
    
    if (x[1] %in% mutchars && x[2] != "") result[1] <- result[1] + 1
    if (x[2] %in% mutchars && x[1] != "") result[2] <- result[2] + 1
    if (x[1] %in% mutchars && x[2] %in% mutchars && x[1] == x[2]) result[3] <- result[3] + 1
    
    if (x[1] %in% delchars && x[2] != "") result[1] <- result[1] + 1
    if (x[2] %in% delchars && x[1] != "") result[2] <- result[2] + 1
    if ( (grepl(delchars[1],x[1]) && grepl(delchars[1],x[2])) ||
           (grepl(delchars[2],x[1]) && grepl(delchars[2],x[2])) ) result[3] <- result[3] + 1
    
    if (grepl("[ACGTN]",x[3]) && x[4] != "") result[1] <- result[1] + 1
    if (grepl("[ACGTN]",x[4]) && x[3] != "") result[2] <- result[2] + 1
    if (grepl("[ACGTN]",x[3]) && grepl("[ACGTN]",x[4]) && adist(x[3],x[4]) < 2) result[3] <- result[3] + 1
    return(result)
  })
}


calculateProfile <- function(mutmat,refseq) {
  profile <- data.frame(Pos=1:nchar(as.character(refseq)))
  profile$Base <- unlist(strsplit(as.character(refseq),""))
  profile$Reads <- apply(mutmat,2,function(x) {sum(grepl("[.ACGTN]",x))})
  profile$Subs <- apply(mutmat,2,function(x) {sum(grepl("[ACGT]",x))})
  profile$Y <- ifelse(profile$Reads > 0, profile$Subs/profile$Reads, 0)
  return(profile)
}

calculateDeletionProfile <- function(mutmat,refseq) {
  profile <- data.frame(Pos=1:nchar(as.character(refseq)))
  profile$Base <- unlist(strsplit(as.character(refseq),""))
  profile$Reads <- colSums(mutmat != "")
  profile$Dels <- apply(mutmat,2,function(x) {sum(grepl("[<->]",x))})
  profile$Y <- ifelse(profile$Reads > 0, profile$Dels/profile$Reads, 0)
  return(profile)
}

tictactoePlot <- function (subs, dels, blocks, ref, tstart, tend, plotrows, cloneIDs) {
  
  bases <- getBases()
  ascii <- getAscii()
  basecolors <- getBasecolors()
  refseq <- paste(ref$Base,collapse="")
  par(mai=c(0.2,0.75,0.2,0.75),omi=c(0.5,0,0,0))
  
  layout(as.matrix(1:plotrows,ncol=1,nrow=plotrows))
  
  
  ymax <- max(5.5,length(cloneIDs))
  
  
  dels$y <- match(dels$Clone,cloneIDs)
  
  
  subs$y <- match(subs$Clone,cloneIDs)
  subs$color <- match(subs$To,bases)
  subs$pch <- match(subs$To,bases)
  
  
  
  

  
  agct <- unlist(gregexpr("AGCT",as.character(refseq)))
  rgyw <- unique(c(unlist(gregexpr("[AG]G[CT][AT]",as.character(refseq))),unlist(gregexpr("[AT][AG]C[CT]",as.character(refseq)))))
  
  rowwidth <- ceiling((tend-tstart+1)/plotrows)
  
  tstarts <- tstart+rowwidth*0:(plotrows-1)
  tends <- tstarts+rowwidth-1
  
  
  
  for (i in 1:plotrows) {
    
    plot(c(),c(),ylab="",xlab="",xaxt="n",xlim=c(max(1,tstarts[i]-2),min(nchar(refseq),tends[i]+2)),ylim=c(0,ymax),xaxs="i",bty="o")
    axis(1,lwd=0,lwd.ticks=1)
    
    rect(xleft=rgyw-0.5,ybottom=-1,xright=rgyw+3.5,ytop=ymax,col=rgb(254,217,142,max=255),border=F)
    rect(xleft=agct-0.5,ybottom=-1,xright=agct+3.5,ytop=ymax,col=rgb(254,153,41,max=255),border=F)
    
    
    grid(ny=0,col=grey(0.5),lty=3)
    points(1:nrow(ref),rep(0,nrow(ref)),col=basecolors[ref$color],pch=ascii[ref$pch],cex=0.6)
    if (length(cloneIDs) > 0) {
      rect(xleft=blocks$Start-0.5,ybottom=blocks$Clone-0.5,xright=blocks$End+0.5,ytop=blocks$Clone+0.5,col=grey(0.1,0.25),border=F)
      points(subs$Pos,subs$y,col=basecolors[subs$color],pch=ascii[subs$pch],cex=0.5)
      segments(x0=dels$Pos-0.5,y0=dels$y,x1=dels$End+0.5,col=basecolors[7],lwd=1)
    } else {
      text(x=(tstarts[i]+tends[i])/2,y=ymax/2-0.5,"No Mutations To Display")
    }
  }
  
}

connectfourSubPlot <- function (subs, blocks, ref, tstart, tend, plotrows, cloneIDs) {
  
  bases <- getBases()
  ascii <- getAscii()
  basecolors <- getBasecolors()
  refseq <- paste(ref$Base,collapse="")
  par(mai=c(0.2,0.75,0.2,0.75),omi=c(0.5,0,0,0))
  
  layout(as.matrix(1:plotrows,ncol=1,nrow=plotrows))
  
  agct <- unlist(gregexpr("AGCT",as.character(refseq)))
  rgyw <- unique(c(unlist(gregexpr("[AG]G[CT][AT]",as.character(refseq))),unlist(gregexpr("[AT][AG]C[CT]",as.character(refseq)))))
  rowwidth <- ceiling((tend-tstart+1)/plotrows)
  
  tstarts <- tstart+rowwidth*0:(plotrows-1)
  tends <- tstarts+rowwidth-1
  
  
  ref$dens_denom <- 0
  ref$dens_numer <- 0
  ref$dens <- 0
  
  if (nrow(subs) > 0) {
    ref$dens_numer <- hist(subs$Pos,0:nrow(ref),plot=F)$counts
  }
  ref$dens_denom <- length(cloneIDs) - unlist(lapply(1:nrow(ref),function(x) {
    return(nrow(blocks[blocks$Start <= ref$Pos[x] & blocks$End >= ref$Pos[x],]))
  }))
  
  ref$dens <- ifelse(ref$dens_denom==0,0,ref$dens_numer/ref$dens_denom)
  
  dens <- data.frame(x=c(ref$Pos-0.49,ref$Pos+0.49),y=c(ref$dens,ref$dens))
  dens <- dens[order(dens$x),]
  ymax <- max(0.05,1.25*dens$y)
  
  refy <- -ymax/20
  
  #   rgyw_tops <- ymax
  #   agct_tops <- ymax
  #   rgyw_bottoms <- 0.82*ymax
  #   agct_bottoms <- 0.82*ymax
  
  #   Bottom of the plot - below the sequence  
  rgyw_bottoms <- 2*refy
  agct_bottoms <- 2*refy
  
  #   Adjust to height of the peak
  rgyw_tops <- unlist(lapply(rgyw,function(x,dens) {max(dens[x:(x+3)])},ref$dens)) + 0.025*ymax
  agct_tops <- unlist(lapply(agct,function(x,dens) {max(dens[x:(x+3)])},ref$dens)) + 0.025*ymax
  
  
  
  for (i in 1:plotrows) {
    plot(c(),c(),ylab="",xaxt="n",xlab="",xlim=c(max(1,tstarts[i]-2),min(nchar(refseq),tends[i]+2)),ylim=c(refy,ymax),xaxs="i",bty="o")
    axis(1,lwd=0,lwd.ticks=1)
    rect(xleft=rgyw-0.5,ybottom=rgyw_bottoms,xright=rgyw+3.5,ytop=rgyw_tops,col=rgb(254,217,142,max=255),border=rgb(254,217,142,max=255),lwd=2)
    rect(xleft=agct-0.5,ybottom=agct_bottoms,xright=agct+3.5,ytop=agct_tops,col=rgb(254,153,41,max=255),border=rgb(254,153,41,max=255),lwd=2)
    
    grid(col=grey(0.5))
    points(1:nrow(ref),rep(refy,nrow(ref)),col=basecolors[ref$color],pch=ascii[ref$pch],cex=0.6)
    lines(dens$x,dens$y)
  }

  return(ref)
}

connectfourDelPlot <- function (dels, blocks, ref, tstart, tend, plotrows, cloneIDs) {
  
  bases <- getBases()
  ascii <- getAscii()
  basecolors <- getBasecolors()
  refseq <- paste(ref$Base,collapse="")
  par(mai=c(0.2,0.75,0.2,0.75),omi=c(0.5,0,0,0))
  
  layout(as.matrix(1:plotrows,ncol=1,nrow=plotrows))
  
  agct <- unlist(gregexpr("AGCT",as.character(refseq)))
  rgyw <- unique(c(unlist(gregexpr("[AG]G[CT][AT]",as.character(refseq))),unlist(gregexpr("[AT][AG]C[CT]",as.character(refseq)))))
  rowwidth <- ceiling((tend-tstart+1)/plotrows)
  
  tstarts <- tstart+rowwidth*0:(plotrows-1)
  tends <- tstarts+rowwidth-1
  
  
  ref$dens_denom <- 0
  ref$dens_numer <- 0
  ref$dens <- 0
  
  if (nrow(dels) > 0) {
    dels_expand <- c()
    for (i in 1:nrow(dels)) {
      dels_expand <- c(dels_expand,seq(dels$Pos[i],dels$Pos[i]+dels$Size[i]-1))
    }
    ref$dens_numer <- hist(dels_expand,c(0,1:nrow(ref)),plot=F)$counts
  }
  
  
  ref$dens_denom <- length(cloneIDs) - unlist(lapply(1:nrow(ref),function(x) {
    return(nrow(blocks[blocks$Start <= ref$Pos[x] & blocks$End >= ref$Pos[x],]))
  }))
  
  ref$dens <- ifelse(ref$dens_denom==0,0,ref$dens_numer/ref$dens_denom)
  
  dens <- data.frame(x=c(ref$Pos-0.49,ref$Pos+0.49),y=c(ref$dens,ref$dens))
  dens <- dens[order(dens$x),]
  ymax <- max(0.01,1.25*dens$y)
  
  refy <- -ymax/20
  
  #   rgyw_tops <- ymax
  #   agct_tops <- ymax
  #   rgyw_bottoms <- 0.82*ymax
  #   agct_bottoms <- 0.82*ymax
  
  #   Bottom of the plot - below the sequence  
  rgyw_bottoms <- 2*refy
  agct_bottoms <- 2*refy
  
  #   Adjust to height of the peak
  rgyw_tops <- unlist(lapply(rgyw,function(x,dens) {max(dens[x:(x+3)])},ref$dens)) + 0.025*ymax
  agct_tops <- unlist(lapply(agct,function(x,dens) {max(dens[x:(x+3)])},ref$dens)) + 0.025*ymax
  
  
  
  for (i in 1:plotrows) {
    plot(c(),c(),ylab="",xaxt="n",xlab="",xlim=c(max(1,tstarts[i]-2),min(nchar(refseq),tends[i]+2)),ylim=c(refy,ymax),xaxs="i",bty="o")
    axis(1,lwd=0,lwd.ticks=1)
    rect(xleft=rgyw-0.5,ybottom=rgyw_bottoms,xright=rgyw+3.5,ytop=rgyw_tops,col=rgb(254,217,142,max=255),border=rgb(254,217,142,max=255),lwd=2)
    rect(xleft=agct-0.5,ybottom=agct_bottoms,xright=agct+3.5,ytop=agct_tops,col=rgb(254,153,41,max=255),border=rgb(254,153,41,max=255),lwd=2)
    
    grid(col=grey(0.5))
    points(1:nrow(ref),rep(refy,nrow(ref)),col=basecolors[ref$color],pch=ascii[ref$pch],cex=0.6)
    lines(dens$x,dens$y)
  }
  
  
}