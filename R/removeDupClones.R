ARGS <- c(
  "mutfile","character","file path of mutation file",
  "clonefile","character","file path of summary file"
  
)

OPTS <- c(
  "perc_similar","numeric",75,"percent similarity threshold for clones to be repeats"
)

source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local("Rsub.R")

parseArgs("removeDupClones.R", ARGS, OPTS)


muts <- read.delim(mutfile,header=F,as.is=T,col.names=c("Expt","ID","Pos","Type","From","To","Size","End","Ins"))
clones <- read.delim(clonefile,header=T,as.is=T)

clones$Index <- clones$Subs + clones$DelBp + clones$InsBp

clones <- clones[rev(order(clones$Index)),]
clones$Dup <- ""

for (i in 2:nrow(clones)) {
  if (is.na(clones$Index[i]) || clones$Index[i] < 1) next
  clone1 <- clones[i,]
  muts1 <- muts[muts$ID == clone1$ID,]
  coords1 <- data.frame(start=integer(),end=integer())
  
  for (x in unlist(strsplit(clone1$Coords,","))) {
    coords1[nrow(coords1)+1,] <- as.integer(unlist(strsplit(x,"-")))
  }
      
  for (j in 1:(i-1)) {
    if (is.na(clones$Index[j]) || clones$Index[j] < 1) next
    clone2 <- clones[j,]
    muts2 <- muts[muts$ID == clone2$ID,]
    coords2 <- data.frame(start=integer(),end=integer())
    
    for (x in unlist(strsplit(clone2$Coords,","))) {
      coords2[nrow(coords2)+1,] <- as.integer(unlist(strsplit(x,"-")))
      
    }
    
    pts_possible <- 0
    pts_awarded <- 0
    
    for (k in 1:nrow(muts2)) {
      
      pos <- muts2$Pos[k]
      type <- muts2$Type[k]
      
      if (type == "sub") {
        if (nrow(coords1[ coords1$start <= pos & coords1$end >= pos,]) > 0) {
          pts_possible <- pts_possible + 1
          if (nrow( muts1[ muts1$Type == type & muts1$Pos > pos-2 & muts1$Pos < pos+2 & muts1$To == muts2$To[k],] ) > 0) {
            pts_awarded <- pts_awarded + 1
          }
        }
      } else if (type == "del") {
        end <- muts2$End[k]
        size <- muts2$Size[k]
        if (nrow(coords1[ (coords1$start <= pos & coords1$end >= pos) | (coords1$start <= end & coords1$end >= end) ,]) > 0) {
          pts_possible <- pts_possible + size
          if (nrow( muts1[ muts1$Type == type & muts1$Pos > pos-2 & muts1$Pos < pos+2 & muts1$End > end-2 & muts1$End < end+2,] ) > 0) {
            pts_awarded <- pts_awarded + size
          }
        }
      } else if (type == "ins") {
        size <- muts2$Size[k]
        insert <- muts2$Ins[k]
        if (nrow(coords1[ coords1$start <= pos & coords1$end >= pos,]) > 0) {
          pts_possible <- pts_possible + size
          if (nrow( muts1[ muts1$Type == type & muts1$Pos > pos-2 & muts1$Pos < pos+2 & as.vector(adist(muts1$Ins,insert)) < 3,] ) > 0) {
            pts_awarded <- pts_awarded + size
          }
        }
      }
    }
    
    if (pts_possible > 0 && pts_awarded/pts_possible >= perc_similar/100) {
      #cat(clones$ID[i],"is a clone of",clones$ID[j],"with",pts_awarded,"out of",pts_possible,"\n")
      clones$Dup[i] <- clones$ID[j]
      clones$Bp[i] <- 0
      break
    } else {
      #cat(clones$ID[i],"is not a clone of",clones$ID[j],"with",pts_awarded,"out of",pts_possible,"\n")   
    }
  }
  
  if (clones$Dup[i] != "") {
    while (clones$Dup[clones$ID==clones$Dup[i]] != "") {
      clones$Dup[i] <- clones$Dup[clones$ID==clones$Dup[i]]
    }
  }
  
}

clones$Index <- NULL

write.table(clones,clonefile,quote=F,sep="\t",na="",row.names=F,col.names=T)