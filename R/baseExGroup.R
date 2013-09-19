ARGS <- c(
  "exptfile","character","file path of expt file",
  "metafile","character","file path of meta file",
  "statsfile","character","output"
  
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

parseArgs("mutationStats.R", ARGS, OPTS)

grouping <- unlist(strsplit(grouping,","))

meta <- read.delim(metafile,header=T,as.is=T)
expts <- read.delim(exptfile,header=T,as.is=T)

meta <- meta[order(meta$experiment),]
expts <- expts[order(expts$Expt),]

if (!all.equal(meta$experiment,expts$Expt)) { stop("experiment names do not match")}

expts <- cbind(expts,meta[,grouping])



# agg_stats <- aggregate(cbind(Clones,Bp,Subs,Del,DelBp,Ins,InsBp,RefA,RefC,RefG,RefT,RefN) ~ genotype + allele + mouse + tissue + pna, expts, sum)
# agg_stats <- aggregate(eval(as.formula(paste( "cbind(",paste(colnames(clones[,sapply(clones, is.numeric)]),collapse="," ),") ~ ", paste(grouping,collapse=" + "),sep=""))), clones, sum)
agg_stats <- aggregate(eval(as.formula(paste( "cbind(",paste(colnames(expts[,sapply(expts, is.numeric)]),collapse="," ),") ~ ", paste(grouping,collapse=" + "),sep=""))), expts, sum)
write.table(agg_stats,statsfile,quote=F,sep="\t",na="",row.names=F,col.names=T)