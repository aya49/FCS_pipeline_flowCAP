## Input: original normalized count matrix --> Output: count statistic plots
# aya43@sfu.ca 20161220

#Directory
root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

#Options
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

#Input
matrix_dir = paste(result_dir, "/matrix", sep="")

#Output
png_dir = paste(result_dir, "/count_stats.png", sep="")

#Libraries/Functions
library(foreach)
library(doMC)
library(colorspace)
source("~/projects/IMPC/code/_funcAlice.R")

#Setup Cores
no_cores = detectCores()-1
registerDoMC(no_cores)











#Options for script
countThres = seq(20,2000,20)
matrix_count = c("CountAdj") #specify countmatrix to use



#Prepare data
mm = get(load(paste0(matrix_dir, matrix_count, ".Rdata")))
markers = unlist(strsplit(colnames(mm)[which.max(nchar(colnames(mm)))],"[+-]"))
phenolevel = getPhen(colnames(mm), phenotype=F, markers=markers)$phenoLevel

k0 = seq(1,max(phenolevel)) #layers to plot











start = Sys.time()

start1 = Sys.time()

#trim low count phenotypes
underCountlist = foreach(c = 1:length(countThres)) %dopar% { return(which(apply(mm,2,function(x) all(x<=countThres[c])))) }

#trim layers
underklist = foreach(k = 1:length(k0)) %dopar% { return(which(phenolevel<=k0[k])) }

TimeOutput(start1)
start1 = Sys.time()

## collate data according to count and layer thresholds
underBoth = foreach(c = 1:length(countThres)) %dopar% {
  ub = rep(0,length(k0))
  for (k in 1:length(k0)) {
    l = length(intersect(underCountlist[[c]],underklist[[k]]))
    if (!(length(l)>0)) l = 0
    ub[k] = l
  }
  return(ub)
}
underBoth = do.call(rbind, underBoth)

TimeOutput(start1)

## plot
png(png_dir, width=400, height=400)
colour = rainbow_hcl(ncol(underBoth))
plot(countThres, underBoth[,ncol(underBoth)], type="l", col=colour[ncol(underBoth)], 
     xlab="Count threshold", ylab="# of cell populations with count <= Count threshold, for all samples")
for(i in (ncol(underBoth)-1):1) { lines(countThres, underBoth[,i], col=colour[i]) }
legend("topleft",legend=paste0("<= phenolevel ", c(1:ncol(underBoth))), fill=colour)
graphics.off()

TimeOutput(start)





