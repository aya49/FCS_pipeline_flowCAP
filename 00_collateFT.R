## Input: flowtype file and meta data paths --> Output: matrixCount, sampleMeta, phenoMeta
#aya43@sfu.ca 20161220

#Directory
root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

#Options
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

#Input
data_dir = "/mnt/f/Brinkman group/current/Alice/flowCAP-II/data" #main data directory
ft_dir = paste0(data_dir,"/FT") #flowtype file directory
csv_dir = paste0(data_dir,"/AML.csv") #meta file directory

#Output
phenoMeta_dir = paste(result_dir, "/phenoMeta", sep="")
sampleMeta_dir = paste(result_dir, "/sampleMeta", sep="")
markers_dir = paste(result_dir, "/markers", sep="")
matrixCount_dir = paste(result_dir, "/matrixCount", sep="")
matrixProp_dir = paste(result_dir, "/matrixProp", sep="")

#Libraries/Functions
library(flowCore)
library(flowType)
library(foreach)
library(doMC)
source("~/projects/IMPC/code/_funcAlice.R")

#Setup Cores
no_cores = detectCores()-1
registerDoMC(no_cores)










#Prepare directories
ftFile_dir = sort( dir(ft_dir, pattern=".Rda", all.files=T, full.names=T, recursive=T) )
ftGT = folderNames(ftFile_dir)
ftFileNames = fileNames(ftFile_dir, "Rda")










start = Sys.time()

#get meta data on files for sampleMeta
sampleMetatemp = data.frame(read.csv(csv_dir))
colnames(sampleMetatemp) = c("fileName", "tube", "specimen", "aml") #rename columns

#get meta data on phenotypes (cell populations) for phenoMeta
ft = get(load(ftFile_dir[1]))
markers = ft@MarkerNames
save(markers, file=paste0(markers_dir,".Rdata"))

phenotype = rownames(ft@MFIs)
phenocode = unlist(lapply(phenotype, function(x){return( encodePhenotype(x, markers) )}))
phenolevel = unlist(lapply(phenocode, function(x){return(length(markers) - charOccurences("0", x)) } ))
phenoMeta = data.frame(phenotype, phenocode, phenolevel, stringsAsFactors=F)

#load flowtype files for matrixCount & fill in sampleMeta
matrixCount = NULL
sampleMeta = NULL

result = foreach (i=1:length(ftFile_dir), .combine="rbind") %dopar% {
  ft = get(load(ftFile_dir[i]))
  print(ft@MarkerNames)
  npheno = length(ft@CellFreqs)
  # matrixCount = rbind(matrixCount, ft@CellFreqs)
  cat("\n", i, "Loaded file: ", ftGT[i], "/", ftFileNames[i], "; # of phenotypes x samples: ", npheno, " x ", length(ft), "", sep="")
  
  fn = ftFileNames[i]
  fn = strsplit(fn, "\\.")[[1]][1]
  fn = strsplit(fn, "S")[[1]]
  tube = substr(fn[1],2,nchar(fn[1]))
  specimen = gsub("FT","",fn[2])
  # sampleMeta = rbind(sampleMeta, sampleMetatemp[which(as.numeric(sampleMetatemp$specimen)==as.numeric(specimen) & as.numeric(sampleMetatemp$tube)==as.numeric(tube)),])
  return(c(sampleMetatemp[which(as.numeric(sampleMetatemp$specimen)==as.numeric(specimen) & as.numeric(sampleMetatemp$tube)==as.numeric(tube)),], ft@CellFreqs))
  
  rm(ft) #save memory
}

#collate matrixCount
matrixCount = as.matrix(apply(result[,(ncol(sampleMetatemp)+1):ncol(result)], 2, as.numeric))
sm = result[,1:(ncol(sampleMetatemp))]

#delete phenotypes with no cells for all files
all0Phen = apply(matrixCount[-1,], 2, function(x) all(x==0))
if (sum(all0Phen)>0) {
  matrixCount = matrixCount[,!all0Phen]
  phenoMeta[!all0Phen,]
}

#order phenoMeta and matrixCount
phenoorder = order(phenoMeta$phenocode)
phenoMeta = phenoMeta[phenoorder,]
matrixCount = matrixCount[,phenoorder]

#make sampleMeta
rownames(matrixCount) = ftFileNames
sampleMeta = as.data.frame(ftFileNames)
sampleMeta$tube = unlist(sm[,2])
sampleMeta$specimen = unlist(sm[,3])
sampleMeta$aml = unlist(sm[,4])
colnames(sampleMeta) = colnames(sampleMetatemp)

rownames(matrixCount) = sampleMeta[,1] = ftFileNames
colnames(matrixCount) = phenoMeta$phenotype
matrixProp = matrixCount/matrixCount[,1]
dimnames(matrixProp) = dimnames(matrixCount)

#save
save(phenoMeta, file=paste0(phenoMeta_dir,".Rdata")); write.csv(phenoMeta, file=paste0(phenoMeta_dir,".csv"), row.names=F)
save(sampleMeta, file=paste0(sampleMeta_dir,".Rdata")); write.csv(sampleMeta, file=paste0(sampleMeta_dir,".csv"), row.names=F)
save(matrixCount, file=paste0(matrixCount_dir,".Rdata")); write.csv(matrixCount, file=paste0(matrixCount_dir,".csv"), row.names=F)
save(matrixProp, file=paste0(matrixProp_dir,".Rdata")); write.csv(matrixProp, file=paste0(matrixProp_dir,".csv"), row.names=F)

TimeOutput(start)






