# aya43@sfu.ca 20161220
# Reads in flowCAP-II AML flowtype files and collates/save cell counts into matrix

root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

options(stringsAsFactors=FALSE)
options(device="cairo")
options(na.rm=T)

#Input
ft_dir = "data/FT"
csv_dir = "data/AML.csv"

#Output
phenoMeta_dir = paste(result_dir, "/phenoMeta.Rdata", sep="")
phenoMeta_dircsv = paste(result_dir, "/phenoMeta.csv", sep="")
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="")
sampleMeta_dircsv = paste(result_dir, "/sampleMeta.csv", sep="")
markers_dir = paste(result_dir, "/markers.Rdata", sep="")
matrixCount_dir = paste(result_dir, "/matrixCount.Rdata", sep="")
matrixCount_dircsv = paste(result_dir, "/matrixCount.csv", sep="")
matrixProp_dir = paste(result_dir, "/matrixProp.Rdata", sep="")
matrixProp_dircsv = paste(result_dir, "/matrixProp.csv", sep="")

library(flowCore)
library(flowType)
source("code/_funcAlice.R")



start = Sys.time()

ftFile_dir = sort( dir(ft_dir, pattern<-".Rda", all.files<-TRUE, full.names<-TRUE, recursive<-TRUE) )
ftGT = folderNames(ftFile_dir)
ftFileNames = fileNames(ftFile_dir, "Rda")

sampleMetatemp = data.frame(read.csv(csv_dir))
colnames(sampleMetatemp) = c("filename", "tube", "specimen", "aml")

ft = get(load(ftFile_dir[1]))
markers = ft@MarkerNames
save(markers, file=markers_dir)

phenotype = rownames(ft@MFIs)
phenocode = unlist(lapply(phenotype, function(x){return( encodePhenotype(x, markers) )}))
phenoLevel = unlist(lapply(phenocode, function(x){return(length(markers) - charOccurences("0", x)) } ))
phenoMeta = data.frame(phenotype, phenocode, phenoLevel, stringsAsFactors=F)

matrixCount = NULL
sampleMeta = NULL

for (i in 1:length(ftFile_dir)) {
  ft = get(load(ftFile_dir[i]))
  npheno = length(ft@CellFreqs)
  matrixCount = rbind(matrixCount, ft@CellFreqs)
  cat("\n", i, "Loaded file: ", ftGT[i], "/", ftFileNames[i], "; # of phenotypes x samples: ", npheno, " x ", length(ft), "", sep <-"")
  
  fn = ftFileNames[i]
  fn = strsplit(fn, "\\.")[[1]][1]
  fn = strsplit(fn, "S")[[1]]
  tube = substr(fn[1],2,nchar(fn[1]))
  specimen = gsub("FT","",fn[2])
  sampleMeta = rbind(sampleMeta, sampleMetatemp[which(as.numeric(sampleMetatemp$specimen)==as.numeric(specimen) & as.numeric(sampleMetatemp$tube)==as.numeric(tube)),])
  
  rm(ft)
}

#delete phenotypes with no cells
all0Phen = which(apply(matrixCount[,-1], 1, function(x) all(x==0))==T)
if (length(all0Phen)>0) {
  matrixCount = matrixCount[,-all0Phen]
  phenoMeta[-all0Phen,]
}

#order matrix phenotypes
phenoorder = order(phenoMeta$phenocode)
phenoMeta = phenoMeta[phenoorder,]
matrixCount = matrixCount[,phenoorder]

rownames(matrixCount) = ftFileNames
colnames(matrixCount) = phenoMeta$phenotype
matrixProp = matrixCount%*%diag(1/matrixCount[1,])

save(phenoMeta, file=phenoMeta_dir)
write.csv(phenoMeta, file=phenoMeta_dircsv, row.names=F)
save(sampleMeta, file=sampleMeta_dir)
write.csv(sampleMeta, file=sampleMeta_dircsv, row.names=F)
save(matrixCount, file=matrixCount_dir)
write.csv(matrixCount, file=matrixCount_dircsv, row.names=F)
save(matrixProp, file=matrixProp_dir)
write.csv(matrixProp, file=matrixProp_dircsv, row.names=F)

TimeOutput(start)






