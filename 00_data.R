## Input: flowtype file and meta data paths --> Output: feat_file_cell_count, meta_file, meta_cell
#aya43@sfu.ca 20161220

## root directory
root = "~/projects/flowCAP-II"
setwd(root)

result_dir = "result"; suppressWarnings(dir.create (result_dir))
data_dir = "/mnt/f/Brinkman group/current/Alice/flowCAP-II/data" #main data directory



## input directories
ft_dir = paste0(data_dir,"/FT") #flowtype file directory
csv_dir = paste0(data_dir,"/AML.csv") #meta file directory

## output directories
meta_dir = paste0(result_dir,"/meta"); dir.create(meta_dir, showWarnings=F)
meta_cell_dir = paste(meta_dir, "/cell", sep="")
meta_file_dir = paste(meta_dir, "/file", sep="")
markers_dir = paste(meta_dir, "/cell_markers", sep="")

feat_dir = paste(result_dir, "/feat", sep=""); dir.create(feat_dir, showWarnings=F)
feat_file_cell_count_dir = paste(feat_dir, "/file-cell-count", sep="")
feat_file_cell_prop_dir = paste(feat_dir, "/file-cell-prop", sep="")



## libraries
source("~/projects/IMPC/code/_funcAlice.R")
libr("flowCore")
libr("flowType")
libr("foreach")
libr("doMC")

## cores
no_cores = detectCores()-1
registerDoMC(no_cores)

## options
options(stringsAsFactors=F)
options(device="cairo")
countThres = 0 #delete columns/rows where all values equal or below countThres
levelThres = 8 #Inf if delete no layers; >levelcutoff are deleted











start = Sys.time()

## prepare directories
ftFile_dir = sort( dir(ft_dir, pattern=".Rda", all.files=T, full.names=T, recursive=T) )
ftGT = folderNames(ftFile_dir)
ftFileNames = fileNames(ftFile_dir, "Rda")

## get meta data on files for meta_file
meta_filetemp = data.frame(read.csv(csv_dir))
colnames(meta_filetemp) = c("fileName", "tube", "specimen", "aml") #rename columns

## get meta data on phenotypes (cell populations) for meta_cell
ft = get(load(ftFile_dir[1]))
markers = ft@MarkerNames
save(markers, file=paste0(markers_dir,".Rdata"))

meta_cell = getPhen(rownames(ft@MFIs))

## load flowtype files for feat_file_cell_count & fill in meta_file
feat_file_cell_count = NULL
meta_file = NULL

result = foreach (i=1:length(ftFile_dir), .combine="rbind") %dopar% {
  ft = get(load(ftFile_dir[i]))
  print(ft@MarkerNames)
  npheno = length(ft@CellFreqs)
  # feat_file_cell_count = rbind(feat_file_cell_count, ft@CellFreqs)
  cat("\n", i, "Loaded file: ", ftGT[i], "/", ftFileNames[i], "; # of phenotypes x samples: ", npheno, " x ", length(ft), "", sep="")
  
  fn = ftFileNames[i]
  fn = strsplit(fn, "\\.")[[1]][1]
  fn = strsplit(fn, "S")[[1]]
  tube = substr(fn[1],2,nchar(fn[1]))
  specimen = gsub("FT","",fn[2])
  # meta_file = rbind(meta_file, meta_filetemp[which(as.numeric(meta_filetemp$specimen)==as.numeric(specimen) & as.numeric(meta_filetemp$tube)==as.numeric(tube)),])
  return(c(meta_filetemp[which(as.numeric(meta_filetemp$specimen)==as.numeric(specimen) & as.numeric(meta_filetemp$tube)==as.numeric(tube)),], ft@CellFreqs))
  
  rm(ft) #save memory
}

#collate feat_file_cell_count
feat_file_cell_count = as.matrix(apply(result[,(ncol(meta_filetemp)+1):ncol(result)], 2, as.numeric))
sm = result[,1:(ncol(meta_filetemp))]

#order meta_cell and feat_file_cell_count
pheno_order = order(meta_cell$phenocode)
meta_cell = meta_cell[pheno_order,]
feat_file_cell_count = feat_file_cell_count[,pheno_order]

#make meta_file
rownames(feat_file_cell_count) = ftFileNames
meta_file = as.data.frame(ftFileNames)
meta_file$tube = unlist(sm[,2])
meta_file$specimen = unlist(sm[,3])
meta_file$aml = unlist(sm[,4])
colnames(meta_file) = colnames(meta_filetemp)

rownames(feat_file_cell_count) = meta_file[,1] = ftFileNames
colnames(feat_file_cell_count) = meta_cell$phenotype
feat_file_cell_prop = feat_file_cell_count/feat_file_cell_count[,1]
dimnames(feat_file_cell_prop) = dimnames(feat_file_cell_count)



## trim matrix/meta_cell ------------------------------
rowIndex = apply(feat_file_cell_count, 1, function(x) any(x > countThres)) #delete rows of all 0 or too little count
colIndex1 = apply(feat_file_cell_count, 2, function(x) any(x > countThres)) #delete cols of all 0 or too little count
colIndex2 = meta_cell$phenolevel <= levelThres #delete cols of too high level
colIndex = colIndex1 & colIndex2

feat_file_cell_count <- feat_file_cell_count[rowIndex,colIndex]
feat_file_cell_prop <- feat_file_cell_prop[rowIndex,colIndex]
meta_cell <- meta_cell[colIndex,]
meta_file <- meta_file[rowIndex,]



#save
save(meta_cell, file=paste0(meta_cell_dir,".Rdata")); write.csv(meta_cell, file=paste0(meta_cell_dir,".csv"), row.names=F)
save(meta_file, file=paste0(meta_file_dir,".Rdata")); write.csv(meta_file, file=paste0(meta_file_dir,".csv"), row.names=F)
save(feat_file_cell_count, file=paste0(feat_file_cell_count_dir,".Rdata"))
write.csv(feat_file_cell_count, file=paste0(feat_file_cell_count_dir,".csv"), row.names=T)
save(feat_file_cell_prop, file=paste0(feat_file_cell_prop_dir,".Rdata"))
write.csv(feat_file_cell_prop, file=paste0(feat_file_cell_prop_dir,".csv"), row.names=T)

TimeOutput(start)






