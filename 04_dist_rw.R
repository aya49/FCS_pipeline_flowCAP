## Input: random walk paths --> min jaccard distance
# aya43@sfu.ca 20180521

#Directory
root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

## input directories
meta_dir = paste0(result_dir,"/meta") # meta files directory
meta_file_dir = paste(meta_dir, "/file", sep="") #meta for rows (sample)
meta_cell_dir = paste(meta_dir, "/cell", sep="") #meta for rows (sample)
feat_dir = paste(result_dir, "/feat", sep="") #feature files directory
meta_cell_child_names_dir = paste(meta_dir, "/cell_child_names",sep="") #specifies a phenotypes children

## output directories
dist_dir = paste0(result_dir,"/dist"); suppressWarnings(dir.create (dist_dir))

## libraries
library(Matrix)
library(foreach)
library(doMC)
library(stringr)
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")

#Setup Cores
no_cores = 15#detectCores()-6
registerDoMC(no_cores)








## options for script
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

overwrite = F
writecsv = T

cellCountThres = c(200) #(needed if sample x cell pop matrices are used) insignificant if count under
good_sample = c(3)
good_count = c(3)
id_col = "fileName"
target_col = "aml"
order_cols = c("aml","tube")
split_col = "tube"

#data paths
feat_type_paths = list.dirs(path=feat_dir, full.names=F)
feat_type_paths = feat_type_paths[grepl("rw",feat_type_paths)]

feat_count = c("file-cell-countAdj") #(needed if sample x cell pop matrices are used) count matrix, to get rid of low cell count cells
























start1 = Sys.time()

meta_file = get(load(paste0(meta_file_dir,".Rdata")))

#load different features matrix and calculate distance matrix
# for (feat_type_path in feat_type_paths_) {
a = foreach(feat_type_path=feat_type_paths) %dopar% {
  tryCatch({
    cat("\n", feat_type_path, " ",sep="")
    start2 = Sys.time()
    
    feat_type_path_files = list.files(path=paste0(feat_dir,"/",feat_type_path), full.names=F)
    feat_type_path_filenames = gsub(".Rdata","",feat_type_path_files)

    #load feature matrix
    mp = random_paths_all_files = foreach(x=feat_type_path_files) %dopar% { 
      return( get(load(paste0(feat_dir,"/",feat_type_path,"/",x))) )
    }
    path_names = sapply(mp, function(x) names(x))
    path_sum = sapply(mp, function(x) sum(x))
    path_ratio = path_sum/min(path_sum)
    random_paths_all_files1 = foreach(x=random_paths_all_files, .combine=rbind) %dopar% {
      x0 = x[match(path_names, names(x))]
      x0[is.na(x0)] = 0
      return(x0*path_ratio)
    }
    # random_paths_all_files1 = Matrix(random_paths_all_files1,sparse=T)
    # rownames(random_paths_all_files1) = feat_type_path_files
    # colnames(random_paths_all_files1) = path_names
    # mp = random_paths_all_files1
    
    # dist = matrix(0,nrow=nrow(mp),ncol=nrow(mp))
    dist = foreach (i = 1:(length(mp)-1), .combine=rbind) %dopar% {
      dtemp = rep(0,length(mp))
      pathsi = mp[[i]]
      for (j in (i+1):length(mp)) {
        pathsj = mp[[j]]
        dtemp[j] = sum(pmin(pathsi,pathsj,na.rm=T))
        # dist[i,j] = dist[j,i] = sum(pmin(pathsi,pathsj,na.rm=T))
      }
    }
    dist = as.matrix(as.dist(rbind(rep(0,nrow(mp)),dist)))
    colnames(dist) = rownames(dist) = feat_type_path_files
    
    save(dist, dist_dir,"/",feat_type_path, "_layer", str_pad(k, 2, pad = "0"), "_countThres-", countThres,
         "_dist-jaccardmin.Rdata")

    TimeOutput(start2)
  }, error = function(err) { cat(paste("ERROR:  ",err)) })
}
TimeOutput(start1)





