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

## output directories
stat_dir = paste(result_dir, "/stat", sep=""); suppressWarnings(dir.create (stat_dir))
dist_dir = paste0(result_dir,"/dist"); suppressWarnings(dir.create (dist_dir))

## libraries
libr(Matrix)
libr(foreach)
libr(doMC)
libr(plyr)
libr(stringr)
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")

#Setup Cores
no_cores = 14#detectCores()-6
registerDoMC(no_cores)








## options for script
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

overwrite = F
writecsv = T

minwalks = 5 #if number of walks on a path is less than minwalks, set to 0

cellCountThres = c(200) #(needed if sample x cell pop matrices are used) insignificant if count under
good_sample = c(3)
good_count = c(3)
id_col = "fileName"
target_col = "aml"
order_cols = c("aml","tube")
split_col = "tube"

#data paths
feat_type_paths = list.dirs(path=feat_dir, full.names=F)
feat_type_paths = feat_type_paths[grepl("path-rw",feat_type_paths)]

feat_count = c("file-cell-countAdj") #(needed if sample x cell pop matrices are used) count matrix, to get rid of low cell count cells
























start = Sys.time()

meta_file = get(load(paste0(meta_file_dir,".Rdata")))

#load different features matrix and calculate distance matrix
# for (feat_type_path in feat_type_paths_) {
for (feat_type_path in feat_type_paths) {
  tryCatch({
    cat("\n", feat_type_path, " ",sep="")
    start2 = Sys.time()
    
    feat_type_path_files = list.files(path=paste0(feat_dir,"/",feat_type_path), full.names=F)
    feat_type_path_filenames = gsub(".Rdata","",feat_type_path_files)
    
    #load feature matrix
    # mp = random_paths_all_files = foreach(x=feat_type_path_files) %dopar% { 
    #   a = get(load(paste0(feat_dir,"/",feat_type_path,"/",x)))
    #   an = names(a)
    #   as = sum(a)
    # }
    
    # start1 = Sys.time()
    # mp0 = lapply(feat_type_path_files, function(x)
    #   # data.frame(as.list(get(load(paste0(feat_dir,"/",feat_type_path,"/",x))))) ) #30sec for flowcap
    #   get(load(paste0(feat_dir,"/",feat_type_path,"/",x))) ) #30sec for flowcap
    # TimeOutput(start1)
    
    if (file.exists(paste0(feat_dir,"/",feat_type_path,".Rdata"))) {
      mp2 = get(load(paste0(feat_dir,"/",feat_type_path,".Rdata")))
    } else {
      
      start1 = Sys.time()
      mp0 = foreach(x=feat_type_path_files) %dopar% {
        # data.frame(as.list(get(load(paste0(feat_dir,"/",feat_type_path,"/",x))))) ) #30sec for flowcap
        return( get(load(paste0(feat_dir,"/",feat_type_path,"/",x))) ) } #30sec for flowcap
      TimeOutput(start1)
      
      # start1 = Sys.time()
      # mp0 = foreach(x=feat_type_path_files) %dopar% {
      #   return( data.frame(as.list(get(load(paste0(feat_dir,"/",feat_type_path,"/",x))))) ) #8 min 14 cores for flowcap
      # }
      # TimeOutput(start1)
      
      path_sum = sapply(mp0, function(x) sum(x))
      
      path_ratio = path_sum/min(path_sum)
      mp = lapply(1:length(mp0), function(xi) mp0[[xi]]*path_ratio[xi])
      if (minwalks>0) mp = lapply(mp, function(x) x[x>=minwalks])
      
      start1 = Sys.time()
      path_names = lapply(mp0, function(x) names(x))
      path_name = unique(unlist(path_names))
      # path_order = lapply(path_names, function(x) match(path_name,path_names))
      
      # mp2 = Matrix(0,ncol=length(path_names),nrow=length(mp), 
      #              dimnames=list(gsub(".Rdata","",feat_type_path_files),path_name), sparse=T)
      
      loop.inds = loopInd(1:length(mp),no_cores)
      mp2 = foreach (ii = loop.inds, .combine=rbind) %dopar% {
        mtemp = Matrix(0,nrow=length(ii),ncol=length(path_name),sparse=T)
        for(i in 1:length(ii)) {
          mpi = mp[[ii[i]]]
          mpiorder = match(path_name,names(mpi))
          mpi2 = mpi[mpiorder]
          mpi2[is.na(mpiorder)] = 0
          # mtempi = c.sparseVector(mpi2)
          mtemp[i,] = c.sparseVector(mpi2)
          # return(mtempi)
        }
        # mtemp = lapply(mtempt, as, "sparseMatrix")
        # mtemp = Reduce(cbind, mtemp)
        return(mtemp)
      }
      # mp1 = foreach(x=mp) %dopar% {data.frame(as.list(x))}
      # mp2 = Reduce(rbind.fill, mp1)
      # mp2[is.na(mp2)] = 0
      # mp2 = Matrix(mp2, sparse=T)
      
      save(mp2, file=paste0(feat_dir,"/",feat_type_path,".Rdata"))
      
      TimeOutput(start1)
      
      # random_paths_all_files1 = Matrix(random_paths_all_files1,sparse=T)
      # rownames(random_paths_all_files1) = feat_type_path_files
      # colnames(random_paths_all_files1) = path_names
      # mp = random_paths_all_files1
      
    }
    
    png(paste0(stat_dir,"/",feat_type_path,"_pathfreqdensity.png"), width=900, height=400)
    par(mfrow=c(1,2))
    plot(density(mp2),main="file vs path feature (val = freq of walks per path)",
         xlab="frequency of walks on a path", ylab="density (smoothed histogram) of paths with that frequency")
    hist(mp2,main="same plot; zoomed in", xlab="", ylab="histogram",
         breaks=1000,xlim=c(0,3000000), ylim=c(0,500))
    graphics.off()
    
    
    start2 = Sys.time()
    
    dname = paste0(dist_dir,"/",feat_type_path, "_layer", str_pad(0, 2, pad = "0"), "_countThres-", 200,
                   "_dist-jaccardmin")
    dir.create(dname, showWarnings=F)
    
    # dist = matrix(0,nrow=nrow(mp),ncol=nrow(mp))
    # dist = foreach (i = 1:(nrow(mp2)-1), .combine=rbind) %dopar% {
    a = foreach (i = 1:(nrow(mp2)-1)) %dopar% {
      dtemp = rep(0,nrow(mp2))
      pathsi = mp2[i,]
      for (j in (i+1):nrow(mp2)) {
        pathsj = mp2[j,]
        dtemp[j] = sum(pmin(pathsi,pathsj,na.rm=T))
        # dist[i,j] = dist[j,i] = sum(pmin(pathsi,pathsj,na.rm=T))
      }
      save(dtemp, file=paste0(dname,"/",i,".Rdata"))
      # return(dtemp)
    }
    dist = lapply(list.files(dname,full.names=T), function(x) get(load(x)))
    dist = Reduce("rbind",dist)
    
    dist1 = as.matrix(as.dist(rbind(rep(0,nrow(mp2)),dist)))
    colnames(dist1) = rownames(dist1) = gsub(".Rdata","",feat_type_path_files)
    
    dist2 = -(dist1-max(dist1))
    
    save(dist2, file=paste0(dname,".Rdata"))
    
    png(paste0(stat_dir,"/",feat_type_path,"_rwminjaccarddist.png"), width=450, height=400)
    # plot(density(dist),main="kernel density vs random walk min jaccard distances")
    hist(dist2,main="file vs file distance (values = -(max - all jaccard min sim)", breaks=200,
         xlab="distance", ylab="histogram of each distance value")
    graphics.off()
    
    
    TimeOutput(start2)
  }, error = function(err) { cat(paste("ERROR:  ",err)) })
}
TimeOutput(start)





