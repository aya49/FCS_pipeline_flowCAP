## Input: original features --> Output: distance matrices
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
phenoMeta_dir = paste(result_dir,  "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir,  "/sampleMeta.Rdata", sep="")
rw_dir = paste(result_dir,  "/rw", sep="")

#Output
dist_dir = paste(result_dir,  "/dist", sep=""); dir.create (dist_dir, showWarnings = F)
dist_type_dir = NULL
phenoChild_dir = paste(result_dir,  "/phenoChild.Rdata",sep="")
phenoChild_ind_dir = paste(result_dir,  "/phenoChild_ind.Rdata",sep="")
phenoChildpn_dir = paste(result_dir,  "/phenoChildpn.Rdata",sep="")
phenoChildpn_ind_dir = paste(result_dir,  "/phenoChildpn_ind.Rdata",sep="")

#Libraries/Functions
library(stringr)
library(colorspace)
library(vegan) # library(proxy)
library(foreach)
library(doMC)
library(lubridate) #if there are date variables
library(kernlab)
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")

#Setup Cores
no_cores = detectCores()-1
registerDoMC(no_cores)








#Options for script
# sampletimes = 5
overwrite = F #overwrite distances?
writecsv = F

cellCountThres = c(200) #insignificant if count under

dis = c("manhattan", "euclidean") #distances metrics to use #, "binomial", "cao", "jaccard","mahalanobis", "gower","altGower","morisita","horn", "manhattan", "maximum", "binary", "minkowski", "raup", "mountford"
disnoavg = c("euclidean", "manhattan", "canberra", "bray", "kulczynski", "morisita", "horn", "binomial")  #dis measures that don't average over number of features
disindist = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski") #use dist function when possible, instead of vegdist
disnoneg = c("canberra") #dis measures that can't handle negative values
# disinkernel = c("rbf","poly","vanilla","tanh","laplace","bessel","anova","spline") #kernels

normalize = c("none","cellpop", "layer") # by none (for child matrices only), cell pop, layer

matrix_count = c("CountAdj")
























start1 = Sys.time()


#Prepare data
matrix_type = gsub(".Rdata","",list.files(path=rw_dir))

# matrix_type = c("CountAdj", "Parent_entropy", "Child_entropy", 
#                 paste0("Child_entropyTRIM_CountAdj_",c(1:sampletimes)), "Child_entropyTRIM_CountAdj", "Child_entropyTRIM_Prop", 
#                 paste0("Parent_entropyTRIM_CountAdj_",c(1:sampletimes)), "Parent_entropyTRIM_CountAdj", "Parent_entropyTRIM_Prop", 
#                 paste0("LogFoldTRIM_CountAdj_",c(1:sampletimes)), "LogFoldTRIM_CountAdj", "LogFoldTRIM_Prop", 
#                 paste0("PvalTRIM_CountAdj_",c(1:sampletimes)), "PvalTRIM_CountAdj", "PvalTRIM_Prop",
#                 paste0("LogFold_CountAdj_",c(1:sampletimes)), "LogFold_CountAdj", "LogFold_Prop", 
#                 paste0("Pval_CountAdj_",c(1:sampletimes)), "Pval_CountAdj", "Pval_Prop",
#                 "Child_pnratio", "Child_prop",
#                 paste0("Parent_effortTRIM_CountAdj_",c(1:sampletimes)), "Parent_effortTRIM_CountAdj", "Parent_effortTRIM_Prop", 
#                 paste0("Parent_contribTRIM_CountAdj_",c(1:sampletimes)), "Parent_contribTRIM_CountAdj", "Parent_contribTRIM_Prop", 
#                 paste0("Child_pnratioTRIM_CountAdj_",c(1:sampletimes)), "Child_pnratioTRIM_CountAdj", "Child_pnratioTRIM_Prop", 
#                 paste0("Child_propTRIM_CountAdj_",c(1:sampletimes)), "Child_propTRIM_CountAdj", "Child_propTRIM_Prop",
#                 paste0("Parent_effort_CountAdj_",c(1:sampletimes)), "Parent_effort_CountAdj", "Parent_effort_Prop", 
#                 paste0("Parent_contrib_CountAdj_",c(1:sampletimes)), "Parent_contrib_CountAdj", "Parent_contrib_Prop")







sampleMeta = get(load(sampleMeta_dir))

















#load different features matrix and calculate distance matrix
# for (mcp in matrix_type_) {
a = foreach(mcp=matrix_type) %dopar% {
  tryCatch({
    cat("\n", mcp, " ",sep="")
    start2 = Sys.time()
    
    #start a list of phenotypes included in each distance matrix calculation such that none are repeated
    leavePhenotype = list()
    doneAll = F
    
    #load feature matrix
    m0 = get(load(paste0(rw_dir,"/", mcp,".Rdata")))
    
    phenotype = sapply(str_split(colnames(m0),"_"), function(x) x[[2]])
    colsplitlen = cell_type_layers(phenotype); k0 = 0; if (!is.null(colsplitlen)) k0 = c(1,2,4,max(colsplitlen)) #1,4, # how many layers to consider i.e. k=max(phenolevel) only
    
    
    #get to-delete high no of marker phenotypes
    for (k in k0) { cat(" level",k," ",sep="")
      
      #trim feature matrix
      m = m0[,colsplitlen<=k]
      
      if (is.null(m)) next
      if (all(m==0)) next
      
      #for every distance type
      a = match(disnoneg,dis)
      loop.ind = 1:length(dis); if (sum(!is.na(a))>0) loop.ind = loop.ind[-a]
      # foreach(i=loop.ind) %dopar% { #for each phenotype
      for (i in loop.ind) {
        cat(", ", length(dis)-i+1, ":", dis[i], " ", sep="")
        
        #assume after splitting dist filename by "_", distance is second element
        dname = paste(dist_dir, "/", mcp, "_", dis[i], "_layer", str_pad(k, 2, pad = "0"), sep = "" )
        
        ## calculate distances
        dname1 = paste0(dname, "_cellpop.Rdata")
        if (overwrite | !file.exists(dname1)) {
          if (dis[i]%in%disindist) { d = dist(m, method=dis[i])
          } else { d = vegdist(m, method=dis[i]) }
          
          save(d, file=dname1)
          if (writecsv) write.csv(as.matrix(d), file=gsub(".Rdata",".csv",checkm(d,dname1)))
        }
        
      } #dis
    } #layer
    
    TimeOutput(start2)
  }, error = function(err) { cat(paste("ERROR:  ",err)) })
}
TimeOutput(start1)





