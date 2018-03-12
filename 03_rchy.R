## Input: original/trimmed features --> Output: trimmed features + rchyoptimyx additional nodes
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
phenoMeta_dir = paste(result_dir, "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="")
matrix_dir = paste(result_dir, "/matrix", sep="")

#Output
rchy_dir = paste(result_dir, "/rchy", sep="")


#Libraries/Functions
library(Matrix)
library(stringr)
library(foreach)
library(doMC)
library(RchyOptimyx)
source("~/projects/IMPC/code/_funcAlice.R")

#Setup Cores
no_cores = detectCores()-1
registerDoMC(no_cores)








#Options for script
cellCountThres = c(1200) #insignificant if count under

kpaths = 5

matrix_count = c("CountAdj")
matrix_TRIM = c("matrixPvalTRIM_CountAdj.Rdata") #matrix with 0s for insignificant nodes
matrix_rchy = c("matrixPval_CountAdj.Rdata") #matrix as input into rchy

weight_matrix = c("MaxCountAdj_CountAdj") #assume by cell pop

matrix_type = c("PvalTRIM_CountAdj")
matrix_all_type = c("Pval_CountAdj")
matrix_edge_type = c("Parent_contrib_CountAdj") #plot only


#Prepare data
m0 = get(load(paste0(matrix_dir, matrix_count,".Rdata")))
phenoMeta = get(load(phenoMeta_dir))

k0 = c(4)#1,4,max(phenoMeta$phenolevel)) # how many markers to consider i.e. k=max(phenolevel) only











start = Sys.time()

for (mti in 1:length(matrix_type)) {
  tryCatch({
    mcp = matrix_type[mti]
    map = matrix_all_type[mti]
    mep = matrix_edge_type[mti]
    cat("\n", mcp, " ",sep="")
    start2 = Sys.time()
    
    #start a list of phenotypes included in each distance matrix calculation such that none are repeated
    leavePhenotype = list()
    doneAll = F
    
    #load feature matrix
    mresult = Loadintermatrices(c(paste0(matrix_dir, mcp,".Rdata"),paste0(matrix_dir, map,".Rdata"),paste0(matrix_dir, mcp,".Rdata")))
    mml0 = mresult$mml
    mmlname = names(mml0)
    pt = mpi = mresult$pt
    gt = mresult$gt
    
    #get to-delete low count phenotype indices; CountAdj should be first one
    for (countThres in cellCountThres) { cat("\ncountThres: ",countThres," >",sep="")
      
      #get to-delete high no of marker phenotypes
      for (k in k0) { cat(" level",k," ",sep="")
        
        #trim feature matrix
        mmlresult = trimMatrices(mml0,m0,pt,gt,phenoMeta,leavePhenotype,doneAll, countThres,k)
        if (is.null(mmlresult)) return(NULL)
        m = mmlresult$mml[[1]]
        ma = mmlresult$mml[[2]]
        e = mmlresult$mml[[3]]
        pm = mmlresult$pm
        leavePhenotype = mmlresult$leavePhenotype
        doneAll = mmlresult$doneAll
        
        if (is.null(m)) return(NULL)
        if (is.null(dim(m))) {
          a = Reduce('cbind',m); if (all(a==0)) next
        } else { if (all(m==0)) next }
        
        fname = paste0(rchy_dir,"/dir_k=",kpaths,"_layer", str_pad(k, 2, pad = "0"), "_countThres-", countThres); suppressWarnings(dir.create(fname))
        fname0 = paste0(rchy_dir,"/all_k=",kpaths,"_layer", str_pad(k, 2, pad = "0"), "_countThres-", countThres); suppressWarnings(dir.create(fname0))
        
        a = foreach(ii=1:nrow(m)) %dopar% {
          i = rownames(m)[ii]
          phenocodes = phenoMeta$phenocode[match(colnames(ma),phenoMeta$phenotype)]
          
          
          # either change pvalues
          rchy0 = NULL
          if (sum(m[ii,]!=0)>0) {
            phenoscore = abs(ma[i,])
            startpheno = phenoMeta$phenocode[match(colnames(m)[m[ii,]!=0],phenoMeta$phenotype)]
            startpheno = setdiff(startpheno,paste0(rep(0,nchar(phenoMeta$phenocode[1]))))
            rchy0 = lapply(startpheno, function(sp) {
              rr = RchyOptimyx(pheno.codes=phenocodes, phenotypeScores=phenoscore, startPhenotype=sp, pathCount=kpaths, trimPaths=F)
              rr@nodes[2,] = ma[ii,match(rr@nodes[1,],phenoMeta$phenocode[match(colnames(ma),phenoMeta$phenotype)])]
              rownames(rr@nodes) = c("phenocode",matrix_type[mti],"colour")
              return(rr)
            })
            names(rchy0) = startpheno
            # rchy = rchy1[[1]]; if (length(rchy1)>1) for (ri in 2:length(rchy1)) { rchy = merge(rchy,rchy1[[ri]]) }
          }
          
          rchy0_m = NULL
          if (!is.null(rchy0)) {
            rchy0_m = rchy0[[1]]
            if (length(rchy0)>1) {
              for (ri in 2:length(rchy0)) {
                rchy0_m = merge(rchy0_m,rchy0[[ri]])
              }
            }
          }
          
          
          
          rchy = NULL
          # positive change pvalues
          rchy1 = NULL
          if (sum(m[ii,]>0)>0) {
            phenoscore = ma[ii,]
            if (min(phenoscore)<0) { phenoscore = phenoscore+min(-ma[ii,]);  }
            startpheno = phenoMeta$phenocode[match(colnames(m)[m[i,]>0],phenoMeta$phenotype)]
            startpheno = setdiff(startpheno,paste0(rep(0,nchar(phenoMeta$phenocode[1]))))
            rchy1 = lapply(startpheno, function(sp) {
              rr = RchyOptimyx(pheno.codes=phenocodes, phenotypeScores=phenoscore, startPhenotype=sp, pathCount=kpaths, trimPaths=F)
              rr@nodes[2,] = ma[ii,match(rr@nodes[1,],phenoMeta$phenocode[match(colnames(ma),phenoMeta$phenotype)])]
              rownames(rr@nodes) = c("phenocode",matrix_type[mti],"colour")
              return(rr)
            })
            names(rchy1) = startpheno
            # rchy = rchy1[[1]]; if (length(rchy1)>1) for (ri in 2:length(rchy1)) { rchy = merge(rchy,rchy1[[ri]]) }
          }
          # negative change pvalues
          rchy2 = NULL
          if (sum(m[ii,]<0)>0) {
            phenoscore = -ma[ii,]
            if (min(phenoscore)<0) phenoscore = phenoscore+min(-ma[ii,])
            startpheno = phenoMeta$phenocode[match(colnames(m)[m[i,]<0],phenoMeta$phenotype)]
            startpheno = setdiff(startpheno,paste0(rep(0,nchar(phenoMeta$phenocode[1])),collapse=""))
            rchy2 = lapply(startpheno, function(sp) {
              rr = RchyOptimyx(pheno.codes=phenocodes, phenotypeScores=phenoscore, startPhenotype=sp, pathCount=kpaths, trimPaths=F)
              rr@nodes[2,] = ma[ii,match(rr@nodes[1,],phenoMeta$phenocode[match(colnames(ma),phenoMeta$phenotype)])]
              rownames(rr@nodes) = c("phenocode",matrix_type[mti],"colour")
              return(rr)
            })
            names(rchy2) = startpheno
            # if (!sum(m[i,]>0)>0) rchy = rchy2[[1]]
            # rchy = merge(rchy,rchy2[[1]]); if (length(rchy2)>1) for (ri in 2:length(rchy2)) { rchy = merge(rchy,rchy2[[ri]]) }
          }
          rchy = append(rchy1,rchy2)
          
          rchy_m = NULL
          if (!is.null(rchy)) {
            rchy_m = rchy[[1]]
            if (length(rchy)>1) {
              for (ri in 2:length(rchy)) {
                rchy_m = merge(rchy_m,rchy[[ri]])
              }
            }
          }
          
          if (!is.null(rchy0)) save(rchy0,file=paste0(fname0,"/",i,".Rdata"))
          if (!is.null(rchy)) save(rchy,file=paste0(fname,"/",i,".Rdata"))
          if (!is.null(rchy0_m)) save(rchy0,file=paste0(fname0,"/",i,"_merged.Rdata"))
          if (!is.null(rchy_m)) save(rchy,file=paste0(fname,"/",i,"_merged.Rdata"))
          
        }
        
      } #layer
    } #countThres
    
    TimeOutput(start2)
  }, error = function(err) { cat(paste("ERROR:  ",err)) })
}

TimeOutput(start)












#delete distance matrices with all 0
distmfile = list.files(dist_dir, recursive=F, full.names=T, pattern=".Rdata$")
a = foreach (i = 1:length(distmfile),.combine="c") %dopar% {
  a = F
  d = get(load(distmfile[i]))
  if (any(is.na(as.matrix(d)))) { rm(d); return(T) }
  if (all(as.matrix(d)==as.matrix(d)[1,1])) { rm(d); return(T) }
  rm(d)
  return(a)
}
file.remove(distmfile[a])










