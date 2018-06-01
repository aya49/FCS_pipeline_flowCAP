## Input: original features --> Output: prints whether there are irregular values in features
# aya43@sfu.ca 20170316

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

#Libraries/Functions
libr(Matrix)
libr(foreach)
libr(doMC)
source("~/projects/IMPC/code/_funcAlice.R")

#Setup Cores
no_cores = detectCores()-1
registerDoMC(no_cores)











#Options for script
sampletimes = 5



#Prepare data
matrix_type = c("CountAdj", "Parent_entropy", "Child_entropy", 
                paste0("Child_entropyTRIM_CountAdj_",c(1:sampletimes)), "Child_entropyTRIM_CountAdj", "Child_entropyTRIM_Prop", 
                paste0("Parent_entropyTRIM_CountAdj_",c(1:sampletimes)), "Parent_entropyTRIM_CountAdj", "Parent_entropyTRIM_Prop", 
                paste0("LogFoldTRIM_CountAdj_",c(1:sampletimes)), "LogFoldTRIM_CountAdj", "LogFoldTRIM_Prop", 
                paste0("PvalTRIM_CountAdj_",c(1:sampletimes)), "PvalTRIM_CountAdj", "PvalTRIM_Prop",
                paste0("LogFold_CountAdj_",c(1:sampletimes)), "LogFold_CountAdj", "LogFold_Prop", 
                paste0("Pval_CountAdj_",c(1:sampletimes)), "Pval_CountAdj", "Pval_Prop",
                "Child_pnratio", "Child_prop",
                paste0("Parent_effortTRIM_CountAdj_",c(1:sampletimes)), "Parent_effortTRIM_CountAdj", "Parent_effortTRIM_Prop", 
                paste0("Parent_contribTRIM_CountAdj_",c(1:sampletimes)), "Parent_contribTRIM_CountAdj", "Parent_contribTRIM_Prop", 
                paste0("Child_pnratioTRIM_CountAdj_",c(1:sampletimes)), "Child_pnratioTRIM_CountAdj", "Child_pnratioTRIM_Prop", 
                paste0("Child_propTRIM_CountAdj_",c(1:sampletimes)), "Child_propTRIM_CountAdj", "Child_propTRIM_Prop",
                paste0("Parent_effort_CountAdj_",c(1:sampletimes)), "Parent_effort_CountAdj", "Parent_effort_Prop", 
                paste0("Parent_contrib_CountAdj_",c(1:sampletimes)), "Parent_contrib_CountAdj", "Parent_contrib_Prop")










start = Sys.time()

for (mcp in sort(matrix_type)) {
  cat("\n", mcp, " ",sep="")
  start2 = Sys.time()
  
  #load feature matrix
  if (!file.exists(paste0(matrix_dir, mcp,".Rdata"))) {cat("doesn't exist"); next}
  mm = get(load(paste0(matrix_dir, mcp,".Rdata")))
  
  ## check for irregular values in matrix
  m = mm
  naed = F
  Infed = F
  mcpneg = F
  if (!is.null(dim(m))) {
    if (sum(m[!is.na(m)]<0)>1) mcpneg = T
    cat(dim(m))
    cat(" ",checkm(m," "))
    if(mcpneg) cat(" neg ")
    if (is.null(rownames(m))) cat(" norownames; ")
  } else {
    cat(dim(m[[length(m)]])," ")
    nainfneg = foreach (mind = 1:length(m),.combine="rbind") %dopar% {
      nai = F
      infi = F
      negi = F
      checkmed = checkm(m[[mind]],paste(" ",mind))
      if (grepl("_NA",checkmed,ignore.case=T)) {nai=T}
      if (grepl("_Inf",checkmed,ignore.case=T)) {infi=T}
      if (sum(m[[mind]][!is.na(m[[mind]])]<0)>1) { negi=T }
      return(c(nai,infi,negi))
    }
    if (sum(nainfneg[,1])>0) cat (min(which(nainfneg[,1]==T)),".na ",sep="")
    if (sum(nainfneg[,2])>0) cat (min(which(nainfneg[,2]==T)),".inf ",sep="")
    if (sum(nainfneg[,3])>0) cat (min(which(nainfneg[,3]==T)),".neg ",sep="")
    if (is.null(rownames(m[[1]]))) cat(" norownames; ")
  }
}

TimeOutput(start)

