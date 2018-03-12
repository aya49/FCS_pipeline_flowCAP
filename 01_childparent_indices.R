## Input: phenotype list --> Output: lists of children/parent for each non-leaf/root node
#aya43@sfu.ca 20170131

#Directory
root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

#Options
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/phenoMeta", sep="")

#Output
phenoChild_dir = paste(result_dir, "/phenoChild",sep="") #specifies a phenotypes children
phenoChild_ind_dir = paste(result_dir, "/phenoChild_ind",sep="")
phenoChildpn_dir = paste(result_dir, "/phenoChildpn",sep="") #specifies a phenotypes children and splits them into +/- (only for when both -/+ exists)
phenoChildpn_ind_dir = paste(result_dir, "/phenoChildpn_ind",sep="")
phenoParent_dir = paste(result_dir, "/phenoParent",sep="") #specifies a phenotypes parents
phenoParent_ind_dir = paste(result_dir, "/phenoParent_ind",sep="")
phenoParentpn_dir = paste(result_dir, "/phenoParentpn",sep="") #specifies a phenotypes parents and splits them into +/- (only for when both -/+ exists)
phenoParentpn_ind_dir = paste(result_dir, "/phenoParentpn_ind",sep="")

#Libraries/Functions
library(foreach)
library(doMC)
source("~/projects/IMPC/code/_funcAlice.R")

#Setup Cores
no_cores = detectCores()-1
registerDoMC(no_cores)










#Prepare Data
phenoMeta = get(load(paste0(phenoMeta_dir, ".Rdata")))
# colnames(phenoMeta) = c("phenotype","phenocode","phenolevel")
# save(phenoMeta,file=paste0(phenoMeta_dir, pcp,".Rdata"))










start = Sys.time()

start1 = Sys.time()

cat("\ncreating children indices ")
pc = getphenoChild(phenoMeta, no_cores=no_cores)

phenoChild = pc$phenoChild
phenoChild_ind = pc$phenoChild_ind
phenoChildpn = pc$phenoChildpn
phenoChildpn_ind = pc$phenoChildpn_ind

#save
save(phenoChild, file=paste0(phenoChild_dir, ".Rdata"))
save(phenoChild_ind, file=paste0(phenoChild_ind_dir, ".Rdata"))
save(phenoChildpn, file=paste0(phenoChildpn_dir, ".Rdata"))
save(phenoChildpn_ind, file=paste0(phenoChildpn_ind_dir, ".Rdata"))

TimeOutput(start1)










cat("\ncreating parents indices ")
pp = getphenoParent(phenoMeta, phenoChildpn, phenoChildpn_ind, no_cores)

phenoParent = pp$phenoParent
phenoParent_ind = pp$phenoParent_ind
phenoParentpn = pp$phenoParentpn
phenoParentpn_ind = pp$phenoParentpn_ind

#save
save(phenoParent, file=paste0(phenoParent_dir, ".Rdata"))
save(phenoParent_ind, file=paste0(phenoParent_ind_dir, ".Rdata"))
save(phenoParentpn, file=paste0(phenoParentpn_dir, ".Rdata")) 
save(phenoParentpn_ind, file=paste0(phenoParentpn_ind_dir, ".Rdata"))

TimeOutput(start1)

TimeOutput(start)




