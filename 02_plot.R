# Plot: plot samples & divide dates up for p value calculation
# aya43@sfu.ca 20161220

root = "~/projects/flowCAP-II"
result_dir = "result"
setwd(root)


options(stringsAsFactors=FALSE)
options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="")
matrix_dir = paste(result_dir, "/matrix", sep="")
matrix_type = c("CountAdj","Prop")
# matrix_type = c("Child_entropyTRIM_CountAdj", "Child_entropyTRIM_Prop", "Parent_entropyTRIM_CountAdj", "Parent_entropyTRIM_Prop", 
#                 "LogFoldTRIM_CountAdj", "LogFoldTRIM_Prop", "PvalTRIM_CountAdj", "PvalTRIM_Prop",
#                 "LogFold_CountAdj", "LogFold_Prop", "Pval_CountAdj", "Pval_Prop",
#                 "Child_entropy", "Parent_entropy")
dist_dir = paste(result_dir, "/dist", sep="")
sampleCountThres = 1 #only compare if >=3 samples available

#Output
plot_dir = paste(result_dir, "/plots", sep=""); for(i in 1:length(plot_dir)) { suppressWarnings(dir.create(plot_dir[i])) }
cp_dir = paste(plot_dir, "/changepoint", sep=""); for(i in 1:length(cp_dir)) { suppressWarnings(dir.create(cp_dir[i])) }
single_dir = paste(plot_dir, "/singlephen", sep=""); for(i in 1:length(single_dir)) { suppressWarnings(dir.create(single_dir[i])) }


library(stringr)
library(colorspace)
library(changepoint) # library(proxy)
library(FKF)
library(fastcluster)
library(dendextend)
library(circlize)
library(Rtsne)
library(MASS)
library(RDRToolbox)
library(scater, quietly = TRUE)
library(foreach)
library(doMC)
library(lubridate) #if there are date variables
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")


countThres = 20
plotpc = 5 #number of pca pc to plot
doHC = F
doISO = T
doTsne = T
theta=.5 #for Tsne
dofullPCA = F #do PCA for all cell popoulations not just first level

methods = c("BinSeg","AMOC","PELT") #Changepoint analysis; AMOC at most one change, PELT is poisson (quick exact), BinSeg (quick, approx), SegNeigh (slow, exact)
usemethod="AMOC" #use this method to determine changepoints for pvalue calculation
mds_type = c("iso", "mds")
link = c("ward.D", "ward.D2", "mcquitty") #, "single", "complete") # , "median", "average", "centroid"

interested = c("aml", "tube", "specimen") #sampleMeta columns to plot
interestedcont = "" #continuous variables
splitby = c("tube") #"none" if split by nothing

layers = c(1,3,5)



start = Sys.time()

no_cores = 8#detectCores() - 1
registerDoMC(no_cores)

sampleMeta = get(load(sampleMeta_dir))
phenoMeta = get(load(phenoMeta_dir))
#interestedCols0 = which(colnames(sampleMeta)%in%interested)

#order samples by date, exclude genotypes with less than 3 samples
sampleMeta2 <- sampleMeta[order(sampleMeta$specimen, sampleMeta$aml, sampleMeta$tube),] #sanger centre

ftWTGT = "control"

a = foreach (mcp=matrix_type) %dopar% { cat("\n  ", mcp, ": ")
  start2 = Sys.time()
  
  mm = as.matrix(get(load(paste0(matrix_dir, mcp,".Rdata"))))
  mind = match(sampleMeta2[,1], rownames(mm))
  mind = mind[!is.na(mind)]
  m0 = mm[match(sampleMeta2[,1], rownames(mm)),]
  sm0 = as.data.frame(sapply(sampleMeta2, function(x) as.numeric(factor(x, ordered=T))))
  sm1 = sm0[match(rownames(m0),sampleMeta2[,1]),]
  
  tubeind = list()
  if (splitby=="none") { tubeind[["all"]] = 1:nrow(sm1)
  } else {
    for (tube in unique(sm1[,splitby])) { if (!is.na(tube)) tubeind[[as.character(tube)]] = which(sm1[,splitby]==tube) }
  }
  
  ## get interested columns
  interestedCols = which(colnames(sm1)%in%interested)
  
  for (layer in layers) {
    pl = colnames(mm)%in%phenoMeta$phenotype[phenoMeta$phenolevel==layer]
    if (!sum(pl)>0) next
    m1 = m0[,pl]
    
    for (tube in names(tubeind)) {
      m = m1[tubeind[[tube]],]
      if (!sum(m!=0)>0) next
      sm = sm1[tubeind[[tube]],]
      
      ## pca analysis ------------------------------------------------
      iso <- 0
      if (doISO) { cat("iso; ")
        iso <- 1
        fit <- Isomap(m,k=2)
      }
      cat("pca; ")
      pc <- prcomp(m)
      
      #pca scatterplot
      pngname <- paste0(plot_dir, "/pca_iso_", mcp, "_",splitby,"-",tube,"_layer", str_pad(layer, 2, pad = "0"), ".png")
      png(filename=pngname, width=length(interestedCols)*400, height=(plotpc+1+iso)*400)
      layout(matrix(c(rep(1,length(interestedCols)), 2:(((2*plotpc)+iso)*length(interestedCols)+1)),ncol=length(interestedCols),byrow=T))
      
      plot(pc$sdev^2/sum(pc$sdev^2), xlim = c(0, 15), type = "b", pch = 16, xlab = "principal components", ylab = "variance explained")
      for (i in 1:plotpc) {
        for (col in interestedCols) {
          coloursm <- sm[,col]
          plot(pc$x[,i+1], pc$x[,i], col = coloursm, main = paste0("PCA ", colnames(sm)[col]), xlab = paste0("PC_",i+1), ylab = paste0("PC_",i))
          points(0, 0, pch = 3, cex = 4, lwd = 4)
        }
        for (col in interestedCols) {
          x_split = as.integer(factor(sm[,col]))
          cont = F; if (colnames(sm)[col]%in%interestedcont) cont = T
          x_splitnames = sort(unique(sm[,col]))
          
          cor = cor(x_split, pc$x[,i])
          
          coloursm <- sm[,col]
          if (cont) {
            plot(sm[,col], pc$x[,i], col = coloursm, main = paste0("PCA ", colnames(sm)[col]," Pearson Cor = ", cor), xlab = colnames(sm)[col], ylab = paste0("PC_",i))
          } else {
            xy_boxplot = lapply(sort(unique(x_split)), function(x) pc$x[x_split==x,i])
            boxplot(xy_boxplot, lwd = 1, outline=F, ylim=c(min(pc$x[,i]),max(pc$x[,i])),
                    main = paste0("PCA ", colnames(sm)[col]," Pearson Cor (ok if 2 var values) = ", cor),
                    xaxt = 'n', xlab = colnames(sm)[col], ylab = paste0("PC_",i)) #,xaxt ='n', xlab=testcol
            axis(1, at=1:length(x_splitnames), labels=x_splitnames)
            points(jitter(x_split, factor=1), pc$x[,i], col = coloursm)
          }
        }
      }
      if (doISO) {
        for (col in interestedCols) {
          coloursm <- sm[,col]
          plot(fit$dim2, col = coloursm, main = paste0("ISO ", colnames(sm)[col]))
          points(0, 0, pch = 3, cex = 4, lwd = 4)
        }
      }
      graphics.off()
    }
  }  
}

TimeOutput(start)
