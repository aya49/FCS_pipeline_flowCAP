# aya43@sfu.ca 20161220
# Normalizes cell count matrix

root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

options(stringsAsFactors=FALSE)
options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="")
matrixCount_dir = paste(result_dir, "/matrixCount.Rdata", sep="")

#Output
matrixCountAdj_dir = paste(result_dir, "/matrixCountAdj.Rdata",sep="")
matrixCountAdj_dircsv = paste(result_dir, "/matrixCountAdj.csv",sep="")
norm_dir = paste(result_dir, "/cellCountNormFactor",sep="")
normFactor_dir = paste(result_dir, "/normFactor.Rdata", sep="")
normFactorDiffLog_dir = paste(result_dir, "/normFactorDiffDensLogged.Rdata", sep="")

library(stringr)
library(pracma)
source("code/_funcAlice.R")



start = Sys.time()

suppressWarnings(dir.create(norm_dir))
sampleMeta = get(load(sampleMeta_dir))
matrixCount = get(load(matrixCount_dir))
#phenotype on rows
if (ncol(matrixCount)!=nrow(sampleMeta)) { matrixCount = t(matrixCount) }


## Taken from TMM ------------------------------------------------------
x = as.matrix(matrixCount)[-1,] #take out total cell count

p=0.75
lib.size = matrixCount[1,]
refColumn = which( matrixCount[1,]==median(matrixCount[1,which(sampleMeta$aml=="normal")]) ) #Get reference column: median normal patient
#  f75 = apply(t(t(x)/lib.size),2,function(data) quantile(data,p=p))
#  refColumn = which.min(abs(f75-mean(f75)))

f = rep(NA,ncol(x))
fdiff = rep(NA,ncol(x)) #diff between density peak and value (note: logged)
ref = x[,refColumn]
nR = lib.size[refColumn]
doWeighting = T
Acutoff = -1e10
logratioTrim = .3
sumTrim = 0.05
for(i in ncol(x):1) { cat(i," ",sep="")
  obs = x[,i]
  nO = lib.size[i]
  logR = log2((obs/nO)/(ref/nR))			#log ratio of expression, accounting for library size
  absE = (log2(obs/nO) + log2(ref/nR))/2	#absolute expression
  v = (nO-obs)/nO/obs + (nR-ref)/nR/ref	 #estimated asymptotic variance
  
  #remove infinite values, cutoff based on A
  fin = is.finite(logR) & is.finite(absE) & (absE > Acutoff)
  logR = logR[fin]
  absE = absE[fin]
  v = v[fin]
  
  if(max(abs(logR)) < 1e-6) { f[i] = 1; next }
  
  #taken from the original mean() function
  n = length(logR)
  loL = floor(n * logratioTrim) + 1
  hiL = n + 1 - loL
  loS = floor(n * sumTrim) + 1
  hiS = n + 1 - loS
  
  #keep = (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
  #a fix from leonardo ivan almonacid cardenas, since rank() can return
  #non-integer values when there are a lot of ties
  keep = (rank(logR)>=loL & rank(logR)<=hiL) & (rank(absE)>=loS & rank(absE)<=hiS)
  
  if(doWeighting) {
    f[i] = sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE)
  } else { f[i] = mean(logR[keep], na.rm=TRUE) }
  
  #Results will be missing if the two libraries share no features with positive counts
  #In this case, return unity
  if(is.na(f[i])) { f[i] = 0 }
  
  #check if close to peak; if not, switch to peak
  d = density(log2((matrixCount[,i])/matrixCount[,refColumn]), na.rm=T)
  p = as.matrix(findpeaks(d$y)); if(ncol(p)==1) { p = t(p) }
  p1 = d$x[p[which.max(p[,1]),2]]
  fdiff[i] = p1-f[i]
  
  pngname = paste(norm_dir, "/cellCountNormFactor_",str_pad(i, 4, pad = "0"), "_",sampleMeta$gene[i], ".png", sep = "" )
  png (file=pngname , width=700, height=1800)
  par(mfrow=c(3,1), mar=(c(5, 5, 4, 2) + 0.1))
  
  plot(d); abline(v=f[i], col="red"); abline(v=p1, col="blue"); 
  
  plot((matrixCount[,i]+matrixCount[,refColumn])/2, log2(matrixCount[,i]/matrixCount[,refColumn]), cex=.5, main=paste("mean count vs. log2 fold change: ", sampleMeta$gene[i]," over refColumn ", sampleMeta$gene[refColumn],": f=",f[i], sep=""))
  abline(h=f[i], col="red")
  
  #if f[i] too far from peak
  if (abs(f[i]-p1)>.5) {
    abline(h=p1, col="blue")
    #f[i] = p1
  }
  
  f[i] = 1/2^f[i]
  
  plot((matrixCount[,i]+matrixCount[,refColumn])/2, log2((matrixCount[,i]*f[i])/matrixCount[,refColumn]), cex=.5, main=paste("AFTER CHANGE: mean count vs. log2 fold change: ", sampleMeta$gene[i]," over refColumn ", sampleMeta$gene[refColumn],": f=",f[i], sep=""))
  abline(h=0, col="red")
  dev.off()
}
#multiple of 1
rm(x)

matrixCountAdj = sapply(c(1:ncol(matrixCount)), function(x) {matrixCount[,x]*f[x]})
colnames(matrixCountAdj) = colnames(matrixCount)
#phenotype on cols
matrixCountAdj = t(matrixCountAdj)

save(f, file=normFactor_dir)
save(fdiff, file=normFactorDiffLog_dir)
save(matrixCountAdj, file=matrixCountAdj_dir)
write.csv(matrixCountAdj, file=matrixCountAdj_dircsv, row.names=F)

TimeOutput(start)







