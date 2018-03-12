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
matrix_dir = paste(result_dir,  "/matrix", sep="")

#Output
biclust_dir = paste(result_dir,  "/biclust_201803", sep=""); suppressWarnings(dir.create (biclust_dir))

#Libraries/Functions
library(biclust)
library(foreach)
library(doMC)
library(stringr)
source("~/projects/IMPC/code/_funcAlice.R")

#Setup Cores
no_cores = 2#detectCores()-1
registerDoMC(no_cores)








#Options for script
# sampletimes = 5
overwrite = F #overwrite biclust?
writecsv = F

cellCountThres = c(2000) #insignificant if count under
matrix_count = c("CountAdj") #...

plot_size = c(500,500)
plot_size_bar = c(1300,2000)

#for every distance type
bcmethods = c("plaid","CC","Xmotifs","spectral","bimax","quest")    # have to change this manually in function...


tube = 4














start = Sys.time()



#Prepare data
matrix_type = list.files(path=result_dir,pattern=glob2rx("matrix*.Rdata"))
matrix_type = matrix_type[!matrix_type=="matrixCount.Rdata" & !matrix_type=="matrixProp.Rdata"]
matrix_type = gsub("matrix|.Rdata","",matrix_type)
matrix_type = matrix_type[!grepl("leaf",matrix_type)]
matrix_type = matrix_type[!grepl("KO|Max",matrix_type)]
matrix_type = matrix_type[grepl("Pval",matrix_type)]
matrix_type = matrix_type[grepl("CountAdj",matrix_type)]
matrix_type = matrix_type[!grepl("_[0-9]$",matrix_type)]

m0 = get(load(paste0(matrix_dir, matrix_count,".Rdata")))
phenoMeta = get(load(phenoMeta_dir))
sampleMeta = get(load(sampleMeta_dir))
k0 = c(1,2,4,max(phenoMeta$phenolevel)) #1,4, # how many layers to consider i.e. k=max(phenolevel) only

fe = foreach(mcp=matrix_type) %dopar% {
  tryCatch({
    cat("\n", mcp, " ",sep="")
    start2 = Sys.time()
    
    #start a list of phenotypes included in each distance matrix calculation such that none are repeated
    leavePhenotype = list()
    doneAll = F
    
    #load feature matrix
    mresult = Loadintermatrices(paste0(matrix_dir, mcp,".Rdata"))
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
        pm = mmlresult$pm
        leavePhenotype = mmlresult$leavePhenotype
        doneAll = mmlresult$doneAll
        
        m_rownames = rownames(m)
        rownames(m) = gsub(".labelled.fcs","",rownames(m))
        
        
        for (bcmethod  in bcmethods) {
          
          #assume after splitting dist filename by "_", distance is second element
          dname = paste(biclust_dir, "/", mcp, "_", bcmethod, "_layer", str_pad(k, 2, pad = "0"), "_countThres-", countThres, sep = "" )
          
          ## calculate distances
          if (is.null(m)) return(NULL)
          dname0 = paste0(dname, ".Rdata", sep="")
          if (overwrite | !file.exists(dname0)) {
            if (is.null(dim(m))) { #edge feature
              m = Reduce('cbind',m); if (all(a==0)) next
            } else { #node feature
              if (all(m==0)) next 
            }
            
            m = m[sampleMeta$tube[match(rownames(m),sampleMeta$fileName)],]
            
            # rownames(m) = sampleMeta$gene[match(rownames(m),sampleMeta$fileName)]
            
            if (bcmethod == "plaid") d = biclust(as.matrix(m), method=BCPlaid())
            if (bcmethod == "CC") d = biclust(as.matrix(m), method=BCCC())
            if (bcmethod == "Xmotifs") d = biclust(as.matrix(m), method=BCXmotifs())
            if (bcmethod == "spectral") d = biclust(as.matrix(m), method=BCSpectral())
            if (bcmethod == "bimax") d = biclust(as.matrix(m), method=BCBimax())
            if (bcmethod == "quest") d = biclust(as.matrix(m), method=BCQuest())
            
            save(d, file=dname0)
            if (writecsv) write.csv(as.matrix(d), file=gsub(".Rdata",".csv",checkm(d,dname0)))
            
          } else {
            d = get(load(dname0))
          }
          
          if (d@Number > 0) {
            png(paste0(dname, "_bar.png", sep=""), height=plot_size_bar[1], width=plot_size_bar[2])
            par(mar=c(2,6,3,2))
            plotclust(d,as.matrix(m))
            graphics.off()
            
            png(paste0(dname, "_bubble.png", sep=""), height=plot_size_bar[1], width=plot_size_bar[2])
            par(mar=c(20,20,20,20))
            bubbleplot(as.matrix(m),d)
            graphics.off()
            
            png(paste0(dname, "_heatmap0.png", sep=""), height=plot_size_bar[1], width=plot_size_bar[2])
            par(mar=c(5,3,6,5))
            try({
              heatmapBC(as.matrix(m),d)
            })
            graphics.off()
            
            rowcolno = ceiling(sqrt(d@Number))
            png(paste0(dname, "_heatmap.png", sep=""), height=plot_size[1]*rowcolno, width=plot_size[2]*rowcolno)
            par(mar=c(25,50,20,50))
            par(mfrow=rep(rowcolno,2))
            for (BCi in 1:d@Number) {
              try({
                drawHeatmap(as.matrix(m),d,BCi,plotAll=T)
              })
            }
            graphics.off()
            
            for (attri in c("aml")) {
              png(paste0(dname, "_stats_",attri,".png", sep=""), height=plot_size[1]*rowcolno, width=plot_size[2]*rowcolno)
              par(mar=c(10,5,3,2), mfrow=rep(rowcolno,2))
              for (BCi in 1:d@Number) {
                filenames = m_rownames[d@RowxNumber[,BCi]]
                attri_values = sampleMeta[match(filenames,sampleMeta$fileName),attri]
                # try({
                barplot(table(attri_values), las=2)
                # })
              }
              graphics.off()
            }
            
          }
          
          
          
        }
      } #layer
    } #countThres
    
    TimeOutput(start2)
  }, error = function(err) { cat(paste("ERROR:  ",err)); return(T) })
  return(F)
}

TimeOutput(start)




