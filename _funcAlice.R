

### Basic Functions ######################################################################



# Prints out the time since start_time. Used for optimizing code.
TimeOutput = function(start_time) {
  start_time = as.POSIXct(start_time)
  dt = difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  cat(format(.POSIXct(dt,tz="GMT"), "%H:%M:%S"))
}



'%between%'=function(x,rng) (x>rng[1] & x<rng[2]) | (x<rng[1] & x>rng[2])


## Input: matrix, threshold
## Output: columns where all values are below= a threshold
colIndBelow = function(m,thres) {
  return(which(apply(m,2,function(x) all(x<=thres))))
}

## Output: Random matrix
randomMatrix = function(nrow,ncol) {
  return(matrix(rexp(ncol*nrow, rate=.1), ncol=ncol))
}

## Input: package
## Output: TRUE if package is installed, FALSE otherwise
is.installed = function(pkg){
  is.element(pkg, installed.packages()[,1])
} 

charOccurences = function(char, str) {
  str2 = gsub(char,"",str)
  return (nchar(str) - nchar(str2))
}

## Input: file paths of matrices with column names you want to compare
## Output: list of unique column names; last vector in list corresponds to which filepath the unique column names belong to
colNamesSame = function(pathList) {
  count = 1
  colNameSame = NULL
  colNameSame[[1]] = colnames(read.csv(pathList[1]))
  colNameSameIndex = c(1)
  for (i in 2:length(pathList)) {
    s = colnames(read.csv(pathList[i]))
    if(!TRUE %in% (colNameSame%in%list(s)) ) {
      count = count+1
      colNameSameIndex[count] = i
      colNameSame[[count]] = s
    }
  }
  colNameSame[[count+1]] = colNameSameIndex
  return(colNameSame)
}

## Input: file path and a file extension 
## Output: List of all file names in given path with the given extension
fileNames = function(pathList, ext="fcs") {
  for (i in 1:length(pathList)) {
    temp.str = strsplit(pathList[i], split="/")
    pathList[i] = sub(paste(".", ext, sep=""), "", temp.str[[1]][ length(temp.str[[1]]) ], ignore.case=TRUE)
  }
  return(pathList)
}

## Input: file path
## Output: List of all last folder names in given path (Genotype)
folderNames = function(pathList) {
  folders = c()
  for (i in 1:length(pathList)) {
    temp.str = strsplit(pathList[i],split="/")
    folders[i] = temp.str[[1]][length(temp.str[[1]])-1]
  }
  return(folders)
}




getGTindex = function(attlist,control,sampleCountThres,explist=NULL, all=F) {
  if (is.null(explist)) explist = attlist
  if (length(control)==1) {
    controlIndex = grep(control, attlist)
  } else {
    controlIndex = which(attlist%in%control)
  }
  if (all) {
    exp = unique(explist[])
  } else {
    exp = unique(explist[-controlIndex]) #Unique KO Genotypes
  }
  expIndex = list() #Index of unique KO Genotypes (one genotype per vector in the list)
  for (i in length(exp):1) { expIndex[[i]] = which(explist==exp[i]) }
  goodexpIndex = which( unlist(lapply(expIndex, length)) >=sampleCountThres ) #Index of unique KO genotypes that has 3+ samples available expIndex(expIndexIndex)
  return(list(attlist=attlist, control=control, controlIndex=controlIndex, exp=exp, expIndex=expIndex, goodexpIndex=goodexpIndex))
}







## Input: Clt, a matrix
## Output: Outputs column inds of duplicate columns
duplicateindM = function(clt) {
  delcol_dup = duplicated(t(clt))
  col_dup = 1:ncol(clt); col_dup[delcol_dup] = sapply(which(delcol_dup),function(x) {
    for (xi in which(!delcol_dup)) {
      if (identical(clt[,xi],clt[,x]) ) { a = xi; break }
    }
    return(a)
  })
  return(col_dup)
}


## Input: vector
## Output: vector with elements being indices of who it replicates; unique values will simply be its own index
duplicateind = function(cl) {
  delcol_dup = duplicated(cl)
  col_dup = 1:length(cl); col_dup[delcol_dup] = sapply(which(delcol_dup),function(x) {
    for (xi in which(!delcol_dup)) {
      if (cl[xi]==cl[x]) { a = xi; break }
    }
    return(a)
  })
  return(col_dup)
}



## Input: string
## Output: last n characters
substrRight = function(x, n) substr(x, nchar(x)-n+1, nchar(x))





### Flow Functions ######################################################################


## Input:path list
## Output:list of frames from paths
loadFrames = function(pathList) {
  myFrames = new.env()
  cat("loading", length(pathList), "files:")
  for (i in 1:length(pathList)) {
    myFrames[[as.character(i)]] = get(load(pathList[i]))
    cat(" ", i, sep="")
  }
  return(myFrames)
}

## Input: a phenotype with all markers
## Output: vector of markers
getMarkers = function(phenotype) {
  return (unlist(strsplit(as.character(phenotype), '[-+]')))
}



## Input: number of markers
## Output: number of nodes, edges
getGraphInfo = function(m) { 
  npl = epl = rep(0,m)
  for (i in 1:m) {
    npl[i] = choose(m,i)*(2^i)
    epl[i] = npl[i]*i
  }
  n = 3^m
  e = n*2*m/3
  return (list(n=n, npl=npl, e=e, epl=epl))
}





## Input: PhenoMeta
## Output: Index of leaves (phenotypes with no children)
getleaves = function (phenoMeta, no_cores) {
  require(foreach)
  require(doMC)
  
  registerDoMC(no_cores)
  
  finalLevel = which(phenoMeta$phenolevel==max(phenoMeta$phenolevel))
  notFinalLevel = setdiff(1:nrow(phenoMeta), finalLevel)
  nflleaves = foreach (i=1:length(notFinalLevel), .combine="c") %dopar% {
    pheni = unlist(strsplit(phenoMeta$phenocode[i],""))
    phenind = which(pheni!="0")
    zeroind = setdiff(1:length(pheni),phenind) #must have zeros because notFinalLevel
    childphenocode = as.vector(sapply(zeroind, function(x) { pi1=pi2=pheni; pi1[x]="1"; pi2[x] = "2"; return(c(paste(pi1,collapse=""),paste(pi2,collapse=""))) } ))
    childrenind = match(childphenocode, phenoMeta$phenocode)
    if (all(is.na(childrenind))) return(i)
    return(NA)
  }
  leaves = notFinalLevel[nflleaves[!is.na(nflleaves)]]
  return(list(leaves=leaves, finalLevel=finalLevel))
}










## Input: Phenotype
## Output: Phenotype, PhenoCode, Phenolevel (number of markers per phenotype)
getPhen = function(phen, markers=NULL, phenotype=T) {
  require(flowType)
  require(stringr)
  if (length(grep("+", phen))>0) { #if phen is phenotype, build phenocde
    if (is.null(markers)) markers = unlist(strsplit(phen[which.max(nchar(phen))],"[+-]"))
    phenoCode = unlist(lapply(phen, function(x){return( encodePhenotype(x, markers) )}))
  } else if (phenotype) { #if phen is phenocode
    if (is.null(markers)) {
      cat("input markers")
      return(NULL)
    } else {
      phenoCode = phen
      phen = unlist(lapply(pheno, function(x){return( decodePhenotype(x, markers, rep(2,length(markers))) )}))
    }
  }
  require(stringr)
  phenolevel = str_count(phen, "[+-]")
  
  if (phenotype) return(list(phenotype=phen,phenoCode=phenoCode,phenoLevel=phenolevel))
  return(list(phenoCode=phenoCode,phenoLevel=phenolevel))
}






## Jensen-shannon-divergence (half of kullback divergence both ways)
# http://stackoverflow.com/questions/11226627/jensen-shannon-divergence-in-r

jsd = function(p,q, list=T) { #p,q are two distributions, two lists of distributions if to avg over
  if (list) {
    JS = 0
    for (i in 1:length(p)) {
      m = (p[[i]] + q[[i]])/2
      JS = ((sum(p[[i]]*log(p[[i]] / m)) + sum(q[[i]]*log(q[[i]] / m))) /2) + JS
    }
    JS = JS/length(p)
  } else {
    m = (p + q)/2
    JS = ((sum(p * log(p / m)) + sum(q * log(q / m))) /2) + JS
  }
  return(JS)
}


### FPM (set aside for now) ######################################################################

## Input: 2 rows to compare, distance measure
## Output: Ranking of farthes to nearest features
featDist = function(m, dis) {
  require(vegan)
  m[1,which(m[1,]==0 & m[2,]!=0)] = 1
  m[2,which(m[1,]!=0 & m[2,]==0)] = 1
  
  return(unlist( lapply(c(1:ncol(m)), function(x) { return(vegdist(as.matrix(m[,x]), method=dis))}) ))
}

featDistList = function(m1, m2=NULL, dis) {
  d = NULL
  if (is.null(m2)) {
    for (i in 1:nrow(m1)) {
      for (j in (i+1):nrow(m1)) {
        d[[i]][[j]] = featDist(m1[c(i,j),],dis)
      }
    }
  } else {
    for (i in 1:nrow(m1)) {
      for (j in 1:nrow(m2)) {
        d[[i]][[j]] = featDist(rbind(m1[i,],m2[j,]),dis)
      }
    }
  }
  return(d)
}

createTransactionMatrix = function(phenoCode, markers) {
  m1 = matrix(0,nrow=length(phenoCode),ncol=length(markers))
  m2 = matrix(0,nrow=length(phenoCode),ncol=length(markers))
  for (i in 1:length(phenoCode)) {
    pc = as.numeric(strsplit(as.character(phenoCode[i]),"")[[1]])
    pc1 = pc2 = pc
    pc1[which(pc1==2)] = 0
    pc2[which(pc2==1)] = 0
    pc2[which(pc2==2)] = 1
    m1[i,] = pc1
    m2[i,] = pc2
  }
  colnames(m1) = paste0(markers,"-",sep="")
  colnames(m2) = paste0(markers,"+",sep="")
  return(cbind(m1,m2))
}




### ETC ######################################################################



# from https://aurelienmadouasse.wordpress.com/2012/01/13/legend-for-a-continuous-color-scale-in-r/
legend.col = function(col, lev){
  opar = par
  n = length(col)
  bx = par("usr")
  box.cx = c(bx[2] + (bx[2] - bx[1]) / 1000, bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
  box.cy = c(bx[3], bx[3])
  box.sy = (bx[4] - bx[3]) / n
  
  xx = rep(box.cx, each = 2)
  
  par(xpd = TRUE)
  for(i in 1:n){
    yy = c(box.cy[1] + (box.sy * (i - 1)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i - 1)))
    polygon(xx, yy, col = col[i], border = col[i])
  }
  par(new = TRUE)
  plot(0, 0, type = "n",
       ylim = c(min(lev), max(lev)),
       yaxt = "n", ylab = "",
       xaxt = "n", xlab = "",
       frame.plot = FALSE)
  axis(side = 4, las = 2, tick = FALSE, line = .25)
  par = opar
}












# http://stackoverflow.com/questions/4787332/how-to-remove-outliers-from-a-dataset
remove_outliers = function(x, na.rm = TRUE, ...) {
  qnt = quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H = 1.5 * IQR(x, na.rm = na.rm)
  y = x
  y[x < (qnt[1] - H)] = NA
  y[x > (qnt[2] + H)] = NA
  y
}





































### Plotting Functions ######################################################################


plotTsne = function(x, continuous=F, main, colBarWidth=.08, colBarTitle="", leg=T, text=F) {
  require(lubridate)
  x1 = rownames(x)
  if (is.Date(x1) | is.timepoint(x1)) {
    x1 = as.numeric(x1-min(x1))
  }
  x = x[order(x1),]
  x1 = sort(x1)
  if (continuous) {
    ts = as.numeric(factor(x1))
    colour = heat.colors(length(unique(ts)))
    plot(x, col=colour[ts], pch=20, main=main)
    legend.col(col = colour, lev = ts)
  } else {
    c1 = 7
    if (length(unique(x1))/7>25) {
      c1 = ceiling(length(unique(x1))/25)
    }
    colour = rainbow(c1)
    cp = expand.grid(c(1:c1),c(20,1:19,21:25))
    plot(x, t='n', main=main)
    if (text) {
      text(x,rownames(x))
    } else {
      for (g in 1:length(unique(x1))) {
        points(x[which(x1%in%unique(x1)[g]),1],x[which(x1%in%unique(x1)[g]),2], col=colour[cp[g%%nrow(cp),1]], pch=cp[g%%nrow(cp),2])
      }
    }
    if (leg) legend("topleft", legend=unique(x1), fill=colour[rep(c(1:c1),g)], pch=rep(c(20,1:19,21:25),each=c1)[c(1:g)])
  }
}








### Matrix Functions ######################################################################


## Input: matrices
## Output: check matrix for anything weird, add onto name and return modified filenames
checkm = function(m,fname) {
  d = as.matrix(m)
  if (sum(is.na(d))>0) fname = paste0(fname,"_NA")
  if (sum(is.nan(d))>0) fname = paste0(fname,"_NaN")
  isinf = sum(d==Inf)
  if (!is.na(isinf)) { if (isinf>0) fname = paste0(fname,"_Inf")}
  return(fname)
}

## return coordinates of na in a matrix http://stackoverflow.com/questions/4833300/get-positions-for-nas-only-in-the-middle-of-a-matrix-column
# x1 is an indexing sequence which increments at non-NA entries, counting top-to-bottom. x2 is the same, counting backward (bottom-to-top). They are only both nonzero at internal entries enclosed by non-NAs on both top and bottom, hence their non-NA indices counting in both directions are >0. Finally gate that with a & to filter out just the internal NAs. x1 = nonNA_idx.tb, x2 = nonNA_idx.bt
checkmna = function(m) {
  z = as.matrix(m)
  isNA = is.na(z)
  # Vertical index which increments at non-NA entries, counting top-to-bottom:
  nonNA_idx.tb = apply(!isNA, 2, cumsum)
  # Vertical index which increments at non-NA entries, counting bottom-to-top:
  nonNA_idx.bt = apply(!isNA, 2, function(x) { rev(cumsum(rev(x))) })
  which(isNA & nonNA_idx.tb>0 & nonNA_idx.bt>0, arr.ind=TRUE)
}


## trim a list of matrices according to reference matrix, colnames and rownames must match. 
trimlistmatrix = function(refmatrix,refdel=0, targetlist,targetdel=0, ncore=1) {
  require(foreach)
  require(doMC)
  registerDoMC(ncore)
  
  pt0 = union(names(targetlist), Reduce('union',lapply(targetlist,function(y) return(colnames(y)))) )
  lpi = intersect(pt0, colnames(refmatrix))
  targetlistTRIM = targetlist[names(targetlist)%in%colnames(refmatrix)]
  nrows = nrow(targetlist[[1]])#original
  targetlistTRIM = foreach (i = 1:length(targetlistTRIM)) %dopar% {
    a = as.matrix(targetlistTRIM[[i]],nrow=nrows)
    dimnames(a) = dimnames(targetlistTRIM[[i]])
    row = !is.na(match(rownames(a),rownames(refmatrix)))
    col = !is.na(match(colnames(a),colnames(refmatrix)))
    if (sum(col)>0) {
      a0 = a
      a = as.matrix(a0[row,col],ncol=sum(col))
      colnames(a) = colnames(a0)[col]
      rownames(a) = rownames(a0)[row]
      a[ matrix(refmatrix[,!is.na(match(colnames(refmatrix),colnames(a)))],ncol=sum(col)) ==refdel ] = targetdel
      # if (sum(col)>0) a[ refmatrix[!is.na(match(rownames(refmatrix),rownames(a)))m, !is.na(match(colnames(refmatrix),colnames(a)))]==0 ] = targetdel
    } else { return(NULL) }
    return(Matrix(a,sparse=T))
  }
  names(targetlistTRIM) = names(targetlist)[names(targetlist)%in%colnames(refmatrix)]
  targetlistTRIM = targetlistTRIM[!sapply(targetlistTRIM, is.null)] 
  return(targetlistTRIM)
}





## Input : list of matrix/matrixList with phenotype on colummn, samples on rows (filenames)
## Output: Trimmed version of matrix
Loadintermatrices = function(mfilenames,verbose=T) {
  require(Matrix)
  mind = 1
  mml = list()
  mmlname = c()
  if (verbose) cat("loading feature matrices, ",sep="")
  for (mind in 1:length(mfilenames)) {
    mcp = mfilenames[mind]
    if (!file.exists(mcp)) {cat("doesn't exist"); break}
    mml[[mind]] = get(load(mcp))
    if (data.class(mml[[mind]])=="dist") mml[[mind]] = as.matrix(mml[[mind]])
    mmlname[mind] = mcp
    if (mind==1) {
      if (is.null(dim(mml[[1]]))) {
        pt = union(names(mml[[1]]), Reduce('union',lapply(mml[[1]],function(y) return(colnames(y)))) )
        gt = rownames(mml[[1]][[1]])
      } else {
        pt = colnames(mml[[1]])
        gt = rownames(mml[[1]])
      }
    } else {
      if (is.null(dim(mml[[mind]]))) {
        pt0 = union(names(mml[[mind]]), Reduce('union',lapply(mml[[mind]],function(y) return(colnames(y)))) )
        gt0 = rownames(mml[[mind]][[1]])
      } else {
        pt0 = colnames(mml[[mind]])
        gt0 = rownames(mml[[mind]])
      }
      pt = intersect(pt, pt0)
      gt = intersect(gt, gt0)
    }
  }
  
  #trim matrices (sample) based on samples/phenotypes in common
  if (length(mml)>1) {
    if (verbose) cat("trimming, ",sep="")
    for (i in 1:length(mml)) {
      b = b0 = mml[[i]]
      if (is.null(dim(b))) {
        gt0 = rownames(b[[1]])
        pt0 = union(names(b), Reduce('union',lapply(b,function(y) return(colnames(y)))) )
        ind = match(gt,gt0)
        b = foreach (j = 1:length(b)) %dopar% {
          ptind = match(pt,colnames(b[[j]]))
          if(sum(!is.na(ptind))>0){
            ptind = ptind[!is.na(ptind)]
            # if (is.null(dim(b[[j]]))) {
            #   a = as.matrix(b[[j]],ncol=1)[ind,ptind]
            #   dimnames(a) = dimnames(b[[j]][ind,ptind])
            # } else { 
            a = b[[j]][ind,ptind]
            if (is.null(dim(a))) a = matrix(a,nrow=length(ind))
            rownames(a) = rownames(b[[j]])[ind]
            colnames(a) = colnames(b[[j]])[ptind]
            # }
          } else { a = NULL }
          return(a)
        }
        names(b) = names(b0)
        map = match(pt,names(b0))
        b1 = b[map[!is.na(map)]]
        b1ind = !sapply(b1, is.null)
        if (sum(b1ind)>0) {
          mml[[i]] = b1[which(b1ind)]
        } else {
          mml[[i]] = NULL
        }
      } else {
        gt0 = rownames(b)
        pt0 = colnames(b)
        mml[[i]] = b[match(gt,gt0),match(pt,pt0)]
      }
    }
    if (length(mml)<length(mmlname)) { names(mml) = mmlname[1:length(mml)]
    } else { names(mml) = mmlname }
    mml = mml[!sapply(mml, is.null)] 
  } else { names(mml) = mfilenames }
  if (verbose) cat("done")
  return(list(mml=mml,pt=pt,gt=gt))
}



trimMatrices = function(mml0,m0,pt,gt,phenoMeta,leavePhenotype,doneAll, countThres,k) {
  lowCountpt = colnames(m0)[apply(m0,2,function(x) all(x<=countThres))]
  if (!length(lowCountpt)>0) lowCountpt = c()
  
  if (!is.null(dim(phenoMeta))) phenolevel = phenoMeta$phenolevel
  
  highLevelpt = colnames(m0)[phenolevel>k]
  if (!length(highLevelpt)>0) highLevelpt = c()
  
  # cat("length of countthres: ", length(lowCountInd))
  # cat("\nlength of highlevelind: ", length(highLevelInd))
  
  ## Load & fix cell count/countAdj/proportion matrix 
  dpi0 = union(lowCountpt,highLevelpt)
  # cat("\nlength of dpi: ", length(dpi0), "\n")
  
  lpi = setdiff(pt,dpi0)
  
  #check if indices already calculated for on this matrix
  if (Position(function(x) identical(x, lpi), leavePhenotype, nomatch = 0) > 0) {cat("-skipped ", sep=""); return(0)}
  leavePhenotype[[paste0(k="k.",k,"_countThres.",countThres)]] = lpi
  
  if (!doneAll & length(lpi)==length(pt)) doneAll = T
  
  
  #trim matrices (phenotype)
  mmlname = names(mml0)
  mml = mml0
  for (i in 1:length(mml)) {
    b0 = b = mml[[i]]
    if (is.null(dim(b))) {
      pt0 = union(names(b), Reduce('union',lapply(b,function(y) return(colnames(y)))) )
      b = foreach (j = 1:length(b)) %dopar% {
        ptind = match(lpi,colnames(b[[j]]))
        if(sum(!is.na(ptind))>0){
          ptind = ptind[!is.na(ptind)]
          # if (is.null(dim(b[[j]]))) {
          #   a = as.matrix(b[[j]],ncol=1)[ind,ptind]
          #   dimnames(a) = dimnames(b[[j]][ind,ptind])
          # } else { 
          a = b[[j]][,ptind]
          if (is.null(dim(a))) {
            a = matrix(a,nrow=nrow(b[[j]]))
            rownames(a) = rownames(b[[j]])
            colnames(a) = colnames(b[[j]])[ptind]
          }
        } else { a = NULL }
        return(a)
      }
      names(b) = names(b0)
      map = match(lpi,names(b0))
      b1 = b[map[!is.na(map)]]
      b1ind = !sapply(b1, is.null)
      if (sum(b1ind)>0) {
        mml[[i]] = b1[which(b1ind)]
      } else {
        mml[[i]] = NULL
      }
    } else {
      pt0 = colnames(b)
      mml[[i]] = b[,match(lpi,pt0)]
    }
  }
  if (!length(mml)>0) return(NULL)
  if (length(mml)<length(mmlname)) { names(mml) = mmlname[1:length(mml)]
  } else { names(mml) = mmlname }
  mml = mml[!sapply(mml, is.null)] 
  
  if (!is.null(dim(phenoMeta))) {
    pm = phenoMeta[match(lpi,phenoMeta$phenotype),]
    return(list(mml=mml,pm=pm,leavePhenotype=leavePhenotype,doneAll=doneAll))
  } else {
    return(list(mml=mml,leavePhenotype=leavePhenotype,doneAll=doneAll))
  }
}



## Input: vector
## Output: mode element
Mode = function(x) {
  ux = unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}



## Input: dissimilarity Matrix
## Output: similarity matrix
get_graph = function(d0) {
  d0 = as.matrix(d0)
  d1 = d0
  rownames(d1) = colnames(d0) = rownames(d0)
  d1 = 1/d0
  d1[d0==0] = 1/min(d0[d0>0])
  diag(d1) = 0
  return(d1)
}

get_graphd = function(s0) {
  df = 1/s0
  df[s0==0] = 1/min(s0[s0>0])
  diag(df) = 0
  df = as.dist(df)
  return(df)
}







## Input: matrices with items labelled on row + variables corresponding to those rows to split by
## Output: list of original matrix
splitmatrix = function(m,s=NULL,split=T,defaultlabel="all",missym=T) { #only choose one of avgall=F,includeall=F, else will do includeall=T only
  mm = list()
  if (length(unique(s))>1 & split) {
    for (tubei in 1:length(unique(s))) {
      index = which(s==unique(s)[tubei])
      if (missym) { mm[[tubei]] = m[index,index]
      } else { mm[[tubei]] = m[index,] }
    }
    names(mm) = unique(s)
    # if (avgall) {
    #   mm1 = Reduce("+",mm)/length(mm)
    #   mm[[defaultlabel]] = mm1
    # }
    # if (includeall) {
    #   mm[[defaultlabel]] = m
    # }
  } else {
    mm[[1]] = m
    names(mm) = defaultlabel
  }
  return(mm)
}





## Input: distance matrix filenames
## Output: distMeta
distMetafun = function(distmfile, dis,features=NULL) {
  distmfilenames = fileNames(distmfile)
  dmfns = strsplit(distmfilenames,"_")
  dmfns = lapply(1:length(dmfns),function(y) {
    x = dmfns[[y]]
    disind = which(x%in%dis)
    if (suppressWarnings(is.na(as.numeric(x[disind-1]))) ) { typeind = 1:(disind-1); rand = 0
    } else { typeind = 1:(disind-2); rand = x[disind-1] }
    return(c(type=paste(x[typeind],collapse="_"), rand=rand, x[(disind):length(x)]))
  })
  path = distmfile
  distMeta = as.data.frame(path)
  distMeta$type = sapply(dmfns,function(x) x[1])
  distMeta$layer = as.numeric(gsub("layer","",str_extract(distmfilenames, "layer[0-9]+")))
  distMeta$norm = gsub("normalize-","",str_extract(distmfilenames, "normalize-[a-z]+"))
  distMeta$count = as.numeric(gsub("countThres-","",str_extract(distmfilenames, "countThres-[0-9]+")))
  distMeta$weighted = sapply(distmfilenames,function(x) return(grepl("weighted",x)))
  distMeta$weightedorig = sapply(distmfilenames,function(x) return(grepl("orig-weighted",x)))
  distMeta$dist = sapply(dmfns,function(x) return(x[x%in%dis]))
  distMeta$feature = sapply(dmfns,function(x) x[1])
  distMeta$rand = sapply(dmfns,function(x) as.numeric(x[2]))
  distMeta$rw = sapply (dmfns, function(x) any(grepl("_rw",x)))
  distMeta$sim = grepl("simmatrix",distmfilenames)
  if(!is.null(features)) {
    feature = paste0("^",features)
    distMeta$feature = sapply(distMeta$type, function(x) features[sapply(feature,function(y) grepl(y,x,ignore.case=T))] )
    distMeta$featend = sapply(1:nrow(distMeta), function(i) gsub(distMeta$feature[i],"",distMeta$type[i]) )
  } 
  
  return(distMeta)
}







