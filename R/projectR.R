#' @importFrom stats hclust kmeans prcomp
setOldClass("kmeans")
setOldClass("hclust")
setOldClass("prcomp")

#######################################################################################################################################
#' @import limma
#' @importFrom stats model.matrix
#' @param NP vector of integers indicating which columns of loadings object to use. The default of NP=NA will use entire matrix.
#' @param full logical indicating whether to return the full model solution. By default only the new pattern object is returned.
#' @param model Optional arguements to choose method for projection
#' @param family VGAM family function for model fitting (default: "gaussianff")
#' @param bootstrapPval logical to indicate whether to generate p-values using bootstrap, not available for prcomp and rotatoR objects
#' @param bootIter number of bootstrap iterations, default = 1000
#' @rdname projectR-methods
#' @aliases projectR
setMethod("projectR",signature(data="matrix",loadings="matrix"),function(
  data, # a dataset to be projected onto
  loadings, # a matrix of continous values to be projected with unique rownames
  dataNames = NULL, # a vector with names of data rows
  loadingsNames = NULL, # a vector with names of loadings rows
  NP=NA, # vector of integers indicating which columns of loadings object to use. The default of NP=NA will use entire matrix.
  full=FALSE, # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  family="gaussianff",  # VGAM family function (default: "gaussianff")
  bootstrapPval=FALSE, # logical to indicate whether to generate p-values using bootstrap
  bootIter=1e3 # No of bootstrap iterations
  ){

  ifelse(!is.na(NP),loadings<-loadings[,NP],loadings<-loadings)
  #if(!is.na(NP)){loadings<-loadings[,NP]} was giving warning with subset of patterns
  #match genes in data sets
  if(is.null(dataNames)){
    dataNames <- rownames(data)
  }
  if(is.null(loadingsNames)){
    loadingsNames <- rownames(loadings)
  }
  dataM<-geneMatchR(data1=data, data2=loadings, data1Names=dataNames, data2Names=loadingsNames, merge=FALSE)
  print(paste(as.character(dim(dataM[[2]])[1]),'row names matched between data and loadings'))
  print(paste('Updated dimension of data:',as.character(paste(dim(dataM[[2]]), collapse = ' '))))
  # do projection
  Design <- model.matrix(~0 + dataM[[1]])
  colnames(Design) <- colnames(dataM[[1]])
  projection <- lmFit(as.matrix(t(dataM[[2]])),Design)
  projectionPatterns <- t(projection$coefficients)
  projection.ts<-t(projection$coefficients/projection$stdev.unscaled/projection$sigma)

  #projection<-vglm(dataM$data2 ~ 0 + dataM$data1,family=family)
  #projectionPatterns<-coefvlm(projection,matrix.out=TRUE)

  #For VGAM
  #pval.matrix<-matrix(2*pnorm(-abs(summary(projection)@coef3[,3])),nrow=5,byrow=TRUE)

  #For limma
  pval.matrix<-2*pnorm(-abs(projection.ts))
  #colnames(pval.matrix)<-colnames(projectionPatterns)
  #rownames(pval.matrix)<-rownames(projectionPatterns)

  if(bootstrapPval){
  boots <- lapply(1:bootIter,function(x){
  rows <- sample(nrow(Design),nrow(Design),replace = T)
  projection <- lmFit(as.matrix(t(dataM[[2]][rows,])),Design[rows,])
  return(projection$coefficients)
    })
  bootPval <- compareBoots(projection$coefficients,boots)
  }

  if(full & bootstrapPval){
      #projectionFit <- list('projection'=projectionPatterns, 'fit'=projection,'pval'=pval.matrix)
      projectionFit <- list('projection'=projectionPatterns, 'pval'=pval.matrix, 'bootstrapPval' = bootPval)
      return(projectionFit)
  } else if(full){
      projectionFit <- list('projection'=projectionPatterns, 'pval'=pval.matrix)
      return(projectionFit)
  }
  else{return(projectionPatterns)}
})

#######################################################################################################################################
#' @import MatrixModels
#' @importFrom stats model.matrix
#' @param NP vector of integers indicating which columns of loadings object to use. The default of NP=NA will use entire matrix.
#' @param full logical indicating whether to return the full model solution. By default only the new pattern object is returned.
#' @param model Optional arguements to choose method for projection
#' @param family VGAM family function for model fitting (default: "gaussianff")
#' @rdname projectR-methods
#' @aliases projectR
setMethod("projectR",signature(data="dgCMatrix",loadings="matrix"),function(
  data, # a dataset to be projected onto
  loadings, # a matrix of continous values to be projected with unique rownames
  dataNames = NULL, # a vector with names of data rows
  loadingsNames = NULL, # a vector with names of loadings rows
  NP=NA, # vector of integers indicating which columns of loadings object to use. The default of NP=NA will use entire matrix.
  full=FALSE, # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  family="gaussianff"  # VGAM family function (default: "gaussianff")
  ){

  ifelse(!is.na(NP),loadings<-loadings[,NP],loadings<-loadings)
  #if(!is.na(NP)){loadings<-loadings[,NP]} was giving warning with subset of patterns
  #match genes in data sets
  if(is.null(dataNames)){
    dataNames <- rownames(data)
  }
  if(is.null(loadingsNames)){
    loadingsNames <- rownames(loadings)
  }
  dataM<-geneMatchR(data1=data, data2=loadings, data1Names=dataNames, data2Names=loadingsNames, merge=FALSE)
  print(paste(as.character(dim(dataM[[2]])[1]),'row names matched between data and loadings'))
  print(paste('Updated dimension of data:',as.character(paste(dim(dataM[[2]]), collapse = ' '))))
  # do projection
  Design <- model.matrix(~0 + dataM[[1]])
  colnames(Design) <- colnames(dataM[[1]])
  projection <- MatrixModels:::lm.fit.sparse(t(dataM[[2]]),dataM[[1]])
  projectionPatterns <- t(projection$coefficients)
  projection.ts<-t(projection$coefficients/projection$stdev.unscaled/projection$sigma)

  #projection<-vglm(dataM$data2 ~ 0 + dataM$data1,family=family)
  #projectionPatterns<-coefvlm(projection,matrix.out=TRUE)

  #For VGAM
  #pval.matrix<-matrix(2*pnorm(-abs(summary(projection)@coef3[,3])),nrow=5,byrow=TRUE)

  #For limma
  pval.matrix<-2*pnorm(-abs(projection.ts))
  #colnames(pval.matrix)<-colnames(projectionPatterns)
  #rownames(pval.matrix)<-rownames(projectionPatterns)

  if(full==TRUE){
      #projectionFit <- list('projection'=projectionPatterns, 'fit'=projection,'pval'=pval.matrix)
      projectionFit <- list('projection'=projectionPatterns, 'pval'=pval.matrix)
      return(projectionFit)
  }
  else{return(projectionPatterns)}
})


#######################################################################################################################################
#' @import limma
#' @importFrom NMF fcnnls
#' @examples
#' library("CoGAPS")
#' CR.RNAseq6l3c3t <- CoGAPS(p.RNAseq6l3c3t, params = new("CogapsParams",
#' nPatterns=5))
#' projectR(data=p.ESepiGen4c1l$mRNA.Seq,loadings=CR.RNAseq6l3c3t,
#' dataNames = map.ESepiGen4c1l[["GeneSymbols"]])
#'
#' @rdname projectR-methods
#' @aliases projectR
setMethod("projectR",signature(data="matrix",loadings="LinearEmbeddingMatrix"),function(
  data, # a dataset to be projected onto
  loadings, # a matrix of continous values to be projected with unique rownames
  dataNames = NULL, # a vector with names of data rows
  loadingsNames = NULL, # a vector with names of loadings rows
  NP=NA, # vector of integers indicating which columns of loadings object to use. The default of NP=NA will use entire matrix.
  full=FALSE, # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  model=NA, # optional arguements to choose method for projection
  family="gaussianff", # VGAM family function (default: "gaussianff")
  bootstrapPval=FALSE, # logical to indicate whether to generate p-values using bootstrap
  bootIter=1e3 # No of bootstrap iterations
  ){

  loadings<-loadings@featureLoadings
  ifelse(!is.na(NP),loadings<-loadings[,NP],loadings<-loadings)
  return(projectR(data,loadings = loadings,dataNames = dataNames, loadingsNames = loadingsNames,NP,full,bootstrapPval=bootstrapPval,bootIter=bootIter))

})

#######################################################################################################################################

#' @import limma
#' @importFrom stats var
#' @examples
#' pca.RNAseq6l3c3t<-prcomp(t(p.RNAseq6l3c3t))
#' pca.ESepiGen4c1l<-projectR(data=p.ESepiGen4c1l$mRNA.Seq,
#' loadings=pca.RNAseq6l3c3t, dataNames = map.ESepiGen4c1l[["GeneSymbols"]])
#'
#' @rdname projectR-methods
#' @aliases projectR
setMethod("projectR",signature(data="matrix",loadings="prcomp"),function(
  data, # a dataset to be projected onto
  loadings, # a matrix of continous values to be projected with unique rownames
  dataNames = NULL, # a vector with names of data rows
  loadingsNames = NULL, # a vector with names of loadings rows
  NP=NA, # vector of integers indicating which columns of loadings object to use. The default of NP=NA will use entire matrix.
  full=FALSE # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  ){

  loadings<-loadings$rotation
  ifelse(!is.na(NP),loadings<-loadings[,NP],loadings<-loadings)

  #match genes in data sets
  if(is.null(dataNames)){
    dataNames <- rownames(data)
  }
  if(is.null(loadingsNames)){
    loadingsNames <- rownames(loadings)
  }
  dataM<-geneMatchR(data1=data, data2=loadings, data1Names=dataNames, data2Names=loadingsNames, merge=FALSE)
  print(paste(as.character(dim(dataM[[2]])[1]),'row names matched between data and loadings'))
  print(paste('Updated dimension of data:',as.character(paste(dim(dataM[[2]]), collapse = ' '))))
  # do projection
  dat2P<-apply(dataM[[2]],1,function(x) x-mean(x))
  projectionPatterns<- dat2P %*% dataM[[1]] #head(X %*% PCA$rotation)

  if(full==TRUE){
  #calculate percent varience accoutned for by each PC in newdata
  #Eigenvalues<-eigen(cov(projectionPatterns))$values
  #PercentVariance<-round(Eigenvalues/sum(Eigenvalues) * 100, digits = 2)

  PercentVariance<-apply(projectionPatterns,2, function(x) 100*var(x)/sum(apply(projectionPatterns,2,var)))

    projectionFit <- list(t(projectionPatterns), PercentVariance) #also need to change this to transpose
    return(projectionFit)
  }
  else{return(t(projectionPatterns))}

})
#######################################################################################################################################

#' @examples
#' pca.RNAseq6l3c3t<-prcomp(t(p.RNAseq6l3c3t))
#' r.RNAseq6l3c3t<-rotatoR(1,1,-1,-1,pca.RNAseq6l3c3t$rotation[,1:2])
#' pca.ESepiGen4c1l<-projectR(data=p.ESepiGen4c1l$mRNA.Seq,
#' loadings=r.RNAseq6l3c3t, dataNames = map.ESepiGen4c1l[["GeneSymbols"]])
#'
#' @rdname projectR-methods
#' @aliases projectR

setMethod("projectR",signature(data="matrix",loadings="rotatoR"),function(
  data, # a dataset to be projected onto
  loadings, # a matrix of continous values to be projected with unique rownames
  dataNames = NULL, # a vector with names of data rows
  loadingsNames = NULL, # a vector with names of loadings rows
  NP=NA, # vector of integers indicating which columns of loadings object to use. The default of NP=NA will use entire matrix.
  full=FALSE # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  ){

  loadings <- loadings@rotatedM
  ifelse(!is.na(NP),loadings<-loadings[,NP],loadings<-loadings)

  #match genes in data sets
  if(is.null(dataNames)){
    dataNames <- rownames(data)
  }
  if(is.null(loadingsNames)){
    loadingsNames <- rownames(loadings)
  }
  dataM<-geneMatchR(data1=data, data2=loadings, data1Names=dataNames, data2Names=loadingsNames, merge=FALSE)
  print(paste(as.character(dim(dataM[[2]])[1]),'row names matched between data and loadings'))
  print(paste('Updated dimension of data:',as.character(paste(dim(dataM[[2]]), collapse = ' '))))
  # do projection
  dat2P<-apply(dataM[[2]],1,function(x) x-mean(x))
  projectionPatterns<- dat2P %*% dataM[[1]] #head(X %*% PCA$rotation)

  if(full==TRUE){
  #calculate percent varience accounted for by each PC in newdata
  Eigenvalues<-eigen(cov(projectionPatterns))$values
  PercentVariance<-round(Eigenvalues/sum(Eigenvalues) * 100, digits = 2)

  #PercentVariance<-apply(projectionPatterns,2, function(x) 100*var(x)/sum(apply(p2P,2,var)))

    projectionFit <- list(t(projectionPatterns), PercentVariance)
    return(projectionFit)
  }
  else{return(t(projectionPatterns))}

})

#######################################################################################################################################

#' @import limma
#' @examples
#' c.RNAseq6l3c3t<-correlateR(genes="T", dat=p.RNAseq6l3c3t, threshtype="N",
#' threshold=10, absR=TRUE)
#' cor.ESepiGen4c1l<-projectR(data=p.ESepiGen4c1l$mRNA.Seq, loadings=c.RNAseq6l3c3t,
#' NP="PositiveCOR", dataNames = map.ESepiGen4c1l[["GeneSymbols"]])
#'
#' @rdname projectR-methods
#' @aliases projectR
setMethod("projectR",signature(data="matrix",loadings="correlateR"),function(
  data, # a dataset to be projected onto
  loadings, # a matrix of continous values to be projected with unique rownames
  dataNames = NULL, # a vector with names of data rows
  loadingsNames = NULL, # a vector with names of loadings rows
  NP=NA, #can be used to select for "NegativeCOR" or "PositiveCOR" list from correlateR class obj containing both. By default is NA
  full=FALSE, # logical indicating whether to return the percent variance accounted for by each projected PC. By default only the new pattern object is returned.
  bootstrapPval=FALSE, # logical to indicate whether to generate p-values using bootstrap
  bootIter=1e3 # No of bootstrap iterations
  ){

  patterns <- loadings@corM
  if(!is.na(NP)){
    patterns<-as.matrix(patterns[[NP]])
    colnames(patterns) <- NP
  }
  else {
  patterns<-loadings
}
  #check length of patterns "PositiveCOR" and "NegativeCOR" or just positive
  if(length(patterns)==2){
    patterns <- do.call(rbind,patterns)
  }
  else{
    patterns <- as.matrix(patterns)
}
  return(projectR(data = data, loadings = patterns,dataNames = dataNames, loadingsNames = loadingsNames,  full = full,
    bootstrapPval = bootstrapPval, bootIter = bootIter))

})

#######################################################################################################################################

#' @param targetNumPatterns desired number of patterns with hclust
#' @param sourceData data used to create cluster object
#' @import limma
#' @import cluster
#' @importFrom stats cutree
#' @rdname projectR-methods
#' @aliases projectR
setMethod("projectR", signature(data="matrix", loadings="hclust"),
function(data, loadings, dataNames=NULL, loadingsNames=NULL, full=FALSE,
targetNumPatterns, sourceData,bootstrapPval=FALSE,bootIter=1000)
{
  cut <- cutree(loadings, k=targetNumPatterns)
  patterns <- matrix(0, nrow=nrow(sourceData), ncol=targetNumPatterns)
  rownames(patterns) <- rownames(sourceData)
    for(x in 1:targetNumPatterns)
    {
      patterns[cut==x,x] <- apply(sourceData[cut==x,], 1, cor, y=colMeans(sourceData[cut==x,]))
    }
  return(projectR(data, loadings=patterns, dataNames, loadingsNames, full = full))
})

#' @rdname projectR-methods
#' @aliases projectR
setMethod("projectR", signature(data="matrix", loadings="kmeans"),
function(data, loadings, dataNames=NULL, loadingsNames=NULL, full=FALSE, sourceData,bootstrapPval=FALSE,bootIter=1000)
{
  patterns <- matrix(0, nrow=nrow(sourceData), ncol=length(loadings$size))
  rownames(patterns) <- rownames(sourceData)
  for(x in 1:length(loadings$size))
  {
    patterns[loadings$cluster==x,x] <- apply(sourceData[loadings$cluster==x,], 1, cor, y=colMeans(sourceData[loadings$cluster==x,]))
  }
  return(projectR(data, loadings=patterns, dataNames= dataNames, full = full,bootstrapPval=bootstrapPval,bootIter=bootIter))
})

#########################################################################

#' @examples
#' library("projectR")
#' data(p.RNAseq6l3c3t)
#' nP<-3
#' kClust<-kmeans(t(p.RNAseq6l3c3t),centers=nP)
#' kpattern<-cluster2pattern(clusters = kClust, NP = nP, data = p.RNAseq6l3c3t)
#' p<-as.matrix(p.RNAseq6l3c3t)
#' projectR(p,kpattern)
#'
#' @rdname projectR-methods
#' @aliases projectR

setMethod("projectR", signature(data="matrix", loadings="cluster2pattern"),
function(data, loadings, dataNames=NULL, loadingsNames=NULL, full=FALSE, sourceData,bootstrapPval=FALSE,bootIter=1000)
{
  loadings = loadings@clusterMatrix
  # NA results from cor when sd is zero in some of the groups
  loadings[is.na(loadings)] <- 0
  return(projectR(data, loadings=loadings, dataNames= dataNames, full = full,bootstrapPval=bootstrapPval,bootIter=bootIter))
})

#########################################################################

compareBoots <- function(projection,boots){
mat <- sapply(1:nrow(projection),function(i){
  sapply(1:ncol(projection),function(j){
    val <- sapply(1:length(boots),function(x){
      return(boots[[x]][i,j])
    })
    valD <- ecdf(val)
    qt0 <- valD(0)
    if(qt0 < 0.5){
        return(2*qt0)
      } else {
        return(2*(1-qt0))
      }
  })
})
return(mat)
}
