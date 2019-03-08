#' @importFrom stats hclust kmeans prcomp
setOldClass("kmeans")
setOldClass("hclust")
setOldClass("prcomp")
#' @importFrom CoGAPS CoGAPS
setOldClass("CoGAPS")
#' @importFrom limma lmFit


#######################################################################################################################################
#' @import limma
#' @import stats


projectR.default <- function(
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
  dataM<-geneMatchR(data1=data, data2=loadings, data1Names=dataNames, data2Names=data2Names, merge=FALSE)
  print(dim(dataM[[2]]))
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

  if(full==TRUE){
      #projectionFit <- list('projection'=projectionPatterns, 'fit'=projection,'pval'=pval.matrix)
      projectionFit <- list('projection'=projectionPatterns, 'pval'=pval.matrix)
      return(projectionFit)
  }
  else{return(projectionPatterns)}
}

#' @param AnnotationObj an annotation object for data. If NA (default) the rownames of data will be used.
#' @param IDcol the column of AnnotationData object corresponding to identifiers matching the type used for GeneWeights
#' @param NP vector of integers indicating which columns of loadings object to use. The default of NP = NA will use entire matrix.
#' @param full logical indicating whether to return the full model solution. By default only the new pattern object is returned.
#' @param model Optional arguements to choose method for projection
#' @param family VGAM family function for model fitting (default: "gaussianff")
#' @rdname projectR-methods
#' @aliases projectR
setMethod("projectR",signature(data="matrix",loadings="matrix"),projectR.default)


#######################################################################################################################################
#' @import limma
#' @import stats
#' @importFrom NMF fcnnls

projectR.LEM <- function(
  data, # a dataset to be projected onto
  loadings, # a matrix of continous values to be projected with unique rownames
  dataNames, # a vector with names of data rows
  loadingsNames, # a vector with names of loadings rows
  NP=NA, # vector of integers indicating which columns of loadings object to use. The default of NP=NA will use entire matrix.
  full=FALSE, # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  model=NA, # optional arguements to choose method for projection
  family="gaussianff" # VGAM family function (default: "gaussianff")
  ){

  loadings<-loadings@featureLoadings
  ifelse(!is.na(NP),loadings<-loadings[,NP],loadings<-loadings)
  return(projectR(data,loadings = loadings,dataNames = dataNames, loadingsNames = loadingsNames,NP,full))

}

#' @examples
#' library("CoGAPS")
#' CR.RNAseq6l3c3t <- CoGAPS(p.RNAseq6l3c3t, params = new("CogapsParams",
#' nPatterns=5))
#' projectR(data=p.ESepiGen4c1l$mRNA.Seq,loadings=CR.RNAseq6l3c3t,
#' AnnotationObj=map.ESepiGen4c1l,IDcol="GeneSymbols")
#'
#' @rdname projectR-methods
#' @aliases projectR
setMethod("projectR",signature(data="matrix",loadings="LinearEmbeddingMatrix"),projectR.LEM)

#######################################################################################################################################

#' @import limma
#' @import cluster
#' @import stats

projectR.pclust <- function(
  data, # a dataset to be projected onto
  loadings, # a matrix of continous values to be projected with unique rownames
  dataNames, # a vector with names of data rows
  loadingsNames, # a vector with names of loadings rows
  NP=NA, # number of desired patterns
  full=FALSE # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  ){

  loadings <- loadings@patterns
  ifelse(!is.na(NP),loadings<-loadings[,NP],loadings<-loadings)
  return(projectR(data,loadings = loadings,dataNames = dataNames, loadingsNames = loadingsNames,NP,full))

}

#' @examples
#' k.RNAseq6l3c3t<-kmeans(p.RNAseq6l3c3t,22)
#' k.RNAseq6l3c3t<-cluster2pattern (clusters=k.RNAseq6l3c3t, NP=22, Data=p.RNAseq6l3c3t)
#' k.ESepiGen4c1l<-projectR(data=p.ESepiGen4c1l$mRNA.Seq, loadings=k.RNAseq6l3c3t,
#' AnnotationObj=map.ESepiGen4c1l,IDcol="GeneSymbols")
#'
#' @rdname projectR-methods
#' @aliases projectR

setMethod("projectR",signature(data="matrix",loadings="pclust"),projectR.pclust)
#######################################################################################################################################

#' @import limma
#' @import stats



projectR.prcomp <- function(
  data, # a dataset to be projected onto
  loadings, # a matrix of continous values to be projected with unique rownames
  dataNames, # a vector with names of data rows
  loadingsNames, # a vector with names of loadings rows
  NP=NA, # range of PCs to project. The default of NP=NA will use entire matrix.
  full=FALSE # logical indicating whether to return the percent variance accounted for by each projected PC. By default only the new pattern object is returned.
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
  print(dim(dataM[[2]]))

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

}

#' @examples
#' pca.RNAseq6l3c3t<-prcomp(t(p.RNAseq6l3c3t))
#' pca.ESepiGen4c1l<-projectR(data=p.ESepiGen4c1l$mRNA.Seq, 
#' loadings=pca.RNAseq6l3c3t,AnnotationObj=map.ESepiGen4c1l,IDcol="GeneSymbols")
#'
#' @rdname projectR-methods
#' @aliases projectR
setMethod("projectR",signature(data="matrix",loadings="prcomp"),projectR.prcomp)
#######################################################################################################################################

#' @import stats


projectR.rotatoR <- function(
  data, # a dataset to be projected onto
  loadings, # a matrix of continous values to be projected with unique rownames
  dataNames, # a vector with names of data rows
  loadingsNames, # a vector with names of loadings rows
  NP=NA, # range of PCs to project. The default of NP=NA will use entire matrix.
  full=FALSE # logical indicating whether to return the percent variance accounted for by each projected PC. By default only the new pattern object is returned.
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
  print(dim(dataM[[2]]))

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

}

#' @examples
#' pca.RNAseq6l3c3t<-prcomp(t(p.RNAseq6l3c3t))
#' r.RNAseq6l3c3t<-rotatoR(1,1,-1,-1,pca.RNAseq6l3c3t$rotation[,1:2])
#' pca.ESepiGen4c1l<-projectR(data=p.ESepiGen4c1l$mRNA.Seq, 
#' loadings=r.RNAseq6l3c3t, AnnotationObj=map.ESepiGen4c1l, IDcol="GeneSymbols")
#'
#' @rdname projectR-methods
#' @aliases projectR

setMethod("projectR",signature(data="matrix",loadings="rotatoR"),projectR.rotatoR)

#######################################################################################################################################

#' @import limma
#' @import stats

projectR.correlateR <- function(
  data, # a dataset to be projected onto
  loadings, # a matrix of continous values to be projected with unique rownames
  dataNames, # a vector with names of data rows
  loadingsNames, # a vector with names of loadings rows
  NP=NA, #can be used to select for "NegativeCOR" or "PositiveCOR" list from correlateR class obj containing both. By default is NA
  full=FALSE # logical indicating whether to return the percent variance accounted for by each projected PC. By default only the new pattern object is returned.
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
  print(patterns)
  print(class(patterns))
  return(projectR(data = data, loadings = patterns,dataNames = dataNames, loadingsNames = loadingsNames, IDcol = IDcol, full = full ))
 
}
#' @examples
#' c.RNAseq6l3c3t<-correlateR(genes="T", dat=p.RNAseq6l3c3t, threshtype="N", 
#' threshold=10, absR=TRUE)
#' cor.ESepiGen4c1l<-projectR(data=p.ESepiGen4c1l$mRNA.Seq, loadings=c.RNAseq6l3c3t, 
#' NP="PositiveCOR", AnnotationObj=map.ESepiGen4c1l, IDcol="GeneSymbols")
#'
#' @rdname projectR-methods
#' @aliases projectR
setMethod("projectR",signature(data="matrix",loadings="correlateR"),projectR.correlateR)

#######################################################################################################################################

#' @import limma
#' @import cluster
#' @import stats
#' @rdname projectR-methods
#' @aliases projectR
setMethod("projectR", signature(data="matrix", loadings="hclust"),
function(data, loadings, dataNames=NULL, loadingsNames=NULL, full=FALSE,
targetNumPatterns, sourceData)
{
  cut <- cutree(loadings, k=targetNumPatterns)
  patterns <- matrix(0, nrow=nrow(sourceData), ncol=targetNumPatterns)
  for(x in 1:targetNumPatterns)
  {
    patterns[cut==x,x] <- apply(Data[cut==x,], 1, cor, y=colMeans(sourceData[cut==x,]))
  }
  return(projectR(data, loadings=patterns, dataNames, loadingsNames, full))
})

#' @rdname projectR-methods
#' @aliases projectR
setMethod("projectR", signature(data="matrix", loadings="kmeans"),
function(data, loadings, dataNames=NULL, loadingsNames=NULL, full=FALSE,
targetNumPatterns, sourceData)
{
  patterns <- matrix(0, nrow=nrow(sourceData), ncol=length(loadings$size))
  for(x in 1:length(loadings$size))
  {
    patterns[loadings$cluster==x,x] <- apply(sourceData[loadings$cluster==x,], 1, cor, y=colMeans(sourceData[loadings$cluster==x,]))
  }
  return(projectR(data, loadings=patterns, dataNames, loadingsNames, full))
})
