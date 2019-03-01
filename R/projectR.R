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
  Patterns, # a matrix of continous values to be projected with unique rownames
  AnnotationObj=NA, # an annotation object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol", # the column of AnnotationData object corresponding to identifiers matching the type used for GeneWeights
  NP=NA, # vector of integers indicating which columns of Patterns object to use. The default of NP=NA will use entire matrix.
  full=FALSE, # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  family="gaussianff"  # VGAM family function (default: "gaussianff")
  ){

  ifelse(!is.na(NP),Patterns<-Patterns[,NP],Patterns<-Patterns)
  #if(!is.na(NP)){Patterns<-Patterns[,NP]} was giving warning with subset of patterns
  #match genes in data sets
  dataM<-geneMatchR(data1=data, AnnotationObj=AnnotationObj, IDcol=IDcol, data2=Patterns, merge=FALSE)
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
#' @param NP vector of integers indicating which columns of Patterns object to use. The default of NP = NA will use entire matrix.
#' @param full logical indicating whether to return the full model solution. By default only the new pattern object is returned.
#' @param model Optional arguements to choose method for projection
#' @param family VGAM family function for model fitting (default: "gaussianff")
#' @rdname projectR-methods
#' @aliases projectR
setMethod("projectR",signature(data="matrix",Patterns="matrix"),projectR.default)


#######################################################################################################################################
#' @import limma
#' @import stats
#' @importFrom NMF fcnnls

projectR.CogapsResult <- function(
  data, # a dataset to be projected onto
  AnnotationObj, # an annotation object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol", # the column of AnnotationData object corresponding to identifiers matching the type used for GeneWeights
  Patterns=NA, # a CogapsResult object
  NP=NA, # vector of integers indicating which columns of Patterns object to use. The default of NP=NA will use entire matrix.
  full=FALSE, # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  model=NA, # optional arguements to choose method for projection
  family="gaussianff" # VGAM family function (default: "gaussianff")
  ){

  Patterns<-Patterns@featureLoadings
  ifelse(!is.na(NP),Patterns<-Patterns[,NP],Patterns<-Patterns)

  #match genes in data sets
  dataM<-geneMatchR(data1=data, AnnotationObj=AnnotationObj, IDcol=IDcol, data2=Patterns, merge=FALSE)
  print(dim(dataM[[2]]))
  colnames(dataM[[1]]) <- paste('Pattern ',1:dim(dataM[[1]])[2],sep='') #make option to imput vector or change label

  # do projection
  Design <- model.matrix(~0 + dataM[[1]])
  colnames(Design) <- colnames(dataM[[1]])

  if(!is.na(model) && model=="NonNegative"){
    Projection <- fcnnls(Design,as.matrix(t(dataM[[2]])))
  }else{
    projection <- lmFit(as.matrix(t(dataM[[2]])),Design)
  }

  #projectionPatterns<-coefvlm(Projection$coefficients,matrix.out=TRUE)
  projectionPatterns <- t(projection$coefficients)
  projection.ts<-t(projection$coefficients/projection$stdev.unscaled/projection$sigma)
  pval.matrix<-2*pnorm(-abs(projection.ts))

  if(full==TRUE){
      projectionFit <- list('projection' = projectionPatterns, 'pval' = pval.matrix)
      return(projectionFit)
  }
  else{return(projectionPatterns)}
}

#' @examples
#' library("CoGAPS")
#' CR.RNAseq6l3c3t <- CoGAPS(p.RNAseq6l3c3t, params = new("CogapsParams",
#' nPatterns=5))
#' projectR(data=p.ESepiGen4c1l$mRNA.Seq,Patterns=CR.RNAseq6l3c3t,
#' AnnotationObj=map.ESepiGen4c1l,IDcol="GeneSymbols")
#'
#' @rdname projectR-methods
#' @aliases projectR
setMethod("projectR",signature(data="matrix",Patterns="CogapsResult"),projectR.CogapsResult)

#######################################################################################################################################

#' @import limma
#' @import stats
#' @importFrom NMF fcnnls

projectR.CoGAPS <- function(
  data, # a dataset to be projected onto
  Patterns, # a CoGAPS object
  AnnotationObj=NA, # an annotation object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol", # the column of AnnotationData object corresponding to identifiers matching the type used for GeneWeights
  NP=NA, # vector of integers indicating which columns of Patterns object to use. The default of NP=NA will use entire matrix.
  full=FALSE, # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  family="gaussianff" # VGAM family function (default: "gaussianff")
  ){

  if(is.null(dim(Patterns))){Patterns<-Patterns$Amean}
  ifelse(!is.na(NP),Patterns<-Patterns[,NP],Patterns<-Patterns)
  #match genes in data sets
  dataM<-geneMatchR(data1=data, AnnotationObj=AnnotationObj, IDcol=IDcol, data2=Patterns, merge=FALSE)
  print(dim(dataM[[2]]))
  colnames(dataM[[1]]) <- paste('Pattern ',1:dim(dataM[[1]])[2],sep='') #make option to imput vector or change label

  # do projection
  Design <- model.matrix(~0 + dataM[[1]])
  colnames(Design) <- colnames(dataM[[1]])

    Projection <- lmFit(as.matrix(t(dataM[[2]])),Design)
  #projectionPatterns<-coefvlm(Projection$coefficients,matrix.out=TRUE)
  projectionPatterns <- t(Projection$coefficients)
  projection.ts<-t(Projection$coefficients/Projection$stdev.unscaled/Projection$sigma)

  if(full==TRUE){
      projectionFit <- list(projectionPatterns, Projection)
      return(projectionFit)
  }
  else{return(projectionPatterns)}
}

#' @rdname projectR-methods
#' @aliases projectR
setMethod("projectR",signature(data="matrix",Patterns="CoGAPS"),projectR.CoGAPS)

#######################################################################################################################################

#' @import limma
#' @import cluster
#' @import stats

projectR.pclust <- function(
  data, # a dataset to be projected onto
  Patterns, # an Pclust object from the cluster2pattern function
  AnnotationObj=NA, # an annotation object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol", # the column of AnnotationData object corresponding to identifiers matching the type used for GeneWeights
  NP=NA, # number of desired patterns
  full=FALSE # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  ){

  Patterns <- Patterns@patterns
  ifelse(!is.na(NP),Patterns<-Patterns[,NP],Patterns<-Patterns)

  #match genes in data sets
  dataM<-geneMatchR(data1=data, AnnotationObj=AnnotationObj, IDcol=IDcol, data2=Patterns, merge=FALSE)
  print(dim(dataM[[2]]))
  colnames(dataM[[1]]) <- paste('Pattern ',1:dim(dataM[[1]])[2],sep='') #make option to imput vector or change label

  # do projection
  Design <- model.matrix(~0 + dataM[[1]])
  colnames(Design) <- colnames(dataM[[1]])
  Projection <- lmFit(as.matrix(t(dataM[[2]])),Design)
  projectionPatterns <- t(Projection$coefficients)
  projection.ts<-t(Projection$coefficients/Projection$stdev.unscaled/Projection$sigma)
  pval.matrix<-2*pnorm(-abs(projection.ts))
  if(full==TRUE){
      projectionFit <- list('projection'=projectionPatterns, 'pval'=pval.matrix)
      return(projectionFit)
  }
  else{return(projectionPatterns)}
}

#' @examples
#' k.RNAseq6l3c3t<-kmeans(p.RNAseq6l3c3t,22)
#' k.RNAseq6l3c3t<-cluster2pattern (clusters=k.RNAseq6l3c3t, NP=22, Data=p.RNAseq6l3c3t)
#' k.ESepiGen4c1l<-projectR(data=p.ESepiGen4c1l$mRNA.Seq, Patterns=k.RNAseq6l3c3t,
#' AnnotationObj=map.ESepiGen4c1l,IDcol="GeneSymbols")
#'
#' @rdname projectR-methods
#' @aliases projectR

setMethod("projectR",signature(data="matrix",Patterns="pclust"),projectR.pclust)
#######################################################################################################################################

#' @import limma
#' @import stats



projectR.prcomp <- function(
  data, # a dataset to be projected onto
  Patterns, # an prcomp object with a rotation matrix of genes by PCs
  AnnotationObj=NA, # an annotation object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol", # the column of AnnotationData object corresponding to identifiers matching the type used for GeneWeights
  NP=NA, # range of PCs to project. The default of NP=NA will use entire matrix.
  full=FALSE # logical indicating whether to return the percent variance accounted for by each projected PC. By default only the new pattern object is returned.
  ){

  Patterns<-Patterns$rotation
  ifelse(!is.na(NP),Patterns<-Patterns[,NP],Patterns<-Patterns)

  #match genes in data sets
  dataM<-geneMatchR(data1=data, AnnotationObj=AnnotationObj, IDcol=IDcol, data2=Patterns, merge=FALSE)
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
#' Patterns=pca.RNAseq6l3c3t,AnnotationObj=map.ESepiGen4c1l,IDcol="GeneSymbols")
#'
#' @rdname projectR-methods
#' @aliases projectR
setMethod("projectR",signature(data="matrix",Patterns="prcomp"),projectR.prcomp)
#######################################################################################################################################

#' @import stats


projectR.rotatoR <- function(
  data, # a dataset to be projected onto
  Patterns, # an prcomp object with a rotation matrix of genes by PCs
  AnnotationObj=NA, # an annotation object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol", # the column of AnnotationData object corresponding to identifiers matching the type used for GeneWeights
  NP=NA, # range of PCs to project. The default of NP=NA will use entire matrix.
  full=FALSE # logical indicating whether to return the percent variance accounted for by each projected PC. By default only the new pattern object is returned.
  ){

  Patterns <- Patterns@rotatedM
  ifelse(!is.na(NP),Patterns<-Patterns[,NP],Patterns<-Patterns)

  #match genes in data sets
  dataM<-geneMatchR(data1=data, AnnotationObj=AnnotationObj, IDcol=IDcol, data2=Patterns, merge=FALSE)
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
#' Patterns=r.RNAseq6l3c3t, AnnotationObj=map.ESepiGen4c1l, IDcol="GeneSymbols")
#'
#' @rdname projectR-methods
#' @aliases projectR

setMethod("projectR",signature(data="matrix",Patterns="rotatoR"),projectR.rotatoR)

#######################################################################################################################################

#' @import limma
#' @import stats

projectR.correlateR <- function(
  data, # a dataset to be projected onto
  Patterns, # an correlateR object of correlated genes
  AnnotationObj=NA, # an annotation object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol", # the column of AnnotationData object corresponding to identifiers matching the type used for GeneWeights
  NP=NA, #can be used to select for "NegativeCOR" or "PositiveCOR" list from correlateR class obj containing both. By default is NA
  full=FALSE # logical indicating whether to return the percent variance accounted for by each projected PC. By default only the new pattern object is returned.
  ){

  Patterns <- Patterns@corM
  if(!is.na(NP)){
    Patterns<-as.matrix(Patterns[[NP]])
    colnames(Patterns) <- NP
  }
  else {
  Patterns<-Patterns
}
  #check length of patterns "PositiveCOR" and "NegativeCOR" or just positive
  if(length(Patterns)==2){
    Patterns <- do.call(rbind,Patterns)
  }

  #match genes in data sets
  dataM<-geneMatchR(data1=data, AnnotationObj=AnnotationObj, IDcol=IDcol, data2=Patterns, merge=FALSE)
  print(dim(dataM[[2]]))

  # do projection
  Design <- model.matrix(~0 + dataM[[1]])
  colnames(Design) <- colnames(dataM[[1]])
  Projection <- lmFit(as.matrix(t(dataM[[2]])),Design)
  projectionPatterns <- t(Projection$coefficients)
  projection.ts<-t(Projection$coefficients/Projection$stdev.unscaled/Projection$sigma)
  pval.matrix<-2*pnorm(-abs(projection.ts))
    if(full==TRUE){
      projectionFit <- list('projection'=projectionPatterns, 'pval'=pval.matrix)
      return(projectionFit)
    }
    else{return(projectionPatterns)}
}
#' @examples
#' c.RNAseq6l3c3t<-correlateR(genes="T", dat=p.RNAseq6l3c3t, threshtype="N", 
#' threshold=10, absR=TRUE)
#' cor.ESepiGen4c1l<-projectR(data=p.ESepiGen4c1l$mRNA.Seq, Patterns=c.RNAseq6l3c3t, 
#' NP="PositiveCOR", AnnotationObj=map.ESepiGen4c1l, IDcol="GeneSymbols")
#'
#' @rdname projectR-methods
#' @aliases projectR
setMethod("projectR",signature(data="matrix",Patterns="correlateR"),projectR.correlateR)
#######################################################################################################################################

#' @import limma
#' @import stats

projectR.list <- function(
  data, # a dataset to be projected onto
  Patterns, # a CoGAPS object
  AnnotationObj=NA, # an annotation object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol", # the column of AnnotationData object corresponding to identifiers matching the type used for GeneWeights
  NP=NA, # vector of integers indicating which columns of Patterns object to use. The default of NP=NA will use entire matrix.
  full=FALSE # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  ){
  print(class(Patterns))
  print(Patterns)
  if("CoGAPS" %in% class(Patterns)){
    return(projectR.CoGAPS(data = data, AnnotationObj = AnnotationObj, IDcol = IDcol, Patterns = Patterns, NP = NP, full = full))
  }
  else{
    stop("Invalid object type Patterns. Patterns be from class matrix, pclust, CogapsResult, CoGAPS, correlateR, rotatoR or prcomp")
  }
}

#' @rdname projectR-methods
#' @aliases projectR
setMethod("projectR",signature(data="matrix",Patterns="list"),projectR.list)