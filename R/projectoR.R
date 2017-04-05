
#' @title Projection function (Base)
#'
#' @description a function for the projection of new data into a previously defined feature space
#' @param data a dataset to be projected into the pattern space
#' @param AnnotionObj an annotion object for data. If NA the rownames of data will be used.
#' @param IDcol the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
#' @param Patterns a matrix of continous values with unique rownames to be projected
#' @param NP vector of integers indicating which columns of Patterns object to use. The default of NP = NA will use entire matrix.
#' @param full logical indicating whether to return the full model solution. By default only the new pattern object is returned.
#' @examples \dontrun{
#'    projectoR(data=D,Patterns=AP)
#'}
#' @import limma
#' @importFrom limma lmFit
#' @import stats
#' @export


projectoR <- function(
  data=NA,#a dataset to be projected onto
  AnnotionObj=NA,#an annotion object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol",#the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
  Patterns=NA,#a matrix of continous values with unique rownames to be projected
  NP=NA,#vector of integers indicating which columns of Patterns object to use. The default of NP=NA will use entire matrix.
  full=FALSE # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  ){
  UseMethod("projectoR",Patterns)
}


#######################################################################################################################################

#' @title Projection function (default)
#'
#' @description default version
#' @param data a dataset to be projected into the pattern space
#' @param AnnotionObj an annotion object for data. If NA the rownames of data will be used.
#' @param IDcol the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
#' @param Patterns a matrix of continous values to be projected with unique rownames
#' @param NP vector of integers indicating which columns of Patterns object to use. The default of NP = NA will use entire matrix.
#' @param full logical indicating whether to return the full model solution. By default only the new pattern object is returned.
#' @examples \dontrun{
#'    projectoR(data=D,Patterns=AP)
#'}

#' @import limma
#' @import stats
#' @export


projectoR.default <- function(
  data=NA, # a dataset to be projected onto
  AnnotionObj=NA, # an annotion object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol", # the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
  Patterns=NA, # a matrix of continous values to be projected with unique rownames
  NP=NA, # vector of integers indicating which columns of Patterns object to use. The default of NP=NA will use entire matrix.
  full=FALSE # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  ){

  if(!is.na(NP)){Patterns<-Patterns[,NP]}
  #match genes in data sets
  dataM<-geneMatchR(data1=data, AnnotionObj=AnnotionObj, IDcol=IDcol, data2=Patterns, merge=FALSE)
  print(dim(dataM[[2]]))
  # do projection
  Design <- model.matrix(~0 + dataM[[1]])
  colnames(Design) <- colnames(dataM[[1]])
  Projection <- lmFit(t(dataM[[2]]),Design)
  projectionPatterns <- t(Projection$coefficients)
  if(full==TRUE){
      projectionFit <- list(projectionPatterns, Projection)
      return(projectionFit)
  }
  else{return(projectionPatterns)}
}

#######################################################################################################################################

#' @title Projection function (CoGAPS)
#'
#' @description for use with object of class CoGAPS
#' @param data a dataset to be projected into the pattern space
#' @param AnnotionObj an annotion object for data. If NA the rownames of data will be used.
#' @param IDcol the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
#' @param Patterns a CoGAPS object
#' @param NP vector of integers indicating which columns of Patterns object to use. The default of NP = NA will use entire matrix.
#' @param full logical indicating whether to return the full model solution. By default only the new pattern object is returned.
#' @examples \dontrun{
#'    projectoR(data=D,Patterns=AP,PatternData=D)
#'}
#' @import limma
#' @import stats
#' @export

projectoR.CoGAPS <- function(
  data=NA, # a dataset to be projected onto
  AnnotionObj=NA, # an annotion object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol", # the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
  Patterns=NA, # a CoGAPS object
  NP=NA, # vector of integers indicating which columns of Patterns object to use. The default of NP=NA will use entire matrix.
  full=FALSE # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  ){

  if(is.null(dim(Patterns))){Patterns<-Patterns$Amean}
  if(!is.na(NP)){Patterns<-Patterns[,NP]}

  #match genes in data sets
  dataM<-geneMatchR(data1=data, AnnotionObj=AnnotionObj, IDcol=IDcol, data2=Patterns, merge=FALSE)
  print(dim(dataM[[2]]))
  colnames(dataM[[1]]) <- paste('Pattern ',1:dim(dataM[[1]])[2],sep='') #make option to imput vector or change label

  # do projection
  Design <- model.matrix(~0 + dataM[[1]])
  colnames(Design) <- colnames(dataM[[1]])
  Projection <- lmFit(t(dataM[[2]]),Design)
  projectionPatterns <- t(Projection$coefficients)
  if(full==TRUE){
      projectionFit <- list(projectionPatterns, Projection)
      return(projectionFit)
  }
  else{return(projectionPatterns)}
}


#######################################################################################################################################

#' @title Projection function (clustering)
#'
#' @description for use with object of class Pclust
#' @param data a dataset to be projected into the pattern space
#' @param AnnotionObj an annotion object for data. If NA the rownames of data will be used.
#' @param IDcol the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
#' @param Patterns an Pclust object from the cluster2pattern function
#' @param NP number of desired patterns
#' @param full logical indicating whether to return the full model solution. By default only the new pattern object is returned.
#' @examples \dontrun{
#'    projectoR(data=D,Patterns=cls,PatternData=D)
#'}
#' @import limma
#' @import stats
#' @export


projectoR.pclust <- function(
  data=NA, # a dataset to be projected onto
  AnnotionObj=NA, # an annotion object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol", # the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
  Patterns=NA, # an Pclust object from the cluster2pattern function
  NP=NA, # number of desired patterns
  full=FALSE # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  ){

  if(!is.na(NP)){Patterns<-Patterns[,NP]}

  #match genes in data sets
  dataM<-geneMatchR(data1=data, AnnotionObj=AnnotionObj, IDcol=IDcol, data2=Patterns, merge=FALSE)
  print(dim(dataM[[2]]))
  colnames(dataM[[1]]) <- paste('Pattern ',1:dim(dataM[[1]])[2],sep='') #make option to imput vector or change label

  # do projection
  Design <- model.matrix(~0 + dataM[[1]])
  colnames(Design) <- colnames(dataM[[1]])
  Projection <- lmFit(t(dataM[[2]]),Design)
  projectionPatterns <- t(Projection$coefficients)
  if(full==TRUE){
      projectionFit <- list(projectionPatterns, Projection)
      return(projectionFit)
  }
  else{return(projectionPatterns)}
}

#######################################################################################################################################

#' @title Projection function (PCA)
#'
#' @description for use with object of class prcomp
#' @param data a dataset to be projected into the pattern space
#' @param AnnotionObj an annotion object for data. If NA the rownames of data will be used.
#' @param IDcol the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
#' @param Patterns an prcomp object with a rotation matrix of genes by PCs
#' @param NP range of PCs to project. The default of NP = NA will use entire matrix.
#' @param full logical indicating whether to return the percent variance accounted for by each projected PC. By default only the new pattern object is returned.
#' @examples \dontrun{
#'   projectoR(data=D,Patterns=PCA,full=TRUE)
#'}
#' @import limma
#' @import stats
#' @export


projectoR.prcomp <- function(
  data=NA, # a dataset to be projected onto
  AnnotionObj=NA, # an annotion object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol", # the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
  Patterns=NA, # an prcomp object with a rotation matrix of genes by PCs
  NP=NA, # range of PCs to project. The default of NP=NA will use entire matrix.
  full=FALSE # logical indicating whether to return the percent variance accounted for by each projected PC. By default only the new pattern object is returned.
  ){

  Patterns<-Patterns$rotation
  if(!is.na(NP)){Patterns<-Patterns[,NP]}

  #match genes in data sets
  dataM<-geneMatchR(data1=data, AnnotionObj=AnnotionObj, IDcol=IDcol, data2=Patterns, merge=FALSE)
  print(dim(dataM[[2]]))

  # do projection
  dat2P<-apply(dataM[[2]],1,function(x) x-mean(x))
  projectionPatterns<- dat2P %*% dataM[[1]] #head(X %*% PCA$rotation)

  if(full==TRUE){
  #calculate percent varience accoutned for by each PC in newdata
  #Eigenvalues<-eigen(cov(projectionPatterns))$values
  #PercentVariance<-round(Eigenvalues/sum(Eigenvalues) * 100, digits = 2)
  
  PercentVariance<-apply(projectionPatterns,2, function(x) 100*var(x)/sum(apply(projectionPatterns,2,var)))  

    projectionFit <- list(projectionPatterns, PercentVariance)
    return(projectionFit)
  }
  else{return(projectionPatterns)}

}

#######################################################################################################################################

#' @title Projection function (rotatoR objects)
#'
#' @description for use with object of class rotatoR
#' @param data a dataset to be projected into the pattern space
#' @param AnnotionObj an annotion object for data. If NA the rownames of data will be used.
#' @param IDcol the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
#' @param Patterns an rotatoR object with a rotation matrix of genes by new PCs
#' @param NP range of PCs to project. The default of NP = NA will use entire matrix.
#' @param full logical indicating whether to return the percent variance accounted for by each projected PC. By default only the new pattern object is returned.
#' @examples \dontrun{
#'   projectoR(data=D,Patterns=rPCA,full=TRUE)
#'}
#' @import stats
#' @export


projectoR.rotatoR <- function(
  data=NA, # a dataset to be projected onto
  AnnotionObj=NA, # an annotion object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol", # the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
  Patterns=NA, # an prcomp object with a rotation matrix of genes by PCs
  NP=NA, # range of PCs to project. The default of NP=NA will use entire matrix.
  full=FALSE # logical indicating whether to return the percent variance accounted for by each projected PC. By default only the new pattern object is returned.
  ){

  if(!is.na(NP)){Patterns<-Patterns[,NP]}

  #match genes in data sets
  dataM<-geneMatchR(data1=data, AnnotionObj=AnnotionObj, IDcol=IDcol, data2=Patterns, merge=FALSE)
  print(dim(dataM[[2]]))

  # do projection
  dat2P<-apply(dataM[[2]],1,function(x) x-mean(x))
  projectionPatterns<- dat2P %*% dataM[[1]] #head(X %*% PCA$rotation)

  if(full==TRUE){
  #calculate percent varience accoutned for by each PC in newdata
  Eigenvalues<-eigen(cov(projectionPatterns))$values
  PercentVariance<-round(Eigenvalues/sum(Eigenvalues) * 100, digits = 2)
  
  #PercentVariance<-apply(projectionPatterns,2, function(x) 100*var(x)/sum(apply(p2P,2,var)))  

    projectionFit <- list(projectionPatterns, PercentVariance)
    return(projectionFit)
  }
  else{return(projectionPatterns)}

}



#######################################################################################################################################

#' @title Projection function (correlateR)
#'
#' @description for use with object of class corR
#' @param data a dataset to be projected into the pattern space
#' @param AnnotionObj an annotion object for data. If NA the rownames of data will be used.
#' @param IDcol the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
#' @param Patterns an correlateR object
#' @param NP the number of clusters
#' @param full logical indicating whether to return the full clustering information.  By default only the new pattern object is returned.
#' @examples \dontrun{
#'   projectoR(data=D,Patterns=PCA,full=TRUE)
#'}
#' @import limma
#' @import stats
#' @export


projectoR.correlateR <- function(
  data=NA, # a dataset to be projected onto
  AnnotionObj=NA, # an annotion object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol", # the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
  Patterns=NA, # an prcomp object with a rotation matrix of genes by PCs
  NP=NA, # range of PCs to project. The default of NP=NA will use entire matrix.
  full=FALSE # logical indicating whether to return the percent variance accounted for by each projected PC. By default only the new pattern object is returned.
  ){

  Patterns<-Patterns$rotation
  if(!is.na(NP)){Patterns<-Patterns[,NP]}

  #match genes in data sets
  dataM<-geneMatchR(data1=data, AnnotionObj=AnnotionObj, IDcol=IDcol, data2=Patterns, merge=FALSE)
  print(dim(dataM[[2]]))

  # do projection
  p2P<-apply(dataM[[2]],1,function(x) x-mean(x))
  projectionPatterns<- p2P %*% dataM[[1]] #head(X %*% PCA$rotation)

  #calculate percent varience accoutned for by each PC in newdata
  #Eigenvalues<-eigen(cov(projectionPatterns))$values
  #PercentVariance<-round(Eigenvalues/sum(Eigenvalues) * 100, digits = 2)
  PercentVariance<-apply(projectionPatterns,2, function(x) 100*var(x)/sum(apply(p2P,2,var)))

  if(full==TRUE){
      projectionFit <- list(projectionPatterns, PercentVariance)
      return(projectionFit)
  }
  else{return(projectionPatterns)}

}
