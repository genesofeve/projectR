
#' @title <Projection function (Base)>
#'
#' @description <full description>
#' @param data a dataset to be projected onto
#' @param AnnotionObj an annotion object for data. If NA the rownames of data will be used.
#' @param IDcol the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
#' @param Weights a matrix of continous values with unique rownames to be projected
#' @param NP vector of integers indicating which columns of Patterns object to use. The default of NP = NA will use entire matrix.
#' @param PatternData data used to make Patterns
#' @param full logical indicating whether to return the full model solution. By default only the new pattern object is returned.
#' @param ...
#' @export
#' @seealso
#' @return
#' @examples \dontrun{
#'    projectR(data=D,Patterns=AP)
#'}
#' @import limma


projectR <- function(
  data=NA,#a dataset to be projected onto
  AnnotionObj=NA,#an annotion object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol",#the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
  Patterns=NA,#a matrix of continous values with unique rownames to be projected
  NP=NA,#vector of integers indicating which columns of Patterns object to use. The default of NP=NA will use entire matrix.
  PatternData=NA, # data used to make Patterns
  full=FALSE,# logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  ...){
  UseMethod("projectR",Patterns)
}


#######################################################################################################################################

#' @title <Projection function (default)>
#'
#' @description <default version>
#' @param data a dataset to be projected onto
#' @param AnnotionObj an annotion object for data. If NA the rownames of data will be used.
#' @param IDcol the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
#' @param Patterns a matrix of continous values to be projected with unique rownames
#' @param NP vector of integers indicating which columns of Patterns object to use. The default of NP = NA will use entire matrix.
#' @param full
#' @export
#' @seealso
#' @examples \dontrun{
#'    projectR(data=D,Patterns=AP)
#'}
#' @return
#' @import limma


projectR.default <- function(
  data=NA, # a dataset to be projected onto
  AnnotionObj=NA, # an annotion object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol", # the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
  Patterns=NA, # a matrix of continous values to be projected with unique rownames
  NP=NA, # vector of integers indicating which columns of Patterns object to use. The default of NP=NA will use entire matrix.
  full=FALSE, # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  ...){

  if(!is.na(NP)){Patterns<-Patterns[,NP]}

  if(!is.na(AnnotionObj)){
    uniEGids=unique(AnnotionObj[,IDcol][AnnotionObj[,IDcol]%in%rownames(Patterns)])
    rows1=match(uniEGids,AnnotionObj[,IDcol])
    rnP<-AnnotionObj[rows1,IDcol]
  } else {
    uniEGids=unique(rownames(data)[rownames(data)%in%rownames(Patterns)])
    rows1=match(uniEGids,rownames(data))
    rnP<-rownames(data[rows1,])
  }

  rows2=match(uniEGids,rownames(Patterns))
  data <- as.matrix(data)
  p2P <- as.matrix(data[rows1,])
  rownames(p2P) <- rnP
  As4P <- Patterns[rows2,]
  colnames(As4P) <- paste('Pattern ',1:dim(As4P)[2],sep='') #make option to imput vector or change label
  p4P <- p2P[match(rownames(p2P),rownames(As4P)),]
  print(dim(p4P))
  Design <- model.matrix(~0 + As4P)
  colnames(Design) <- colnames(As4P)
  Projection <- lmFit(t(p4P),Design)
  projectionPatterns <- t(Projection$coefficients)
  if(full==TRUE){
      projectionFit <- list(projectionPatterns, Projection)
      return(projectionFit)
  }
  else{return(projectionPatterns)}
}

#######################################################################################################################################

#' @title <Projection function (CoGAPS NMF)>
#'
#' @description <for use with object of class CoGAPS>
#' @param data a dataset to be projected onto
#' @param AnnotionObj an annotion object for data. If NA the rownames of data will be used.
#' @param IDcol the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
#' @param Patterns a CoGAPS object
#' @param NP vector of integers indicating which columns of Patterns object to use. The default of NP = NA will use entire matrix.
#' @param full logical indicating whether to return the full model solution. By default only the new pattern object is returned.
#' @param ...
#' @export
#' @seealso
#' @return
#' @examples \dontrun{
#'    projectR(data=D,Patterns=AP,PatternData=D)
#'}
#' @import limma

projectR.CoGAPS <- function(
  data=NA, # a dataset to be projected onto
  AnnotionObj=NA, # an annotion object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol", # the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
  Patterns=NA, # a CoGAPS object
  NP=NA, # vector of integers indicating which columns of Patterns object to use. The default of NP=NA will use entire matrix.
  full=FALSE, # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  ...){

  if(is.null(dim(Patterns))){Patterns<-Patterns$Amean}
  if(!is.na(NP)){Patterns<-Patterns[,NP]}

  if(!is.na(AnnotionObj)){
    uniEGids=unique(AnnotionObj[,IDcol][AnnotionObj[,IDcol]%in%rownames(Patterns)])
    rows1=match(uniEGids,AnnotionObj[,IDcol])
    rnP<-AnnotionObj[rows1,IDcol]
  } else {
  uniEGids=unique(rownames(data)[rownames(data)%in%rownames(Patterns)])
  rows1=match(uniEGids,rownames(data))
  rnP<-rownames(data[rows1,])
  }

  rows2=match(uniEGids,rownames(Patterns))
  data <- as.matrix(data)
  p2P <- as.matrix(data[rows1,])
  rownames(p2P) <- rnP
  As4P <- Patterns[rows2,]
  colnames(As4P) <- paste('Pattern ',1:dim(As4P)[2],sep='') #make option to imput vector or change label
  p4P <- p2P[match(rownames(p2P),rownames(As4P)),]
  print(dim(p4P))
  Design <- model.matrix(~0 + As4P)
  colnames(Design) <- colnames(As4P)
  Projection <- lmFit(t(p4P),Design)
  projectionPatterns <- t(Projection$coefficients)
  if(full==TRUE){
      projectionFit <- list(projectionPatterns, Projection)
      return(projectionFit)
  }
  else{return(projectionPatterns)}
}

#######################################################################################################################################

#' @title <Projection function (kmeans clustering)>
#'
#' @description <for use with object of class kmeans>
#' @param data a dataset to be projected onto
#' @param AnnotionObj an annotion object for data. If NA the rownames of data will be used.
#' @param IDcol the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
#' @param Patterns a kmeans object
#' @param NP vector of integers indicating which columns of Patterns object to use. The default of NP = NA will use entire matrix.
#' @param PatternData data used to make kmeans object
#' @param full logical indicating whether to return the full model solution. By default only the new pattern object is returned.
#' @param ...
#' @export
#' @seealso
#' @return
#' @examples \dontrun{
#'    projectR(data=D,Patterns=cls,PatternData=D)
#'}
#' @import limma

projectR.kmeans <- function(
  data=NA, # a dataset to be projected onto
  AnnotionObj=NA, # an annotion object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol", # the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
  Patterns=NA, # a kmeans object
  NP=NA, # vector of integers indicating which columns of Patterns object to use. The default of NP=NA will use entire matrix.
  PatternData=NA, # data used to make kmeans object
  full=FALSE, # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  ...){

  nD<-length(Patterns$size)
  nG<-dim(PatternData)[1]
  tempP<-matrix(data=rep(0,nD*nG),nrow = nG,ncol =nD)
  rownames(tempP)<-rownames(PatternData)
  #for(x in 1:nD) {tempP[Patterns$cluster==x,x]<-rowMeans(PatternData[Patterns$cluster==x,])}
  for(x in 1:nD) {tempP[Patterns$cluster==x,x]<-apply(D[Patterns$cluster==x,],1,cor,y=colMeans(D[Patterns$cluster==x,]))}
  Patterns<-tempP

  if(!is.na(NP)){Patterns<-Patterns[,NP]}

  if(!is.na(AnnotionObj)){
    uniEGids=unique(AnnotionObj[,IDcol][AnnotionObj[,IDcol]%in%rownames(Patterns)])
    rows1=match(uniEGids,AnnotionObj[,IDcol])
    rnP<-AnnotionObj[rows1,IDcol]
  } else {
    uniEGids=unique(rownames(data)[rownames(data)%in%rownames(Patterns)])
    rows1=match(uniEGids,rownames(data))
    rnP<-rownames(data[rows1,])
  }

  rows2=match(uniEGids,rownames(Patterns))
  data <- as.matrix(data)
  p2P <- as.matrix(data[rows1,])
  rownames(p2P) <- rnP
  As4P <- Patterns[rows2,]
  colnames(As4P) <- paste('Pattern ',1:dim(As4P)[2],sep='') #make option to imput vector or change label
  p4P <- p2P[match(rownames(p2P),rownames(As4P)),]
  print(dim(p4P))
  Design <- model.matrix(~0 + As4P)
  colnames(Design) <- colnames(As4P)
  Projection <- lmFit(t(p4P),Design)
  projectionPatterns <- t(Projection$coefficients)
  if(full==TRUE){
      projectionFit <- list(projectionPatterns, Projection)
      return(projectionFit)
  }
  else{return(projectionPatterns)}
}

#######################################################################################################################################

#' @title <Projection function (hierachical clustering)>
#'
#' @description <for use with object of class hclust>
#' @param data a dataset to be projected onto
#' @param AnnotionObj an annotion object for data. If NA the rownames of data will be used.
#' @param IDcol the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
#' @param Patterns an hclust object
#' @param NP number of desired patterns
#' @param PatternData data used to make hclust object
#' @param full logical indicating whether to return the full model solution. By default only the new pattern object is returned.
#' @param ...
#' @export
#' @seealso
#' @return
#' @examples \dontrun{
#'    projectR(data=D,Patterns=cls,PatternData=D)
#'}
#' @import limma


projectR.hclust <- function(
  data=NA, # a dataset to be projected onto
  AnnotionObj=NA, # an annotion object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol", # the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
  Patterns=NA, # an hclust object
  NP=NA, # number of desired patterns
  PatternData=NA, # data used to make hclust object
  full=FALSE, # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  ...){

#  if(is.na(PatternData)){stop("Data used to make hclust object must also be provided.")}
  cut=cutree(Patterns,k=NP)
  nG<-dim(PatternData)[1]
  tempP<-matrix(data=rep(0,NP*nG),nrow = nG,ncol=NP)
  rownames(tempP)<-rownames(PatternData)
  #for(x in 1:NP) {tempP[cut==x,x]<-rowMeans(PatternData[cut==x,])}
  for(x in 1:NP) {tempP[cut==x,x]<-apply(D[cut==x,],1,cor,y=colMeans(D[cut==x,]))}
  Patterns<-tempP

  if(!is.na(AnnotionObj)){
    uniEGids=unique(AnnotionObj[,IDcol][AnnotionObj[,IDcol]%in%rownames(Patterns)])
    rows1=match(uniEGids,AnnotionObj[,IDcol])
    rnP<-AnnotionObj[rows1,IDcol]
  } else {
    uniEGids=unique(rownames(data)[rownames(data)%in%rownames(Patterns)])
    rows1=match(uniEGids,rownames(data))
    rnP<-rownames(data[rows1,])
  }

  rows2=match(uniEGids,rownames(Patterns))
  data <- as.matrix(data)
  p2P <- as.matrix(data[rows1,])
  rownames(p2P) <- rnP
  As4P <- Patterns[rows2,]
  colnames(As4P) <- paste('Pattern ',1:dim(As4P)[2],sep='') #make option to imput vector or change label
  p4P <- p2P[match(rownames(p2P),rownames(As4P)),]
  print(dim(p4P))
  Design <- model.matrix(~0 + As4P)
  colnames(Design) <- colnames(As4P)
  Projection <- lmFit(t(p4P),Design)
  projectionPatterns <- t(Projection$coefficients)
  if(full==TRUE){
      projectionFit <- list(projectionPatterns, Projection)
      return(projectionFit)
  }
  else{return(projectionPatterns)}
}

#######################################################################################################################################

#' @title <Projection function (PCA)>
#'
#' @description <for use with object of class prcomp>
#' @param data a dataset to be projected onto
#' @param AnnotionObj an annotion object for data. If NA the rownames of data will be used.
#' @param IDcol the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
#' @param Patterns an prcomp object with a rotation matrix of genes by PCs
#' @param NP range of PCs to project. The default of NP = NA will use entire matrix.
#' @param full logical indicating whether to return the percent variance accounted for by each projected PC. By default only the new pattern object is returned.
#' @param ...
#' @export
#' @seealso
#' @return
#' @examples \dontrun{
#'   projectR(data=D,Patterns=PCA,full=TRUE)
#'}
#' @import limma


projectR.prcomp <- function(
  data=NA, # a dataset to be projected onto
  AnnotionObj=NA, # an annotion object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol", # the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
  Patterns=NA, # an prcomp object with a rotation matrix of genes by PCs
  NP=NA, # range of PCs to project. The default of NP=NA will use entire matrix.
  full=FALSE, # logical indicating whether to return the percent variance accounted for by each projected PC. By default only the new pattern object is returned.
  ...){

  # to use in Carlo Version, make flag?
  #old.centers<-Patterns$center

  Patterns<-Patterns$rotation
  if(!is.na(NP)){Patterns<-Patterns[,NP]}

  if(!is.na(AnnotionObj)){
    uniEGids=unique(AnnotionObj[,IDcol][AnnotionObj[,IDcol]%in%rownames(Patterns)])
    rows1=match(uniEGids,AnnotionObj[,IDcol])
    rnP<-AnnotionObj[rows1,IDcol]
  } else {
    uniEGids=unique(rownames(data)[rownames(data)%in%rownames(Patterns)])
    rows1=match(uniEGids,rownames(data))
    rnP<-rownames(data[rows1,])
  }

  rows2=match(uniEGids,rownames(Patterns))
  p2P <- as.matrix(data[rows1,])
  rownames(p2P) <- rnP
  As4P <- Patterns[rows2,]
  p2P<-apply(p2P,1,function(x) x-mean(x))

  projectionPatterns<- p2P %*% As4P #head(X %*% PCA$rotation)

  #calculate percent varience accoutned for by each PC in newdata
  Eigenvalues<-eigen(cov(projectionPatterns),only.values=TRUE)
  PercentVariance<-round(Eigenvalues/sum(Eigenvalues) * 100, digits = 2)
  #PercentVariance<-apply(projectionPatterns,2, function(x) 100*var(x)/sum(apply(p2P,2,var)))

  if(full==TRUE){
      projectionFit <- list(projectionPatterns, PercentVariance)
      return(projectionFit)
  }
  else{return(projectionPatterns)}

}
