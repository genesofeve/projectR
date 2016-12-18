
#' @title <Projection driveRs function (Base)>
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
#' @import limma
#' @examples \dontrun{
#'    project(data=D,Patterns=AP)
#'}

driveRs <- function(
  data=NA,#a dataset to be projected onto
  AnnotionObj=NA,#an annotion object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol",#the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
  Patterns=NA,#a matrix of continous values with unique rownames to be projected
  NP=NA,#vector of integers indicating which columns of Patterns object to use. The default of NP=NA will use entire matrix.
  PatternData=NA, # data used to make Patterns
  full=FALSE,# logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  ...){
  UseMethod("driveRs",Patterns)
}


#######################################################################################################################################

#' @title <Projection driveRs function (default)>
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
#' @return
#' @import limma

driveRs.default <- function(
  data=NA, # the projected data set
  AnnotionObj=NA, # an annotion object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol", # the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
  Patterns=NA, # a projectoR object
  NP=NA, # vector of integers indicating which columns of Patterns object to use. The default of NP=NA will use entire matrix.
  full=FALSE, # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  ...){

  if(length(Patterns)>1Patterns<-Patterns$projectionPatterns
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

#' @title <Projection driveRs function (CoGAPS NMF)>
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
#' @import limma
#' @examples \dontrun{
#'    projectR(data=D,Patterns=AP)
#'}

driveRs.CoGAPS <- function(
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

#' @title <Projection driveRs function (kmeans clustering)>
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
#' @import limma
#' @examples \dontrun{
#'    project(data=D,Patterns=cls,PatternData=D)
#'}

driveRs.kmeans <- function(
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

#' @title <Projection driveRs function (hierachical clustering)>
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
#' @import limma

driveRs.hclust <- function(
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

#' @title <Projection driveRs function (PCA)>
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
#' @import limma
#' @import MASS
#' @examples \dontrun{
#'   driveRs(data=D,Patterns=PCA,full=TRUE)
#'}

p.BMP4mesoTcdx2.log2

driveRs.prcomp <- function(
  data=NA, # a dataset to be projected onto
  projectedP=NA, # an prcomp object with a rotation matrix of genes by PCs
  design=NA,
  constrasts=NA,
  NP=NA, # range of PCs to project. The default of NP=NA will use entire matrix.
  full=FALSE, # logical indicating whether to return the percent variance accounted for by each projected PC. By default only the new pattern object is returned.
  ...){

  require(MASS)
  if(!is.na(NP)){x<-projectedP[,NP]}
  EY<- x %*% ginv(t(x)%*%x) %*% t(x) %*% t(data)

  # dif express Exp
  require(limma)
  fit <- lmFit(t(EY), design)
  fit <- eBayes(fit)

  # makecontrast
  cont.dif <- makeContrasts(constrasts,levels=design)
  fit2 <- contrasts.fit(fit, cont.dif)
  fit2 <- eBayes(fit2)
  driverRs<-topTableF(fit2, n=dim(data)[1], ...)
  driverRs<-driverRs[order(driverRs[,1],decreasing=TRUE),]

  if(full==TRUE){
      driverRs <- list("driverRsFit"=fit2, "driverRs"=driverRs)
      return(driverRs)
  }
  else{return(driverRs)}

}
