<<<<<<< HEAD
#' @importFrom stats hclust kmeans prcomp
setOldClass("kmeans")
setOldClass("hclust")
setOldClass("prcomp")

#' @importFrom CoGAPS CoGAPS
setOldClass("CoGAPS")

#' @importFrom limma lmFit


#' @title Projection function (Base)
#'
#' @description a function for the projection of new data into a previously defined feature space
#' @param data a dataset to be projected into the pattern space
#' @param AnnotionObj an annotion object for data. If NA the rownames of data will be used.
#' @param IDcol the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
#' @param Patterns a matrix of continous values with unique rownames to be projected
#' @param NP vector of integers indicating which columns of Patterns object to use. The default of NP = NA will use entire matrix.
#' @param full logical indicating whether to return the full model solution. By default only the new pattern object is returned.
#' @param model  # optional arguements to choose method for projection
#' @return A matrix of sample weights for each input pattern. (if full=TRUE, full model solution is returned)
#' @examples
#' #
#' # do not run, this will break the build.  Fix it instead.
#' # 
#' # projectR(data=p.RNAseq6l3c3t,Patterns=AP.RNAseq6l3c3t)
#'
#' @export
setGeneric("projectR", 
           function(data, AnnotionObj, IDcol, Patterns, NP, full, model=NA) 
           standardGeneric("projectR"))

=======
>>>>>>> 68e7dc7b88caed3ca6bf64468c75fe11aef8577f
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
#' @param model  # optional arguements to choose method for projection
#' @param family # VGAM family function for model fitting (default: "gaussianff")
#' @return A matrix of sample weights for each input pattern. (if full=TRUE, full model solution is returned)
#' @examples
<<<<<<< HEAD
#'    projectR(data=p.ESepiGen4c1l$mRNA.Seq,Patterns=AP.RNAseq6l3c3t$Amean,
#'                AnnotionObj=map.ESepiGen4c1l,IDcol="GeneSymbols")
<<<<<<< HEAD
=======
#' # 
#' # Do not run, will break the build.  Fix and then uncomment. 
#' # 
#' #   projectR(data=p.ESepiGen4c1l$mRNA.Seq,Patterns=AP.RNAseq6l3c3t,
#' #            AnnotionObj=map.ESepiGen4c1l,IDcol="GeneSymbols")
>>>>>>> master
#' @export

=======
# @import VGAM
#' @import limma
#' @import stats
>>>>>>> 68e7dc7b88caed3ca6bf64468c75fe11aef8577f

projectR.default <- function(
  data=NA, # a dataset to be projected onto
  AnnotionObj=NA, # an annotion object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol", # the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
  Patterns=NA, # a matrix of continous values to be projected with unique rownames
  NP=NA, # vector of integers indicating which columns of Patterns object to use. The default of NP=NA will use entire matrix.
  full=FALSE, # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  model=NA,
  family="gaussianff"
  ){

  if(!is.na(NP)){Patterns<-Patterns[,NP]}
  #match genes in data sets
  dataM<-geneMatchR(data1=data, AnnotionObj=AnnotionObj, IDcol=IDcol, data2=Patterns, merge=FALSE)
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

setMethod("projectR",signature(data="matrix",Patterns="matrix"),projectR.default)


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
#' @param model  # optional arguements to choose method for projection
#' @param family # VGAM family function for model fitting (default: "gaussianff")
#' @return A matrix of sample weights for each input pattern. (if full=TRUE, full model solution is returned)
#' @examples
<<<<<<< HEAD
#'    projectR(data=p.ESepiGen4c1l$mRNA.Seq,Patterns=AP.RNAseq6l3c3t,
#'                AnnotionObj=map.ESepiGen4c1l,IDcol="GeneSymbols",model="NonNegative")
#' @import VGAM
=======
#' # 
#' #   Do not run, this will break the build.  Fix the generics instead.
#' #  
#' #   projectR(data=p.ESepiGen4c1l$mRNA.Seq,Patterns=AP.RNAseq6l3c3t,
#' #               AnnotionObj=map.ESepiGen4c1l,IDcol="GeneSymbols")
#' @import limma
>>>>>>> master
#' @import stats
#' @import NMF

projectR.CoGAPS <- function(
  data=NA, # a dataset to be projected onto
  AnnotionObj=NA, # an annotion object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol", # the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
  Patterns=NA, # a CoGAPS object
  NP=NA, # vector of integers indicating which columns of Patterns object to use. The default of NP=NA will use entire matrix.
  full=FALSE, # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  model=NA, # optional arguements to choose method for projection
  family="gaussianff" # VGAM family function (default: "gaussianff")
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

  if(model=="NonNegative"){
<<<<<<< HEAD
    Projection <- fcnnls(Design,as.matrix(t(dataM[[2]])))
  }else{
    Projection <- lmFit(as.matrix(t(dataM[[2]])),Design)
=======
    projection <- fcnnls(Design,as.matrix(t(dataM[[2]])))
  } else{
    projection<-vglm(dataM$data2 ~ 0 + dataM$data1,family=family)
>>>>>>> 68e7dc7b88caed3ca6bf64468c75fe11aef8577f
  }
  projectionPatterns<-coefvlm(projection,matrix.out=TRUE)
  if(full==TRUE){
      projectionFit <- list(projectionPatterns, projection)
      return(projectionFit)
  }
  else{return(projectionPatterns)}
}

setMethod("projectR",signature(data="matrix",Patterns="list"),projectR.CoGAPS)

#######################################################################################################################################
#' @title Projection function (LEM)
#'
#' @description for use with object of class LinearEmbeddingMatrix (from 'SingleCellExperiment' package)
#' @param data a dataset to be projected into the pattern space
#' @param AnnotionObj an annotion object for data. If NA the rownames of data will be used.
#' @param IDcol the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
#' @param Patterns a LinearEmbeddingMatrix object
#' @param NP vector of integers indicating which row(s) of  object to use. The default of NP = NA will use entire matrix.
#' @param full logical indicating whether to return the full model solution. By default only the new pattern object is returned.
#' @param model  # optional arguements to choose method for projection
#' @return A matrix of sample weights for each input pattern. (if full=TRUE, full model solution is returned)
#' @examples
#'    projectR(data=p.ESepiGen4c1l$mRNA.Seq,Patterns=AP.RNAseq6l3c3t,
#'                AnnotionObj=map.ESepiGen4c1l,IDcol="GeneSymbols",model="NonNegative") #Update
#' @import limma
#' @import SingleCellExperiment
#' @importClassesFrom SingleCellExperiment LinearEmbeddingMatrix
#' @import stats
#' @import NMF

projectR.LEM<- function(
  data=NA, # a dataset to be projected onto
  AnnotionObj=NA, # an annotion object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol", # the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
  Patterns=NA, # a LinearEmbeddingMatrix object
  NP=NA, # vector of integers indicating which columns of Patterns object to use. The default of NP=NA will use entire matrix.
  full=FALSE, # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  model=NA # optional arguements to choose method for projection
){

  #if(is.null(dim(Patterns))){Patterns<-Patterns$Amean}
  #if(!is.na(NP)){Patterns<-Patterns[,NP]}

  #match genes in data sets
  dataM<-geneMatchR(data1=data, AnnotionObj=AnnotionObj, IDcol=IDcol, data2=featureLoadings(Patterns)[,NP], merge=FALSE)
  print(dim(dataM[[2]]))
  colnames(dataM[[1]]) <- paste('Pattern ',1:dim(dataM[[1]])[2],sep='') #make option to imput vector or change label

  # do projection
  Design <- model.matrix(~0 + dataM[[1]])
  colnames(Design) <- colnames(dataM[[1]])

  if(model=="NonNegative"){
    projection <- fcnnls(Design,as.matrix(t(dataM[[2]])))
  } else{
    projection<-vglm(dataM$data2 ~ 0 + dataM$data1,family=family)
    projectionPatterns<-coefvlm(projection,matrix.out=TRUE)
  }

  if(full==TRUE){
    projectionFit <- list(projectionPatterns, projection)
    return(projectionFit)
  }
  else{return(projectionPatterns)}
}

setMethod("projectR",signature(data="matrix",Patterns="list"),projectR.CoGAPS)



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
#' @param model  # optional arguements to choose method for projection
#' @return A matrix of sample weights for each input pattern. (if full=TRUE, full model solution is returned)
#' @examples
#' # 
#' #   Do not run, this will break the build.  Fix the generics instead.
#' #  
#' # k.RNAseq6l3c3t<-kmeans(p.RNAseq6l3c3t,22)
#' # k.RNAseq6l3c3t<-cluster2pattern(clusters=k.RNAseq6l3c3t,NP=22,Data=p.RNAseq6l3c3t)
#' # k.ESepiGen4c1l<-projectR(data=p.ESepiGen4c1l$mRNA.Seq,Patterns=k.RNAseq6l3c3t,
#' #                             AnnotionObj=map.ESepiGen4c1l,IDcol="GeneSymbols")
#'
#' @import limma
#' @import cluster
#' @import stats



projectR.pclust <- function(
  data=NA, # a dataset to be projected onto
  AnnotionObj=NA, # an annotion object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol", # the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
  Patterns=NA, # an Pclust object from the cluster2pattern function
  NP=NA, # number of desired patterns
  full=FALSE, # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  model=NA
  ){

  if(!is.na(NP)){Patterns<-Patterns[,NP]}

  #match genes in data sets
  dataM<-geneMatchR(data1=data, AnnotionObj=AnnotionObj, IDcol=IDcol, data2=Patterns, merge=FALSE)
  print(dim(dataM[[2]]))
  colnames(dataM[[1]]) <- paste('Pattern ',1:dim(dataM[[1]])[2],sep='') #make option to imput vector or change label

  # do projection
  Design <- model.matrix(~0 + dataM[[1]])
  colnames(Design) <- colnames(dataM[[1]])
  Projection <- lmFit(as.matrix(t(dataM[[2]])),Design)
  projectionPatterns <- t(Projection$coefficients)
  if(full==TRUE){
      projectionFit <- list(projectionPatterns, Projection)
      return(projectionFit)
  }
  else{return(projectionPatterns)}
}

setMethod("projectR",signature(data="matrix",Patterns="pclust"),projectR.pclust)
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
#' @param model  # optional arguements to choose method for projection
#' @return A matrix of sample weights for each input pattern. (if full=TRUE, full model solution is returned)
#' @examples
#' # 
#' #   Do not run, this will break the build.  Fix the generics instead.
#' #  
#' #  pca.RNAseq6l3c3t<-prcomp(t(p.RNAseq6l3c3t))
#' #  pca.ESepiGen4c1l<-projectR(data=p.ESepiGen4c1l$mRNA.Seq,Patterns=pca.RNAseq6l3c3t,
#' #                                AnnotionObj=map.ESepiGen4c1l,IDcol="GeneSymbols")
#'
#' @import limma
#' @import stats
#' @import MASS



projectR.prcomp <- function(
  data=NA, # a dataset to be projected onto
  AnnotionObj=NA, # an annotion object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol", # the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
  Patterns=NA, # an prcomp object with a rotation matrix of genes by PCs
  NP=NA, # range of PCs to project. The default of NP=NA will use entire matrix.
  full=FALSE, # logical indicating whether to return the percent variance accounted for by each projected PC. By default only the new pattern object is returned.
  model=NA
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

setMethod("projectR",signature(data="matrix",Patterns="prcomp"),projectR.prcomp)
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
#' @param model  # optional arguements to choose method for projection
#' @return A matrix of sample weights for each input pattern. (if full=TRUE, full model solution is returned)
#' @examples
#' # 
#' #   Do not run, this will break the build.  Fix the generics instead.
#' #  
#' # pca.RNAseq6l3c3t<-prcomp(t(p.RNAseq6l3c3t))
#' # r.RNAseq6l3c3t<-rotatoR(1,1,-1,-1,pca.RNAseq6l3c3t$x[,1:2])
#' #  pca.ESepiGen4c1l<-projectR(data=p.ESepiGen4c1l$mRNA.Seq,Patterns=r.RNAseq6l3c3t,
#' #                         AnnotionObj=map.ESepiGen4c1l,IDcol="GeneSymbols")
#'
#' @import stats


projectR.rotatoR <- function(
  data=NA, # a dataset to be projected onto
  AnnotionObj=NA, # an annotion object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol", # the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
  Patterns=NA, # an prcomp object with a rotation matrix of genes by PCs
  NP=NA, # range of PCs to project. The default of NP=NA will use entire matrix.
  full=FALSE, # logical indicating whether to return the percent variance accounted for by each projected PC. By default only the new pattern object is returned.
  model=NA
  ){

  if(!is.na(NP)){Patterns<-Patterns[,NP]}

  #match genes in data sets
  dataM<-geneMatchR(data1=data, AnnotionObj=AnnotionObj, IDcol=IDcol, data2=Patterns, merge=FALSE)
  print(dim(dataM[[2]]))

  # do projection
  dat2P<-apply(dataM[[2]],1,function(x) x-mean(x))
  projectionPatterns<- dat2P %*% dataM[[1]] #head(X %*% PCA$rotation)

  if(full==TRUE){
  #calculate percent varience accounted for by each PC in newdata
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
#' @param NP can be used to select for "NegativeCOR" or "PositiveCOR" list from correlateR class obj containing both. By default is NA
#' @param full logical indicating whether to return the full clustering information.  By default only the new pattern object is returned.
#' @param model  # optional arguements to choose method for projection
#' @return A matrix of sample weights for each input pattern. (if full=TRUE, full model solution is returned)
#' @examples
#' # 
#' #   Do not run, this will break the build.  Fix the generics instead.
#' #  
#' #  c.RNAseq6l3c3t<-correlateR(genes="T", dat=p.RNAseq6l3c3t, threshtype="N", threshold=10, absR=TRUE)
#' #  cor.ESepiGen4c1l<-projectR(data=p.ESepiGen4c1l$mRNA.Seq,Patterns=c.RNAseq6l3c3t,NP="PositiveCOR",
#' #                                   AnnotionObj=map.ESepiGen4c1l,IDcol="GeneSymbols")
#'
#' @import limma
#' @import stats

projectR.correlateR <- function(
  data=NA, # a dataset to be projected onto
  AnnotionObj=NA, # an annotion object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol", # the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
  Patterns=NA, # an correlateR object of correlated genes
  NP=NA, #can be used to select for "NegativeCOR" or "PositiveCOR" list from correlateR class obj containing both. By default is NA
  full=FALSE, # logical indicating whether to return the percent variance accounted for by each projected PC. By default only the new pattern object is returned.
  model=NA
  ){

  if(!is.na(NP)){Patterns<-Patterns[[NP]]}

  #check length of patterns "PositiveCOR" and "NegativeCOR" or just positive
  if(length(Patterns)==2){
    do.call(rbind,Patterns)
  }

  #match genes in data sets
  dataM<-geneMatchR(data1=data, AnnotionObj=AnnotionObj, IDcol=IDcol, data2=Patterns, merge=FALSE)
  print(dim(dataM[[2]]))

  # do projection
  Design <- model.matrix(~0 + dataM[[1]])
  colnames(Design) <- colnames(dataM[[1]])
  Projection <- lmFit(as.matrix(t(dataM[[2]])),Design)
  projectionPatterns <- t(Projection$coefficients)
    if(full==TRUE){
      projectionFit <- list(projectionPatterns, Projection)
      return(projectionFit)
    }
    else{return(projectionPatterns)}
}


