#' @title Projection function - VGAM (Base)
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
#'    projectR(data=p.RNAseq6l3c3t,Patterns=AP.RNAseq6l3c3t)
#'
#' @import VGAM
#' @import stats
#' @import grDevices
#' @import methods
#' @import utils
#' @export


projectR <- function(
  data=NA,#a dataset to be projected onto
  AnnotionObj=NA,#an annotion object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol",#the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
  Patterns=NA,#a matrix of continous values with unique rownames to be projected
  NP=NA,#vector of integers indicating which columns of Patterns object to use. The default of NP=NA will use entire matrix.
  full=FALSE, # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  model=NA 
  ){
  UseMethod("projectR",Patterns)
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
#' @param model  # optional arguements to choose method for projection
#' @return A matrix of sample weights for each input pattern. (if full=TRUE, full model solution is returned)
#' @examples
#'    projectR(data=p.ESepiGen4c1l$mRNA.Seq,Patterns=AP.RNAseq6l3c3t,
#'                AnnotionObj=map.ESepiGen4c1l,IDcol="GeneSymbols")
#' @import VGAM
#' @import stats
#' @export


projectR.VGAM <- function(
  data=NA, # a dataset to be projected onto
  AnnotionObj=NA, # an annotion object for data. If NA, the rownames of data will be used.
  IDcol="GeneSymbol", # the column of AnnotionData object corresponding to identifiers matching the type used for GeneWeights
  Patterns=NA, # a matrix of continous values to be projected with unique rownames
  NP=NA, # vector of integers indicating which columns of Patterns object to use. The default of NP=NA will use entire matrix.
  full=FALSE, # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  model=NA 
  ){

  if(!is.na(NP)){Patterns<-Patterns[,NP]}
  #match genes in data sets
  dataM<-geneMatchR(data1=data, AnnotionObj=AnnotionObj, IDcol=IDcol, data2=Patterns, merge=FALSE)
  print(dim(dataM[[2]]))
  # do projection
  Design <- model.matrix(~0 + dataM[[1]])
  colnames(Design) <- colnames(dataM[[1]])
  Projection <- VGAM::vglm(as.matrix(t(dataM[[2]])),Design,family="gaussianff")
  projectionPatterns <- t(Projection$coefficients)
  if(full==TRUE){
      projectionFit <- list(projectionPatterns, Projection)
      return(projectionFit)
  }
  else{return(projectionPatterns)}
}

#######################################################################################################################################
