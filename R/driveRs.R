
#' @title <Projection driveRs function >
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
#' @import MASS
#' @examples \dontrun{
#'    project(data=D,Patterns=AP)
#'}

driveRs <- function(
  data=NA, # a dataset to be projected onto
  projectedP=NA, # an prcomp object with a rotation matrix of genes by PCs
  design=NA,
  contrasts=NA,
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
  cont.dif <- makeContrasts(contrasts=contrasts,levels=design)
  fit2 <- contrasts.fit(fit, cont.dif)
  fit2 <- eBayes(fit2)

  driveRs<-topTableF(fit2, n=dim(data)[1], ...)
  driveRs<-driveRs[order(driveRs[,1],decreasing=TRUE),]

  if(full==TRUE){
      driveRs <- list("driveRsFit"=fit2, "driveRs"=driveRs)
      return(driveRs)
  } else{return(driveRs)}
}
