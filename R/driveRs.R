
#' @title <Projection driveRs function >
#'
#' @description <full description>
#' @param data a dataset to be projected onto
#' @param projectedP an prcomp object with a rotation matrix of genes by PCs
#' @param design design matrix describing experiment
#' @param contrasts a matrix of continous values with unique rownames to be projected
#' @param NP vector of integers indicating which columns of Patterns object to use. The default of NP = NA will use entire matrix.
#' @param full logical indicating whether to return the full model solution. By default only the new pattern object is returned.
#' @param ... additional arguements to topTable
#' @export
#' @seealso topTable, lmFit, eBayes
#' @import limma
#' @importFrom MASS ginv
#' @examples
#'  driveRs(data=NAprojectedP=NA, design=NA)
#'

driveRs <- function(
  data=NA, # a dataset to be projected onto
  projectedP=NA, # an prcomp object with a rotation matrix of genes by PCs
  design=NA,
  contrasts=NA,
  NP=NA, # range of PCs to project. The default of NP=NA will use entire matrix.
  full=FALSE, # logical indicating whether to return the percent variance accounted for by each projected PC. By default only the new pattern object is returned.
  ...){

  if(!is.na(NP)){
    x<-projectedP[,NP]
  }else{x<-projectedP}

  EY<- x %*% ginv(t(x)%*%x) %*% t(x) %*% t(data)

  # dif express Exp
  fit <- lmFit(t(EY), design)
  fit <- eBayes(fit)

  driveRs<-topTableF(fit, number=dim(data)[1], ...)
  driveRs<-driveRs[order(driveRs[,1],decreasing=TRUE),]

  if(full==TRUE){
      driveRs <- list("driveRsFit"=fit, "driveRs"=driveRs)
      return(driveRs)
  } else{return(driveRs)}
}
