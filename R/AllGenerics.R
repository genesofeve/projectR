###############
# AllGenerics.R
###############


#' Generic projectR function
#'
#' @param data Target dataset into which you will prject
#' @param Patterns Patterns learned from source dataset
#' @param ...
#'
#' @return A matrix of sample weights for each input pattern. (if full=TRUE, full model solution is returned)
#' @export
#'
#' @examples
#'         projectR(data=p.RNAseq6l3c3t,Patterns=AP.RNAseq6l3c3t)
#'
setGeneric("projectR",function(data,Patterns,...) standardGeneric("projectR"))
