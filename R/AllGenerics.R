###############
# AllGenerics.R
###############


#' Generic projectR function
#'
#' @param data Target dataset into which you will prject
#' @param Patterns Patterns learned from source dataset
#' @param ... Additional arguments to projectR
#'
#' @return A matrix of sample weights for each input pattern. (if full=TRUE, full model solution is returned)
#' @export
#'
#' @examples
#'         projectR(data=p.RNAseq6l3c3t,Patterns=AP.RNAseq6l3c3t)
#'
setGeneric("projectR",function(data,Patterns,...) standardGeneric("projectR"))


#' Generic geneMatchR function
#'
#' @param data1 a dataset of genes by samples
#' @param AnnotionObj an annotion object for data1. If NA, the rownames of data will be used.
#' @param IDcol the column of AnnotionData object corresponding to identifiers matching the rownames of data2
#' @param data2 dataset with rownames to be matched
#' @param ... Additional arguments to geneMatchR
#'
#' @return A list of genes (intersection) in both datasets. (if merge = TRUE, Also returns merged data.)
#' @export
#'
#' @examples
#'       geneMatchR(data1=p.RNAseq6l3c3t,AnnotionObj=map.ESepiGen4c1l,
#'                  IDcol="GeneSymbols",data2=p.ESepiGen4c1l$mRNA.Seq)
setGeneric("geneMatchR",function(data1,AnnotionObj,IDcol="GeneSymbol",data2,...) standardGeneric("geneMatchR"))

#TODO: pclust is not an S4 class
#' @title cluster2pattern
#'
#' @description Function to make patterns of continuous weights from clusters.
#'
#' @param clusters an cluster object
#' @param NP number of desired patterns
#' @param Data data used to make clusters object
#' @param ... Additional arguments to geneMatchR
#' @return An object of class 'pclust' containing pattern weights corresponding for each cluster.
#' @export
#' @examples
#'  k.RNAseq6l3c3t<-kmeans(p.RNAseq6l3c3t,22)
#'  cluster2pattern(clusters=k.RNAseq6l3c3t,NP=22,Data=p.RNAseq6l3c3t)
setGeneric("cluster2pattern",function(clusters,NP,Data,...) standardGeneric("cluster2pattern"))
