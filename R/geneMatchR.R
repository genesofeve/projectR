#' @title geneMatchR
#'
#' @description Matches genes accross datasets.
#' @param data1 a dataset of genes by samples
#' @param AnnotionObj an annotion object for data1. If NA, the rownames of data will be used.
#' @param IDcol the column of AnnotionData object corresponding to identifiers matching the rownames of data2
#' @param data2 dataset with rownames to be matched
#' @param merge logical indicating wether or not to merged data sets
#' @return A list of genes (intersection) in both datasets. (if merge = TRUE, Also returns merged data.)
#' @examples
#'  geneMatchR(data1=p.RNAseq6l3c3t,AnnotionObj=map.ESepiGen4c1l,
#'                  IDcol="GeneSymbols",data2=p.ESepiGen4c1l$mRNA.Seq)

geneMatchR.default<-function(
  data1=NA,# a dataset of genes by samples
  AnnotionObj=NA,#an annotion object for data1. If NA, the rownames of data will be used.
  IDcol="GeneSymbol",#the column of AnnotionData object corresponding to identifiers matching the rownames of data2
  data2=NA, # dataset with rownames to be matched
  merge=FALSE # logical indicating wether or not to merged data sets
){
  if(length(AnnotionObj)!=1){
    uniEGids=unique(AnnotionObj[,IDcol][AnnotionObj[,IDcol]%in%rownames(data2)])
    rows1=match(uniEGids,AnnotionObj[,IDcol])
    rnP<-AnnotionObj[rows1,IDcol]
  } else {
    uniEGids=unique(rownames(data1)[rownames(data1)%in%rownames(data2)])
    rows1=match(uniEGids,rownames(data1))
    rnP<-rownames(data1[rows1,])
  }
  rows2=match(uniEGids,rownames(data2))
  data <- data1
  p2P <- data[rows1,]
  rownames(p2P) <- rnP
  As4P <- data2[rows2,]
  p4P <- p2P[match(rownames(p2P),rownames(As4P)),]
  dataM<-list("data1"=As4P,"data2"=p4P)
  if(merge){
    dataME<- do.call(cbind,dataM)
    return(dataME)
  } else(return(dataM))
}

#setMethod("geneMatchR",signature(data1="data.frame",AnnotionObj="data.frame",IDcol="character",data2="matrix"),geneMatchR.default)
setMethod("geneMatchR",signature(data1="ANY",AnnotionObj="ANY",IDcol="ANY",data2="ANY"),geneMatchR.default)
