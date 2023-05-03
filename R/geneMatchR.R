
#' Generic geneMatchR function
#' @rdname geneMatchR
#'  
#' @description Matches genes accross datasets
#' @param data1 a data matrix, typically genes by samples
#' @param data2 an amplitude matrix, typically genes by factors
#' @param data1Names rownames of data matrix, for eg genenames
#' @param data2Names rownames of amplitude matrix to be matched to rownames of datamatrix
#' @param merge logical indicating wether or not to merged data sets
#' @param ... Additional arguments to geneMatchR
#'
#' @return A list of genes (intersection) in both datasets. (if merge = TRUE, Also returns merged data.)
#' @export
#'
#' @examples
#' geneMatchR(data1=p.ESepiGen4c1l$mRNA.Seq,data2=p.RNAseq6l3c3t,
#' data1Names=map.ESepiGen4c1l[["GeneSymbols"]])
#' 

geneMatchR <- function(
  data1,# a dataset of genes by samples
  data2, # dataset with rownames to be matched
  data1Names = NULL, # rownames of data1
  data2Names = NULL, # rownames of data2
  merge=FALSE, # logical indicating wether or not to merged data sets
  ...
){

  if(!is.null(data1Names)){
    if(length(data1Names) != dim(data1)[1]){
      stop("Length of data1Names should be equal to number of rows of data1")
    }
    rownames(data1) <- data1Names
  }
  if(!is.null(data2Names)){
    if(length(data2Names) != dim(data2)[1]){
      stop("Length of data2Names should be equal to number of rows of data2")
    }
    rownames(data2) <- data2Names
  }
  data1Names <- rownames(data1)
  data2Names <- rownames(data2)
  if(nrow(data1) != 1){
    uniEGids = unique(data1Names[data1Names %in% data2Names])
    rows1 = match(uniEGids,data1Names)
    rnP <- rownames(data1[rows1,])
  }
  else{
    uniEGids=unique(rownames(data1)[rownames(data1)%in%rownames(data2)])
    rows1=match(uniEGids,rownames(data1))
    rnP<-rownames(data1[rows1,])  
  }
  rows2=match(uniEGids,rownames(data2))
  data <- data1
  p2P <- data[rows1,]
  rownames(p2P) <- rnP
  As4P <- data2[rows2,  ,drop = F]
  if(!is.matrix(As4P)){As4P <- as.matrix(As4P)}
  p4P <- p2P[match(rownames(p2P),rownames(As4P)),]
  dataM<-list("data1"=As4P,"data2"=p4P)
  if(merge){
    dataME<- do.call(cbind,dataM)
    return(dataME)
  } else(return(dataM))
}
