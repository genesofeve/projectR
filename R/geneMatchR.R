
geneMatchR.default<-function(
  data1=NA,# a dataset of genes by samples
  AnnotationObj=NA,#an Annotation object for data1. If NA, the rownames of data will be used.
  IDcol="GeneSymbol",#the column of AnnotationData object corresponding to identifiers matching the rownames of data2
  data2=NA, # dataset with rownames to be matched
  merge=FALSE # logical indicating wether or not to merged data sets
){
  if(length(AnnotationObj)!=1){
    uniEGids=unique(AnnotationObj[,IDcol][AnnotationObj[,IDcol]%in%rownames(data2)])
    rows1=match(uniEGids,AnnotationObj[,IDcol])
    rnP<-AnnotationObj[rows1,IDcol]
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
  if(class(As4P) != "matrix"){As4P <- as.matrix(As4P)}
  p4P <- p2P[match(rownames(p2P),rownames(As4P)),]
  dataM<-list("data1"=As4P,"data2"=p4P)
  if(merge){
    dataME<- do.call(cbind,dataM)
    return(dataME)
  } else(return(dataM))
}

#' @rdname geneMatchR-methods
#' @aliases geneMatchR
setMethod("geneMatchR",signature(data1="ANY",AnnotationObj="ANY",IDcol="ANY",data2="ANY"),geneMatchR.default)
