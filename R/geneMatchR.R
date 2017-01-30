
geneMatchR<-function(
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
    uniEGids=unique(rownames(data1)[rownames(data)%in%rownames(data2)])
    rows1=match(uniEGids,rownames(data1))
    rnP<-rownames(data[rows1,])
  }
  rows2=match(uniEGids,rownames(data2))
  data <- as.matrix(data1)
  p2P <- as.matrix(data[rows1,])
  rownames(p2P) <- rnP
  As4P <- data2[rows2,]
  p4P <- p2P[match(rownames(p2P),rownames(As4P)),]
  dataM<-list("data1"=As4P,"data2"=p4P)
  if(merge){
    do.call(cbind,dataM)
    return()
  } else(return(dataM))
}
