
#' @importFrom stats kmeans hclust
setOldClass("kmeans")
setOldClass("hclust")

#######################################################################################################################################

cluster2pattern.kmeans<- function(
  clusters, # a kmeans object
  Data=NA # data used to make clusters object
  ){

  nD<-length(clusters$size)
  nG<-dim(Data)[1]
  tempP<-matrix(data=rep(0,nD*nG),nrow = nG,ncol =nD)
  rownames(tempP)<-rownames(Data)
  #for(x in 1:nD) {tempP[Patterns$cluster==x,x]<-rowMeans(Data[Patterns$cluster==x,])}
  for(x in 1:nD) {tempP[clusters$cluster==x,x]<-apply(Data[clusters$cluster==x,],1,cor,y=colMeans(Data[clusters$cluster==x,]))}
  Patterns<-tempP
  pclustObj <- new("pclustKmeans",patterns = Patterns,kmeans = clusters)
  return(pclustObj)
}

#' @rdname cluster2pattern-methods
#' @aliases cluster2pattern
setMethod("cluster2pattern",signature(clusters="kmeans"),cluster2pattern.kmeans)

#######################################################################################################################################

cluster2pattern.hclust<-function(
  clusters, # an hclust object
  NP=NA, # number of desired patterns
  Data=NA # data used to make hclust object
  ){

#  if(is.na(PatternData)){stop("Data used to make hclust object must also be provided.")}
  cut=cutree(clusters,k=NP)
  nG<-dim(Data)[1]
  tempP<-matrix(data=rep(0,NP*nG),nrow = nG,ncol=NP)
  rownames(tempP)<-rownames(Data)
  #for(x in 1:NP) {tempP[cut==x,x]<-rowMeans(Data[cut==x,])}
  for(x in 1:NP) {tempP[cut==x,x]<-apply(Data[cut==x,],1,cor,y=colMeans(Data[cut==x,]))}
  Patterns<-tempP
  pclustObj <- new("pclustHclust",patterns = Patterns,hclust = clusters)
  return(pclustObj)
}

#' @rdname cluster2pattern-methods
#' @aliases cluster2pattern
setMethod("cluster2pattern",signature(clusters="hclust"),cluster2pattern.hclust)

