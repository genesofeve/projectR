




#######################################################################################################################################
#' @importFrom stats kmeans
setOldClass("kmeans")
.cluster2pattern_kmeans<- function(
  clusters, # a kmeans object
  data # data used to make clusters object
  ){

  nD<-length(clusters$size)
  nG<-dim(data)[1]
  tempP<-matrix(data=rep(0,nD*nG),nrow = nG,ncol =nD)
  rownames(tempP)<-rownames(data)
  #for(x in 1:nD) {tempP[Patterns$cluster==x,x]<-rowMeans(data[Patterns$cluster==x,])}
  for(x in 1:nD) {tempP[clusters$cluster==x,x]<-apply(data[clusters$cluster==x,],1,cor,y=colMeans(data[clusters$cluster==x,]))}
  Patterns<-tempP
  Patterns <- new("cluster2pattern",clusterMatrix = Patterns)
  return(Patterns)
}

#' @rdname cluster2pattern-methods
#' @aliases cluster2pattern
setMethod("cluster2pattern",signature(clusters="kmeans"),.cluster2pattern_kmeans)

#######################################################################################################################################
#' @importFrom stats hclust

setOldClass("hclust")
.cluster2pattern_hclust<-function(
  clusters, # an hclust object
  NP, # number of desired patterns
  data=NA # data used to make hclust object
  ){

#  if(is.na(Patterndata)){stop("data used to make hclust object must also be provided.")}
  cut=cutree(clusters,k=NP)
  nG<-dim(data)[1]
  tempP<-matrix(data=rep(0,NP*nG),nrow = nG,ncol=NP)
  rownames(tempP)<-rownames(data)
  #for(x in 1:NP) {tempP[cut==x,x]<-rowMeans(data[cut==x,])}
  for(x in 1:NP) {tempP[cut==x,x]<-apply(data[cut==x,],1,cor,y=colMeans(data[cut==x,]))}
  Patterns<-tempP
  Patterns <- new("cluster2pattern",clusterMatrix = Patterns)
  return(Patterns)
}

#' @rdname cluster2pattern-methods
#' @aliases cluster2pattern
setMethod("cluster2pattern",signature(clusters="hclust"),.cluster2pattern_hclust)

#######################################################################################################################################

