<<<<<<< HEAD
#' @importFrom stats hclust kmeans
setOldClass("kmeans")
setOldClass("hclust")

#' @title cluster2pattern
#'
#' @description Function to make patterns of continuous weights from clusters.
#' @param clusters an cluster object
#' @param NP number of desired patterns
#' @param Data data used to make clusters object
#' @return An object of class 'pclust' containing pattern weights corresponding for each cluster.
#' @export
#' @examples 
#'  k.RNAseq6l3c3t<-kmeans(p.RNAseq6l3c3t,22)
#'  cluster2pattern(clusters=k.RNAseq6l3c3t,NP=22,Data=p.RNAseq6l3c3t)  
#'

setGeneric("cluster2pattern", function(clusters, NP, Data) standardGeneric("cluster2pattern"))

=======
>>>>>>> 68e7dc7b88caed3ca6bf64468c75fe11aef8577f
#' @title cluster2pattern (kmeans)
#'
#' @description Function to make patterns of continuous weights from kmeans clusters.
#' @param clusters an kmeans cluster object
#' @param NP number of desired patterns
#' @param Data data used to make clusters object
#' @return An object of class 'pclust' containing pattern weights corresponding for each cluster.
#'
#' @examples
#'  k.RNAseq6l3c3t<-kmeans(p.RNAseq6l3c3t,22)
#'  cluster2pattern(clusters=k.RNAseq6l3c3t,NP=22,Data=p.RNAseq6l3c3t)
#'

setMethod("cluster2pattern", signature(clusters="kmeans"), function(
  clusters, # a kmeans object
  NP=NA, # number of desired patterns
  Data=NA # data used to make clusters object
  ){

  nD<-length(clusters$size)
  nG<-dim(Data)[1]
  tempP<-matrix(data=rep(0,nD*nG),nrow = nG,ncol =nD)
  rownames(tempP)<-rownames(Data)
  #for(x in 1:nD) {tempP[Patterns$cluster==x,x]<-rowMeans(Data[Patterns$cluster==x,])}
  for(x in 1:nD) {tempP[clusters$cluster==x,x]<-apply(Data[clusters$cluster==x,],1,cor,y=colMeans(Data[clusters$cluster==x,]))}
  Patterns<-tempP
  class(Patterns)<-append(class(Patterns),"pclust") # Can't/shouldn't do this in S4
  return(Patterns)
})

setMethod("cluster2pattern",signature(clusters="kmeans"),cluster2pattern.kmeans)

####################################
#' @title cluster2pattern (hclust)
#'
#' @description Function to make patterns of continuous weights from hierarchical clusters.
#' @param clusters an hclust object
#' @param NP number of desired patterns
#' @param Data data used to make clusters object
#' @return An object of class 'pclust' containing pattern weights corresponding for each cluster.
#' @examples
#'  h.RNAseq6l3c3t<-hclust(as.dist(1-(cor(t(p.RNAseq6l3c3t),use="pairwise.complete.obs"))))
#'  cluster2pattern(clusters=h.RNAseq6l3c3t,NP=22,Data=p.RNAseq6l3c3t)
#'

setMethod("cluster2pattern", signature(clusters="hclust"), function(
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
  class(Patterns)<-append(class(Patterns),"pclust") # Can't/shouldn't do this in S4
  return(Patterns)
<<<<<<< HEAD
})
=======
}

setMethod("cluster2pattern",signature(clusters="hclust"),cluster2pattern.hclust)
>>>>>>> 68e7dc7b88caed3ca6bf64468c75fe11aef8577f
