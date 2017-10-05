#' @title clusterPlotR (base)
#'
#' @description plotting function for clustering objects
#' @param cData data used to get clusters
#' @param cls  an clustering object
#' @param x a vector of length equal to number of samples to use for plotting
#' @param NC number of clusters to cut dendrogram into 
#' @param annoIndx vector indxing into subsets for plotting
#' @param ... additional parameters for plotting. ex. pch,cex,col,labels, xlab, etc.
#' @export
#' @examples \dontrun{
#'  clusterPlotR(cData, cls, x)
#'}
clusterPlotR <- function(
	cData=NA, # data used to get clusters
	cls=NA, # a cluster object
  	x=NA, # a vector of length equal to number of samples to use for plotting
  	NC=NA,# vector of integers indicating which clusters to use 
	annoIndx=NA, #vector indxing into subsets for plotting#vector of integers indicating which columns of Patterns object to use. The default of NP=NA will use entire matrix.
  	... #additional parameters for plotting 
  ){
  UseMethod("clusterPlotR",cls)
}

#' @title clusterPlotR (kmeans)
#'
#' @description plotting function for kmeans clusters
#' @param cData data used to get clusters
#' @param cls  a kmeans object
#' @param x a vector of length equal to number of samples to use for plotting
#' @param NC vector of integers indicating which clusters to use
#' @param annoIndx vector indxing into subsets for plotting
#' @param label character vector to use for plotting text, defaults is NULL
#' @param ... additional parameters for plotting. ex. pch,cex,col,labels, xlab, etc.
#' @export
#' @examples \dontrun{
#'  clusterPlotR(cData=p, cls=pk, x=jitter(pd$days), col=pd$colors)
#'}

clusterPlotR.kmeans <- function(
	cData=NA, # data used to get clusters
	cls=NA, # a kmeans object
	x=NA, # a vector of length equal to number of samples to use for plotting
  	NC=NA,# vector of integers indicating which clusters to use 
	annoIndx=NA, #vector indxing into subsets for plotting
	label=NULL, #character vector to use for plotting text
  	... #additional parameters for plotting 
  ){ 	
if(is.na(NC)){
	cls1<-sort(unique(cls$cluster))
} else(cls1<-NC)
cMNs1=matrix(ncol=dim(cData)[2],nrow=length(cls1))
meanRRs1=vector(length=length(cls1))
for(i in cls1){
	if(sum(cls$cluster==i)>1)
	{
	p1Kc=cData[cls$cluster==i,]
	p1KcMN=colMeans(p1Kc)
	cMNs1[i,]=p1KcMN
	meanRRs1[i]=mean(apply(p1Kc,1,cor,y=p1KcMN))
	}
	if(sum(cls$cluster==i)==1){cMNs1[i,]=pcData[cls$cluster==i,];meanRRs1[i]=1}
	if(sum(cls$cluster==i)==0){print("cluster error !")}
	}	
for(i in cls1){
	plot(x,cMNs1[i,],type="n",main=paste("\nCluster ",i,", N = ",
		sum(cls$cluster==i)," of ",length(cls$cluster)," total genes (",
		round(100*((sum(cls$cluster==i))/(length(cls$cluster))),1),"%)",sep=""),
		ylab=paste("Cluster ",i," : Avg. w/i cluster corr(r) to mean = ",
		round(meanRRs1[i],3),sep=""),...)
	if(length(annoIndx)>0){for(j in unique(annoIndx)){
		lines(x[annoIndx==j],cMNs1[i,annoIndx==j],...)}
		}
	if(is.null(label)){
		points(x,cMNs1[i,],...)
	} else (text(x,cMNs1[i,],labels=label, ...))
	}
}


#' @title clusterPlotR (hclust)
#'
#' @description plotting function for hclust clusters
#' @param cData data used to get clusters
#' @param cls  an hclust object
#' @param x a vector of length equal to number of samples to use for plotting
#' @param NC number of clusters to cut dendrogram into 
#' @param annoIndx vector indxing into subsets for plotting
#' @param ... additional parameters for plotting. ex. pch,cex,col,labels, xlab, etc.
#' @export
#' @examples \dontrun{
#'  clusterPlotR(cData=p, cls=pk, x=jitter(pd$days), col=pd$colors)
#'}

clusterPlotR.hclust <- function(
	cData=NA, # data used to get clusters
	cls=NA, # an hclust object
	x=NA, # a vector of length equal to number of samples to use for plotting
  	NC=NA,  # number of clusters to cut dendrogram into 
  	annoIndx=NA, #vector indxing into subsets for plotting
  	... #additional parameters for plotting. ex. pch,cex,col,labels, xlab, etc. 
  ){
cut1=cutree(cls,k=NC)
cls1=sort(unique(cut1))
cMNs1=matrix(ncol=dim(p1K)[2],nrow=length(cls1))
meanRRs1=vector(length=length(cls1))
for(i in cls1){
	if(sum(cut1==i)>1){
	p1Kc=cData[cut1==i,]
	p1KcMN=colMeans(p1Kc)
	cMNs1[i,]=p1KcMN
	meanRRs1[i]=mean(apply(p1Kc,1,cor,y=p1KcMN))
	}
	if(sum(cut1==i)==1){cMNs1[i,]=p1K[cut1==i,];meanRRs1[i]=1}
	if(sum(cut1==i)==0){print("cluster error !")}
}
for(i in cls1){
	plot(x,cMNs1[i,],type="n",main=paste("\nCluster ",i,", N = ",sum(cut1==i),
		" of ",length(cut1)," total genes (",round(100*((sum(cut1==i))/(length(cut1))),1),
		"%)",sep=""),ylab=paste("Cluster ",i," : Avg. w/i cluster corr(r) to mean = ",
		round(meanRRs1[i],3),sep=""),...)
	if(length(annoIndx)>0){
		for(j in unique(annoIndx)){
			lines(x[annoIndx==j],cMNs1[i,annoIndx==j],...)}}
	if(is.null(label)){
		points(x,cMNs1[i,],...)
	} else (text(x,cMNs1[i,],labels=label, ...))
	}
}