#######################################################################################################################################


# intersectoR.default <- function(
# 	pSet1=NA, #a list for a set of patterns where each entry is a set of genes associated with a single pattern
# 	pSet2=NA, #a list for a second set of patterns where each entry is a set of genes associated with a single pattern
# 	pval=.05, # the maximum p-value considered significant
# 	full=FALSE, #logical indicating whether to return full data frame of signigicantly overlapping sets. Default is false will return summary matrix.
# 	k=NULL #cut height for hclust objects, not used for default
# ){
# 	gns<-lapply(pSet1,function(x) unique(names(x)%in%names(unlist(pSet2))))
# 	indx1<-lapply(pSet1,function(x) match(gns,names(x)))
# 	indx2<-lapply(pSet2,function(x)match(gns,names(x)))
# 	pSet1<-sapply(1:length(pSet1),function(x) pSet1[[x]][indx1[[x]]])
# 	pSet2<-sapply(1:length(pSet2),function(x) pSet2[[x]][indx2[[x]]])

# 	overLPmtx=matrix(nrow=0,ncol=9) #intialize matrix
# 	colnames(overLPmtx)=c("pSet1","NpSet1","pSet2","NpSet2","NoverLP",
# 		"OverLap%1","OverLap%2","pval","pval.REV")

# #calculate overlaps and stats
# 	for(i in 1:length(pSet1)){
# 		for(j in 1:length(pSet2)){
# 	# overlap between genes in p1K cluster i and p2K cluster j
# 		pvalOLP=phyper(
# 			q=length(which(pSet1[[i]] %in% pSet2[[j]])), # q: # white balls drawn without replacement from an urn with both black and white balls.
# 			m=length(pSet1[[i]]), # m: the number of white balls in the urn.
# 			n=length(unlist(pSet1))-length(pSet1[[i]]), # n: the number of black balls in the urn.
# 			k=length(pSet2[[j]]), # k: the number of balls drawn from the urn.
# 			lower.tail = FALSE, log.p = FALSE) # lower.tail: logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
# 		pvalOLP.rev=phyper(q=length(which(pSet2[[j]] %in% pSet1[[i]])), m=length(pSet2[[j]]),
# 			n=length(unlist(pSet2))-length(pSet2[[j]]), k=length(pSet1[[i]]), lower.tail = FALSE, log.p = FALSE)

# 		if(pvalOLP<=pval){overLPmtx=rbind(overLPmtx,c(i,length(pSet1[[i]]),j,length(pSet2[[j]]),
# 			length(which(pSet1[[i]] %in% pSet2[[j]])),
# 			round(100*length(which(pSet1[[i]] %in% pSet2[[j]]))/length(pSet1[[i]]),2),
# 			round(100*length(which(pSet1[[j]] %in% pSet1[[i]]))/length(pSet2[[j]]),2),
# 			pvalOLP,pvalOLP.rev))}
# 		}
# 	}
# 	if(full==FALSE){
# 		return(overLPmtx) #return summary matrix
# 	} else if(full){
# 		overLPindx<-overLPmtx[,c("pSet1","pSet2")] #indx of significantly overlapping sets
# 		overLPsets<-sapply(1:dim(overLPmtx)[1],function(x)
# 			cbind("pSet1"=sort(names(pSet1[pSet1==overLPindx[x,1]&pSet2==overLPindx[x,2]])),
# 				"pSet2"=sort(names(pSet2[pSet2==overLPindx[x,2]&pSet1==overLPindx[x,1]]))))
# 		overLPall=cbind(pSet1,pSet2)
# 		colnames(overLPall)=c("pSet1","pSet2")
# 		rownames(overLPall)=names(pSet1)
# 		ord03=order(pSet1,pSet2)
# 		overLPall=overLPall[ord03,]
# 		return(list(overLPmtx=overLPmtx,overLPindx=overLPindx,overLPsets=overLPsets,overLPall=overLPall))
# 	}
# }

#######################################################################################################################################

intersectoR.kmeans <- function(
	pSet1=NA, #a list for a set of patterns where each entry is a set of genes associated with a single pattern
	pSet2=NA, #a list for a second set of patterns where each entry is a set of genes associated with a single pattern
	pval=.05, #the maximum p-value considered significant
	full=FALSE #logical indicating whether to return full data frame of signigicantly overlapping sets. Default is false will return summary matrix.
  ){
	

	#match gene order
	gns=unique(names(pSet1$cluster)[names(pSet1$cluster)%in%names(pSet2$cluster)])
	ord01=match(gns,names(pSet1$cluster))
	ord02=match(gns,names(pSet2$cluster))
	pSet1$cluster<-pSet1$cluster[ord01]
	pSet2$cluster<-pSet2$cluster[ord02]

	overLPmtx=matrix(nrow=0,ncol=9)
	colnames(overLPmtx)=c("pSet1","NpSet1","pSet2","NpSet2","NoverLP",
		"OverLap%1","OverLap%2","pval","pval.REV")

	for(i in sort(unique(pSet1$cluster))){
		for(j in sort(unique(pSet2$cluster))){
			pvalOLP=phyper(q=sum(pSet1$cluster==i&pSet2$cluster==j),m=sum(pSet1$cluster==i),
				n=sum(pSet1$cluster!=i), k=sum(pSet2$cluster==j),lower.tail = FALSE, log.p = FALSE)
			pvalOLP.rev=phyper(q=sum(pSet2$cluster==j&pSet1$cluster==i), m=sum(pSet2$cluster==j),
				n=sum(pSet2$cluster!=j), k=sum(pSet1$cluster==i), lower.tail = FALSE, log.p = FALSE)
			if(pvalOLP<=pval){overLPmtx=rbind(overLPmtx,c(i,sum(pSet1$cluster==i),j,
				sum(pSet2$cluster==j),sum(pSet1$cluster==i&pSet2$cluster==j),
				round(100*sum(pSet1$cluster==i&pSet2$cluster==j)/sum(pSet1$cluster==i),2),
				round(100*sum(pSet1$cluster==i&pSet2$cluster==j)/sum(pSet2$cluster==j),2),
				pvalOLP,pvalOLP.rev))}
		}
	}
	message(paste(dim(overLPmtx)[1]," cluster pairs have overlap with p<",pval,":",sep=""))
	if(full==FALSE){
		return(overLPmtx)
	} else if(full){
		overLPindx<-overLPmtx[,c("pSet1","pSet2")]
		overLPsets<-sapply(1:dim(overLPmtx)[1],function(x)
			cbind("pSet1"=sort(names(pSet1$cluster[pSet1$cluster==overLPindx[x,1]&pSet2$cluster==overLPindx[x,2]])),
				"pSet2"=sort(names(pSet2$cluster[pSet2$cluster==overLPindx[x,2]&pSet1$cluster==overLPindx[x,1]]))))
		overLPall=cbind(pSet1$cluster,pSet2$cluster)
		colnames(overLPall)=c("pSet1","pSet2")
		rownames(overLPall)=names(pSet1$cluster)
		ord03=order(pSet1$cluster,pSet2$cluster)
		overLPall=overLPall[ord03,]
		return(list(overLPmtx=overLPmtx,overLPindx=overLPindx,overLPsets=overLPsets,overLPall=overLPall))
	}
}


#' @param full logical indicating whether to return full data frame of signigicantly overlapping sets. Default is false will return summary matrix.
#' @examples
#' ESepiGen4c1lmRNASeq <- p.ESepiGen4c1l$mRNA.Seq
#' rownames(ESepiGen4c1lmRNASeq) <- map.ESepiGen4c1l$GeneSymbols
#'
#' k.RNAseq6l3c3t<-kmeans(p.RNAseq6l3c3t,22)
#' k.ESepiGen4c1l<-kmeans(ESepiGen4c1lmRNASeq,10)
#' intersectoR(k.RNAseq6l3c3t, k.ESepiGen4c1l, pval=.05)
#' @rdname intersectoR-methods
#' @aliases intersectoR

setMethod("intersectoR",signature(pSet1 = "kmeans",pSet2 = "kmeans"),intersectoR.kmeans)
#######################################################################################################################################


intersectoR.hclust <- function(
	pSet1=NA, #a hclust obj
	pSet2=NA, #a second set hclust obj
	pval=.05, # the maximum p-value considered significant
	full=FALSE, #logical indicating whether to return full data frame of signigicantly overlapping sets. Default is false will return summary matrix.
	k=NULL #numeric giving cut height for hclust objects, if vector arguments will be applied to pSet1 and pSet2 in that order
  ){
  	overLPmtx=matrix(nrow=0,ncol=9)
  	colnames(overLPmtx)=c("pSet1","NpSet1","pSet2","NpSet2","NoverLP",
		"OverLap%1","OverLap%2","pval","pval.REV")

	if(length(k)==1){cut1<-cutree(pSet1,k=k) ; cut2<-cutree(pSet2,k=k)}
	if(length(k)==2){cut1<-cutree(pSet1,k=k[1]) ; cut2<-cutree(pSet2,k=k[2])}

	#cut1<-cut1[names(cut1)%in%names(cut2)]
	#cut2<-cut2[names(cut2)%in%names(cut1)]
	gns=unique(names(cut1)[names(cut1)%in%names(cut2)])
	ord01=match(gns,names(cut1))
	ord02=match(gns,names(cut2))
	cut1<-cut1[ord01]
	cut2<-cut2[ord02]

	for(i in sort(unique(cut1))){
		for(j in sort(unique(cut2))){
			pvalOLP=phyper(q=sum(cut1==i&cut2==j),m=sum(cut1==i),n=sum(cut1!=i),
				k=sum(cut2==j), lower.tail = FALSE, log.p = FALSE)
			pvalOLP.rev=phyper(q=sum(cut2==j&cut1==i), m=sum(cut2==j), n=sum(cut2!=j),
				k=sum(cut1==i), lower.tail = FALSE, log.p = FALSE)
			if(pvalOLP<=pval){overLPmtx=rbind(overLPmtx,c(i,sum(cut1==i),j,
				sum(cut2==j),sum(cut1==i&cut2==j),sum(cut1==i&cut2==j)/sum(cut1==i),
				sum(cut1==i&cut2==j)/sum(cut2==j),pvalOLP,pvalOLP.rev))}
		}
	}
	print(paste(dim(overLPmtx)[1]," cluster pairs have overlap with p<",pval,":",sep=""))
	if(full==FALSE){
		return(overLPmtx) #return summary matrix
	} else if(full){
		overLPindx<-overLPmtx[,c("pSet1","pSet2")]
		overLPsets<-sapply(1:dim(overLPmtx)[1],function(x)
			cbind("pSet1"=sort(names(cut1[cut1==overLPindx[x,1]&cut1==overLPindx[x,2]])),
				"pSet2"=sort(names(cut2[cut2==overLPindx[x,2]&cut2==overLPindx[x,1]]))))
		overLPall=cbind(cut1,cut2)
		colnames(overLPall)=c("pSet1","pSet2")
		rownames(overLPall)=names(cut1)
		ord03=order(cut1,cut2)
		overLPall=overLPall[ord03,]
		return(list(overLPmtx=overLPmtx,overLPindx=overLPindx,overLPsets=overLPsets,overLPall=overLPall))

	}
}

#' @param k Numeric giving cut height for hclust objects, if a vector is given arguments will be applied to pSet1 and pSet2 in that order
#' @examples 
#'	
#'	h.RNAseq6l3c3t<-hclust(as.dist(1-(cor(t(p.RNAseq6l3c3t)))))
#'	h.ESepiGen4c1l<-hclust(as.dist(1-(cor(t(ESepiGen4c1lmRNASeq)))))
#'	intersectoR(pSet1=h.ESepiGen4c1l, pSet2=h.RNAseq6l3c3t, pval=.05, k=c(3,4))
#'
#' @rdname intersectoR-methods
#' @aliases intersectoR

setMethod("intersectoR",signature(pSet1 = "hclust",pSet2 = "hclust"),intersectoR.hclust)
