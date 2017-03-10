

intersectoR<-function(
	pSet1=NA, #a set of patterns 
	pSet2=NA,
	plot=FALSE
  ){
  UseMethod("intersectoR",Patterns)
}

intersectoR.default <- function(
		pSet1=NA, #a set of patterns 
	pSet2=NA,
	plot=FALSE
  ){
overLPmtx=matrix(nrow=0,ncol=9)
colnames(overLPmtx)=c("pSet1","NpSet1","pSet2","NpSet2","NoverLP","OverLap%1","OverLap%2","pval","pval.REV")


for(i in length(pSet1)){
	for(j in length(pSet2)){
# overlap between genes in p1K cluster i and p2K cluster j
	pvalOLP=phyper(
		q=length(which(pSet1[[i]] %in% pSet2[[j]])), # q: # white balls drawn without replacement from an urn with both black and white balls.
		m=length(pSet1[[i]]), # m: the number of white balls in the urn.
		n=length(unlist(pSet1))-length(pSet1[[i]]), # n: the number of black balls in the urn.
		k=length(pSet2[[j]]), # k: the number of balls drawn from the urn.
		lower.tail = FALSE, log.p = FALSE) # lower.tail: logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
	pvalOLP.rev=phyper(q=length(which(pSet2[[j]] %in% pSet1[[i]])), m=length(Set2[[j]]), 
		n=length(unlist(pSet2))-length(pSet2[[j]]), k=length(pSet1[[i]]), lower.tail = FALSE, log.p = FALSE)

	if(pvalOLP<=KoverLPpval){overLPmtx=rbind(overLPmtx,c(i,sum(pSet1$cluster==i),j,
		sum(pSet2$cluster==j),sum(pSet1$cluster==i&pSet2$cluster==j),NA,NA,pvalOLP,pvalOLP.rev))}
	}
}

if(dim(overLPmtx)[1]>=1){
	overLPmtx[,"OverLap%1"]=round(100*((overLPmtx[,"NoverLP"])/(overLPmtx[,"Ncls1"])),1)
	overLPmtx[,"OverLap%2"]=round(100*((overLPmtx[,"NoverLP"])/(overLPmtx[,"Ncls2"])),1)
	print(paste(dim(overLPmtx)[1]," cluster pairs have overlap with p<",KoverLPpval,":",sep=""))
	# export signif cluster pair info
	overLPmtx

	overlaps<-cbind(pSet1$cluster[indxKoverLP],pSet2$cluster[indxKoverLP])
}

}  	


intersectoR.kmeans <- function(
	pSet1=NA, #a set of patterns 
	pSet2=NA,
	plot=FALSE
  ){

# find all "significant" cluster pair overlaps + then cycle thru them in a loop inside a pdf call:
overLPmtx=matrix(nrow=0,ncol=9)
colnames(overLPmtx)=c("cls1","Ncls1","cls2","Ncls2","NoverLP","OverLap%1","OverLap%2","pval","pval.REV")
cls1<-length(pSet1$cluster)
cls2<-length(pSet2$cluster)

cls2=sort(unique(KB2$cluster))
cMNs2=matrix(ncol=dim(p2K)[2],nrow=length(cls2))
meanRRs2=vector(length=length(cls2))
for(i in cls2)
{
	if(sum(KB2$cluster==i)>1)
	{
	p2Kc=p2K[KB2$cluster==i,]
	p2KcMN=colMeans(p2Kc)
	cMNs2[i,]=p2KcMN
	meanRRs2[i]=mean(apply(p2Kc,1,cor,y=p2KcMN))
	}
	if(sum(KB2$cluster==i)==1){cMNs2[i,]=p2K[KB2$cluster==i,];meanRRs2[i]=1}
	if(sum(KB2$cluster==i)==0){print("cluster error !")}
}


for(i in cls1){
	for(j in cls2){
# overlap between genes in p1K cluster i and p2K cluster j
	pvalOLP=phyper(
		q=sum($cluster==i&pSet2$cluster==j), # q: # white balls drawn without replacement from an urn with both black and white balls.
		m=sum(pSet1$cluster==i), # m: the number of white balls in the urn.
		n=sum(pSet1$cluster!=i), # n: the number of black balls in the urn.
		k=sum(pSet2$cluster==j), # k: the number of balls drawn from the urn.
		lower.tail = FALSE, log.p = FALSE) # lower.tail: logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
	pvalOLP.rev=phyper(q=sum(pSet2$cluster==j&pSet1$cluster==i), m=sum(pSet2$cluster==j), 
		n=sum(pSet2$cluster!=j), k=sum(pSet1$cluster==i), lower.tail = FALSE, log.p = FALSE)

	if(pvalOLP<=KoverLPpval){overLPmtx=rbind(overLPmtx,c(i,sum(pSet1$cluster==i),j,
		sum(pSet2$cluster==j),sum(pSet1$cluster==i&pSet2$cluster==j),NA,NA,pvalOLP,pvalOLP.rev))}
	}
}

if(dim(overLPmtx)[1]>=1){
	overLPmtx[,"OverLap%1"]=round(100*((overLPmtx[,"NoverLP"])/(overLPmtx[,"Ncls1"])),1)
	overLPmtx[,"OverLap%2"]=round(100*((overLPmtx[,"NoverLP"])/(overLPmtx[,"Ncls2"])),1)
	print(paste(dim(overLPmtx)[1]," cluster pairs have overlap with p<",KoverLPpval,":",sep=""))
	# export signif cluster pair info
	overLPmtx

	overlaps<-cbind(pSet1$cluster[indxKoverLP],pSet2$cluster[indxKoverLP])
}

}


intersectoR.hclust <- function(
	pSet1=NA, #a set of patterns 
	pSet2=NA,
	plot=FALSE
  ){


 }