

cor2T<-correlateR(genes="T", dat=p, threshtype="N", threshold=10, absR=TRUE) 
lapply(cor2T,head)

correlateR<-function(genes=NA, #gene or character vector of genes for reference expression pattern
	dat=NA,
	threshtype="R", #Default "R" indicates thresholding by R value or equivalent. Alternatively, "N" indicates a numerical cut off. 
	threshold=.7, # numeric indicating value at which to make threshold
	absR=FALSE, # logical indicating where to include both positive and negatively correlated genes  
	... # addtion imputs to cor, such as method 
	){
	if(threshtype=="N" & threshold<1){
		stop('Threshold must be integer greater than 1 for threshold type "N"')
	} else if(threshtype=="R" & threshold>1){
		stop('Threshold must be between -1 and 1 for threshold type "R"')
	}
	cor2gene<-apply(dat, 1, function(G) cor(t(p[genes,]), G, ...))
	if(absR){
		if(threshtype=="R"){
			corGS<-list("PositiveCOR"=sort(cor2gene,decreasing=TRUE)[sort(cor2gene,decreasing=TRUE)>=threshold],
					"NegativeCOR"=rev(sort(cor2gene,decreasing=TRUE)[sort(cor2gene,decreasing=TRUE) <= -threshold]))
		} else if(threshtype=="N"){
			corGS<-list("PositiveCOR"=sort(cor2gene,decreasing=TRUE)[1:threshold],
					"NegativeCOR"=sort(cor2gene,decreasing=TRUE)[dim(dat)[1]:(dim(dat)[1]-threshold)])
		}
	} else {
		if(threshtype=="R"){
			corGS<-sort(cor2gene,decreasing=TRUE)[which(sort(cor2gene,decreasing=TRUE)>=threshold)]
		} else if(threshtype=="N"){
			corGS<-sort(cor2gene,decreasing=TRUE)[1:threshold]
		}
	}
	return(corGS)
}



