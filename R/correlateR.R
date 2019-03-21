#' @title correlateR
#'
#' @description Function to extract genes highly correlated with a gene or reference expression pattern.
#' @param 	genes  gene or character vector of genes for reference expression pattern
#' @param 	dat matrix or data frame with  genes to be used for to calculate correlation
#' @param 	threshtype Default "R" indicates thresholding by R value or equivalent. Alternatively, "N" indicates a numerical cut off.
#' @param 	threshold numeric indicating value at which to make threshold.
#' @param 	absR logical indicating where to include both positive and negatively correlated genes
#' @param 	...  addtion inputs to cor, such as method
#' @return 	A correlation matrix
#' @export
#' @import stats
#' @examples
#' cor2T<-correlateR(genes="T", dat=p.RNAseq6l3c3t, threshtype="N", threshold=10, absR=TRUE)
#'
#' @details 
#' If threshtype is "R" than threshold must be between -1 and 1. Otherwise if top N correlated genes are required, set \code{threshtype}
#' as "N" and set \code{threshold} = N, i.e, the number of correlated genes required.

correlateR<-function(genes, #gene or character vector of genes for reference expression pattern
	dat, #matrix or data frame with  genes to be used for to calculate correlation
	threshtype="R", #Default "R" indicates thresholding by R value or equivalent. Alternatively, "N" indicates a numerical cut off.
	threshold=.7, # numeric indicating value at which to make threshold
	absR=FALSE, # logical indicating where to include both positive and negatively correlated genes
	... # addtion inputs to cor, such as method
	){
	if(threshtype=="N" & threshold<1){
		stop('Threshold must be integer greater than 1 for threshold type "N"')
	} else if(threshtype=="R" & threshold>1){
		stop('Threshold must be between -1 and 1 for threshold type "R"')
	}
	cor2gene<-apply(dat, 1, function(G) cor(t(dat[genes,]), G, ...))
	if(absR){
		if(threshtype=="R"){
			corGS<-list("PositiveCOR"=as.matrix(sort(cor2gene,decreasing=TRUE)[sort(cor2gene,decreasing=TRUE)>=threshold]),
					"NegativeCOR"=as.matrix(rev(sort(cor2gene,decreasing=TRUE)[sort(cor2gene,decreasing=TRUE) <= -threshold])))
		} else if(threshtype=="N"){
			corGS<-list("PositiveCOR"=as.matrix(sort(cor2gene,decreasing=TRUE)[1:threshold]),
					"NegativeCOR"=as.matrix(sort(cor2gene,decreasing=TRUE)[dim(dat)[1]:(dim(dat)[1]- (threshold-1))]))
		}
	} else {
		if(threshtype=="R"){
			corGS<-as.matrix(sort(cor2gene,decreasing=TRUE)[which(sort(cor2gene,decreasing=TRUE)>=threshold)])
		} else if(threshtype=="N"){
			corGS<-as.matrix(sort(cor2gene,decreasing=TRUE)[1:threshold])
		}
	}
	corR <- new("correlateR",corM = corGS)
	return(corR)
}



