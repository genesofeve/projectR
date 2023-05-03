#' @title getUMAP
#' @import umap
#' @description Function to provide umap of projection
#' @param   projection  martrix, a projection generated from projectR
#' @param   axis		integer, either 1 umap of projection or 2 for umap of transpose of projection
#' @param   umapMethod  character, implementation. Available methods are 'naive' (an implementation written in pure R) and 'umap-learn' (requires python package 'umap-learn')
#' @param 	umapConfig  umap.config, a list of parameters to customize umap embedding
#' @return  A umap of projection
#' @examples
#' projection <- projectR(data=p.ESepiGen4c1l$mRNA.Seq,loadings=AP.RNAseq6l3c3t$Amean, 
#' dataNames = map.ESepiGen4c1l[["GeneSymbols"]], full = TRUE)
#' umapConfig = umap.defaults
#' umapConfig$n_neighbors = 3
#' projectionUMAP <- getUMAP(projection,umapConfig = umapConfig)
#' @export

getUMAP <- function(projection,axis=2,umapMethod="naive",umapConfig=umap.defaults){
	if(inherits(projection,"list")){
		dat <- projection[[1]]
		if(axis == 2) {
			dat <- t(dat)
		}
		return(umap(dat,method=umapMethod,umapConfig))
	} else {
		dat <- projection
		if(axis == 2) {
			dat <- t(dat)
		}
		return(umap(dat,method=umapMethod,umapConfig))
	}
}

#' @title getTSNE
#' @import tsne
#' @description Function to provide tSNE of projection
#' @param   projection  martrix, a projection generated from projectR
#' @param   axis		integer, either 1 umap of projection or 2 for umap of transpose of projection
#' @param   ...			addtional arguments passed to tsne
#' @examples
#' projection <- projectR(data=p.ESepiGen4c1l$mRNA.Seq,loadings=AP.RNAseq6l3c3t$Amean, 
#' dataNames = map.ESepiGen4c1l[["GeneSymbols"]], full = TRUE)
#' projectionTSNE <- getTSNE(projection)
#' @export

getTSNE <- function(projection,axis=2,...){
	if(inherits(projection,"list")){
		dat <- projection[[1]]
		if(axis == 2) {
			dat <- t(dat)
		}
		return(tsne(dat,...))
	} else {
		dat <- projection
		if(axis == 2) {
			dat <- t(dat)
		}
		return(tsne(dat,...))
	}
}