#' @importFrom stats hclust kmeans
setOldClass("kmeans")
setOldClass("hclust")

####################
# Classes
# Author
####################

# setClass("ProjectionSet"),
# 	representation( targetData="matrix",
# 					sourcePatterns="matrix",
# 					AnnotionObj="data.frame",
# 					IDcol="character"
# 					)
# 	)
#
# setMethod("initialize","ProjectionSet",
# 			function(.Object,
# 					targetData,
# 					sourcePatterns,
# 					AnnotationObj=NA,
# 					IDcol=NA,
# 					idField,
# 					... ){
# 				.Object<-callNextMethod(.Object,
# 						targetData = targetData,
# 						sourcePatterns = sourcePatterns,
# 						AnnotationObj = AnnotationObj,
# 						IDcol = IDcol,
# 						...)
# 		}
# )
#
# setValidity("ProjectionSet",function(object){
# 		TRUE
# 		}
# )
#
# ################
# #Class Methods
# ################
# setMethod("show","ProjectionSet",
# 		function(object){
# 			#######
# 		}
# )
#
# setMethod("dim","ProjectionSet",
# 		function(x){
# 			#######
# 		}
# )


#' pclust
#' @export
#'
#' @slot pattern pattern found from clusters using cluster2pattern
#' @description parent class of plcustKmeans and plcusltHclust

setClass("pclust", slots=c(
	patterns = "matrix"      
))

#' Constructor for pclust
#' @param .Object pclust object
#' @param hclust hclust object passed to cluster2pattern function
#' @param kmeans kmeans object passed to cluster2pattern function
#' @return initialized plclust object

#' @importFrom methods callNextMethod

setMethod("initialize", "pclust",
function(.Object, patterns, ... )
{
	.Object <- callNextMethod(.Object, ...)
	.Object@patterns <- patterns
    .Object
})

# setValidity("hclust",
#     function(object)
#     {    }
# )    

#' pclustKmeans
#' @export
#' @slot kmeans kmeans input to cluster2pattern
#' @description defines class of cluster2pattern with kmeans input
setClass("pclustKmeans", slots=c(
	kmeans = "kmeans"	      
),contains = "pclust")

#' Constructor for pclustKmeans
#' @param .Object pclustKmeans object
#' @return initialized plclustKmeans object
#' @importFrom methods callNextMethod
setMethod("initialize", "pclustKmeans",
function(.Object, kmeans, ... )
{
	.Object <- callNextMethod(.Object, ...)
	.Object@kmeans <- kmeans
    .Object
})

#' pclustHclust
#' @export
#' @slot hclust hclust input to cluster2pattern
#' @description defines class of cluster2pattern with hclust input
setClass("pclustHclust", slots=c(
	hclust = "hclust"	      
),contains = "pclust")

#' Constructor for pclust
#' @param .Object pclust object
#' @param hclust hclust object passed to cluster2pattern function
#' @return initialized plclustHclust object
#' @importFrom methods callNextMethod
setMethod("initialize", "pclustHclust",
function(.Object, hclust, ... )
{
	.Object <- callNextMethod(.Object, ...)
	.Object@hclust <- hclust
    .Object
})

#' rotatoR
#' @export
#'
#' @slot pattern pattern found from clusters using cluster2pattern
#' @description class of function roatoR's output

setClass("rotatoR", slots=c(
	rotatedM = "matrix"      
))

#' Constructor for rotatoR
#' @param .Object rotatoR object
#' @param rotatedM rotated matrix from rotatoR function
#' @return initialized rotatoR object

#' @importFrom methods callNextMethod

setMethod("initialize", "rotatoR",
function(.Object, rotatedM, ... )
{
	.Object <- callNextMethod(.Object, ...)
	.Object@rotatedM <- rotatedM
    .Object
})