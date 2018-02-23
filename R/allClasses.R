####################
# Classes
#
# Author:
#
####################

setClass("ProjectionSet"),
	representation( targetData="matrix", 
					sourcePatterns="matrix",
					AnnotionObj="data.frame",
					IDcol="character"
					)
	)

setMethod("initialize","ProjectionSet",
			function(.Object,
					targetData,
					sourcePatterns,
					AnnotationObj=NA,
					IDcol=NA,
					idField,
					... ){
				.Object<-callNextMethod(.Object,
						targetData = targetData,
						sourcePatterns = sourcePatterns,
						AnnotationObj = AnnotationObj,
						IDcol = IDcol,
						...)				
		}
)

setValidity("ProjectionSet",function(object){
		TRUE
		}
)			

################
#Class Methods
################
setMethod("show","ProjectionSet",
		function(object){
			#######
		}
)

setMethod("dim","ProjectionSet",
		function(x){
			#######
		}
)

#Below borrowed from CummeRbund as an example method.
# .addFeatures<-function(object,features,...){
# 	if(!is.data.frame(features)){
# 		stop("features must be a data.frame")
# 	}
# 	colnames(features)[1]<-object@idField
# 	colnames(features)<-make.db.names(object@DB,colnames(features),unique=T)
# 	dbWriteTable(object@DB,object@tables$featureTable,features,row.names=F,overwrite=T)
# 	indexQuery<-paste("CREATE INDEX ",object@idField," ON ", object@tables$featureTable," (",object@idField,")",sep="")
# 	res<-dbGetQuery(object@DB,indexQuery)
# }

# setMethod("addFeatures",signature="ProjectionSet",.addFeatures)
