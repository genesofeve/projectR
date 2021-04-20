#######################################################################################################################################
#' @import
#' @importFrom 
#' @param 
#' @rdname 
#' @aliases 
#' 
.drivers_matrix<-function(
  data, # a dataset to be projected onto
  loadings, # a matrix of continous values to be projected with unique rownames
  dataNames = NULL, # a vector with names of data rows
  loadingsNames = NULL, # a vector with names of loadings rows
  NP=NA, # vector of integers indicating which columns of loadings object to use. The default of NP=NA will use entire matrix.
){
  
  ifelse(!is.na(NP),loadings<-loadings[,NP],loadings<-loadings)
  #if(!is.na(NP)){loadings<-loadings[,NP]} was giving warning with subset of patterns
  #match genes in data sets
  if(is.null(dataNames)){
    dataNames <- rownames(data)
  }
  if(is.null(loadingsNames)){
    loadingsNames <- rownames(loadings)
  }
  dataM<-geneMatchR(data1=data, data2=loadings, data1Names=dataNames, data2Names=loadingsNames, merge=FALSE)
  print(paste(as.character(dim(dataM[[2]])[1]),'row names matched between data and loadings'))
  print(paste('Updated dimension of data:',as.character(paste(dim(dataM[[2]]), collapse = ' '))))
  
setMethod("geneDriveR",signature(data="matrix",loadings="matrix"),.drivers_matrix)

