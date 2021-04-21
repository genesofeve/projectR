#######################################################################################################################################
#' @import
#' @importFrom 
#' @param 
###STATS FUNCTIONS HERE


#######################################################################################################################################
#' @import
#' @importFrom 
#' @param cellgroup1 gene x cell count matrix for cell group 1
#' @param cellgroup2 gene x cell count matrix for cell group 2
#' @param loadings A matrix of continuous values defining the features
#' @param feature_name column of loadings for which drivers will be calculated.
#' @param alpha confidence value for the bonferroni confidence intervals
#' @param loadingsNames a vector with names of loading rows. Defaults to rownames.
#' @rdname 
#' 
#' 
.drivers_matrix<-function(
  cellgroup1, #gene x cell count matrix for cell group 1
  cellgroup2, #gene x cell count matrix for cell group 2
  loadings, # a matrix of continous values to be projected with unique rownames
  loadingsNames = NULL # a vector with names of loadings rows
  
){
  
  #check that alpha significance level is appropriate
  if(alpha < 0 | alpha > 1){
    stop("alpha must be numeric between 0 and 1")
  }
  
  #select specified feature to calculate drivers for
  if(feature_name %in% colnames(loadings)){
    feature <- loadings[,feature_name, drop = F] #data.frame
  } else  {
    stop(paste0(feature_name, " is not a column in provided loadings"))
  }
  
  # #match genes in data sets
  # if(is.null(dataNames)){
  #   dataNames <- rownames(data)
  # }
  if(is.null(loadingsNames)){
    loadingsNames <- rownames(loadings)
  }
  
  #shared rows in two data matrices
  filtered_data <-geneMatchR(data1=cellgroup1, data2=cellgroup2, data1Names=NULL, data2Names=NULL, merge=FALSE)
  print(paste(as.character(dim(filtered_data[[2]])[1]),'row names matched between datasets'))
  
  cellgroup1 <- filtered_data[[2]] #TODO: geneMatchR flips the indexes, correct?
  cellgroup2 <- filtered_data[[1]]
  
  #shared rows in data matrices and loadings
  filtered_weights <- geneMatchR(data1 = cellgroup1, data2 = feature, data1Names = NULL, data2Names = NULL, merge = F)
  print(paste(as.character(dim(filtered_weights[[2]])[1]),'row names matched between data and loadings'))
  print(paste('Updated dimension of data:',as.character(paste(dim(filtered_weights[[2]])[1], collapse = ' '))))
  
  feature_filtered <- filtered_weights[[1]]
  cellgroup1_filtered <- filtered_weights[[2]]
  
}  
  
#setMethod("geneDriveR",signature(data="matrix",loadings="matrix"),.drivers_matrix)

