#######################################################################################################################################
#' @import
#' @importFrom 
#' @param 
bonferroniCorrectedDifferences(cellgroup1, cellgroup2, alpha){
  
  #if passed from projectionDrivers, cellgroup1 and cellgroup 1 will have the same rows (genes)
  #TODO: if going to be called in other places, or directly, need to do that filtering prior to call
  
  if(!(dim(cellgroup1)[1] == dim(cellgroup2)[1])){
    stop("Rows of two cell group matrices are not identical")
  }
  
  dimensionality <- dim(cellgroup1)
  
}

#######################################################################################################################################
#' @import stats
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
  
  #TODO: Do sparse and dense matrices need to be handled differently?
  
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
  #do second filtering on other cell group so all genes are consistent
  cellgroup2_filtered <- cellgroup2[rownames(cellgroup1_filtered),]
  
  ##Take means over all genes and calculate differences
  group1_mean <- apply(cellgroup1_filtered, 1, mean)
  group2_mean <- apply(cellgroup2_filtered, 1, mean)
  
  mean_diff <- group1_mean - group2_mean
  
  #normalize feature weights
  weight_norm <- norm(feature_filtered) #square of sums of squares (sum for all positive values)
  feature_normalized <- feature_filtered #TODO: add normalization here
  
  if(!(names(mean_diff) == rownames(feature_normalized) && T)){
    stop("Names of loadings and counts do not match")
  }
  drivers <- mean_diff * feature_normalized
  
  ##Stats and corrections beginning here
  
  dimensionality <- dim(drivers)[1] #number of measurements (genes)
  
  n1_samples <- dim(cellgroup1_filtered)[2] 
  n2_samples <- dim(cellgroup2_filtered)[2]
  pval <- 1 - alpha
  bon_correct <- pval / (2*dimensionality)
  qval <- 1 - bon_correct
  
  tval <- qt(p = qval, df = n1_samples + n2_samples -2) #critical value
  
  ##TODO: Compute covariance and pool here. What exactly are these values?
  pooled <- matrix(rep(1, dimensionality*dimensionality), nrow = dimensionality)
  
  
  plusminus = data.frame(low = rep(NA_integer_, dimensionality), high = rep(NA_integer_, dimensionality))
  #for each gene
  for(i in 1:dimensionality){
      
    scale = tval * sqrt(pooled[i,i] * (1/n1_samples + 1/n2_samples))
      
    plusminus[i, "low"] <- drivers[i, 1]  - scale
    plusminus[i, "high"] <- drivers[i, 1] + scale
      
  }

  ##Plot CI here
}  
  
#setMethod("geneDriveR",signature(data="matrix",loadings="matrix"),.drivers_matrix)

