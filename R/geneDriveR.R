#######################################################################################################################################
#' 
#' Calculate the (weighted) difference in means for each measurement between two groups. 
#' 
#' 
#' @importFrom stats var
#' @param group1
#' @export
bonferroniCorrectedDifferences <- function(
  group1,
  group2,
  diff_weights = NULL,
  alpha){
  
  #if passed from projectionDrivers, cellgroup1 and cellgroup 1 will have the same rows (genes)
  #TODO: if going to be called in other places, or directly, need to do that filtering prior to call
  
  if(!(dim(group1)[1] == dim(group2)[1])){
    stop("Rows of two cell group matrices are not identical")
  }
  
  ##Take means over all genes and calculate differences
  group1_mean <- apply(group1, 1, mean)
  group2_mean <- apply(group2, 1, mean)
  
  mean_diff <- group1_mean - group2_mean
  
  
  

  #if weights are provided, use them to weight the difference in means
  if(!is.null(diff_weights)){
    
    #check that genes exactly match between difference vector and weight vector
    if(!(names(mean_diff) == names(diff_weights) && T)){
      stop("Names of loadings and counts do not match")
    }
    
    mean_diff <- mean_diff * diff_weights
  }
  
  
  
  ##Stats and corrections beginning here
  dimensionality <- length(mean_diff) #number of measurements (genes)
  
  n1_samples <- dim(group1)[2] #number of samples (cells)
  n2_samples <- dim(group2)[2]
  pval <- 1 - alpha
  bon_correct <- pval / (2*dimensionality)
  qval <- 1 - bon_correct
  
  tval <- qt(p = qval, df = n1_samples + n2_samples -2) #critical value
  
  ##TODO: Compute covariance and pool here. What exactly are these values?
  group1_var <- apply(group1, 1, var)
  group2_var <- apply(group2, 1, var)
  
  #TODO: is this right?
  pooled <- ((n1_samples-1)*group1_var + (n2_samples-1)*group2_var)  /  (n1_samples+n2_samples-2)
    
  #establish dataframe to populate in the following for loop
  plusminus = data.frame(low = rep(NA_integer_, dimensionality), high = rep(NA_integer_, dimensionality))
  rownames(plusminus) <- names(mean_diff)
  
  #for each gene
  for(i in 1:dimensionality){
    
    scale = tval * sqrt(pooled[i] * (1/n1_samples + 1/n2_samples))
    
    plusminus[i, "low"] <- mean_diff[i]  - scale
    plusminus[i, "high"] <- mean_diff[i] + scale
    
  }
  
  
  ##Plot CI here
  
  return(plusminus)
}  



#######################################################################################################################################
#' projectionDriveR
#' 
#' Calculate the weighted difference in expression between two groups (group1 - group2)
#' 
#' @param cellgroup1 gene x cell count matrix for cell group 1
#' @param cellgroup2 gene x cell count matrix for cell group 2
#' @param loadings A matrix of continuous values defining the features
#' @param feature_name column of loadings for which drivers will be calculated.
#' @param alpha confidence value for the bonferroni confidence intervals
#' @param loadingsNames a vector with names of loading rows. Defaults to rownames.
#' @param display boolean. Whether or not to plot and display confidence intervals
#' @param normalize_feature Boolean. Whether or not to normalize feature weights.
#' @return A list with weighted mean differences, mean differences, and differential genes that meet the provided signficance threshold.
#' @export
#' 
#' 
projectionDriveR<-function(
  cellgroup1, #gene x cell count matrix for cell group 1
  cellgroup2, #gene x cell count matrix for cell group 2
  loadings, # a matrix of continous values to be projected with unique rownames
  loadingsNames = NULL, # a vector with names of loadings rows
  feature_name,
  alpha,
  display = TRUE,
  normalize_feature = TRUE
){
  
  #TODO: Something isnt right with the documentation, arguments do not autosuggest
  #TODO: Do sparse and dense matrices need to be handled differently?
  #TODO: assert rownames and colnames exist where needed, and that things are matrices (or can be cast to)
  #TODO: Enforce that loadings is one row? or loop over patterns?
  
  #check that alpha significance level is appropriate
  if(alpha < 0 | alpha > 1){
    stop("alpha must be numeric between 0 and 1")
  }
  
  #select specified feature to calculate drivers for
  if(length(feature_name) != 1 && class(feature_name) == "character"){
    stop("provided feature name must not a character vector of length one")
  }
  
  if(is.null(loadingsNames)){
    loadingsNames <- rownames(loadings)
    #TODO: set loadingsnames if provided
  }
  
  #feature weights must be formatted as a matrix for normalization
  if(feature_name %in% colnames(loadings)){
    feature <- loadings[,feature_name, drop = F] #data.frame
    feature <- as.matrix(feature)
  } else  {
    stop(paste0(feature_name, " is not a column in provided loadings"))
  }
  
  # #match genes in data sets
  # if(is.null(dataNames)){
  #   dataNames <- rownames(data)
  # }

  
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
  
  
  #normalize feature weights
  if(normalize_feature){
    weight_norm <- norm(feature_filtered) #square of sums of squares (sum for all positive values)
    num_nonzero <- sum(feature_filtered > 0) #number of nonzero weights
    feature_filtered <- feature_filtered * num_nonzero / weight_norm
  }
  
  #cast feature weights to a named vector
  feature_normalized_vec <- feature_filtered[,1]
  names(feature_normalized_vec) <- rownames(feature_filtered)
  
  #weighted confidence intervals of differences in cluster means
  weighted_drivers_bonferroni <- bonferroniCorrectedDifferences(group1 = cellgroup1_filtered,
                                group2 = cellgroup2_filtered,
                                diff_weights = feature_normalized_vec,
                                alpha = alpha)
  
  
   #unweighted confidence intervals of difference in cluster means
  mean_bonferroni <- bonferroniCorrectedDifferences(group1 = cellgroup1_filtered,
                                                    group2 = cellgroup2_filtered,
                                                    diff_weights = NULL,
                                                    alpha = alpha)
  
  #Determine which genes have significantly non-zero mean difference and weighted mean difference
  #significant
  weighted_sig_idx <- apply(weighted_drivers_bonferroni, 1, function(interval){
    (interval[1] > 0 & interval[2] > 0) | (interval[1] < 0 & interval[2] < 0)
  })
  
  mean_sig_idx <- apply(mean_bonferroni, 1, function(interval){
    (interval[1] > 0 & interval[2] > 0) | (interval[1] < 0 & interval[2] < 0)
  })
  
  #genes that are collectively either up or down
  shared_genes <- base::intersect(
    rownames(weighted_drivers_bonferroni)[weighted_sig_idx],
    rownames(mean_bonferroni)[mean_sig_idx])
  
  
  sorted_conf_intervals <- mean_bonferroni[shared_genes,]
  
  if(display){
    plotConfidenceIntervals(sorted_conf_intervals)
    #TODO: return plot below too?
  }
  
  return(list(
    weighted_mean_differences = weighted_drivers_bonferroni,
    mean_differences = mean_bonferroni,
    significant_genes = shared_genes))
}
#setMethod("projectionDriveR",signature(data="matrix",loadings="matrix"),.drivers_matrix)

