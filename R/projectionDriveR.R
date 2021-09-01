#######################################################################################################################################
#' 
#' Calculate the (weighted) difference in means for each measurement between two groups. 
#' 
#' 
#' @importFrom stats var
bonferroniCorrectedDifferences <- function(
  group1, #count matrix 1
  group2, #count matrix 2
  diff_weights = NULL, #loadings to weight the differential expression between the groups
  pvalue) #signficance value to threshold at
  {
  
  #if passed from projectionDrivers, cellgroup1 and cellgroup 1 will have the same rows (genes)

  if(!(dim(group1)[1] == dim(group2)[1])){
    stop("Rows of two cell group matrices are not identical")
  }
  
  if(any(is.na(group1)) | any(is.na(group2))){
    stop("NA values in count matrices not allowed")
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
  bon_correct <- pvalue / (2*dimensionality) #bonferroni correction
  qval <- 1 - bon_correct
  
  tval <- qt(p = qval, df = n1_samples + n2_samples -2) #critical value
  
  group1_var <- apply(group1, 1, var)
  group2_var <- apply(group2, 1, var)
  
  pooled <- ((n1_samples-1)*group1_var + (n2_samples-1)*group2_var)  /  (n1_samples+n2_samples-2)
    
  #establish dataframe to populate in the following for loop
  plusminus = data.frame(low = rep(NA_integer_, dimensionality), high = rep(NA_integer_, dimensionality))
  rownames(plusminus) <- names(mean_diff)
  
  
  #for each gene, calculate confidence interval around mean
  for(i in 1:dimensionality){
    
    scale = tval * sqrt(pooled[i] * (1/n1_samples + 1/n2_samples))
    
    plusminus[i, "low"] <- mean_diff[i]  - scale #low estimate
    plusminus[i, "high"] <- mean_diff[i] + scale #high estimate
    
  }
  
  
  return(plusminus)
}  



#######################################################################################################################################
#' projectionDriveR
#' 
#' Calculate the weighted difference in expression between two groups (group1 - group2)
#' 
#' @importFrom cowplot plot_grid
#' @param cellgroup1 gene x cell count matrix for cell group 1
#' @param cellgroup2 gene x cell count matrix for cell group 2
#' @param loadings A matrix of continuous values defining the features
#' @param pattern_name column of loadings for which drivers will be calculated.
#' @param pvalue confidence level for the bonferroni confidence intervals. Default 1e-5
#' @param loadingsNames a vector with names of loading rows. Defaults to rownames.
#' @param display boolean. Whether or not to plot and display confidence intervals
#' @param normalize_pattern Boolean. Whether or not to normalize pattern weights.
#' @return A list with weighted mean differences, mean differences, and differential genes that meet the provided signficance threshold.
#' @export
#' 
#' 
projectionDriveR<-function(
  cellgroup1, #gene x cell count matrix for cell group 1
  cellgroup2, #gene x cell count matrix for cell group 2
  loadings, # a matrix of continous values to be projected with unique rownames
  loadingsNames = NULL, # a vector with names of loadings rows
  pattern_name,
  pvalue = 1e-5,
  display = TRUE,
  normalize_pattern = TRUE
){
  
  #Count matrices can be class matrix, data.frame, sparse.matrix, ... anything that is coercible by as.matrix()
  
  #TODO: assert rownames and colnames exist where needed, and that things are matrices (or can be cast to)

 
  #check that alpha significance level is appropriate
  if(pvalue <= 0 | pvalue >= 1){
    stop("pvalue must be numeric between 0 and 1")
  }
  
  #Make sure provided pattern string is a character vector of length one
  if(length(pattern_name) != 1 | !is.character(pattern_name)){
    stop("provided pattern_name must be a character vector of length one")
  }

  #set loadings rownames if provided
  if(!is.null(loadingsNames)){
    rownames(loadings) <- loadingsNames
  }
  
  #pattern weights must be formatted as a matrix for normalization
  if(pattern_name %in% colnames(loadings)){
    pattern <- loadings[,pattern_name, drop = F] #data.frame
    pattern <- as.matrix(pattern)
  } else  {
    stop(paste0("Provided pattern_name ",pattern_name, " is not a column in provided loadings"))
  }
  
  #Filter the two count matrices and the pattern weights to include the intersection of their features
  #shared rows in two data matrices
  filtered_data <-geneMatchR(data1=cellgroup1, data2=cellgroup2, data1Names=NULL, data2Names=NULL, merge=FALSE)
  print(paste(as.character(dim(filtered_data[[2]])[1]),'row names matched between datasets'))
  
  cellgroup1 <- filtered_data[[2]] #geneMatchR flips the indexes
  cellgroup2 <- filtered_data[[1]]
  
  
  #shared rows in data matrices and loadings
  filtered_weights <- geneMatchR(data1 = cellgroup1, data2 = pattern, data1Names = NULL, data2Names = NULL, merge = F)
  dimensionality_final <- dim(filtered_weights[[2]])[1]
  
  print(paste(as.character(dimensionality_final,'row names matched between data and loadings')))
  print(paste('Updated dimension of data:',as.character(paste(dimensionality_final, collapse = ' '))))
  
  if(dimensionality_final == 0){
    stop("No features matched by rownames of count matrix and rownames of loadings")
  }
  
  pattern_filtered <- filtered_weights[[1]]
  
  cellgroup1_filtered <- filtered_weights[[2]]
  #do second filtering on other cell group so all genes are consistent
  cellgroup2_filtered <- cellgroup2[rownames(cellgroup1_filtered),]
  
  
  #normalize pattern weights
  if(normalize_pattern){
    weight_norm <- norm(pattern_filtered) #square of sums of squares (sum for all positive values)
    num_nonzero <- sum(pattern_filtered > 0) #number of nonzero weights
    pattern_filtered <- pattern_filtered * num_nonzero / weight_norm
  }
  
  #cast feature weights to a named vector
  pattern_normalized_vec <- pattern_filtered[,1]
  names(pattern_normalized_vec) <- rownames(pattern_filtered)
  
  
  
  
  #weighted confidence intervals of differences in cluster means
  weighted_drivers_bonferroni <- bonferroniCorrectedDifferences(group1 = cellgroup1_filtered,
                                group2 = cellgroup2_filtered,
                                diff_weights = pattern_normalized_vec,
                                pvalue = pvalue)
  
  #unweighted confidence intervals of difference in cluster means
  mean_bonferroni <- bonferroniCorrectedDifferences(group1 = cellgroup1_filtered,
                                                    group2 = cellgroup2_filtered,
                                                    diff_weights = NULL,
                                                    pvalue = pvalue)
  
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
  
  if(length(shared_genes) == 0){
    #no genes were significant. Return info we have and skip plotting.
    warning("No features (and weighted features) were significantly differentially used between the two groups")
    return(list(
      mean_ci = mean_bonferroni,
      weighted_mean_ci = weighted_drivers_bonferroni,
      normalized_weights = pattern_normalized_vec,
      significant_genes = shared_genes,
      plotted_ci = NULL))
  }
  
  
  conf_intervals <- mean_bonferroni[shared_genes,]
  sig_weights <- pattern_normalized_vec[shared_genes]
    
  #create confidence interval plot
  pl <- plotConfidenceIntervals(conf_intervals,
                                weights = sig_weights,
                                pattern_name = pattern_name)
  
  if(display){
    #print confidence interval pointrange plot
    print(cowplot::plot_grid(pl[["ci_estimates_plot"]],
                             pl[["weights_heatmap"]],
                             ncol = 2,
                             align = "h",
                             rel_widths = c(1,.3)))
  }
  
  return(list(
    mean_ci = mean_bonferroni,
    weighted_mean_ci = weighted_drivers_bonferroni,
    normalized_weights = pattern_normalized_vec,
    significant_genes = shared_genes,
    plotted_ci = pl))
}

