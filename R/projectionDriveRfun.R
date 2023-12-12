#######################################################################################################################################
#' bonferroniCorrectedDifferences
#'
#' Calculate the weighted and unweighted difference in means for each measurement between two groups.
#' @param group1 count matrix 1
#' @param group2 count matrix 2
#' @param diff_weights loadings to weight the differential expression between the groups
#' @param pvalue significance value to threshold at
#' @importFrom stats var
#' @importFrom ggrepel geom_label_repel
#' @import dplyr
bonferroniCorrectedDifferences <- function(
    group1,
    group2, 
    diff_weights = NULL, 
    pvalue) 

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
  mean_diff <- group1_mean - group2_mean #if this is log normalized counts the mean difference is actually log(group1/group2)
  
  
  #if weights are provided, use them to weight the difference in means
  if(!is.null(diff_weights)){
    
    #check that genes exactly match between difference vector and weight vector
    if(!(all(names(mean_diff) == names(diff_weights)))){
      stop("Names of loadings and counts do not match")
    }
    
    mean_diff <- mean_diff * diff_weights
  }
  
  
  
  ##Stats and corrections beginning here
  #calculate confidence intervals
  dimensionality <- length(mean_diff) #number of measurements (genes)
  
  n1_samples <- dim(group1)[2] #number of samples (cells)
  n2_samples <- dim(group2)[2]
  bon_correct <- pvalue / (2*dimensionality) #bonferroni correction
  qval <- 1 - bon_correct
  
  tval <- qt(p = qval, df = n1_samples + n2_samples -2) #critical value
  
  #calculate p values
  group1_var <- apply(group1, 1, var) #variance of genes across group 1 
  group2_var <- apply(group2, 1, var) #variance of genes across group 2
  
  #vartest <- group1_var / group2_var #test to see if variance across groups is equal, often not equal
  
  pooled <- ((n1_samples-1)*group1_var + (n2_samples-1)*group2_var)  /  (n1_samples+n2_samples-2) #pooled standard deviation
  
  #welch t test
  deltaS <- sqrt((group1_var / n1_samples) + (group2_var / n2_samples)) #variance calculation
  
  welch_estimate <- round(mean_diff / deltaS, digits = 10) #welch t statistic, rounded to 10 digits to avoid infinite decimals
  
  df <- (((group1_var / n1_samples) + (group2_var / n2_samples))^2) / ((((group1_var / n1_samples)^2) / (n1_samples - 1)) + (((group2_var / n2_samples)^2) / (n2_samples - 1))) #Welch-Satterthwaite equation for degrees of freedom
  
  welch_pvalue <- 2*pt(-abs(welch_estimate), df=df) #calculate p value from estimate (tvalue)
  
  welch_p_value_boncorrected <- p.adjust(welch_pvalue, method = "bonferroni", n = dimensionality) #bonferroni correction
  
  # replace p values equal to zero with the smallest machine value possible
  if (min(welch_p_value_boncorrected, na.rm=TRUE) <= .Machine$double.xmin) {
    zp <- length(which(welch_p_value_boncorrected <= .Machine$double.xmin))
    warning(paste(zp,"P value(s) equal 0.",
                  "Converting any values less than", .Machine$double.xmin, "to minimum possible value..."),
            call. = FALSE)
    welch_p_value_boncorrected[welch_p_value_boncorrected <= .Machine$double.xmin] <- .Machine$double.xmin
  }
  
  #establish dataframe to populate in the following for loop
  plusminus = data.frame(low = rep(NA_integer_, dimensionality), 
                         high = rep(NA_integer_, dimensionality), 
                         estimate = rep(NA_integer_, dimensionality),
                         mean_diff = rep(NA_integer_, dimensionality),
                         welch_pvalue = rep(NA_integer_, dimensionality),
                         welch_p_value_boncorrected = rep(NA_integer_, dimensionality),
                         ref_mean = rep(NA_integer_, dimensionality),
                         test_mean = rep(NA_integer_, dimensionality),
                         gene = rep(NA_integer_, dimensionality))
  rownames(plusminus) <- names(mean_diff)
  
  
  #for each gene, calculate confidence interval around mean
  for(i in 1:dimensionality){
    
    scale = tval * sqrt(pooled[i] * (1/n1_samples + 1/n2_samples))
    
    
    plusminus[i, "low"] <- mean_diff[i]  - scale #low estimate
    plusminus[i, "high"] <- mean_diff[i] + scale #high estimate
    plusminus[i, "estimate"] <- welch_estimate[i]  
    plusminus[i, "mean_diff"] <- mean_diff[i]
    plusminus[i, "welch_pvalue"] <- welch_pvalue[i]
    plusminus[i, "welch_p_value_boncorrected"] <- welch_p_value_boncorrected[i]
    plusminus[i, "ref_mean"] <- group2_mean[i]
    plusminus[i, "test_mean"] <- group1_mean[i]   
    plusminus[i, "gene"] <- names(mean_diff[i])
    
    
  }
  return(plusminus)
}


#######################################################################################################################################
#' projectionDriveR
#'
#' Calculate the weighted difference in expression between two groups (group1 - group2)
#'
#' @importFrom cowplot plot_grid
#' @importFrom ggpubr ggarrange
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
  
  #extract names of data objects
  group1name <- deparse(substitute(cellgroup1))

  group2name <- deparse(substitute(cellgroup2))

  
  #Filter the two count matrices and the pattern weights to include the intersection of their features
  #shared rows in two data matrices
  filtered_data <-geneMatchR(data1=cellgroup1, data2=cellgroup2, data1Names=NULL, data2Names=NULL, merge=FALSE)
  print(paste(as.character(dim(filtered_data[[2]])[1]),'row names matched between datasets'))
  
  cellgroup1 <- filtered_data[[2]] #geneMatchR flips the indexes
  cellgroup2 <- filtered_data[[1]]
  
  
  #shared rows in data matrices and loadings
  filtered_weights <- geneMatchR(data1 = cellgroup1, data2 = pattern, data1Names = NULL, data2Names = NULL, merge = F)
  dimensionality_final <- dim(filtered_weights[[2]])[1]
  
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
  weighted_welch_sig <- rownames(weighted_drivers_bonferroni[which(weighted_drivers_bonferroni$welch_p_value_boncorrected <= pvalue),])
  
  # Apply the t test function to unweighted expression matrices and call significance
  n_tests <- nrow(cellgroup1_filtered)
  cat("ntests are:", n_tests, "\n")
  
  
  #unweighted confidence intervals of difference in cluster means
  mean_bonferroni <- bonferroniCorrectedDifferences(group1 = cellgroup1_filtered,
                                                    group2 = cellgroup2_filtered,
                                                    diff_weights = NULL,
                                                    pvalue = pvalue)
  welch_sig <- rownames(mean_bonferroni[which(mean_bonferroni$welch_p_value_boncorrected <= pvalue),])
  
  #Determine which genes have significantly non-zero mean difference and weighted mean difference
  #significant
  weighted_sig_idx <- apply(weighted_drivers_bonferroni[,1:2], 1, function(interval){
    (interval[1] > 0 & interval[2] > 0) | (interval[1] < 0 & interval[2] < 0)
  })
  
  weighted_sig_genes <- weighted_drivers_bonferroni[weighted_sig_idx,]
  
  mean_sig_idx <- apply(mean_bonferroni[,1:2], 1, function(interval){
    (interval[1] > 0 & interval[2] > 0) | (interval[1] < 0 & interval[2] < 0)
  })
  
  unweighted_sig_genes <- mean_bonferroni[mean_sig_idx,]
  
  
  #genes that are collectively either up or down
  shared_genes <- base::intersect(
    rownames(weighted_drivers_bonferroni)[weighted_sig_idx],
    rownames(mean_bonferroni)[mean_sig_idx])
  cat("the length of shared genes are:", length(shared_genes), '\n')
  
  shared_genes2 <- base::intersect(
    welch_sig, weighted_welch_sig)

  
  if(length(shared_genes) == 0){
    #no genes were significant. Return info we have and skip plotting.
    warning("No features (and weighted features) were significantly differentially used between the two groups")
    return(list(
      mean_ci = mean_bonferroni,
      weighted_mean_ci = weighted_drivers_bonferroni,
      significant_shared_genes = shared_genes,
      plotted_ci = NULL,
      weighted_sig_genes = rownames(weighted_sig_genes),
      unweighted_sig_genes = rownames(unweighted_sig_genes),
      reference_matrix = paste0(group2name),
      test_matrix = paste0(group1name),
      welch_sig = welch_sig,
      weighted_welch_sig = weighted_welch_sig, 
      welch_significant_shared_genes = shared_genes2,
      pvalue = pvalue))
  }
  
  
  conf_intervals <- mean_bonferroni[shared_genes,]
  sig_weights <- pattern_normalized_vec[shared_genes]
  
  weighted_conf_intervals <- weighted_drivers_bonferroni[shared_genes,]

  #create confidence interval plot (unweighted)
  pl <- plotConfidenceIntervals(conf_intervals,
                                weights = sig_weights,
                                pattern_name = pattern_name,
                                weighted = F)
  #weighted
  pl_w <- plotConfidenceIntervals(weighted_conf_intervals,
                                weights = sig_weights,
                                pattern_name = pattern_name,
                                weighted = T)
  plots <- list(unweighted = pl,weighted = pl_w)
  if(display){
    #print confidence interval pointrange plot
    pl1_u <- (cowplot::plot_grid(pl[["ci_estimates_plot"]],
                             pl[["weights_heatmap"]],
                             ncol = 2,
                             align = "h",
                             rel_widths = c(1,.3)))
    print(pl1_u)
    pl2_w <- (cowplot::plot_grid(pl_w[["ci_estimates_plot"]],
                                 pl_w[["weights_heatmap"]],
                                 ncol = 2,
                                 align = "h",
                                 rel_widths = c(1,.3)))
    print(pl2_w)
    plt <- ggpubr::ggarrange(pl1_u, pl2_w, common.legend = TRUE, legend = "bottom")
    print(plt)
  }
  
  return(list(
    mean_ci = mean_bonferroni,
    weighted_mean_ci = weighted_drivers_bonferroni,
    significant_shared_genes = shared_genes,
    plotted_ci = plots,
    weighted_sig_genes = rownames(weighted_sig_genes),
    unweighted_sig_genes = rownames(unweighted_sig_genes),
    reference_matrix = paste0(group2name),
    test_matrix = paste0(group1name),
    welch_sig = welch_sig,
    weighted_welch_sig = weighted_welch_sig, 
    welch_significant_shared_genes = shared_genes2,
    pvalue = pvalue))
}

