################################################################################
#' bonferroniCorrectedDifferences
#'
#' Calculate weighted/unweighted mean difference for each gene between 2 groups
#' @param group1 count matrix 1
#' @param group2 count matrix 2
#' @param pvalue significance value to threshold
#' @param diff_weights loadings to weight the differential expression
#' @param mode statistical approach, confidence intervals(CI) or pvalues(PV)
#' @importFrom stats var
#' @importFrom ggrepel geom_label_repel
#' @import dplyr
bonferroniCorrectedDifferences <- function(
    group1,
    group2,
    pvalue,
    diff_weights = NULL,
    mode = "CI") {
  if (!(dim(group1)[1L] == dim(group2)[1L])) {
    #if passed from projectionDrivers, cellgroups will have the same rows
    stop("Rows of two cell group matrices are not identical")
  }

  if (anyNA(group1) || anyNA(group2)) {
    stop("NA values in count matrices not allowed")
  }

  ##Take means over all genes and calculate differences
  group1_mean <- apply(group1, 1L, mean)
  group2_mean <- apply(group2, 1L, mean)
  mean_diff <- group1_mean - group2_mean


  #if weights are provided, use them to weight the difference in means
  if (!is.null(diff_weights)) {

    #check that genes exactly match between difference vector and weight vector
    if (!(all(names(mean_diff) == names(diff_weights)))) {
      stop("Names of loadings and counts do not match")
    }

    mean_diff <- mean_diff * diff_weights
  }

  ##Stats and corrections beginning here
  #calculate confidence intervals
  dimensionality <- length(mean_diff) #number of measurements (genes)

  n1_samples <- dim(group1)[2L] #number of samples (cells)
  n2_samples <- dim(group2)[2L]
  bon_correct <- pvalue / (2L * dimensionality) #bonferroni correction
  qval <- 1L - bon_correct

  tval <- qt(p = qval, df = n1_samples + n2_samples - 2L) #critical value

  group1_var <- apply(group1, 1L, var) #variance of genes across group 1
  group2_var <- apply(group2, 1L, var) #variance of genes across group 2

  if (mode == "CI") {
    #pooled standard deviation
    pool <- ((n1_samples - 1L) * group1_var + (n2_samples - 1L) * group2_var) / (n1_samples + n2_samples - 2L)
    plusminus <- data.frame(low = mean_diff - tval * sqrt(pool * (1L / n1_samples + 1L / n2_samples)),
                            high = mean_diff + tval * sqrt(pool * (1L / n1_samples + 1L / n2_samples)),
                            gene = names(mean_diff))
    rownames(plusminus) <- names(mean_diff)

  } else if (mode == "PV") {
    #welch t test
    #variance calculation
    delta_s <- sqrt((group1_var / n1_samples) + (group2_var / n2_samples))
    #welch t statistic, rounded to 10 digits to avoid infinite decimals
    welch_estimate <- round(mean_diff / delta_s, digits = 10L)
    #Welch-Satterthwaite equation for degrees of freedom
    df <- (((group1_var / n1_samples) + (group2_var / n2_samples)) ^ 2L) /
      ((((group1_var / n1_samples) ^ 2L) / (n1_samples - 1L)) +
         (((group2_var / n2_samples) ^ 2L) / (n2_samples - 1L)))
    #calculate p value from estimate/tvalue
    welch_pvalue <- 2L * pt(-abs(welch_estimate), df = df)
    #bonferroni correction
    welch_padj <- p.adjust(welch_pvalue,
                           method = "bonferroni",
                           n = dimensionality)
    #replace p values equal to zero with the smallest machine value possible
    if (min(welch_padj, na.rm=TRUE) <= .Machine$double.xmin) {
      zp <- length(which(welch_padj <= .Machine$double.xmin))
      message(zp, " P value(s) equal 0. Converting values less than ", .Machine$double.xmin, " to minimum possible value...", call. = FALSE)
      welch_padj[welch_padj <= .Machine$double.xmin] <- .Machine$double.xmin
    }
    plusminus <- data.frame(
      ref_mean = group2_mean,
      test_mean = group1_mean,
      mean_diff = mean_diff,
      estimate = welch_estimate,
      welch_pvalue = welch_pvalue,
      welch_padj = welch_padj,
      gene = names(mean_diff)
    )
  } else {
    stop("Invalid mode selection")
  }
  return(plusminus)
}

################################################################################
#' projectionDriveR
#'
#' Calculate weighted expression difference between two groups (group1 - group2)
#'
#' @importFrom cowplot plot_grid
#' @importFrom Matrix as.matrix
#' @param cellgroup1 gene x cell count matrix for cell group 1
#' @param cellgroup2 gene x cell count matrix for cell group 2
#' @param loadings A matrix of continuous values defining the features
#' @param pattern_name column of loadings for which drivers will be calculated
#' @param pvalue confidence level. Default 1e-5
#' @param loadingsNames a vector with names of loading rows defaults to rownames
#' @param display boolean. Whether or not to display confidence intervals
#' @param normalize_pattern Boolean. Whether or not to normalize pattern weights
#' @param mode statistical approach, confidence intervals or pvalues. default CI
#' @return A list with unweighted/weighted mean differences and differential
#' genes that meet the provided signficance threshold.
#' @export
#'
#'
projectionDriveR <- function(
    cellgroup1, #gene x cell count matrix for cell group 1
    cellgroup2, #gene x cell count matrix for cell group 2
    loadings, # a matrix of continous values to be projected with unique rownames
    pattern_name,
    loadingsNames = NULL, # a vector with names of loadings rows
    pvalue = 1e-5,
    display = TRUE,
    normalize_pattern = TRUE,
    mode = "CI"
) {
  message("Mode: ", mode)
  #Count matrices can be anything that is coercible by as.matrix()
  #check that alpha significance level is appropriate
  if (pvalue <= 0L || pvalue >= 1L) {
    stop("pvalue must be numeric between 0 and 1")
  }

  #Make sure provided pattern string is a character vector of length one
  if (length(pattern_name) != 1L || !is.character(pattern_name)) {
    stop("provided pattern_name must be a character vector of length one")
  }

  #set loadings rownames if provided
  if (!is.null(loadingsNames)) {
    rownames(loadings) <- loadingsNames
  }

  #pattern weights must be formatted as a matrix for normalization
  if (pattern_name %in% colnames(loadings)) {
    pattern <- loadings[, pattern_name, drop = FALSE] #data.frame
    pattern <- Matrix::as.matrix(pattern)
  } else {
    stop("Provided pattern_name ", pattern_name, " is not a column in provided loadings")
  }

  #extract names of data objects
  group1name <- deparse(substitute(cellgroup1))

  group2name <- deparse(substitute(cellgroup2))

  #Filter the two count matrices and the pattern weights to include
  #the intersection of their features
  #shared rows in two data matrices
  filtered_data <- geneMatchR(data1 = cellgroup1,
                             data2 = cellgroup2,
                             data1Names = NULL,
                             data2Names = NULL,
                             merge = FALSE)
  message(as.character(dim(filtered_data[[2L]])[1L]),
          " row names matched between datasets")

  cellgroup1 <- filtered_data[[2L]] #geneMatchR flips the indexes
  cellgroup2 <- filtered_data[[1L]]

  #shared rows in data matrices and loadings
  filtered_weights <- geneMatchR(data1 = cellgroup1,
                                 data2 = pattern,
                                 data1Names = NULL,
                                 data2Names = NULL,
                                 merge = FALSE)

  dimensionality_final <- dim(filtered_weights[[2L]])[1L]

  message("Updated dimension of data: ",
          as.character(paste(dimensionality_final, collapse = " ")))

  if (dimensionality_final == 0L) {
    stop("No features matched by rownames of count matrix and rownames of loadings")
  }

  pattern_filtered <- filtered_weights[[1L]]

  cellgroup1_filtered <- filtered_weights[[2L]]
  #do second filtering on other cell group so all genes are consistent
  cellgroup2_filtered <- cellgroup2[rownames(cellgroup1_filtered), ]


  #normalize pattern weights
  if (normalize_pattern) {
    weight_norm <- norm(pattern_filtered) #square of sums of squares
    num_nonzero <- sum(pattern_filtered > 0L) #number of nonzero weights
    pattern_filtered <- pattern_filtered * num_nonzero / weight_norm
  }
  #cast feature weights to a named vector
  pattern_normalized_vec <- pattern_filtered[, 1L]
  names(pattern_normalized_vec) <- rownames(pattern_filtered)

  #weighted confidence intervals of differences in cluster means
  weighted_bonferroni <- bonferroniCorrectedDifferences(
    group1 = cellgroup1_filtered,
    group2 = cellgroup2_filtered,
    diff_weights = pattern_normalized_vec,
    pvalue = pvalue,
    mode = mode)
  #unweighted confidence intervals of difference in cluster means
  mean_bonferroni <- bonferroniCorrectedDifferences(
    group1 = cellgroup1_filtered,
    group2 = cellgroup2_filtered,
    diff_weights = NULL,
    pvalue = pvalue,
    mode = mode)
#generate confidence interval mode
  if (mode == "CI") {
    #Determine which genes have unweighted/ weighted mean difference
    weighted_sig_idx <- apply(weighted_bonferroni[, 1L:2L], 1L, function(interval) {
      (interval[1L] > 0L & interval[2L] > 0L) | (interval[1L] < 0L & interval[2L] < 0L)
    })

    mean_sig_idx <- apply(mean_bonferroni[, 1L:2L], 1L, function(interval) {
      (interval[1L] > 0L & interval[2L] > 0L) | (interval[1] < 0L & interval[2L] < 0L)
    })

    weighted_sig_genes <- weighted_bonferroni[weighted_sig_idx,]
    unweighted_sig_genes <- mean_bonferroni[mean_sig_idx,]
    #genes that are collectively either up or down
    shared_genes <- base::intersect(
      rownames(weighted_bonferroni)[weighted_sig_idx],
      rownames(mean_bonferroni)[mean_sig_idx])
    message("the length of shared genes are: ", length(shared_genes))
    conf_intervals <- mean_bonferroni[shared_genes, ]
    sig_weights <- pattern_normalized_vec[shared_genes]

    weighted_conf_intervals <- weighted_bonferroni[shared_genes, ]
    #create confidence interval plot (unweighted)
    pl <- plotConfidenceIntervals(conf_intervals,
                                  weights = sig_weights,
                                  pattern_name = pattern_name,
                                  weighted = FALSE)
    #weighted
    pl_w <- plotConfidenceIntervals(weighted_conf_intervals,
                                    weights = sig_weights,
                                    pattern_name = pattern_name,
                                    weighted = TRUE)

    plots <- list(unweighted = pl,weighted = pl_w)
    if (display) {
      #print confidence interval pointrange plot
      pl1_u <- (cowplot::plot_grid(pl[["ci_estimates_plot"]],
                                   pl[["weights_heatmap"]],
                                   ncol = 2L,
                                   align = "h",
                                   rel_widths = c(1.0,0.3)))
      pl2_w <- (cowplot::plot_grid(pl_w[["ci_estimates_plot"]],
                                   pl_w[["weights_heatmap"]],
                                   ncol = 2L,
                                   align = "h",
                                   rel_widths = c(1.0,0.3)))
      plt <- cowplot::plot_grid(pl1_u, pl2_w, ncol = 2L, align = "h")
      print(plt)
    }

    if (length(shared_genes) == 0) {
      #no genes were significant. Return info we have and skip plotting.
      warning("No features were significantly differentially used", 
              call. = FALSE)
      
      result <- list(mean_ci = mean_bonferroni,
                     weighted_mean_ci = weighted_bonferroni,
                     normalized_weights = pattern_normalized_vec,
                     significant_shared_genes = shared_genes,
                     plotted_ci = NULL,
                     meta_data = list(reference_matrix = paste0(group2name),
                                      test_matrix = paste0(group1name))
                     )
      return(result)
    }

    result <- list(
      mean_ci = mean_bonferroni,
      weighted_mean_ci = weighted_bonferroni,
      normalized_weights = pattern_normalized_vec,
      sig_genes = list(unweighted_sig_genes = rownames(unweighted_sig_genes),
                       weighted_sig_genes = rownames(weighted_sig_genes),
                       significant_shared_genes = shared_genes),
      plotted_ci = plots,
      meta_data = list(reference_matrix = paste0(group2name),
                       test_matrix = paste0(group1name))
      )
  } else if (mode == "PV") {
    #create vector of significant genes from weighted and unweighted tests
    weighted_PV_sig <- rownames(weighted_bonferroni[which(weighted_bonferroni$welch_padj <= pvalue),])
    PV_sig <- rownames(mean_bonferroni[which(mean_bonferroni$welch_padj <= pvalue),])
    #create vector of significant genes shared between weighted and unweighted tests
    shared_genes_PV <- base::intersect(
      PV_sig, weighted_PV_sig)
    if (length(shared_genes_PV) == 0L){ 
      #no genes were significant. Return info we have and skip plotting.
      warning("No features were significantly differentially used ",
              call. = FALSE)
      result <- list(mean_stats = mean_bonferroni,
                     weighted_mean_stats = weighted_bonferroni,
                     normalized_weights = pattern_normalized_vec,
                     meta_data = list(reference_matrix = paste0(group2name),
                                      test_matrix = paste0(group1name),
                                      pvalue = pvalue)
      )
      return(result)
    }
    result <- list(mean_stats = mean_bonferroni,
                   weighted_mean_stats = weighted_bonferroni,
                   normalized_weights = pattern_normalized_vec,
                   sig_genes = list(PV_sig = PV_sig,
                                    weighted_PV_sig = weighted_PV_sig,
                                    PV_significant_shared_genes = shared_genes_PV),
                   meta_data = list(reference_matrix = paste0(group2name),
                                    test_matrix = paste0(group1name),
                                    pvalue = pvalue)
                   )
    #apply pdVolcano function to result
    result <- pdVolcano(result, display = FALSE)
    if (display) {
      #print volcano plots
      pltgrid <- cowplot::plot_grid(result$plt$differential_expression +
                                      theme(legend.position = "none"),
                                result$plt$weighted_differential_expression +
                                  theme(legend.position = "none"),
                                ncol = 2L, align = "h")
      legend <- cowplot::get_plot_component(result$plt$differential_expression, "guide-box-bottom")
      plt <- cowplot::plot_grid(pltgrid, legend, ncol = 1L, rel_heights = c(1.0, 0.1))
      print(plt)
    }
  } else {
    stop("Invalid mode selection")
  }
  return(result)
}
