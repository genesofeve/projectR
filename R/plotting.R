# Example plots to be functionalized
# devtools::install_github('rlbarter/superheat')
#
# test<-projectR(data=p.ESepiGen4c1l$mRNA.Seq,Patterns=AP.RNAseq6l3c3t$Amean,AnnotionObj=map.ESepiGen4c1l,IDcol="GeneSymbols",full=TRUE)
#
# tmp<-matrix("",nrow = 5,ncol=9)
# tmp[test$pval<0.01]<-"*"
#
# superheat(test$projection,row.dendrogram=TRUE, pretty.order.cols = TRUE,
#           heat.pal.values = c(0, 0.5, 1),yt=colSums(test$projection),yt.plot.type='scatterline',yt.axis.name="Sum of\nProjections",X.text=tmp,X.text.size=8,bottom.label.text.angle = 90)
#
#######################################################################################################################################
#' 
#' plotConfidenceIntervals
#' 
#' Generate point and line confidence intervals from provided estimates. 
#' 
#' @import ggplot2
#' @import viridis
#' @importFrom scales squish
#' @importFrom dplyr %>% mutate dense_rank
#' @importFrom cowplot plot_grid
#' @param confidence_intervals A dataframe of features x estimates. 
#' @param interval_name names of columns that contain the low and high estimates, respectively. Default: c("low","high")
#' @param pattern_name string to use as the title for plots.
#' @param sort Boolean. Whether or not to sort genes by their estimates (default = T)
#' @param genes a vector with names of genes to include in plot. If sort=F, estimates will be plotted in this order.
#' @param weights optional. weights of features to include as annotation. 
#' @param weights_clip optional. quantile of data to clip color scale for improved visualization. Default: 0.99
#' @param weights_vis_norm Which processed version of weights to visualize as a heatmap.
#' @param weighted specifies whether the confidence intervals in use are weighted by the pattern and labels plots accordingly
#' Options are "none" (which uses provided weights) or "quantiles". Default: none
#' @return A list with pointrange estimates and, if requested, a heatmap of pattern weights.
#' @export
plotConfidenceIntervals <- function(
    confidence_intervals,
    interval_name = c("low","high"),
    pattern_name = NULL,
    sort = T,
    genes = NULL,
    weights = NULL,
    weights_clip = 0.99,
    weights_vis_norm = "none",
    weighted = F){
  
  if(weights_clip < 0 | weights_clip > 1){
    stop("weights_clip must be numeric between 0 and 1")
  }
  
  if(!(weights_vis_norm %in% c("none","quantiles"))){
    stop("weights_vis_norm must be either 'none' or 'quantiles'")
  }
  
  if(weighted == F){
    lab = "Unweighted"
  } else{
    lab = "Weighted"
  }
  #gene names were stored as rownames, make sure high and low estimates are stored
  confidence_intervals$gene_names <- rownames(confidence_intervals)
  confidence_intervals$low <- confidence_intervals[,interval_name[1]]
  confidence_intervals$high <- confidence_intervals[,interval_name[2]]
  
  n <- dim(confidence_intervals)[1]
  confidence_intervals <- confidence_intervals %>%
    mutate(
      mid = (high+low)/2, #estimate, used for point position
      positive = mid > 0) #upregulated, used for color scheme
  
  if(!is.null(genes)){
    #select genes provided and get them in that order
    if(!(is.character(genes))){ stop("Genes must be provided as a character vector") }
    n <- length(genes)
    message(paste0("Selecting ", n, " features"))
    confidence_intervals <- confidence_intervals[genes,]
    
  }
  
  if(sort){
    #order in increasing order on estimates, and create index variable
    message("sorting genes in increasing order of estimates...")
    confidence_intervals <- confidence_intervals %>% mutate(idx = dense_rank(mid)) %>%
      arrange(mid)
    
  } else{
    #if not sorted, create index variable for current order
    confidence_intervals <- confidence_intervals %>% mutate(idx = 1:n)
  }
  
  #genereate point range plot
  ci_plot <- ggplot(data = confidence_intervals, aes(y = idx, x = mid)) + geom_pointrange(aes(xmin = low, xmax = high, color = positive)) +
    geom_point(aes(x = mid, y = idx), fill ="black",color = "black") +
    theme_minimal() +
    xlab("Difference in group means") +
    ylab("Genes") +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
    theme(legend.position = "none") +
    ggtitle(lab)
  
  #if provided, create heatmap for pattern weights
  if(!is.null(weights)){
    
    #check that weights are formatted as a named vector
    if(!(is.numeric(weights))){ stop("Weights must be provided as a numeric vector") }
    if(is.null(names(weights))){ stop("Weights must have names that match estimates")}
    
    #either use pattern_name, or if not provided, just label heatmap with "weights"
    hm_name <- ifelse(is.null(pattern_name), "weights", pattern_name)
    
    #maintain established order from the pointrange plot
    ordered_weights <- weights[rownames(confidence_intervals)]
    
    if(weights_vis_norm == "quantiles"){
      #transform to percentiles from 0 to 1
      ordered_weights <- trunc(rank(ordered_weights))/length(ordered_weights)
      hm_name <- paste0(hm_name, " (quantiles)") #append quantile to plot name
    }
    
    confidence_intervals$weights <- ordered_weights
    
    #generate heatmap
    wt_heatmap <- ggplot(data = confidence_intervals) +
      geom_tile(aes(x = 1, y = 1:n, fill = weights)) +
      scale_fill_viridis(limits=c(0, quantile(ordered_weights,weights_clip )),
                         oob=squish,
                         name = hm_name) +
      theme_void()
    
  } else{ wt_heatmap = NULL} #if weights aren't provided, return NULL
  
  return(list("ci_estimates_plot" = ci_plot,
              "feature_order" = rownames(confidence_intervals),
              "weights_heatmap" = wt_heatmap))
}

#######################################################################################################################################
#' pdVolcano
#'
#' Generate volcano plot and gate genes based on fold change and pvalue, includes vectors that can be used with fast gene set enrichment (fgsea)
#' @param result result output from projectionDriveR function with PI method selected
#' @param FC fold change threshold, default at 0.2
#' @param pvalue significance threshold, default set to pvalue stored in projectionDriveR output
#' @param subset vector of gene names to subset the plot by
#' @param filter.inf remove genes that have pvalues below machine double minimum value
#' @param label.num Number of genes to label on either side of the volcano plot, default 5
#' @import ggpubr
#' @importFrom stats var
#' @importFrom ggrepel geom_label_repel
#' @importFrom ggrepel geom_text_repel
#' @import dplyr
#' @return A list with weighted and unweighted differential expression metrics
#' @export
#plot FC, weighted and unweighted. Designed for use with the output of projectionDriveRs
pdVolcano <- function(result, 
                      FC = 0.2,
                      pvalue = NULL,
                      subset = NULL, 
                      filter.inf = FALSE,
                      label.num = 5) {
  
  #if a genelist is provided, use them to subset the output of projectiondrivers
  if(!is.null(subset)){
    #subset the mean_stats object by provided gene list
    result$mean_stats <- result$mean_stats[(which(rownames(result$mean_stats) %in% subset)),]
    #subset the weighted_mean_stats object by provided gene list
    result$weighted_mean_stats <- result$weighted_mean_stats[(which(rownames(result$weighted_mean_stats) %in% subset)),]
    
  }
  
  if(filter.inf == TRUE){
    #remove p values below the machine limit representation for plotting purposes 
    cat("Filtering", length(which(result$mean_stats$welch_padj <= .Machine$double.xmin)),"unweighted genes and", 
        length(which(result$weighted_mean_stats$welch_padj <= .Machine$double.xmin)), "weighted genes", "\n")
    result$mean_stats <- subset(result$mean_stats, welch_padj > .Machine$double.xmin)
    result$weighted_mean_stats <- subset(result$weighted_mean_stats, welch_padj > .Machine$double.xmin)
  }
  
  if(is.numeric(FC) == FALSE){
    stop('FC must be a number')
  }
  
  if(is.null(pvalue) == FALSE) {
    pvalue = pvalue
  } else {
      pvalue <- result$meta_data$pvalue
  }
  
  #extract object meta data 
  metadata <- result$meta_data
  
  #volcano plotting unweighted
  #extract unweighted confidence intervals / statistics
  mean_stats <- result$mean_stats
  #fold change / significance calls
  mean_stats$Color <- paste("NS or FC", FC)
  mean_stats$Color[mean_stats$welch_padj < pvalue & mean_stats$mean_diff > FC] <- paste("Enriched in", metadata$test_matrix)
  mean_stats$Color[mean_stats$welch_padj < pvalue & mean_stats$mean_diff < -FC] <- paste("Enriched in", metadata$reference_matrix)
  mean_stats$Color <- factor(mean_stats$Color,
                             levels = c(paste("NS or FC <", FC), paste("Enriched in", metadata$reference_matrix), paste("Enriched in", metadata$test_matrix)))
  
  #label the most significant genes for enrichment
  mean_stats$invert_P <- (-log10(mean_stats$welch_padj)) * (mean_stats$mean_diff)
  
  top_indices <- order(mean_stats$invert_P, decreasing = TRUE)[1:label.num]
  bottom_indices <- order(mean_stats$invert_P)[1:label.num]
  
  # Add labels to the dataframe
  mean_stats$label <- NA
  mean_stats$label[top_indices] <- paste(rownames(mean_stats)[top_indices])
  mean_stats$label[bottom_indices] <- paste(rownames(mean_stats)[bottom_indices])
  
  #set custom colors
  myColors <- c("gray","red","dodgerblue")
  names(myColors) <- levels(mean_stats$Color)
  custom_colors <- scale_color_manual(values = myColors, drop = FALSE)
  
  #plot
  unweightedvolcano = ggplot(data = mean_stats, aes(x = mean_diff, y = -log10(welch_padj), color = Color, label = label)) + 
    geom_vline(xintercept = c(FC, -FC), lty = "dashed") +
    geom_hline(yintercept = -log10(pvalue), lty = "dashed") +
    geom_point(na.rm = TRUE) + 
    custom_colors + 
    coord_cartesian(ylim = c(0, 350), xlim = c(-2, 2)) +
    ggrepel::geom_text_repel(size = 3, point.padding = 1, color = "black",
                             min.segment.length = .1, box.padding = 0.15,
                             max.overlaps = Inf, na.rm = TRUE) +
    labs(x = "FC",
         y = "Significance (-log10pval)",
         color = NULL) +
    ggtitle("Differential Expression") + 
    theme_bw() +
    theme(plot.title = element_text(size = 16),
          legend.position = "bottom",
          axis.title=element_text(size=14),
          legend.text = element_text(size=12))
  
  #weighted volcano plot
  weighted_mean_stats <- result$weighted_mean_stats
  weighted_mean_stats$Color <- paste("NS or FC <", FC)
  weighted_mean_stats$Color[weighted_mean_stats$welch_padj < pvalue & weighted_mean_stats$mean_diff > FC] <- paste("Enriched in", metadata$test_matrix)
  weighted_mean_stats$Color[weighted_mean_stats$welch_padj < pvalue & weighted_mean_stats$mean_diff < -FC] <- paste("Enriched in", metadata$reference_matrix)
  weighted_mean_stats$Color <- factor(weighted_mean_stats$Color,
                                      levels = c(paste("NS or FC <", FC), paste("Enriched in", metadata$reference_matrix), paste("Enriched in", metadata$test_matrix)))
  
  weighted_mean_stats$invert_P <- (-log10(weighted_mean_stats$welch_padj)) * (weighted_mean_stats$mean_diff)
  
  
  top_indices <- order(weighted_mean_stats$invert_P, decreasing = TRUE)[1:label.num]
  bottom_indices <- order(weighted_mean_stats$invert_P)[1:label.num]
  
  # Add labels to the dataframe
  weighted_mean_stats$label <- NA
  weighted_mean_stats$label[top_indices] <- paste(rownames(weighted_mean_stats)[top_indices])
  weighted_mean_stats$label[bottom_indices] <- paste(rownames(weighted_mean_stats)[bottom_indices])
  
  myColors <- c("gray","red","dodgerblue")
  names(myColors) <- levels(weighted_mean_stats$Color)
  custom_colors <- scale_color_manual(values = myColors, drop = FALSE)
  
  weightedvolcano = ggplot(data = weighted_mean_stats, aes(x = mean_diff, y = -log10(welch_padj), color = Color, label = label)) + 
    geom_vline(xintercept = c(FC, -FC), lty = "dashed") +
    geom_hline(yintercept = -log10(pvalue), lty = "dashed") +
    geom_point(na.rm = TRUE) + 
    custom_colors +
    coord_cartesian(ylim = c(0, 350), xlim = c(-2, 2)) +
    ggrepel::geom_text_repel(size = 3, point.padding = 1, color = "black",
                             min.segment.length = .1, box.padding = 0.15,
                             max.overlaps = Inf, na.rm = TRUE) +
    labs(x = "FC",
         y = "Significance (-log10pval)",
         color = NULL) +
    ggtitle("Weighted Differential Expression") + 
    theme_bw() +
    theme(plot.title = element_text(size = 16),
          legend.position = "bottom",
          axis.title=element_text(size=14),
          legend.text = element_text(size=12))
  
  plt <- ggpubr::ggarrange(unweightedvolcano, weightedvolcano, common.legend = TRUE, legend = "bottom")
  print(plt)
  
  #return a list of genes that can be used as input to fgsea
  difexdf <- subset(mean_stats, Color == paste("Enriched in", metadata$reference_matrix) | Color == paste("Enriched in", metadata$test_matrix))
  vec <- difexdf$estimate
  names(vec) <- rownames(difexdf)
  
  weighted_difexdf <- subset(weighted_mean_stats, Color == paste("Enriched in", metadata$reference_matrix) | Color == paste("Enriched in", metadata$test_matrix))
  weighted_vec <- weighted_difexdf$estimate
  names(weighted_vec) <- rownames(weighted_difexdf)
  names(vec) <- rownames(difexdf)
  vol_result <- list(mean_stats = mean_stats,
                     weighted_mean_stats = weighted_mean_stats,
                     sig_genes = result$sig_genes,
                     difexpgenes = difexdf,
                     weighted_difexpgenes = weighted_difexdf,
                     fgseavecs = list(unweightedvec = vec,
                                      weightedvec = weighted_vec),
                     meta_data = metadata,
                     plt = plt)
  return(vol_result)
}



