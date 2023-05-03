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
#' Options are "none" (which uses provided weights) or "quantiles". Default: none
#' @return A list with pointrange estimates and, if requested, a heatmap of pattern weights.
#' @export
plotConfidenceIntervals <- function(
  confidence_intervals, #confidence_interval is a data.frame or matrix with two columns (low, high). Genes must be rownames
  interval_name = c("low","high"),
  pattern_name = NULL,
  sort = T,
  genes = NULL,
  weights = NULL,
  weights_clip = 0.99,
  weights_vis_norm = "none"){
  
  if(weights_clip < 0 | weights_clip > 1){
    stop("weights_clip must be numeric between 0 and 1")
  }
  
  if(!(weights_vis_norm %in% c("none","quantiles"))){
    stop("weights_vis_norm must be either 'none' or 'quantiles'")
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
    confidence_intervals <- confidence_intervals %>% 
      mutate(
        idx = dense_rank(mid)) %>%
      arrange(mid)
    
  } else{ 
      #if not sorted, create index variable for current order
      confidence_intervals <- confidence_intervals %>% 
        mutate(idx = 1:n)
  }
  
  #genereate point range plot
  ci_plot <- ggplot(data = confidence_intervals, aes(y = idx, x = mid)) + geom_pointrange(aes(xmin = low, xmax = high, color = positive)) +
    geom_point(aes(x = mid, y = idx), fill ="black",color = "black") +
    theme_minimal() + 
    xlab("Difference in group means") + 
    ylab("Genes") + 
    geom_vline(xintercept = 0, color = "black", linetype = "dashed") + 
    theme(legend.position = "none") + 
    ggtitle(pattern_name)
 
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
