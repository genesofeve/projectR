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
#' @importFrom dplyr %>% mutate dense_rank
#' @param confidence_intervals A dataframe of features x estimates. 
#' @param interval_name names of columns that contain the low and high estimates, respectively. Default: c("low","high")
#' @param sort Boolean. Whether or not to sort genes by their estimates (default = T)
#' @param genes a vector with names of genes to include in plot. If sort=F, estimates will be plotted in this order.
#' @param weights optional. weights of features to include as annotation.
#' @return A list with pointrange estimates and, if requested, a heatmap of pattern weights.
#' @export
plotConfidenceIntervals <- function(
  confidence_intervals, #confidence_interval is a data.frame or matrix with two columns (low, high). Genes must be rownames
  interval_name = c("low","high"),
  sort = T,
  genes = NULL,
  weights = NULL){
  
  #gene names were stored as rownames, make sure high and low estimates are stored
  confidence_intervals$gene_names <- rownames(confidence_intervals)
  confidence_intervals$low <- confidence_intervals[,interval_name[1]]
  confidence_intervals$high <- confidence_intervals[,interval_name[2]]
  
  
  n <- dim(confidence_intervals)[1]
  confidence_intervals <- confidence_intervals %>%
    mutate(
      mid = (high+low)/2,
      positive = mid > 0)
  
  if(!is.null(genes)){
    #select genes provided and get them in that order
    if(!(is.character(genes))){ stop("Genes must be provided as a character vector") }
    message(paste0("Selecting ", length(genes), " features"))
    confidence_intervals <- confidence_intervals[genes,]
    
  }
  
  if(sort){
    #order in increasing order on estimates
    confidence_intervals <- confidence_intervals %>% 
      mutate(
        idx = dense_rank(mid)
      )
  } else{ 
      #if not sorted, create index variable for current order
      confidence_intervals <- confidence_intervals %>% 
        mutate(idx = 1:n)
      
  }
  
  ci_plot <- ggplot(data = confidence_intervals, aes(y = idx, x = mid)) + geom_pointrange(aes(xmin = low, xmax = high, color = positive)) +
    geom_point(aes(x = mid, y = idx), fill ="black",color = "black") +
    theme_minimal() + 
    xlab("Difference in group means") + 
    ylab("Genes") + 
    geom_vline(xintercept = 0, color = "black", linetype = "dashed") + 
    theme(legend.position = "none")
 
  if(!is.null(weights)){
    
    
  }
  
  return(ci_plot)
}  
