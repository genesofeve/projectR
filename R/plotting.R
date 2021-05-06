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
#' Generate point and line confidence intervals from provided estimates
#' 
#' @import ggplot2
#' @importFrom dplyr %>% mutate dense_rank
#' @param confidence_intervals a vector with names of loading rows. Defaults to rownames.
#' @param interval_name names of columns that contain the low and high estimates, respectively
#' @param sort Boolean. Whether or not to sort genes by their estimates (default = T)
#' @param feature_weights optional. weights of features to include as annotation.
#' @export
plotConfidenceIntervals <- function(
  confidence_intervals, #confidence_interval is a data.frame or matrix with two columns (low, high). Genes must be rownames
  interval_name = c("low","high"),
  sort = T){
  
  #gene names were stored as rownames, make sure high and low estimates are stored
  confidence_intervals$gene_names <- rownames(confidence_intervals)
  confidence_intervals$low <- confidence_intervals[,interval_name[1]]
  confidence_intervals$high <- confidence_intervals[,interval_name[2]]
  
  
  n <- dim(confidence_intervals)[1]
  confidence_intervals <- confidence_intervals %>%
    mutate(
      mid = (high+low)/2,
      positive = mid > 0)
  
  #order in descending order on estimates
  if(sort){
    confidence_intervals <- confidence_intervals %>% 
      mutate(
        idx = dense_rank(mid)
      )
  }
  
  ggplot(data = confidence_intervals, aes(y = idx, x = mid)) + geom_pointrange(aes(xmin = low, xmax = high, color = positive)) +
    theme_minimal() + 
    xlab("Difference in group means") + 
    ylab("Genes") + 
    geom_vline(xintercept = 0, color = "black", linetype = "dashed") + 
    theme(legend.position = "none")
  #TODO: add gene names to y axis labels
  #TODO: add annotation with pattern weights
  
}  
