################################################################################
#' plotConfidenceIntervals
#'
#' Generate point and line confidence intervals from provided estimates.
#' @import ggplot2
#' @import viridis
#' @importFrom scales squish
#' @importFrom dplyr %>% mutate dense_rank
#' @importFrom cowplot plot_grid
#' @param confidence_intervals A dataframe of features x estimates.
#' @param interval_name Estimate column names. Default: c("low","high")
#' @param pattern_name string to use as the title for plots.
#' @param sort Boolean. Sort genes by their estimates (default = TRUE)
#' @param genes a vector with names of genes to include in plot.
#' If sort=F, estimates will be plotted in this order.
#' @param weights optional. weights of features to include as annotation.
#' @param weights_clip optional. quantile of data to clip color scale for
#' improved visualization. Default: 0.99
#' @param weights_vis_norm Which version of weights to visualize as a heatmap.
#' Options are "none" (uses provided weights) or "quantiles". Default: none
#' @param weighted specifies whether the confidence intervals in use are
#' weighted by the pattern and labels plots accordingly
#' @return A list with pointrange estimates and a heatmap of pattern weights.
#' @export
plotConfidenceIntervals <- function(
    confidence_intervals,
    interval_name = c("low", "high"),
    pattern_name = NULL,
    sort = TRUE,
    genes = NULL,
    weights = NULL,
    weights_clip = 0.99,
    weights_vis_norm = "none",
    weighted = FALSE) {

  #bind variables locally to the function
  mid <- idx <- low <- high <- positive <- NULL
  if (weights_clip < 0L || weights_clip > 1L) {
    stop("weights_clip must be numeric between 0 and 1")
  }

  if (!(weights_vis_norm %in% c("none", "quantiles"))) {
    stop("weights_vis_norm must be either 'none' or 'quantiles'")
  }

  if (weighted) {
    lab <- "Weighted"
  } else {
    lab <- "Unweighted"
  }
  #gene names were stored as rownames, store high and low estimates
  confidence_intervals$gene_names <- rownames(confidence_intervals)
  confidence_intervals$low <- confidence_intervals[, interval_name[1L]]
  confidence_intervals$high <- confidence_intervals[, interval_name[2L]]

  n <- dim(confidence_intervals)[1L]
  confidence_intervals <- confidence_intervals %>%
    mutate(
      mid = (high + low) / 2L, #estimate, used for point position
      positive = mid > 0L) #upregulated, used for color scheme

  if (!is.null(genes)) {
    #select genes provided and get them in that order
    if (!(is.character(genes))) {
      stop("Genes must be provided as a character vector")
      }
    n <- length(genes)
    message("Selecting ", n, " features")
    confidence_intervals <- confidence_intervals[genes, ]

  }

  if (sort) {
    #order in increasing order on estimates, and create index variable
    message("sorting genes in increasing order of estimates...")
    confidence_intervals <- confidence_intervals %>% mutate(
      idx = dense_rank(mid)) %>%
      arrange(mid)

  } else {
    #if not sorted, create index variable for current order
    confidence_intervals <- confidence_intervals %>% mutate(idx = 1L:n)
  }

  #generate point range plot
  ci_plot <- ggplot(data = confidence_intervals,
                    aes(y = idx, x = mid)) +
    geom_pointrange(aes(xmin = low,
                        xmax = high,
                        color = positive)) +
    geom_point(aes(x = mid,
                   y = idx),
               fill = "black",
               color = "black") +
    theme_minimal() +
    xlab("Difference in group means") +
    ylab("Genes") +
    geom_vline(xintercept = 0L, color = "black", linetype = "dashed") +
    theme(legend.position = "none") +
    ggtitle(lab)

  #if provided, create heatmap for pattern weights
  if (!is.null(weights)) {

    #label with pattern name if provided
    hm_name <- ifelse(is.null(pattern_name), "weights", pattern_name)

    #maintain established order from the pointrange plot
    ordered_weights <- weights[rownames(confidence_intervals)]
    confidence_intervals$weights <- ordered_weights
    
    #generate heatmap
    wt_heatmap <- ggplot2::ggplot(data = confidence_intervals) +
      geom_tile(aes(x = 1L, y = 1L:n, fill = weights)) +
      scale_fill_viridis(limits = c(0L, quantile(ordered_weights, weights_clip)),
                         oob = scales::squish,
                         name = hm_name) +
      theme_void()

    #check that weights are formatted as a named vector
    if (!(is.numeric(weights))) {
      stop("Weights must be provided as a numeric vector")
      }
    if (is.null(names(weights))) {
      stop("Weights must have names that match estimates")
      }

    if (weights_vis_norm == "quantiles") {
      #transform to percentiles from 0 to 1
      ordered_weights <- trunc(rank(ordered_weights)) / length(ordered_weights)
      hm_name <- paste0(hm_name, " (quantiles)") #append quantile to plot name
    }

  } else {
    #if weights aren't provided, return NULL
    return(list("ci_estimates_plot" = ci_plot,
                "feature_order" = rownames(confidence_intervals),
                "weights_heatmap" = NULL))
  }

  return(list("ci_estimates_plot" = ci_plot,
              "feature_order" = rownames(confidence_intervals),
              "weights_heatmap" = wt_heatmap))
}
################################################################################
#' plotVolcano
#'
#' Volcano plotting function
#' @param stats data frame with differential expression statistics
#' @param metadata #metadata from pdVolcano
#' @param FC Fold change threshold
#' @param pvalue p value threshold
#' @param title plot title
#' @export
plotVolcano <- function(
    stats,
    metadata,
    FC,
    pvalue,
    title
) {
  #bind variables locally 
  mean_diff <- welch_padj <- Color <- NULL
  #set custom colors
  myColors <- c("gray", "red", "dodgerblue")
  names(myColors) <- levels(stats$Color)
  custom_colors <- scale_color_manual(values = myColors, drop = FALSE)

  #plot
  volcano <- ggplot(data = stats,
                    aes(x = mean_diff, y = -log10(welch_padj),
                        color = Color,
                        label = stats$label)) +
    geom_vline(xintercept = c(FC, -FC), lty = "dashed") +
    geom_hline(yintercept = -log10(pvalue), lty = "dashed") +
    geom_point(na.rm = TRUE) +
    custom_colors +
    coord_cartesian(ylim = c(0L, 350L),
                    xlim = c(min(stats$mean_diff), max(stats$mean_diff))) +
    ggrepel::geom_text_repel(size = 3L, point.padding = 1L, color = "black",
                             min.segment.length = 0.1, box.padding = 0.15,
                             max.overlaps = Inf, na.rm = TRUE) +
    labs(x = "FC",
         y = "Significance (-log10pval)",
         color = NULL) +
    ggtitle(paste(title)) +
    theme_bw() +
    theme(plot.title = element_text(size = 16L),
          legend.position = "bottom",
          axis.title = element_text(size = 14L),
          legend.text = element_text(size = 12L))
  return(volcano)
}

################################################################################
#' pdVolcano
#'
#' Generate volcano plot and gate genes based on fold change and pvalue, 
#' includes vectors that can be used with fast gene set enrichment (fgsea)
#' @param result result output from projectionDriveR function in PV mode 
#' @param FC fold change threshold, default at 0.2
#' @param pvalue significance threshold, default set stored pvalue
#' @param subset vector of gene names to subset the plot by
#' @param filter.inf remove genes that have pvalues below machine double minimum value
#' @param label.num Number of genes to label on either side of the volcano plot, default 5
#' @param display boolean. Whether or not to plot and display volcano plots
#' @importFrom stats var
#' @importFrom ggrepel geom_label_repel
#' @importFrom ggrepel geom_text_repel
#' @import msigdbr
#' @import fgsea
#' @import dplyr
#' @return A list with weighted and unweighted differential expression metrics
#' @export
#plot FC, weighted and unweighted. Designed for use with the output of projectionDriveRs
pdVolcano <- function(
    result,
    FC = 0.2,
    pvalue = NULL,
    subset = NULL,
    filter.inf = FALSE,
    label.num = 5L,
    display = TRUE) {
  
  #bind local variables 
  welch_padj <- Color <- NULL
  #if a genelist is provided, use them to subset the output of projectiondrivers
  if (!is.null(subset)) {
    #subset the mean_stats object by provided gene list
    result$mean_stats <- result$mean_stats[(which(rownames(result$mean_stats) %in% subset)), ]
    #subset the weighted_mean_stats object by provided gene list
    result$weighted_mean_stats <- result$weighted_mean_stats[(which(rownames(result$weighted_mean_stats) %in% subset)), ]

  }

  if (filter.inf) {
    #remove p values below the machine limit representation for plotting purposes 
    cat("Filtering", length(which(result$mean_stats$welch_padj <= .Machine$double.xmin)), "unweighted genes and",
        length(which(result$weighted_mean_stats$welch_padj <= .Machine$double.xmin)), "weighted genes", "\n")
    result$mean_stats <- subset(result$mean_stats, welch_padj > .Machine$double.xmin)
    result$weighted_mean_stats <- subset(result$weighted_mean_stats, welch_padj > .Machine$double.xmin)
  }

  if (!is.numeric(FC)) {
    stop('FC must be a number')
  }

  if (!is.null(pvalue)) {

    if (!is.numeric(pvalue)) {
      stop('p value must be a number')
    }

    message("Updating sig_genes...")
    #update previously stored pvalue
    pvalue <- pvalue
    result$meta_data$pvalue <- pvalue
    #update sig_genes with new pvalue
    #recreate vector of significant genes from weighted and unweighted tests
    weighted_PV_sig <- rownames(result$weighted_mean_stats[which(result$weighted_mean_stats$welch_padj <= pvalue), ])
    PV_sig <- rownames(result$mean_stats[which(result$mean_stats$welch_padj <= pvalue), ])
    #create vector of significant genes shared between weighted and unweighted tests
    shared_genes_PV <- base::intersect(
      PV_sig, weighted_PV_sig)
    result$sig_genes <- list(PV_sig = PV_sig,
                             weighted_PV_sig = weighted_PV_sig,
                             PV_significant_shared_genes = shared_genes_PV)
  } else {
      pvalue <- result$meta_data$pvalue
  }

  #extract object meta data 
  metadata <- result$meta_data

  #volcano plotting unweighted
  #extract unweighted confidence intervals / statistics
  mean_stats <- result$mean_stats
  #fold change / significance calls
  mean_stats$Color <- paste("NS or FC <", FC)
  mean_stats$Color[mean_stats$welch_padj < pvalue & mean_stats$mean_diff > FC] <- paste("Enriched in", metadata$test_matrix)
  mean_stats$Color[mean_stats$welch_padj < pvalue & mean_stats$mean_diff < -FC] <- paste("Enriched in", metadata$reference_matrix)
  mean_stats$Color <- factor(mean_stats$Color,
                        levels = c(paste("NS or FC <", FC),
                                   paste("Enriched in", metadata$reference_matrix),
                                   paste("Enriched in", metadata$test_matrix)))

  #label the most significant genes for enrichment
  mean_stats$invert_P <- (-log10(mean_stats$welch_padj)) * (mean_stats$mean_diff)

  top_indices <- order(mean_stats$invert_P, decreasing = TRUE)[1L:label.num]
  bottom_indices <- order(mean_stats$invert_P)[1L:label.num]

  # Add labels to the dataframe
  mean_stats$label <- NA
  mean_stats$label[top_indices] <- paste(rownames(mean_stats)[top_indices])
  mean_stats$label[bottom_indices] <- paste(rownames(mean_stats)[bottom_indices])
  #unweighted volcano plot
  unweightedvolcano <- plotVolcano(stats = mean_stats,
                                   metadata = metadata,
                                   FC = FC,
                                   pvalue = pvalue,
                                   title = "Differential Enrichment")
  #weighted volcano plot
  weighted_mean_stats <- result$weighted_mean_stats
  weighted_mean_stats$Color <- paste("NS or FC <", FC)
  weighted_mean_stats$Color[weighted_mean_stats$welch_padj < pvalue & weighted_mean_stats$mean_diff > FC] <- paste("Enriched in", metadata$test_matrix)
  weighted_mean_stats$Color[weighted_mean_stats$welch_padj < pvalue & weighted_mean_stats$mean_diff < -FC] <- paste("Enriched in", metadata$reference_matrix)
  weighted_mean_stats$Color <- factor(weighted_mean_stats$Color,
                                      levels = c(paste("NS or FC <", FC),
                                                 paste("Enriched in", metadata$reference_matrix),
                                                 paste("Enriched in", metadata$test_matrix)))

  weighted_mean_stats$invert_P <- (-log10(weighted_mean_stats$welch_padj)) * (weighted_mean_stats$mean_diff)

  top_indices_w <- order(weighted_mean_stats$invert_P, decreasing = TRUE)[1L:label.num]
  bottom_indices_w <- order(weighted_mean_stats$invert_P)[1L:label.num]

  # Add labels to the dataframe
  weighted_mean_stats$label <- NA
  weighted_mean_stats$label[top_indices_w] <- paste(rownames(weighted_mean_stats)[top_indices_w])
  weighted_mean_stats$label[bottom_indices_w] <- paste(rownames(weighted_mean_stats)[bottom_indices_w])

  #weighted volcano plot
  weightedvolcano <- plotVolcano(stats = weighted_mean_stats,
                                 FC = FC,
                                 pvalue = pvalue,
                                 title = "Weighted Differential Enrichment")

  #return a list of genes that can be used as input to fgsea
  difexdf <- subset(mean_stats,
                    Color == paste("Enriched in", metadata$reference_matrix) | Color == paste("Enriched in", metadata$test_matrix))
  vec <- difexdf$estimate
  names(vec) <- rownames(difexdf)

  weighted_difexdf <- subset(weighted_mean_stats,
                             Color == paste("Enriched in", metadata$reference_matrix) | Color == paste("Enriched in", metadata$test_matrix))
  weighted_vec <- weighted_difexdf$estimate
  names(weighted_vec) <- rownames(weighted_difexdf)
  names(vec) <- rownames(difexdf)
  vol_result <- list(mean_stats = mean_stats,
                     weighted_mean_stats = weighted_mean_stats,
                     normalized_weights = result$normalized_weights,
                     sig_genes = result$sig_genes,
                     difexpgenes = difexdf,
                     weighted_difexpgenes = weighted_difexdf,
                     fgseavecs = list(unweightedvec = vec,
                                      weightedvec = weighted_vec),
                     meta_data = metadata,
                     plt = list(differential_expression = unweightedvolcano,
                                weighted_differential_expression = weightedvolcano))
  if (display) {
    #print volcano plots
    pltgrid <- cowplot::plot_grid(vol_result$plt$differential_expression +
                                    theme(legend.position = "none"),
                                  vol_result$plt$weighted_differential_expression +
                                    theme(legend.position = "none"),
                                  ncol = 2L, align = "h")
    legend <- cowplot::get_legend(vol_result$plt$differential_expression +
                                    guides(color = guide_legend(nrow = 1L)) +
                                    theme(legend.position = "bottom"))
    plt <- cowplot::plot_grid(pltgrid,
                              legend,
                              ncol = 1L,
                              rel_heights = c(1.0, 0.1))
    print(plt)
  }
  return(vol_result)
}
