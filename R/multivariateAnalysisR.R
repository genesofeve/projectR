
#' Generic multivariateAnalysisR function
#' @rdname multivariateAnalysisR
#'  
#' @description Performs multivariate analysis across specified clusters in datasets
#' @import ggplot2 RColorBrewer
#' 
#' @param significanceLevel double value for testing significance in ANOVA test
#' @param patternKeys list of strings indicating pattern subsets from seuratobj to be analyzed
#' @param seuratobj Seurat Object Data containing patternKeys in meta.data
#' @param dictionaries list of dictionaries indicating clusters to be compared
#' @param customNames list of custom names for clusters in corresponding order
#' @param exclusive boolean value for determining interpolation between params in clusters
#' @param exportFolder name of folder to store exported graphs and CSV files
#' @param ANOVAwidth width of ANOVA png
#' @param ANOVAheight height of ANOVA png
#' @param CIwidth width of CI png
#' @param CIheight height of CI png
#' @param CIspacing spacing between each CI in CI graph
#'
#' @return none; ANOVA and Confidence Intervals are visualized and exported in both PNG and CSV
#' @export
#'
#' @examples
#' dict1 <- list("myclusters" = "RPCs-1", "Condition" = c("72h_Mtz", "24h_Mtz"))
#' dict2 <- list("myclusters" = "Muller Glia", "Condition" = c("72h_Mtz", "24h_Mtz"))
#' multivariateAnalysisR(significanceLevel = 0.05, patternKeys = list("Muller_Pattern1", "Muller_Pattern2", "Muller_Pattern3"), seuratobj = RNA4time_seurat, dictionaries = list(dict1, dict2), customNames = list("test1", "test2"), exclusive = TRUE, exportFolder = "results", ANOVAwidth = 1000, ANOVAheight = 1000, CIwidth = 1000, CIheight = 1000, CIspacing = 2)
#' 

multivariateAnalysisR <- function(
    significanceLevel = 0.05, # double value for testing significance in ANOVA test
    patternKeys, # list of strings indicating pattern subsets from seuratobj to be analyzed
    seuratobj, # Seurat Object Data containing patternKeys in meta.data
    dictionaries, # list of dictionaries indicating clusters to be compared
    customNames = NULL, # list of custom names for clusters in corresponding order
    exclusive = TRUE, # boolean value for determining interpolation between params in clusters
    exportFolder = "", # name of folder to store exported graphs
    ANOVAwidth = 1000, # width of ANOVA png
    ANOVAheight = 1000, # height of CI png
    CIwidth = 1000, # width of CI png
    CIheight = 1000, # height of CI png
    CIspacing = 1 # spacing between each CI in CI graph
    ){
  
  getSpecificData <- function(key, seuratobj) {
    expr <- parse(text = key)
    result <- eval(expr, envir = seuratobj@meta.data)
    return (result)
  }
  
  getIndividualIndices <- function(key, values, seuratobj) {
    result <- c()
    subset <- getSpecificData(key = key, seuratobj = seuratobj)
    for (value in values) {
      indices <- which(subset == value)
      result <- unique(c(result, indices))
    }
    return (result)
  }
  
  getFinalIndices <- function(exclusive, seuratobj, dictionary) {
    allLists <- list()
    allListsIndices <- 1
    for (key in names(dictionary)) {
      values <- dictionary[[key]]
      individualIndices <- getIndividualIndices(key = key, values = values, seuratobj = seuratobj)
      allLists[[allListsIndices]] <- individualIndices
      allListsIndices <- allListsIndices + 1
    }
    if (exclusive == TRUE) {
      result <- Reduce(intersect, allLists)
    } else {
      result <- unique(unlist(allLists))
    }
    return (result)
  }
  
  getPatternFromIndices <- function(indices, patternKey, seuratobj) {
    result <- getSpecificData(key = patternKey, seuratobj = seuratobj)
    return (result[indices])
  }
  
  ANOVA_find_mean <- function(patternList) {
    combined_vector <- unlist(patternList)
    result <- mean(combined_vector)
    return (result)
  }
  
  ANOVA_find_ssb <- function(patternList, totalMean) {
    numGroups <- length(patternList)
    result <- 0
    for (i in 1 : numGroups) {
      groupMean <- mean(patternList[[i]])
      result <- result + length(patternList[[i]]) * (groupMean - totalMean)^2
    }
    return (result)
  }
  
  ANOVA_find_ssw <- function(patternList) {
    groupMeans <- sapply(patternList, mean)
    result <- 0
    for (i in seq_along(patternList)) {
      groupValues <- patternList[[i]]
      groupMean <- groupMeans[i]
      squaredDiff <- sum((groupValues - groupMean)^2)
      result <- result + squaredDiff
    }
    return(result)
  }
  
  # ANOVA is one-way for this case, as we are representing different dictionaries as one big unidirectional change rather than a change across two independent factors
  ANOVA <- function(significanceLevel = 0.05, patternKey, exclusive, seuratobj, dictionaries) {
    patternList <- list()
    patternListIndices <- 1
    for (dictionary in dictionaries) {
      dictIndices <- getFinalIndices(exclusive = exclusive, seuratobj = seuratobj, dictionary = dictionary)
      dictPatterns <- getPatternFromIndices(indices = dictIndices, patternKey = patternKey, seuratobj = seuratobj)
      patternList[[patternListIndices]] <- dictPatterns
      patternListIndices <- patternListIndices + 1
    }
    # ANOVA analysis begins
    # step 1: find total mean
    totalMean <- ANOVA_find_mean(patternList)
    # step 2: find sum of squares between groups (ssb)
    ssb <- ANOVA_find_ssb(patternList, totalMean)
    # step 3: find sum of squares within groups (ssw)
    ssw <- ANOVA_find_ssw(patternList)
    # step 4: find degrees of freedom (df) for between groups (dfB) and within groups (dfW)
    k <- length(patternList)
    N <- sum(sapply(patternList, length))
    dfB <- k - 1
    dfW <- N - k
    # step 5: find mean square between groups (msb) and the mean square within groups (msw)
    msb <- ssb / dfB
    msw <- ssw / dfW
    # step 6: find F-statistic, p-value, and boolean value for significance
    FStats <- msb / msw
    pValue <- 1 - pf(FStats, dfB, dfW)
    isSignificant <- if (pValue < significanceLevel) TRUE else FALSE 
    # return result
    return (c(FStats, pValue, isSignificant))
  }
  
  CI <- function(patternKey, exclusive, seuratobj, dictionary1, dictionary2) {
    # for dictionary1
    pattern1DictIndices <- getFinalIndices(exclusive = exclusive, seuratobj = seuratobj, dictionary = dictionary1)
    pattern1DictPatterns <- getPatternFromIndices(indices = pattern1DictIndices, patternKey = patternKey, seuratobj = seuratobj)
    # for dictionary2
    pattern2DictIndices <- getFinalIndices(exclusive = exclusive, seuratobj = seuratobj, dictionary = dictionary2)
    pattern2DictPatterns <- getPatternFromIndices(indices = pattern2DictIndices, patternKey = patternKey, seuratobj = seuratobj) 
    # CI analysis begins
    # Welch’s t-test is used as it is designed to be more reliable when the two samples have unequal variances and/or unequal sample sizes
    # step 1: perform Welch's t-test
    welchResults <- t.test(pattern1DictPatterns, pattern2DictPatterns, var.equal=FALSE) # var.equal=TRUE if variances are approximately equal
    # step 2: get confidence interval (95%)
    CILowerBound <- welchResults$conf.int[1]
    CIUpperBound <- welchResults$conf.int[2]
    # flip CI to keep it positively directional
    if (mean(c(CILowerBound, CIUpperBound)) < 0) {
      CILowerBound <- -CILowerBound
      CIUpperBound <- -CIUpperBound
    }
    if (CILowerBound > CIUpperBound) {
      temp <- CIUpperBound
      CIUpperBound <- CILowerBound
      CILowerBound <- temp
    }
    # return result
    return (c(CILowerBound, CIUpperBound))
  }
  
  multivariateNaming <- function(index, customNames) {
    if ((index < 1) || (index > length(customNames))) {
      return (paste("cluster", index, sep = ""))
    } else {
      return (customNames[index])
    }
  }
  
  # check basic conditions
  errorMsg <- ""
  if (significanceLevel < 0) {
    errorMsg <- paste(errorMsg, "significanceLevel must be non-negative.\n", sep="")
  }
  if (length(patternKeys) < 1) {
    errorMsg <- paste(errorMsg, "patternKeys must have at least one pattern.\n", sep="")
  }
  if (length(dictionaries) < 2) {
    errorMsg <- paste(errorMsg, "dictionaries must have at least two clusters to be compared.\n", sep="")
  }
  if (!is.logical(exclusive)) {
    errorMsg <- paste(errorMsg, "exclusive must be a boolean value.\n", sep="")
  }
  
  if (errorMsg != "") {
    stop(cat(errorMsg), "Invalid parameter(s). Please fix and try again.")
  }
  
  # begin
  grossList <- list()
  for (patternKey in patternKeys) {
    # ANOVA
    patternANOVA <- ANOVA(
      significanceLevel = significanceLevel, 
      patternKey = patternKey, 
      exclusive = exclusive, 
      seuratobj = seuratobj, 
      dictionaries = dictionaries
    )
    # CI
    patternCI <- list()
    # looping through all pairs
    dictionariesLength <- length(dictionaries)
    for (firstIndex in 1 : (dictionariesLength - 1)) {
      for (secondIndex in (firstIndex + 1) : dictionariesLength) {
        pair1Name <- multivariateNaming(firstIndex, customNames)
        pair2Name <- multivariateNaming(secondIndex, customNames)
        pairCI <- CI(
          patternKey = patternKey, 
          exclusive = exclusive,
          seuratobj = seuratobj,
          dictionary1 = dictionaries[[firstIndex]],
          dictionary2 = dictionaries[[secondIndex]]
        )
        pairCIresults <- list(
          "pair1Name" = pair1Name,
          "pair2Name" = pair2Name,
          "CI" = pairCI
        )
        patternCI <- c(patternCI, list(pairCIresults))
      }
    }
    patternResults <- list("patternKey" = patternKey, 
                           "ANOVA" = patternANOVA, 
                           "CI" = patternCI)
    grossList <- c(grossList, list(patternResults))
  }
  # reorder in increasing order (why not decreasing: ggplot reverses)
  sorted_indices <- order(sapply(grossList, function(x) x[["ANOVA"]][1]), decreasing = FALSE)
  sortedGrossList <- grossList[sorted_indices]
  
  # ggplotting ANOVA
  # step 1: initialization
  name <- sapply(sortedGrossList, function(x) x$patternKey)
  value <- sapply(sortedGrossList, function(x) x[["ANOVA"]][1])
  pValue <- sapply(sortedGrossList, function(x) x[["ANOVA"]][2]) # for CSV export
  significance <- sapply(sortedGrossList, function(x) x[["ANOVA"]][3])
  color <- ifelse(significance == 1, "blue", "red")
  # Create a data frame
  df <- data.frame(
    name = factor(name, levels = unique(name)),
    value = value,
    color = color
  )
  # Create barplot
  ggplot(df, aes(x = name, y = value)) +
    geom_bar(stat = 'identity', aes(fill = color), position = 'dodge') +
    coord_flip() +
    labs(title = "Multivariate Analysis - ANOVA",
         x = "Patterns",
         y = "F-statistic") +
    scale_fill_manual(name = "Significance", 
                      values = c("red" = "red", "blue" = "blue"),
                      labels = c(
                        "red" = paste("p > ", significanceLevel, sep = ""), 
                        "blue" = paste("p < ", significanceLevel, sep = ""))
    ) + 
    theme_minimal() + 
    theme(
      axis.line = element_line(linewidth = 0.5, linetype = "solid"),
      plot.title = element_text(size = 20, face = "bold"),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    )
  # save ANOVA plot
  ANOVA_file_path <- file.path(getwd(), exportFolder, "multivariateAnalysisR_ANOVA.png")
  ggsave(filename = ANOVA_file_path, plot = last_plot(), width = ANOVAwidth, height = ANOVAheight, units = "px", limitsize = FALSE, device = "png", bg = "white")
  
  # ggplotting CI
  group <- unlist(lapply(sortedGrossList, function(x) rep(x$patternKey, times = length(x$CI))))
  group <- factor(group, levels = unique(group)) # convert 'group' to a factor with original order
  ci_lower <- unlist(
    lapply(sortedGrossList, function(x) {
      unlist(
        lapply(x$CI, function(y) {
          rep(y$CI[1], times = length(y$CI[1]))
        })
      )
    })
  )
  ci_upper <- unlist(
    lapply(sortedGrossList, function(x) {
      unlist(
        lapply(x$CI, function(y) {
          rep(y$CI[2], times = length(y$CI[2]))
        })
      )
    })
  )
  lengthSortedGrossListCI <- length(sortedGrossList[[1]]$CI)
  y_pos <- seq(length(group) + length(sortedGrossList) - 1) * CIspacing # y-axis positions for each group and dividing lines
  dividing_lines_indices <- seq(length(sortedGrossList) - 1) * (lengthSortedGrossListCI + 1)
  dividing_lines <- y_pos[dividing_lines_indices]
  y_pos <- y_pos[-dividing_lines_indices]
  if (lengthSortedGrossListCI < 3) {
    colorPalette <- switch(lengthSortedGrossListCI, "red", c("red", "blue"))
  } else {
    colorPalette <- grDevices::rainbow(lengthSortedGrossListCI) # different color-generating tools may be used for better distinctiveness
  }
  colors <- rep(colorPalette, length(sortedGrossList))
  labels <- group # identify which group (patternKey) the CI refers to
  color_labels <- unlist(
    lapply(sortedGrossList, function(x) {
      unlist(
        lapply(x$CI, function(y) {
          paste(y$pair1Name, y$pair2Name, sep = "-")
        })
      )
    })
  )
  # Create a data frame
  data <- data.frame(
    group, 
    ci_lower, 
    ci_upper, 
    y_pos,
    colors, 
    labels
  )
  # Create the horizontal plot
  ggplot_obj <- ggplot(data, aes(y = y_pos, x = (ci_lower + ci_upper) / 2)) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = colors), height = 0.2, size = 1.5) +
    geom_text(aes(label = labels), vjust = -1.6, size = 4, fontface = "bold") +
    scale_color_manual(name = "Color Labels", 
                       values = colors,
                       labels = color_labels) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") + # Add a reference line at x = 0
    labs(title = "Multivariate Analysis - Confidence Intervals",
         x = "Confidence Intervals",
         y = "Patterns") +
    theme_minimal() +
    theme(
      axis.line = element_line(linewidth = 0.5, linetype = "solid"),
      plot.title = element_text(size = 20, face = "bold"),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    )
  # dividing lines per pattern
  ggplot_obj <- ggplot_obj +
    geom_hline(yintercept = dividing_lines, linetype = "dashed", color = "black", size = 0.5)
  # save CI plot
  CI_file_path <- file.path(getwd(), exportFolder, "multivariateAnalysisR_CI.png")
  ggsave(filename = CI_file_path, plot = last_plot(), width = CIwidth, height = CIheight, units = "px", limitsize = FALSE, device = "png", bg = "white")
  
  # CSV export
  # ANOVA
  ANOVA_significance <- ifelse(significance == 1, "TRUE", "FALSE")
  ANOVA_CSV_df <- data.frame(
    "Pattern" = rev(name),
    "F-statistic value" = rev(value),
    "p-value" = rev(pValue),
    "Statistically Significant" = rev(ANOVA_significance)
  )
  ANOVA_export_path <- file.path(getwd(), exportFolder, "multivariateAnalysisR_ANOVA.csv")
  write.csv(ANOVA_CSV_df, ANOVA_export_path, row.names = FALSE)
  
  # CI
  color_labels_pair1name <- unlist(
    lapply(sortedGrossList, function(x) {
      unlist(
        lapply(x$CI, function(y) {
          y$pair1Name
        })
      )
    })
  )
  color_labels_pair2name <- unlist(
    lapply(sortedGrossList, function(x) {
      unlist(
        lapply(x$CI, function(y) {
          y$pair2Name
        })
      )
    })
  )
  CI_CSV_df <- data.frame(
    "Pattern" = rev(group),
    "CI Lower Value" = rev(ci_lower),
    "CI Higher Value" = rev(ci_upper),
    "CI Cluster 1" = color_labels_pair1name,
    "CI Cluster 2" = color_labels_pair2name
  )
  CI_export_path <- file.path(getwd(), exportFolder, "multivariateAnalysisR_CI.csv")
  write.csv(CI_CSV_df, CI_export_path, row.names = FALSE)

}