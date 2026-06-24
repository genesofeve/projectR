#' @importFrom stats hclust kmeans prcomp
setOldClass("kmeans")
setOldClass("hclust")
setOldClass("prcomp")

###############################################################################
#' @import limma
#' @importFrom stats model.matrix pt
#' @param NP vector of integers indicating which columns of loadings object to
#'   use. The default of NP=NA will use entire matrix.
#' @param full logical indicating whether to return the full model solution. By
#'   default only the new pattern object is returned.
#' @param include_intercept Logical, default \code{FALSE}. If \code{FALSE}
#'   (default), fits \code{~ 0 + A} — the model must explain absolute
#'   expression levels including the global mean, which is absorbed into the
#'   pattern weights. If \code{TRUE}, fits \code{~ 1 + A} — an explicit
#'   intercept absorbs the per-sample global mean shift, leaving pattern
#'   coefficients to describe only pattern-specific variation above that
#'   baseline. When \code{full = TRUE}, the per-sample intercept estimates are
#'   returned in \code{$intercept}. Only applies to the matrix, dgCMatrix, and
#'   LinearEmbeddingMatrix (CoGAPS/NMF) dispatches.
#' @param model Optional arguments to choose method for projection
#' @param bootstrapPval logical to indicate whether to generate p-values using
#'   bootstrap, not available for prcomp and rotatoR objects
#' @param bootIter number of bootstrap iterations, default = 1000
#'
#' @return \subsection{Matrix, NMF, CoGAPS, correlateR, and clustering
#'   dispatch}{
#'   When \code{full = FALSE} (default), a numeric matrix (patterns x samples)
#'   of projected sample weights in the pattern space defined by
#'   \code{loadings}.
#'
#'   When \code{full = TRUE}, a named list with three or four elements (plus
#'   \code{bootstrapPval} if \code{bootstrapPval = TRUE}):
#'   \describe{
#'     \item{\code{projection}}{Numeric matrix (patterns x samples). Projected
#'       weights for each sample in the pattern space defined by
#'       \code{loadings}. Equivalent to the matrix returned when
#'       \code{full = FALSE}.}
#'     \item{\code{pval}}{Numeric matrix (patterns x samples). Two-sided
#'       p-values from t-tests with \code{lmFit} residual degrees of freedom
#'       (\eqn{n_{\text{samples}} - n_{\text{patterns}}}). Uses the
#'       t-distribution rather than the normal approximation, which is exact
#'       for all sample sizes and converges to the normal for large degrees of
#'       freedom.}
#'     \item{\code{r_squared}}{Named numeric vector (length = n samples).
#'       Per-sample coefficient of determination computed as
#'       \eqn{1 - ||r||^2 / ||X||^2} where \eqn{||r||^2} is the sum of
#'       squared residuals and \eqn{||X||^2} is the raw (uncentered) total sum
#'       of squares for the gene-matched data. Using raw \eqn{||X||^2} ensures
#'       comparability with the prcomp method and across basis methods. When
#'       \code{include_intercept = TRUE}, \eqn{R^2} will be higher because the
#'       intercept absorbs mean-shift variance; the difference between the two
#'       values quantifies that contribution.}
#'     \item{\code{intercept}}{(only when \code{include_intercept = TRUE} and
#'       \code{full = TRUE}) Named numeric vector (length = n samples). The
#'       per-sample intercept estimate from \code{lmFit}, representing the
#'       global mean expression shift not explained by the pattern basis.}
#'     \item{\code{bootstrapPval}}{(only when \code{bootstrapPval = TRUE})
#'       Bootstrap p-values for projection weights.}
#'   }
#' }
#'
#' @rdname projectR-methods
#' @aliases projectR
setMethod(
  "projectR", signature(data = "matrix", loadings = "matrix"),
  function(data, # new dataset. genes x samples
           loadings, # basis matrix defining the target feature space (genes x patterns)
           dataNames = NULL, # a vector with names of data rows
           loadingsNames = NULL, # a vector with names of loadings rows
           NP = NA, # vector of integers indicating which columns of loadings to use
           full = FALSE, # logical indicating whether to return the full model solution
           include_intercept = FALSE, # logical to include an intercept term in the design
           bootstrapPval = FALSE, # logical to indicate whether to generate p-values using bootstrap
           bootIter = 1e3 # No of bootstrap iterations
  ) {
    if (!anyNA(NP)) loadings <- loadings[, NP, drop = FALSE]

    # match genes in data sets
    if (is.null(dataNames)) {
      dataNames <- rownames(data)
    }
    if (is.null(loadingsNames)) {
      loadingsNames <- rownames(loadings)
    }
    dataM <- geneMatchR(
      data1 = data, data2 = loadings, data1Names = dataNames,
      data2Names = loadingsNames, merge = FALSE
    )
    message(dim(dataM[[2]])[1], " row names matched between data and loadings")
    message("Updated dimension of data:", as.character(paste(dim(dataM[[2]]), collapse = " ")))

    # dataM[[1]]: genes x patterns (matched loadings / A matrix)
    # dataM[[2]]: genes x samples  (matched data)
    if (is.null(colnames(dataM[[1]]))) {
      colnames(dataM[[1]]) <- paste0("pattern", seq_len(ncol(dataM[[1]])))
    }
    Design <- if (include_intercept) {
      d <- model.matrix(~ 1 + dataM[[1]])
      colnames(d) <- c("intercept", colnames(dataM[[1]]))
      d
    } else {
      d <- model.matrix(~ 0 + dataM[[1]])
      colnames(d) <- colnames(dataM[[1]])
      d
    }
    projection <- limma::lmFit(as.matrix(t(dataM[[2]])), Design)

    # projection$coefficients is n_new_samples x n_coefs; pattern_cols
    # excludes the intercept column (present only when include_intercept = TRUE).
    pattern_cols <- colnames(dataM[[1]])
    projectionPatterns <- t(projection$coefficients[, pattern_cols, drop = FALSE])
    intercept_vec <- if (include_intercept) {
      projection$coefficients[, "intercept"]
    } else {
      NULL
    }

    # Two-sided p-values using the t-distribution.
    # Previously computed with pnorm() (normal approximation), which is
    # anti-conservative when df.residual (n_samples - n_patterns) is small.
    # pt() is exact for all sample sizes and converges to normal for large df.
    # note: this is a change from the original projectR implementation, which
    # used pnorm rather than the t-distribution
    projection.ts <- t(
      projection$coefficients[, pattern_cols, drop = FALSE] /
        projection$stdev.unscaled[, pattern_cols, drop = FALSE] /
        projection$sigma
    )
    pval.matrix <- 2 * pt(abs(projection.ts),
      df         = projection$df.residual,
      lower.tail = FALSE
    )

    if (bootstrapPval) {
      boots <- lapply(1:bootIter, function(x) {
        rows <- sample(nrow(Design), nrow(Design), replace = TRUE)
        proj_boot <- limma::lmFit(as.matrix(t(dataM[[2]][rows, ])), Design[rows, ])
        return(proj_boot$coefficients)
      })
      bootPval <- compareBoots(
        projection$coefficients[, pattern_cols, drop = FALSE],
        lapply(boots, \(b) b[, pattern_cols, drop = FALSE])
      )
    }

    if (full) {
      reconstruction <- dataM[[1]] %*% projectionPatterns # genes x samples
      if (include_intercept) {
        reconstruction <- reconstruction +
          matrix(intercept_vec,
            nrow = nrow(dataM[[1]]), ncol = ncol(dataM[[2]]), byrow = TRUE
          )
      }
      ss_res <- colSums((dataM[[2]] - reconstruction)^2)
      ss_tot <- colSums(dataM[[2]]^2)
      r_squared <- 1 - ss_res / ss_tot

      projectionFit <- list(
        projection = projectionPatterns,
        pval       = pval.matrix,
        r_squared  = r_squared
      )
      if (include_intercept) projectionFit$intercept <- intercept_vec
      if (bootstrapPval) projectionFit$bootstrapPval <- bootPval
      return(projectionFit)
    } else {
      return(projectionPatterns)
    }
  }
)

###############################################################################
#' @import MatrixModels
#' @importFrom stats model.matrix
#' @param NP vector of integers indicating which columns of loadings object to
#'   use. The default of NP=NA will use entire matrix.
#' @param full logical indicating whether to return the full model solution. By
#'   default only the new pattern object is returned.
#' @param model Optional arguments to choose method for projection
#' @param chopBy number of columns to chop the data into (chopping helps running
#'   large datasets)
#' @rdname projectR-methods
#' @aliases projectR
setMethod(
  "projectR", signature(data = "dgCMatrix", loadings = "matrix"),
  function(data, # new dataset. genes x samples
           loadings, # basis matrix defining the target feature space (genes x patterns)
           dataNames = NULL, # a vector with names of data rows
           loadingsNames = NULL, # a vector with names of loadings rows
           NP = NULL, # vector of integers indicating which columns of loadings to use
           full = FALSE, # logical indicating whether to return the full model solution
           include_intercept = FALSE, # logical to include an intercept term in the design
           chopBy = 1000 # number of columns per chunk
  ) {
    if (!is.null(NP)) {
      loadings <- loadings[, NP, drop = FALSE]
    }

    message("dgCMatrix detected, projecting in chunks.")

    chop <- function(sparsematrix) {
      coln <- ncol(sparsematrix)
      bins <- seq(1, coln, by = chopBy)
      lapply(seq_along(bins), function(i) {
        start <- bins[i]
        end <- ifelse(i < length(bins), bins[i + 1] - 1, coln)
        return(start:end)
      })
    }

    # Suppress per-chunk geneMatchR print statements; print once below
    w <- invisible(capture.output(
      projectionList <- lapply(chop(data), function(i) {
        projectR(as.matrix(data[, i]), loadings, full = full,
          include_intercept = include_intercept)
      })
    ))
    message(w[1])

    if (full) {
      if (length(projectionList) == 1) {
        res <- projectionList[[1]]
      } else {
        projections <- do.call(
          cbind,
          lapply(projectionList, function(x) x[["projection"]])
        )
        pvalues <- do.call(
          cbind,
          lapply(projectionList, function(x) x[["pval"]])
        )
        r_squared <- unlist(
          lapply(projectionList, function(x) x[["r_squared"]])
        )

        res <- list(
          projection = projections,
          pval       = pvalues,
          r_squared  = r_squared
        )
        if (include_intercept) {
          res$intercept <- unlist(
            lapply(projectionList, \(x) x[["intercept"]])
          )
        }
      }
    } else {
      res <- do.call(cbind, projectionList)
    }
    return(res)
  }
)


###############################################################################
#' @import limma
#' @import SingleCellExperiment
#' @importFrom NMF fcnnls
#' @examples
#' library("CoGAPS")
#' # CR.RNAseq6l3c3t <- CoGAPS(p.RNAseq6l3c3t, params = new("CogapsParams", nPatterns=5))
#' projectR(
#'   data = p.ESepiGen4c1l$mRNA.Seq, loadings = CR.RNAseq6l3c3t,
#'   dataNames = map.ESepiGen4c1l[["GeneSymbols"]]
#' )
#'
#' @rdname projectR-methods
#' @aliases projectR
setMethod(
  "projectR", signature(data = "matrix", loadings = "LinearEmbeddingMatrix"),
  function(data, # new dataset. genes x samples
           loadings, # a LinearEmbeddingMatrix (e.g. CoGAPS result)
           dataNames = NULL,
           loadingsNames = NULL,
           NP = NA,
           full = FALSE,
           model = NA,
           include_intercept = FALSE,
           bootstrapPval = FALSE,
           bootIter = 1e3) {
    loadings <- loadings@featureLoadings
    if (!anyNA(NP)) loadings <- loadings[, NP, drop = FALSE]
    return(projectR(data,
      loadings          = loadings,
      dataNames         = dataNames,
      loadingsNames     = loadingsNames,
      NP                = NP,
      full              = full,
      include_intercept = include_intercept,
      bootstrapPval     = bootstrapPval,
      bootIter          = bootIter
    ))
  }
)

###############################################################################
#' @import limma
#' @importFrom stats var
#' @param center_by_loadings Logical, default \code{FALSE}. If \code{TRUE},
#'   centers the new data using the gene means stored in \code{loadings$center}
#'   (computed from the original training data) before projecting onto the
#'   loading vectors. If \code{FALSE}, each gene is centered by its own mean in
#'   the new data prior to projection. Only valid when \code{loadings$center}
#'   is not \code{FALSE}, i.e. when \code{prcomp} was called with
#'   \code{center = TRUE}.
#' @param full Logical, default \code{FALSE}. If \code{FALSE}, returns only the
#'   scores matrix. If \code{TRUE}, returns a named list with three elements;
#'   see Value.
#'
#' @return \subsection{prcomp dispatch}{
#'   When \code{full = FALSE} (default), a numeric matrix (PCs x samples) of
#'   projected sample coordinates in the PC space defined by \code{loadings}.
#'   Row names are PC labels (e.g. \code{PC1}, \code{PC2}, ...); column names
#'   are sample names from \code{data}.
#'
#'   When \code{full = TRUE}, a named list with three elements. Note that the
#'   second element is \code{$pvar} (variance), not \code{$pval} (p-values) as
#'   returned by the matrix/NMF/CoGAPS dispatch:
#'   \describe{
#'     \item{\code{projection}}{Numeric matrix (PCs x samples). Projected
#'       coordinates of each sample in the PC space defined by \code{loadings}.
#'       Equivalent to the matrix returned when \code{full = FALSE}.}
#'     \item{\code{pvar}}{Named numeric vector (length = number of PCs).
#'       Percentage of variance in the \emph{projected} data accounted for by
#'       each PC, computed as \eqn{100 \times \sigma^2_k / \sum_j \sigma^2_j}
#'       where \eqn{\sigma^2_k} is the variance of the scores on PC \eqn{k}.
#'       Note that these percentages are relative to the projected data and
#'       will differ from the variance explained reported by the original
#'       \code{prcomp} fit.}
#'     \item{\code{r_squared}}{Named numeric vector (length = n samples).
#'       Per-sample coefficient of determination computed as
#'       \eqn{1 - ||r||^2 / ||X||^2} where \eqn{||r||^2} is the sum of
#'       squared residuals (centered data minus reconstruction) and
#'       \eqn{||X||^2} is the raw (uncentered) total sum of squares for the
#'       gene-matched data. Using raw \eqn{||X||^2} ensures comparability with
#'       the matrix dispatch (NMF/CoGAPS) method across projection methods.
#'       Note that because PCA explicitly removes the mean before fitting, its
#'       R^2 will be lower than NMF/CoGAPS R^2 when a global mean shift exists
#'       between reference and new data — this difference is itself
#'       informative. Comparing \code{center_by_loadings = TRUE} vs
#'       \code{FALSE} R^2 values quantifies the contribution of the mean shift
#'       to total variance.}
#'   }
#' }
#'
#' @examples
#' pca.RNAseq6l3c3t <- prcomp(t(p.RNAseq6l3c3t))
#' pca.ESepiGen4c1l <- projectR(
#'   data = p.ESepiGen4c1l$mRNA.Seq,
#'   loadings = pca.RNAseq6l3c3t, dataNames = map.ESepiGen4c1l[["GeneSymbols"]]
#' )
#'
#' @rdname projectR-methods
#' @aliases projectR
setMethod(
  "projectR", signature(data = "matrix", loadings = "prcomp"),
  function(data,
           loadings,
           dataNames = NULL,
           loadingsNames = NULL,
           NP = NA,
           center_by_loadings = FALSE,
           full = FALSE) {
    if (center_by_loadings) {
      loadings_center <- loadings$center

      if (is.null(loadings_center) || identical(loadings_center, FALSE)) {
        stop(
          "center_by_loadings = TRUE but loadings$center is missing or FALSE. ",
          "Rerun prcomp() with center = TRUE, or set center_by_loadings = FALSE (default) ",
          "to center each new sample by its own gene means."
        )
      }

      if (is.null(names(loadings_center))) {
        stop(
          "loadings$center is unnamed. Cannot gene-match the centering vector. ",
          "Ensure prcomp() was called with a matrix whose columns are named by gene."
        )
      }
    }

    loadings <- loadings$rotation
    if (!anyNA(NP)) loadings <- loadings[, NP, drop = FALSE]

    if (is.null(dataNames)) {
      dataNames <- rownames(data)
    }
    if (is.null(loadingsNames)) {
      loadingsNames <- rownames(loadings)
    }
    dataM <- geneMatchR(
      data1 = data, data2 = loadings, data1Names = dataNames,
      data2Names = loadingsNames, merge = FALSE
    )
    message(dim(dataM[[2]])[1], " row names matched between data and loadings")
    message("Updated dimension of data: ", nrow(dataM[[2]]), " x ", ncol(dataM[[2]]))

    # dat2P is samples x genes
    dat2P <- if (center_by_loadings) {
      loadings_center <- loadings_center[rownames(dataM[[2]])]
      if (any(is.na(loadings_center))) {
        missing <- rownames(dataM[[2]])[is.na(loadings_center)]
        stop(
          length(missing), " gene(s) in the matched data are absent from ",
          "loadings$center. This should not happen if loadings is a valid ",
          "prcomp object"
        )
      }
      t(sweep(dataM[[2]], 1, loadings_center, "-"))
    } else {
      apply(dataM[[2]], 1, function(x) x - mean(x))
    }

    # projectionPatterns is samples x PCs
    projectionPatterns <- dat2P %*% dataM[[1]]

    if (full) {
      # Percent variance accounted for by each PC in the projected data
      PercentVariance <- apply(
        projectionPatterns, 2,
        function(x) 100 * var(x) / sum(apply(projectionPatterns, 2, var))
      )

      # R^2 per sample using raw (uncentered) ||X||^2 as denominator.
      # Raw SS_tot is used for comparability with the matrix dispatch
      # (NMF/CoGAPS). PCA R^2 will be lower than NMF/CoGAPS R^2 when a
      # global mean shift exists, since PCA does not model the mean.
      reconstruction <- projectionPatterns %*% t(dataM[[1]]) # samples x genes
      ss_res <- rowSums((dat2P - reconstruction)^2)
      ss_tot <- colSums(dataM[[2]]^2) # raw ||X||^2, genes x samples → per sample
      r_squared <- 1 - ss_res / ss_tot

      projectionFit <- list(
        projection = t(projectionPatterns), # PCs x samples
        pvar       = PercentVariance,
        r_squared  = r_squared
      )
      return(projectionFit)
    } else {
      return(t(projectionPatterns))
    }
  }
)

##############################################################################
#' @examples
#' pca.RNAseq6l3c3t <- prcomp(t(p.RNAseq6l3c3t))
#' r.RNAseq6l3c3t <- rotatoR(1, 1, -1, -1, pca.RNAseq6l3c3t$rotation[, 1:2])
#' pca.ESepiGen4c1l <- projectR(
#'   data = p.ESepiGen4c1l$mRNA.Seq,
#'   loadings = r.RNAseq6l3c3t, dataNames = map.ESepiGen4c1l[["GeneSymbols"]]
#' )
#'
#' @rdname projectR-methods
#' @aliases projectR
setMethod(
  "projectR", signature(data = "matrix", loadings = "rotatoR"),
  function(data,
           loadings,
           dataNames = NULL,
           loadingsNames = NULL,
           NP = NA,
           full = FALSE) {
    loadings <- loadings@rotatedM
    if (!anyNA(NP)) loadings <- loadings[, NP, drop = FALSE]

    if (is.null(dataNames)) {
      dataNames <- rownames(data)
    }
    if (is.null(loadingsNames)) {
      loadingsNames <- rownames(loadings)
    }
    dataM <- geneMatchR(
      data1 = data, data2 = loadings, data1Names = dataNames,
      data2Names = loadingsNames, merge = FALSE
    )
    message(dim(dataM[[2]])[1], " row names matched between data and loadings")
    message("Updated dimension of data: ", nrow(dataM[[2]]), " x ", ncol(dataM[[2]]))

    dat2P <- apply(dataM[[2]], 1, function(x) x - mean(x))
    projectionPatterns <- dat2P %*% dataM[[1]]

    if (full) {
      PercentVariance <- round(
        eigen(cov(projectionPatterns))$values /
          sum(eigen(cov(projectionPatterns))$values) * 100,
        digits = 2
      )

      # R^2 per sample using raw ||X||^2
      reconstruction <- projectionPatterns %*% t(dataM[[1]]) # samples x genes
      ss_res <- rowSums((dat2P - reconstruction)^2)
      ss_tot <- colSums(dataM[[2]]^2)
      r_squared <- 1 - ss_res / ss_tot

      projectionFit <- list(
        projection = t(projectionPatterns),
        pvar       = PercentVariance,
        r_squared  = r_squared
      )
      return(projectionFit)
    } else {
      return(t(projectionPatterns))
    }
  }
)

##############################################################################
#' @import limma
#' @examples
#' c.RNAseq6l3c3t <- correlateR(
#'   genes = "T", dat = p.RNAseq6l3c3t, threshtype = "N",
#'   threshold = 10, absR = TRUE
#' )
#' cor.ESepiGen4c1l <- projectR(
#'   data = p.ESepiGen4c1l$mRNA.Seq,
#'   loadings = c.RNAseq6l3c3t, NP = "PositiveCOR",
#'   dataNames = map.ESepiGen4c1l[["GeneSymbols"]]
#' )
#'
#' @rdname projectR-methods
#' @aliases projectR
setMethod(
  "projectR", signature(data = "matrix", loadings = "correlateR"),
  function(data,
           loadings,
           dataNames = NULL,
           loadingsNames = NULL,
           NP = NA,
           full = FALSE,
           bootstrapPval = FALSE,
           bootIter = 1e3) {
    patterns <- loadings@corM
    if (!is.na(NP)) {
      patterns <- as.matrix(patterns[[NP]])
      colnames(patterns) <- NP
    } else {
      patterns <- loadings
    }

    if (length(patterns) == 2) {
      patterns <- do.call(rbind, patterns)
    } else {
      patterns <- as.matrix(patterns)
    }

    return(projectR(
      data = data, loadings = patterns, dataNames = dataNames,
      loadingsNames = loadingsNames, full = full,
      bootstrapPval = bootstrapPval, bootIter = bootIter
    ))
  }
)

###############################################################################
#' @param targetNumPatterns desired number of patterns with hclust
#' @param sourceData data used to create cluster object
#' @import limma
#' @import cluster
#' @importFrom stats cutree
#' @rdname projectR-methods
#' @aliases projectR
setMethod(
  "projectR", signature(data = "matrix", loadings = "hclust"),
  function(data, loadings, dataNames = NULL, loadingsNames = NULL, full = FALSE,
           targetNumPatterns, sourceData, bootstrapPval = FALSE, bootIter = 1000) {
    cut <- cutree(loadings, k = targetNumPatterns)
    patterns <- matrix(0, nrow = nrow(sourceData), ncol = targetNumPatterns)
    rownames(patterns) <- rownames(sourceData)
    for (x in 1:targetNumPatterns) {
      patterns[cut == x, x] <- apply(sourceData[cut == x, ], 1, cor,
        y = colMeans(sourceData[cut == x, ])
      )
    }
    return(projectR(data,
      loadings = patterns, dataNames = dataNames,
      loadingsNames = loadingsNames, full = full
    ))
  }
)

#' @rdname projectR-methods
#' @aliases projectR
setMethod(
  "projectR", signature(data = "matrix", loadings = "kmeans"),
  function(data, loadings, dataNames = NULL, loadingsNames = NULL, full = FALSE,
           sourceData, bootstrapPval = FALSE, bootIter = 1000) {
    patterns <- matrix(0, nrow = nrow(sourceData), ncol = length(loadings$size))
    rownames(patterns) <- rownames(sourceData)
    for (x in 1:length(loadings$size)) {
      patterns[loadings$cluster == x, x] <- apply(sourceData[loadings$cluster == x, ], 1,
        cor,
        y = colMeans(sourceData[loadings$cluster == x, ])
      )
    }
    return(projectR(data,
      loadings = patterns, dataNames = dataNames, full = full,
      bootstrapPval = bootstrapPval, bootIter = bootIter
    ))
  }
)

###############################################################################
#' @examples
#' library("projectR")
#' data(p.RNAseq6l3c3t)
#' nP <- 3
#' kClust <- kmeans(t(p.RNAseq6l3c3t), centers = nP)
#' kpattern <- cluster2pattern(clusters = kClust, NP = nP, data = p.RNAseq6l3c3t)
#' p <- as.matrix(p.RNAseq6l3c3t)
#' projectR(p, kpattern)
#'
#' @rdname projectR-methods
#' @aliases projectR
setMethod(
  "projectR", signature(data = "matrix", loadings = "cluster2pattern"),
  function(data, loadings, dataNames = NULL,
           loadingsNames = NULL, full = FALSE,
           sourceData, bootstrapPval = FALSE, bootIter = 1000) {
    loadings <- loadings@clusterMatrix
    loadings[is.na(loadings)] <- 0
    return(projectR(data,
      loadings = loadings, dataNames = dataNames, full = full,
      bootstrapPval = bootstrapPval, bootIter = bootIter
    ))
  }
)

###############################################################################
compareBoots <- function(projection, boots) {
  mat <- sapply(1:nrow(projection), function(i) {
    sapply(1:ncol(projection), function(j) {
      val <- sapply(1:length(boots), function(x) boots[[x]][i, j])
      valD <- ecdf(val)
      qt0 <- valD(0)
      if (qt0 < 0.5) {
        return(2 * qt0)
      } else {
        return(2 * (1 - qt0))
      }
    })
  })
  return(mat)
}
