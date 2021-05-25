#' RNAseqing and ChIPSeq of matched genes in differentiated human iPS cells 
#'
#' map.ESepiGen4c1l contains gene annotations 
#' 
#' @references
#' 1. Gifford, C. A. et al. Transcriptional and epigenetic dynamics 
#' during specification of human embryonic stem cells. Cell 153, 1149-1163 (2013).
#' 
#' 
#' @name map.ESepiGen4c1l 
#' @docType data
#' @format A data frames with 93 rows and 9 variables:
"map.ESepiGen4c1l"

#' RNAseqing from human 3 iPS & 3 ES cell lines
#' in 3 experimental condition at 3 time points
#'
#' map.RNAseq6l3c3 contains gene annotations for polyA bulk sequencing of 6 cell lines
#' in 3 experimental condition at 3 time points. 
#' 
#' @name  map.RNAseq6l3c3t
#' @docType data
#' @format A data frames with 108 rows and 54 variables:
"map.RNAseq6l3c3t"

#' RNAseqing from human 3 iPS & 3 ES cell lines
#' in 3 experimental condition at 3 time points
#' 
#' p.RNAseq6l3c3 contains log2(RPKM + 1) values for polyA bulk sequencing of 6 cell lines
#' in 3 experimental condition at 3 time points. 
#' 
#' @name  p.RNAseq6l3c3t
#' @docType data
#' @format A data frames with 108 rows and 54 variables:
"p.RNAseq6l3c3t"

#' RNAseqing and ChIPSeq of matched genes in differentiated human iPS cells 
#'
#' pd.ESepiGen4c1l.4cond contains sample phenotype and experimental information
#' 
#' @references
#' 1. Gifford, C. A. et al. Transcriptional and epigenetic dynamics 
#' during specification of human embryonic stem cells. Cell 153, 1149-1163 (2013). 
#' 
#' @name pd.ESepiGen4c1l 
#' @docType data
#' @format A data frames with 9 rows and 2 variables:
"pd.ESepiGen4c1l"

#' RNAseqing from human 3 iPS & 3 ES cell lines
#' in 3 experimental condition at 3 time points
#'
#' pd.RNAseq6l3c3t contains sample phenotype and experimental information
#' for polyA bulk sequencing of 6 cell lines in 3 experimental condition at 3 time points. 
#' 
#' @name  pd.RNAseq6l3c3t
#' @docType data
#' @format A data frames with 54 rows and 38 variables:
"pd.RNAseq6l3c3t"

#' RNAseqing and ChIPSeq of matched genes in differentiated human iPS cells 
#'
#' p.ESepiGen4c1l contains log2(RPKM + 1) values for polyA bulk sequencing 
#' and log2 counts of normalized ChIPSeq reads of 1 cell lines with 2
#' replicates in 4 experimental conditions at a single time point. 
#' 
#' 
#' @references
#' 1. Gifford, C. A. et al. Transcriptional and epigenetic dynamics 
#' during specification of human embryonic stem cells. Cell 153, 1149-1163 (2013).
#' 
#' 
#' @name p.ESepiGen4c1l 
#' @docType data
#' @format p.ESepiGen4c1l is a list of 6 data frames each with with 93 rows and between 4 and 9 variables:
"p.ESepiGen4c1l"

#' CoGAPS patterns and genes weights for p.RNAseq6l3c3t
#'
#' AP.RNAseq6l3c3t contains the output of the gapsRun function in the 
#' CoGAPS package for data = p.RNAseq6l3c3t
#'
#' @name AP.RNAseq6l3c3t
#' @docType data
#' @format A list of 12 items 
"AP.RNAseq6l3c3t"

#' CoGAPS patterns learned from the developing mouse retina.
#'
#' @references 
#' 1. Clark, B.S., & Stein-O'Brien G.L., et al. Single-Cell RNA-Seq Analysis of Development Identifies NFI 
#' Factors as Regulating Mitotic Exit and Late-Born Cell Specification. Cell 102, 1111-1126 (2019).
#' @name retinal_patterns
#' @docType data
#' @format A gene (rows) by pattern (column) matrix
"retinal_patterns"

#' log-normalized count data from microglial cells in the p6 mouse cortex.
#' @name microglial_counts
#' @docType data 
#' @format A gene (rows) by cell (column) matrix
"microglial_counts"

#' log-normalized count data from astrocytes and oligodendrocytes in the p6 mouse cortex.
#' @name glial_counts
#' @docType data 
#' @format A gene (rows) by cell (column) matrix
"glial_counts"

