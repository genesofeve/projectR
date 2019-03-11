#' @importFrom ROCR performance prediction
#' @title auc_mat
#'
#' @description Calculates AUC values for each set of weights for each label and outputs the results as a matrix
#' @param labels a vector of labels whose length is equal to the number of columns in the weight matrix
#' @param weights  a matrix of weights from projection analysis
#' @return A matrix of AUC values for each set of weights classifying each label.
#' @examples
#' projectR(data=p.ESepiGen4c1l$mRNA.Seq,loadings=AP.RNAseq6l3c3t$Amean, 
#' dataNames = map.ESepiGen4c1l[["GeneSymbols"]]) -> projection
#' auc_mat(pd.ESepiGen4c1l$Condition,projection)
#' @export

auc_mat<-function(labels=NA, weights=NA){
  results <- model.matrix(~labels-1)
  auc_matrix = matrix(nrow = dim(weights)[1], ncol = dim(results)[2])
  rownames(auc_matrix) = rownames(weights)
  colnames(auc_matrix) = colnames(results)
  for (i in 1:dim(weights)[1]) {
    for (j in 1:dim(results)[2]) {
      auc_matrix[i,j] = performance(prediction(weights[i,], results[,j]), measure='auc')@"y.values"[[1]]
    }
  }
  colnames(auc_matrix) = sort(unique(labels))
  return(auc_matrix)
}
