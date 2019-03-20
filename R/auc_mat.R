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

auc_mat<-function(labels, weights){
results <- model.matrix(~labels-1)
colnames(results) <- unique(labels)
weights1 <- dim(weights)[1]
results2 <- dim(results)[2]
i <- rep(seq_len(weights1), each = results2)
j <- rep(seq_len(results2), times = weights1)
auc_res <- vapply(seq_len(weights1 * results2), function(k) {
    performance(prediction(weights[i[k],], results[,j[k]]), measure='auc')@"y.values"[[1]]
}, numeric(1))
auc_res <- unlist(auc_res)
auc_matrix <- matrix(auc_res, nrow=weights1, ncol=results2, byrow=TRUE)
rownames(auc_matrix) = rownames(weights)
colnames(auc_matrix) = colnames(results)
return(auc_matrix)
}
