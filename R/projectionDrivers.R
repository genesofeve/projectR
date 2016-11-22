

#' @title <brief desc>
#'
#' @description <full description>
#' @param delta a vector of weights describing
#' @param A
#' @param Order
#' @param ... 
#' @export
#' @keywords
#' @seealso
#' @return
#' @aliases
#' @examples \dontrun{
#'
#'}
projectionDrivers <- function(delta=d4, #a vector of weights describing 
	A=Aneu,
	Order="increasing",
	...){
    delta <- delta[rownames(A)]
    A <- A[match(names(delta),rownames(A)),]
    axd <- apply(A,2,function(x) x*delta)
    if(Order=="increasing"){AxDgenes <- apply(axd,2,function(x) names(x[order(x)]))}
    if(Order=="decreasing"){AxDgenes <- apply(axd,2,function(x) names(x[order(-x)]))}
    if(Order=="maxdelta"){AxDgenes <- apply(axd,2,function(x) names(x[order(-abs(x))]))}
    return(AxDgenes)
}

#ProjectionDrivers(delta=d4,A=Aneu,Save=T,Fname="AxDgenes",)




