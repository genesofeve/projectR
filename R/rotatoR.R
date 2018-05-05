#TODO: Define rotatoR S4 class here or make rotatoR a generic function


#' @title rotatoR
#'
#' @description a function for rotating two basis about a point or line in that plain
#' @param x1  a value describing a the coordinate of a point in the first basis. If no values are provided for x2
#' @param y1  a value describing a the coordinate of a point in the second basis
#' @param x2  a value describing a the coordinate of the second point in the second basis
#' @param y2 a value describing a the coordinate of the second point in the second basis
#' @param basisSET the basis to be rotated
#' @return An object of class rotatoR.
#' @examples
#'  pca.RNAseq6l3c3t<-prcomp(t(p.RNAseq6l3c3t))
#'  r.RNAseq6l3c3t<-rotatoR(1,1,-1,-1,pca.RNAseq6l3c3t$x[,1:2])
#' @export


rotatoR <- function(x1,y1,x2,y2,basisSET){

if(dim(basisSET)[2]!=2){print("basisSET must have 2 and only 2 columns in it for this function.");return()}

slp1=(y1-y2)/(x1-x2)
slp2=1/(-slp1)
atan2.mn=atan2(slp2,1)

theta=(pi/2)-atan2.mn

R=rbind(c(cos(theta),-sin(theta)),c(sin(theta),cos(theta)))

rotaNEW=t(R%*%t(basisSET))

class(rotaNEW) <- append(class(rotaNEW),"rotatoR") #Can't do this directly with S4 withouth a class definition.

return(rotaNEW)
}

# rotaNEW=rotaPCArota(PC1.mean.UNDIFF,PC2.mean.UNDIFF,PC1.mean.KSR,PC2.mean.KSR,pca$rota[,1:2])

# to chk rotatin worked:
# dataCEN=sweep(data,1,pca$center,"-")
# rotaXnew=(t(dataCEN))%*%rotaNEW

#R=rbind(c(cos(t),-sin(t)),c(sin(t),cos(t)))

