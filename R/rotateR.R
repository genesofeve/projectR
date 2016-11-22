
#EG:
#x as PC1
#y as PC2
# where x1,y1 could be the center of a datatype cloud and
#       x2,y2 could be the center of a second , the
# line between which you want re-defined as (space rotatd to be) PC1

rotaPCArota=function(x1,y1,x2,y2,rotaSTART)
{

if(dim(rotaSTART)[2]!=2){print("rotaSTART must have 2 and only 2 columns in it for this function.");return()}

slp1=(y1-y2)/(x1-x2)
slp2=1/(-slp1)
atan2.mn=atan2(slp2,1)

theta=(pi/2)-atan2.mn

R=rbind(c(cos(theta),-sin(theta)),c(sin(theta),cos(theta)))

rotaNEW=t(R%*%t(rotaSTART))

return(rotaNEW)
}

# rotaNEW=rotaPCArota(PC1.mean.UNDIFF,PC2.mean.UNDIFF,PC1.mean.KSR,PC2.mean.KSR,pca$rota[,1:2])

# to chk rotatin worked:
# dataCEN=sweep(data,1,pca$center,"-")
# rotaXnew=(t(dataCEN))%*%rotaNEW
