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
