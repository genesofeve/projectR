## ----prcomp, warning=FALSE-------------------------------------------------
# data to define PCs
library(projectR)
data(RNAseq6l3c3t)

# do PCA on RNAseq6l3c3t expression data 
pc.RNAseq6l3c3t<-prcomp(t(p.RNAseq6l3c3t))
pcVAR <- round(((pc.RNAseq6l3c3t$sdev)^2/sum(pc.RNAseq6l3c3t$sdev^2))*100,2)
dPCA <- data.frame(cbind(pc.RNAseq6l3c3t$x,pd.RNAseq6l3c3t))

#plot pca
library(ggplot2)
setCOL <- scale_colour_manual(values = c("blue","black","red"), name="Condition:") 
setFILL <- scale_fill_manual(values = c("blue","black","red"),guide = FALSE) 
setPCH <- scale_shape_manual(values=c(23,22,25,25,21,24),name="Cell Line:")

pPCA <- ggplot(dPCA, aes(x=PC1, y=PC2, colour=ID.cond, shape=ID.line, 
        fill=ID.cond)) +
        geom_point(aes(size=days),alpha=.6)+ 
        setCOL + setPCH  + setFILL +
        scale_size_area(breaks = c(2,4,6), name="Day") + 
        theme(legend.position=c(0,0), legend.justification=c(0,0),
              legend.direction = "horizontal",
              panel.background = element_rect(fill = "white",colour=NA),
              legend.background = element_rect(fill = "transparent",colour=NA),
              plot.title = element_text(vjust = 0,hjust=0,face="bold")) +
        labs(title = "PCA of hPSC PolyA RNAseq",
            x=paste("PC1 (",pcVAR[1],"% of varience)",sep=""),
            y=paste("PC2 (",pcVAR[2],"% of varience)",sep=""))

## ----projectR.prcomp, warning=FALSE----------------------------------------
# data to project into PCs from RNAseq6l3c3t expression data 
data(ESepiGen4c1l4)

library(projectR)
PCA2ESepi <- projectR(p.ESepiGen4c1l$mRNA.Seq,Patterns=pc.RNAseq6l3c3t,full=TRUE, 
    AnnotionObj=map.ESepiGen4c1l, IDcol="GeneSymbols")

pd.ESepiGen4c1l<-data.frame(Condition=sapply(colnames(p.ESepiGen4c1l$mRNA.Seq), function(x) unlist(strsplit(x,'_'))[1]),stringsAsFactors=FALSE)
pd.ESepiGen4c1l$color<-c("red","red","green","green","green","blue","blue","black","black")
names(pd.ESepiGen4c1l$color)<-pd.ESepiGen4c1l$Cond

dPCA2ESepi<- data.frame(cbind(PCA2ESepi,pd.ESepiGen4c1l))

#plot pca
library(ggplot2)
setEpiCOL <- scale_colour_manual(values = c("red","green","blue","black"),                              guide = guide_legend(title="Lineage"))

pPC2ESepiGen4c1l <- ggplot(dPCA2ESepi, aes(x=PC1, y=PC2, colour=Condition)) + 
      geom_point(size=5) + setEpiCOL + 
      theme(legend.position=c(0,0), legend.justification=c(0,0),
        panel.background = element_rect(fill = "white"),
        legend.direction = "horizontal",
        plot.title = element_text(vjust = 0,hjust=0,face="bold")) +
      labs(title = "Encode RNAseq in target PC1 & PC2", 
          x=paste("Projected PC1 (",round(PCA2ESepi[[2]][1],2),"% of varience)",sep=""),
          y=paste("Projected PC2 (",round(PCA2ESepi[[2]][2],2),"% of varience)",sep=""))


## ---- fig.show='hold', fig.width=10, fig.height=5, echo=FALSE--------------
library(gridExtra)
grid.arrange(pPCA,pPC2ESepiGen4c1l,nrow=1)

## ----correlateR-exp--------------------------------------------------------
# data to 
library(projectR)
data(RNAseq6l3c3t)

# get genes correlated to T 
cor2T<-correlateR(genes="T", dat=p.RNAseq6l3c3t, threshtype="N", threshold=10, absR=TRUE)

### heatmap of genes more correlated to T 
indx<-unlist(sapply(cor2T,rownames))
colnames(p.RNAseq6l3c3t)<-pd.RNAseq6l3c3t$sampleX
library(reshape2)
pm.RNAseq6l3c3t<-melt(cbind(p.RNAseq6l3c3t[indx,],indx))

library(gplots)
library(ggplot2)
library(viridis)
pCorT<-ggplot(pm.RNAseq6l3c3t, aes(variable, indx, fill = value)) + 
  geom_tile(colour="gray20", size=1.5, stat="identity") + 
  scale_fill_viridis(option="B") +
  xlab("") +  ylab("") +
  scale_y_discrete(limits=indx) + 
  ggtitle("Ten genes most highly pos & neg correlated with T") +
  theme(
    panel.background = element_rect(fill="gray20"),
    panel.border = element_rect(fill=NA,color="gray20", size=0.5, linetype="solid"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(), 
    axis.text = element_text(size=rel(1),hjust=1),
    axis.text.x = element_text(angle = 90,vjust=.5),
    legend.text = element_text(color="white", size=rel(1)),
    legend.background = element_rect(fill="gray20"),
    legend.position = "bottom",
    legend.title=element_blank()
)


## ---- fig.show='hold', fig.width=10, fig.height=5, echo=FALSE--------------
pCorT

