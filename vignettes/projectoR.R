## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
options(scipen = 1, digits = 2)
set.seed(1234)

## ----install, eval=FALSE-------------------------------------------------
#  source("https://bioconductor.org/biocLite.R")
#  biocLite("projectoR")

## ----git.install, eval=FALSE---------------------------------------------
#  library(devtools)
#  install_github("projectoR", "genesofeve")

## ----base.projectoR, eval=FALSE------------------------------------------
#  library(projectoR)
#  projectoR(data = NA, AnnotionObj = NA, IDcol = "GeneSymbol",
#      Patterns = NA, NP = NA, full = FALSE)

## ----prcomp, echo=FALSE--------------------------------------------------
# data to define PCs
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
              panel.background = element_rect(fill = "transparent",colour=NA),
              legend.background = element_rect(fill = "transparent",colour=NA),
              plot.title = element_text(vjust = 0,hjust=0,face="bold")) +
        labs(title = "PCA of hPSC PolyA RNAseq",
            x=paste("PC1 (",pcVAR[1],"% of varience)",sep=""),
            y=paste("PC2 (",pcVAR[2],"% of varience)",sep=""))

## ----projectoR.prcomp, echo=FALSE----------------------------------------
# data to project into PCs from RNAseq6l3c3t expression data 
data(ESepiGen4c1l4)

library(projectoR)
PCA2ESepi <- projectoR(p.ESepiGen4c1l$mRNA.Seq,Patterns=pc.RNAseq6l3c3t,full=T, 
    AnnotionObj=map.ESepiGen4c1l, IDcol="GeneSymbols")



#dPCA2ESepi<- data.frame(cbind(t(PCA2ESepi[[1]]),pd.ESepiGen4c1l.4cond))
#str(dPCA2ESepi)

#plot pca
#library(ggplot2)
#setEpiCOL <- scale_colour_manual(values = c("red","green","blue","black"),
#                                 guide = guide_legend(title="Lineage"))

#pPC2ESepiGen4c1l <- ggplot(dPCA2ESepi, aes(x=PC1, y=PC2, colour=Condition)) + 
#      geom_point(size=5) + setEpiCOL + 
#      theme(legend.position=c(0,1), legend.justification=c(0,1),
#        panel.background = element_rect(fill = "white"),
#        legend.title = element_blank(),
#        text=element_text(family="Helvetica-Narrow",size=16)) +
#      labs(title = "Encode RNAseq in target PC1 & PC2", 
#          x=paste("Projected PC1 (",round(PCA2ESepi[[2]][1],2),"% of varience)",sep=""),
#          y=paste("Projected PC2 (",round(PCA2ESepi[[2]][2],2),"% of varience)",sep=""))

## ---- fig.show='hold', fig.width=10, fig.height=3, echo=FALSE------------
library(gridExtra)
#grid.arrange(pPCA,pPC2ESepiGen4c1l,nrow=1)

