#' @export

alluvialMat<-function(projection, annotations, annotationName = "celltype", annotationSample = "cell", plot = FALSE, minPropExplained = 0.75){
  require(dplyr)
  require(reshape2)
  require(ggalluvial)
  require(viridis)
  require(RColorBrewer)
  if(!('pval' %in% names(projection))){
    stop("Please set arguemnt full = TRUE in projectR to generate projection with p-values")
  }
  sigPatternIdx<-apply(projection$pval,1,function(x){if(min(x,na.rm=TRUE)<=0.05){return(TRUE)} else{return(FALSE)}})
  projection$qval<-t(apply(projection$pval,1,function(x){p.adjust(x,method="BH")}))
  sigPatternIdx<-apply(projection$qval,1,function(x){if(min(x,na.rm=T)<=0.01){return(TRUE)} else{return(FALSE)}})
  sig<-as.data.frame(t(projection$qval[sigPatternIdx,]<=0.01))
  DM<-as.data.frame(cbind('celltype'=annotations,sig))  #possible issue when the numbe of annotations is less than significant patterns
  celltype_cells<-as.data.frame(table(annotations)) 
  colnames(celltype_cells)<-c('celltype','nCells_per_type')
  colnames(DM)[1] <- 'celltype' 
  pattern_cells<-as.data.frame(colSums(sig*1,na.rm=T))
  colnames(pattern_cells)<-c('nCells_per_pattern')
  DM.summary<- DM %>%
    dplyr::select('celltype',starts_with("Pat")) %>%
    melt(id.vars=c('celltype'))
  DM.summary$value<-as.numeric(DM.summary$value)
  DM.summary<- as_tibble(DM.summary) %>%
    group_by(celltype,variable) %>%
    summarize(nCells=sum(value,na.rm=T))
  DM.summary<-merge(DM.summary,celltype_cells,by.x='celltype',by.y='celltype')
  DM.summary<-mutate(DM.summary,prop=nCells/nCells_per_type)
  DM.summary<-merge(DM.summary,pattern_cells,by.x='variable',by.y=0)
  DM.summary<- DM.summary %>%
    mutate(pattern_prop=nCells/nCells_per_pattern)
  if(plot == TRUE){
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
minProp <- minPropExplained
plot.data<-subset(DM.summary,prop>minProp)
nPatterns<-length(unique(plot.data$variable))
nCelltype<-length(unique(plot.data$celltype))
p<-ggplot(plot.data,aes(y=prop,axis1=celltype,axis2=variable)) +
  geom_alluvium(aes(fill=celltype),color="black",size=0.2) + 
  geom_stratum(width=1/12,fill="grey50",color="black") + 
  geom_label(stat="stratum",label.strata=TRUE) +
  scale_x_continuous(breaks=1:2, labels=c(annotationName, "Pattern")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=getPalette(nCelltype)) + guides(fill=FALSE) + 
  ggtitle(paste0("Pattern explains at least ",minProp*100,"% of ",annotationSample,"s in a given type"))
plot(p)
  }
  return(DM.summary)
}