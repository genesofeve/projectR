#' @title alluvial_mat
#'
#' @description Function to compute alluvial matrix
#' @param   new.projections new projection to visualize
#' @param   ct_anno Cell type annotation 
#' @return  A dataframe that can be used to plot an alluvial
#' 
#' @rawNamespace import(dplyr, except = c(filter,lag))
#' @import reshape2 tidyverse
# #' @examples
#' @export

#plot logical to ask if the plot should be generated or not
alluvial_mat<-function(new.projections, ct_anno){
  sigPatternIdx<-apply(new.projections$pval,1,function(x){if(min(x,na.rm=TRUE)<=0.05){return(TRUE)} else{return(FALSE)}})
  new.projections$qval<-t(apply(new.projections$pval,1,function(x){p.adjust(x,method="BH")}))
  sigPatternIdx<-apply(new.projections$qval,1,function(x){if(min(x,na.rm=T)<=0.01){return(TRUE)} else{return(FALSE)}})
  sig<-as.data.frame(t(new.projections$qval[sigPatternIdx,]<=0.01))
  DM<-as.data.frame(cbind("celltype"=ct_anno,sig))

  celltype_cells<-as.data.frame(table(ct_anno))
  colnames(celltype_cells)<-c("celltype","nCells_per_type")

  pattern_cells<-as.data.frame(colSums(sig*1,na.rm=T))
  colnames(pattern_cells)<-c("nCells_per_pattern")

  DM.summary<- DM %>%
    dplyr::select("celltype",starts_with("Patt")) %>%
    melt(id.vars=c('celltype'))
  DM.summary$value<-as.numeric(DM.summary$value)
  DM.summary<- as_tibble(DM.summary) %>%
    group_by('celltype',variable) %>%
    summarize(nCells=sum(value,na.rm=T))
  DM.summary<-merge(DM.summary,celltype_cells,by.x='celltype',by.y='celltype')
  DM.summary<-mutate(DM.summary,prop=nCells/nCells_per_type)
  DM.summary<-merge(DM.summary,pattern_cells,by.x='variable',by.y=0)
  DM.summary<- DM.summary %>%
    mutate(pattern_prop=nCells/nCells_per_pattern)
  return(DM.summary)
}
