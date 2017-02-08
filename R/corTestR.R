corTestR=function(x,y,r=T,p=T,df=F){
if(sum(!is.na(x))<=3|sum(!is.na(y))<=3){c=c(NA,NA,NA);names(c)=c("r","p-value","df");return(c)}
c=cor.test(x,y,use="pairwise.complete.obs")
c=c(c$estimate,c$p.value,c$parameter)
names(c)=c("r","p-value","df")
if(r&p&df){return(c)}
if(r&!p&!df){c=c["r"];return(c)}
if(!r&p&!df){c=c["p-value"];return(c)}
if(!r&!p&df){c=c["df"];return(c)}
if(r&p&!df){c=c[c("r","p-value")];return(c)}
if(r&!p&df){c=c[c("r","df")];return(c)}
if(!r&p&df){c=c[c("p-value","df")];return(c)}
if(!r&!p&!df){print("You have not defined any info to return (r==F, p==F, df==F).");return()}
}