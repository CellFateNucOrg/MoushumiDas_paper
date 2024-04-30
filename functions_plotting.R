#' Format scientific exponent
#'
#'  rstatix::p_format does good job but still leaves scientific
#'  notation with an "e". This function converts the e to standard
#'  x10 and superscript exponent. html verion not tested.
#'  @param pvals vector of p values
#'  @param html boolean to say whether to convert to html markdown (T) or plotmath notation (F)
#'  @return vector of formatted pvals
#'  @export
prettyExponents<-function(pvals,html=F){
  if(html==T){
    pvals<-gsub("e","x10<sup>",pvals)
    pvals<-gsub("x10<sup>(.*)$","x10<sup>\\1</sup>", pvals)
  } else {
    pvals<-gsub("e","%*%10",pvals)
    pvals<-gsub("10(.*)$","10^\\1", pvals)
    pvals<-gsub("<","p<",pvals,fixed=T)
  }
  return(pvals)
}

# usage:
# st<-compare_means(width~SMC,df,group.by=c("AB","XvA"),p.adjust.method="fdr",
# method="wilcox.test") %>% p_format(new.col=T,accuracy=1e-32)

# st$html.format<-prettyExponents(st$p.adj.format)
#   geom_text(data=st,mapping=aes(x=SMC,y=width,label=html.format),x=1.5,y=97000/1e3,
# parse=T)+
