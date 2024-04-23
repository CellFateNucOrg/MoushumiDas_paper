library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)
library(ggplot2)
library(dplyr)
library(tidyr)
library(zoo)
library(ggpubr)
library(ufs)

projectDir="."
otherDataDir=paste0(projectDir,"/otherData")
publicDataDir=paste0(projectDir,"/publicData")
finalFigDir=paste0(projectDir,"/finalFigures")
if(!dir.exists(finalFigDir)){
  dir.create(finalFigDir)
}


theme_set(
  theme_classic()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          title=element_text(size=9),
          axis.title.y=ggtext::element_markdown(size=9),
          axis.title.x=ggtext::element_markdown(size=9)
    )
)


RNAseqPath=paste0(projectDir,"/tracks/2021_RNAseq_MDas")

##########################-
# domain sizes -------
##########################-

pcas<-data.frame(SMC=c("TEVonly","dpy26"),
                 strain =c("366","382"),
                 prettyNames=c("TEV only", "dpy-26cs"),
                 E1=NA, E2=NA)



E1files=list.files(paste0(projectDir,"/otherData"),
                   pattern="_merge_2000\\.oriented_E1\\.vecs\\.bw")
E2files=list.files(paste0(projectDir,"/otherData"),
                   pattern="_merge_2000\\.oriented_E2\\.vecs\\.bw")
pcas$E1<-E1files[match(pcas$strain,unlist(strsplit(E1files,"_merge_2000\\.oriented_E1\\.vecs\\.bw")))]
pcas$E2<-E2files[match(pcas$strain,unlist(strsplit(E2files,"_merge_2000\\.oriented_E2\\.vecs\\.bw")))]

mergeAdjacentDomains<-function(pca){
  pcaA<-pca[pca$AB=="A"]
  pcaB<-pca[pca$AB=="B"]
  pca0<-pca[pca$AB==0]
  pcaA<-GenomicRanges::reduce(pcaA)
  pcaB<-GenomicRanges::reduce(pcaB)
  pca0<-GenomicRanges::reduce(pca0)
  pcaA$AB<-"A"
  pcaB$AB<-"B"
  pca0$AB<-0
  pcamerged<-sort(c(pcaA,pcaB,pca0))
  return(pcamerged)
}
#gr<-pca1
smootheByChr<-function(gr,winWidth=25){
  grl<-GenomicRanges::split(gr$score,seqnames(gr))
  tmpl<-lapply(grl,rollmean,k=winWidth,fill=0)
  gr$smScore<-unlist(tmpl)
  return(gr)
}


listdf<-list()
grp=pcas$SMC[1]
for (grp in pcas$SMC){
  pca1<-import(paste0(projectDir,"/otherData/",pcas$E1[pcas$SMC==grp]))
  pca2<-import(paste0(projectDir,"/otherData/",pcas$E2[pcas$SMC==grp]))

  pca1$AB<-ifelse(pca1$score>0,"A",0)
  pca1$AB[pca1$score<0]<-"B"
  pca1<-mergeAdjacentDomains(pca1)
  idx<-which(pca1$AB==0)
  idx<-idx[ifelse(idx[1]==1,2,1):length(idx)] # remove first index if it is 1
  sameChr<-as.vector(seqnames(pca1)[idx]==seqnames(pca1)[idx-1])
  pca1$AB[idx[sameChr]]<-pca1$AB[(idx-1)[sameChr]]
  pca1<-mergeAdjacentDomains(pca1)


  pca2$AB<-ifelse(pca2$score>0,"A",0)
  pca2$AB[pca2$score<0]<-"B"
  pca2<-mergeAdjacentDomains(pca2)
  idx<-which(pca2$AB==0)
  idx<-idx[ifelse(idx[1]==1,2,1):length(idx)] # remove first index if it is 1
  sameChr<-as.vector(seqnames(pca2)[idx]==seqnames(pca2)[idx-1])
  pca2$AB[idx[sameChr]]<-pca2$AB[(idx-1)[sameChr]]
  pca2<-mergeAdjacentDomains(pca2)

  pca1$SMC<-grp
  pca2$SMC<-grp
  pca1$eigen<-"E1"
  pca2$eigen<-"E2"
  listdf[[paste0(grp,"_E1")]]<-as.data.frame(pca1)
  listdf[[paste0(grp,"_E2")]]<-as.data.frame(pca2)
}

lapply(listdf,dim)
df<-do.call(rbind,listdf)
row.names(df)<-NULL
df$XvA<-ifelse(df$seqnames=="chrX","chrX","Autosomes")

df<-df[df$eigen=="E2",]

df$AB<-factor(df$AB,levels=c("A","B"))
df<-df[!is.na(df$AB),]
df$SMC<-factor(df$SMC,levels=pcas$SMC,labels=pcas$prettyNames)


med<-df%>%filter(SMC=="TEVonly") %>% group_by(XvA,AB) %>% summarise(med=median(width))

dfsum<-df %>% group_by(SMC,XvA,AB) %>% summarise(mean=mean(width),median=median(width))
dfsum$percentIncrease<-100*dfsum$mean/dfsum$mean[dfsum$SMC=="TEV only"]


st<-compare_means(width~SMC,df,group.by=c("AB","XvA"),p.adjust.method="fdr")
st$padj.format<-ufs::formatPvalue(st$p.adj,digits=3)
p1a<-ggplot(df,aes(x=SMC,y=width)) +
  geom_boxplot(outlier.shape=NA,notch=T,varwidth=T,fill="lightblue",alpha=0.5) +
  scale_fill_manual(values=c("lightblue"))+
  coord_cartesian(ylim = c(0, 100000))+
  facet_grid(rows=vars(XvA),cols=vars(AB)) +
  theme(text = element_text(size = 12))+
  stat_summary(fun="mean",geom="point",shape=4,size=2)+
  geom_segment(x=1,y=91000,xend=2,yend=91000,linewidth=0.01)+
  geom_text(data=st,mapping=aes(x=SMC,y=width,label=padj.format),x=1.5,y=97000)#+

p1a


df1<-df %>% group_by(SMC,XvA,AB,eigen) %>% summarise(averageWidth=mean(width),
                                                     medianWidth=median(width))

df1$pcAvrIncrease<-df1$averageWidth/df1$averageWidth[1:4]
df1

ggpubr::ggexport(p1a,filename=paste0(finalFigDir, "/Fig_3i_compartmentSizes.pdf"),
                 device="pdf",width=3,height=5, units="cm")

