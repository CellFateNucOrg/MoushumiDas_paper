library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)
#library(eulerr)
library(ggplot2)
library(dplyr)
library(tidyr)
#library(ComplexHeatmap)
#library(gridExtra)
#library(gtable)
library(zoo)
library(ggpubr)

workDir<-getwd()
if(!dir.exists(paste0(workDir,"/plots"))) {
  dir.create(paste0(workDir,"/plots"))
}

#hicFeaturePath="/Users/semple/Documents/MeisterLab/otherPeopleProjects/Moushumi/hicFeatures"
#RNAseqPath="/Users/semple/Documents/MeisterLab/otherPeopleProjects/Moushumi/2021_RNAseq_MDas"
hicFeaturePath=workDir
RNAseqPath=paste0(workDir,"/tracks/2021_RNAseq_MDas")

##########################-
# domain sizes -------
##########################-

pcas<-data.frame(SMC=c("TEVonly","dpy26"),
                 strain =c("366","382"),
                 prettyNames=c("TEV only", "dpy-26cs"),
                 E1=NA, E2=NA)



E1files=list.files(paste0(hicFeaturePath,"/otherData"),
                   pattern="_merge_2000\\.oriented_E1\\.vecs\\.bw")
E2files=list.files(paste0(hicFeaturePath,"/otherData"),
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
  pca1<-import(paste0(hicFeaturePath,"/otherData/",pcas$E1[pcas$SMC==grp]))
  pca2<-import(paste0(hicFeaturePath,"/otherData/",pcas$E2[pcas$SMC==grp]))

  #pca1<-smootheByChr(pca1,winWidth=50)
  #pcasave<-pca1
  #plot(1:length(pca1),pca1$score,type="l")
  #lines(1:length(pca1),pca1$smScore,type="l",col="red")
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

saveRDS(df[df$eigen=="E1",],"./otherData/TCcomp_TEVonly&dpy26_removed0s.RDS")
saveRDS(df[df$eigen=="E2",],"./otherData/ABcomp_TEVonly&dpy26_removed0s.RDS")
df<-df[df$eigen=="E2",]

df$AB<-factor(df$AB,levels=c("A","B"))
df<-df[!is.na(df$AB),]
df$SMC<-factor(df$SMC,levels=pcas$SMC,labels=pcas$prettyNames)

#df<-readRDS("./otherData/ABcomp_TEVonly&dpy26.RDS")

med<-df%>%filter(SMC=="TEVonly") %>% group_by(XvA,AB) %>% summarise(med=median(width))

dfsum<-df %>% group_by(SMC,XvA,AB) %>% summarise(mean=mean(width),median=median(width))
dfsum$percentIncrease<-100*dfsum$mean/dfsum$mean[dfsum$SMC=="TEV only"]

# st<-compare_means(width~AB,df,group.by=c("SMC","XvA"),
#                   method="t.test",p.adjust.method="holm")
# st$padj.format<-ufs::formatPvalue(st$p.adj,digits=3)
# p1<-ggplot(df,aes(x=AB,y=width)) +
#   geom_boxplot(aes(fill=AB),outlier.shape=NA,notch=T,varwidth=T) +
#   scale_fill_manual(values=c("red","lightblue"))+
#   coord_cartesian(ylim = c(0, 100000))+
#   facet_grid(rows=vars(XvA),cols=vars(SMC)) + theme_bw() +
#   theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
#         text = element_text(size = 12))+
#   #geom_hline(yintercept=8000) +
#   stat_summary(fun="mean",geom="point",shape=4,size=2)+
#   geom_segment(x=1,y=91000,xend=2,yend=91000,size=0.01)+
#   geom_text(data=st,mapping=aes(x=SMC,y=width,label=padj.format),x=1.5,y=97000)+
#   geom_text(data = df %>% group_by(AB, SMC, XvA) %>%
#              summarize(Count = n(),width=86000),
#            aes(label = paste0("n=", Count)),
#            position = position_dodge(0.85), size=3.5, show.legend = FALSE)
# p1

st<-compare_means(width~SMC,df,group.by=c("AB","XvA"),p.adjust.method="fdr")
st$padj.format<-ufs::formatPvalue(st$p.adj,digits=3)
p1a<-ggplot(df,aes(x=SMC,y=width)) +
  geom_boxplot(outlier.shape=NA,notch=T,varwidth=T,fill="lightblue",alpha=0.5) +
  scale_fill_manual(values=c("lightblue"))+
  coord_cartesian(ylim = c(0, 100000))+
  facet_grid(rows=vars(XvA),cols=vars(AB)) + theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        text = element_text(size = 12))+
  #geom_hline(yintercept=8000) +
  stat_summary(fun="mean",geom="point",shape=4,size=2)+
  geom_segment(x=1,y=91000,xend=2,yend=91000,size=0.01)+
  geom_text(data=st,mapping=aes(x=SMC,y=width,label=padj.format),x=1.5,y=97000)#+
  # geom_text(data = df %>% group_by(AB, SMC, XvA) %>%
  #             summarize(Count = n(),width=86000),
  #           aes(label = paste0("n=", Count)),
  #           position = position_dodge(0.85), size=3.5, show.legend = FALSE)
p1a


df1<-df %>% group_by(SMC,XvA,AB,eigen) %>% summarise(averageWidth=mean(width),
                                                 medianWidth=median(width))

df1$pcAvrIncrease<-df1$averageWidth/df1$averageWidth[1:4]
df1





###########################-
# compartments - digitized: 366tpm-----
###########################-

### Autosomes

pcas<-data.frame(SMC=c("TEVonly"),
                 strain =c("366"),
                 E1=NA, E2=NA)

tpm366<-import(paste0(RNAseqPath,"/PMW366_TPM_avr.bedgraph"),format="bedgraph")
cov366<-coverage(tpm366,weight="score")

E1files=list.files(paste0(hicFeaturePath,"/otherData"),
                   pattern="_merge_2000\\.saddle_trans_A_noX_E1\\.digitized\\.rds")
E2files=list.files(paste0(hicFeaturePath,"/otherData"),
                   pattern="_merge_2000\\.saddle_trans_A_noX_E2\\.digitized\\.rds")
pcas$E1<-E1files[match(pcas$strain,unlist(strsplit(E1files,"_merge_2000\\.saddle_trans_A_noX_E1\\.digitized\\.rds")))]
pcas$E2<-E2files[match(pcas$strain,unlist(strsplit(E2files,"_merge_2000\\.saddle_trans_A_noX_E2\\.digitized\\.rds")))]
listdf<-NULL
for (grp in pcas$SMC){
  pca1<-readRDS(paste0(hicFeaturePath,"/otherData/",pcas$E1[pcas$SMC==grp]))
  pca2<-readRDS(paste0(hicFeaturePath,"/otherData/",pcas$E2[pcas$SMC==grp]))

  pca1<-binnedAverage(pca1,cov366,varname="tpm366")
  pca2<-binnedAverage(pca2,cov366,varname="tpm366")

  df1<-data.frame(pca1)
  df2<-data.frame(pca2)

  df<-rbind(df1,df2)
  #df$compartment<-factor(df$compartment)
  df$SMC<-grp

  listdf[[grp]]<-df
}

df<-do.call(rbind,listdf)
df$SMC<-factor(df$SMC,levels=pcas$SMC)
df$bin<-factor(df$bin,levels=1:50)

subdf<-df[df$SMC %in% c("wt","TEVonly"),]
subdf<-subdf[subdf$bin %in% 1:50,]
p2<-ggplot(subdf,aes(x=bin,y=log2(tpm366),fill=bin)) +
  geom_boxplot(outlier.shape=NA,size=0.1,fill="lightblue") + facet_grid(SMC~pca)+
  coord_cartesian(ylim=c(-14,14)) + theme_bw()+
  #scale_fill_manual(values=c("white","grey70"))+
  geom_hline(yintercept=0,col="black")+
  ggtitle(paste0("TEV control animals")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none",axis.text.x=element_blank(),axis.ticks.x=element_blank())

p2

head(df)


###########################-
# compartments - correlation: 366tpm-----
###########################-

### Autosomes

pcas<-data.frame(SMC=c("TEVonly"),
                 strain =c("366"),
                 E1=NA, E2=NA)

tpm366<-import(paste0(RNAseqPath,"/PMW366_TPM_avr.bedgraph"),format="bedgraph")
cov366<-coverage(tpm366,weight="score")

E1files=list.files(paste0(hicFeaturePath,"/otherData"),
                   pattern="_merge_2000\\.oriented_E1\\.vecs\\.bw")
E2files=list.files(paste0(hicFeaturePath,"/otherData"),
                   pattern="_merge_2000\\.oriented_E2\\.vecs\\.bw")
pcas$E1<-E1files[match(pcas$strain,unlist(strsplit(E1files,"_merge_2000\\.oriented_E1\\.vecs\\.bw")))]
pcas$E2<-E2files[match(pcas$strain,unlist(strsplit(E2files,"_merge_2000\\.oriented_E2\\.vecs\\.bw")))]

listdf<-NULL
for (grp in pcas$SMC){
  pca1<-import(paste0(hicFeaturePath,"/otherData/",pcas$E1[pcas$SMC==grp]))
  pca2<-import(paste0(hicFeaturePath,"/otherData/",pcas$E2[pcas$SMC==grp]))
  seqlevels(pca1)<-seqlevels(Celegans)
  seqlevels(pca2)<-seqlevels(Celegans)
  pca1<-binnedAverage(pca1,cov366,varname="tpm366")
  pca2<-binnedAverage(pca2,cov366,varname="tpm366")

  df1<-data.frame(pca1)
  df1$eigen<-"E1"
  df2<-data.frame(pca2)
  df2$eigen<-"E2"
  df<-rbind(df1,df2)
  #df$compartment<-factor(df$compartment)
  df$SMC<-grp

  listdf[[grp]]<-df
}

df<-do.call(rbind,listdf)
df<-df[df$tpm366>0,]# remove genes with little or no expression
df<-df[df$seqnames!="chrX",]

cor(df[df$eigen=="E1","score"],log2(df[df$eigen=="E1","tpm366"]),method="spearman")
cor(df[df$eigen=="E2","score"],log2(df[df$eigen=="E2","tpm366"]),method="spearman")



############################-
# overlap of states with digitized compartment---------
############################-


states<-import.bed("./publicData/chromStates_L3_Evans2016_ce11.bed")
#domains<-import.bed("./publicData/chromDomains_L3_Evans2016_ce11.bed")
#seqlevels(domains)<-seqlevels(states)

stateClrs<-c("#fe0003","#f59745","#008100","#74943a",
             "#c4d69c","#05ff7f","#ceff65","#fd0082",
             "#ff70d0","#ffb6c4","#7f00ff","#1900fe",
             "#528dd4","#814006","#b8cde4","#808080",
             "#dcdac3","#c4bc97","#938b54","#141414")

getStateOLtable<-function(gr,states){
  ol<-data.frame(findOverlaps(gr,states,ignore.strand=T,minoverlap=10))
  gr$name<-paste0(seqnames(gr),":",start(gr),"-",end(gr))
  ol$name<-gr$name[ol$queryHits]
  ol$XvA<-ifelse(unlist(strsplit(ol$name,":.*$"))=="chrX","X","A")
  ol$state<-factor(states$score[ol$subjectHits],levels=1:20)
  ol$width<-width(states)[ol$subjectHits]
  df<-ol%>%group_by(state,XvA) %>% summarise(stateFrequency=n(),stateWidth=sum(width))
  allStates<-data.frame(state=factor(1:20, levels=1:20))
  df<-left_join(allStates,df)
  df[is.na(df)]<-0
  return(df)
}



############################-
# overlap of states with digitized compartment X in cis---------
############################-


pcas<-data.frame(SMC=c("wt","TEVonly"),
                 strain =c("N2","366"),
                 E1=NA, E2=NA)

E1files=list.files(paste0(hicFeaturePath,"/otherData"),
                   pattern="_merge_2000\\.saddle_cis_X_noA_E1\\.digitized\\.rds")
E2files=list.files(paste0(hicFeaturePath,"/otherData"),
                   pattern="_merge_2000\\.saddle_cis_X_noA_E2\\.digitized\\.rds")
pcas$E1<-E1files[match(pcas$strain,unlist(strsplit(E1files,"_merge_2000\\.saddle_cis_X_noA_E1\\.digitized\\.rds")))]
pcas$E2<-E2files[match(pcas$strain,unlist(strsplit(E2files,"_merge_2000\\.saddle_cis_X_noA_E2\\.digitized\\.rds")))]


listdf<-NULL

#grp=pcas$SMC[1]
for (grp in pcas$SMC){
  pca1<-readRDS(paste0(hicFeaturePath,"/otherData/",pcas$E1[pcas$SMC==grp]))
  pca2<-readRDS(paste0(hicFeaturePath,"/otherData/",pcas$E2[pcas$SMC==grp]))

  for(binNum in 1:50){
    pca1bin<-pca1[pca1$bin==binNum]
    pca2bin<-pca2[pca2$bin==binNum]
    df1<-getStateOLtable(gr=pca1bin,states)
    df2<-getStateOLtable(gr=pca2bin,states)
    df1$bin<-binNum
    df2$bin<-binNum

    df1$SMC<-grp
    df2$SMC<-grp

    df1$compartment<-"E1"
    df2$compartment<-"E2"

    listdf[[paste0(grp,"_",binNum)]]<-rbind(df1,df2)
  }
}


df<-do.call(rbind,listdf)
df$SMC<-factor(df$SMC,levels=pcas$SMC)


df$bin<-factor(df$bin,levels=1:50)

# p3<-ggplot(df[df$XvA=="X",],aes(x=bin,y=stateFrequency,fill=state)) +
#   geom_col(position="stack") + scale_fill_manual(values=stateClrs) +
#   theme_bw() +
#   theme(legend.position = "none") +
#   facet_grid(rows=vars(compartment),cols=vars(SMC)) + ggtitle(paste0("Frequency of chrX chromosome states in digitized eigen vector bins"))

p3<-ggplot(df[df$XvA=="X",],aes(x=bin,y=stateWidth,fill=state)) +
  geom_col(position="stack") + scale_fill_manual(values=stateClrs) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(cols=vars(compartment),rows=vars(SMC))+ ggtitle(paste0("X chromosome"))







############################-
# overlap of states with digitized compartment autosomes in trans---------
############################-


pcas<-data.frame(SMC=c("wt","TEVonly"),
                 strain =c("N2","366"),
                 E1=NA, E2=NA)

E1files=list.files(paste0(hicFeaturePath,"/otherData"),
                   pattern="_merge_2000\\.saddle_trans_A_noX_E1\\.digitized\\.rds")
E2files=list.files(paste0(hicFeaturePath,"/otherData"),
                   pattern="_merge_2000\\.saddle_trans_A_noX_E2\\.digitized\\.rds")
pcas$E1<-E1files[match(pcas$strain,unlist(strsplit(E1files,"_merge_2000\\.saddle_trans_A_noX_E1\\.digitized\\.rds")))]
pcas$E2<-E2files[match(pcas$strain,unlist(strsplit(E2files,"_merge_2000\\.saddle_trans_A_noX_E2\\.digitized\\.rds")))]


listdf<-NULL

#grp=pcas$SMC[1]
for (grp in pcas$SMC){
  pca1<-readRDS(paste0(hicFeaturePath,"/otherData/",pcas$E1[pcas$SMC==grp]))
  pca2<-readRDS(paste0(hicFeaturePath,"/otherData/",pcas$E2[pcas$SMC==grp]))

  for(binNum in 1:50){
    pca1bin<-pca1[pca1$bin==binNum]
    pca2bin<-pca2[pca2$bin==binNum]
    df1<-getStateOLtable(gr=pca1bin,states)
    df2<-getStateOLtable(gr=pca2bin,states)
    df1$bin<-binNum
    df2$bin<-binNum

    df1$SMC<-grp
    df2$SMC<-grp

    df1$compartment<-"large (E1)"
    df2$compartment<-"small (E2)"

    listdf[[paste0(grp,"_",binNum)]]<-rbind(df1,df2)
  }
}


df<-do.call(rbind,listdf)
df$SMC<-factor(df$SMC,levels=pcas$SMC)


df$bin<-factor(df$bin,levels=1:50)


p4<-ggplot(df[df$XvA!="X",],aes(x=bin,y=stateWidth,fill=state)) +
  geom_col(position="stack") + scale_fill_manual(values=stateClrs) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(cols=vars(compartment),rows=vars(SMC))+ ggtitle(paste0("Autosomes"))
p4

p4a<-ggplot(df[df$XvA!="X" & df$SMC=="wt",],aes(x=bin,y=stateWidth,fill=state)) +
  geom_col(position="stack") + scale_fill_manual(values=stateClrs) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(cols=vars(compartment),rows=vars(SMC))+ ggtitle(paste0("Autosomes"))
p4a





# p<-ggarrange(ggarrange(p2,p4a,nrow=2,labels=c("","D")),p1,ncol=2,labels=c("C","H"))
# p
#
# p<-annotate_figure(p, top = text_grob("Das et al., Figure 2", size = 14))
# ggsave(paste0(workDir,"/plots/HiCfeatures_Fig2.pdf"),p,device="pdf",width=21,height=14,units="cm")

p<-ggarrange(ggarrange(p2,p4a,nrow=2,labels=c("","D")),p1a,ncol=2,labels=c("C","H"))
p

p<-annotate_figure(p, top = text_grob("Das et al., Figure 2", size = 14))
ggsave(paste0(workDir,"/plots/HiCfeatures_Fig2a.pdf"),p,device="pdf",width=21,height=14,units="cm")


p<-ggarrange(p4,p3,ncol=2,labels=c("B",""))
p
p<-annotate_figure(p, top = text_grob("Das et al., Figure S2", size = 14))
ggsave(paste0(workDir,"/plots/HiCfeatures_FigS2.pdf"),p,device="pdf",width=21,height=7,units="cm")


