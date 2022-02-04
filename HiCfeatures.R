library(rtracklayer)
#library(BSgenome.Celegans.UCSC.ce11)
#library(eulerr)
library(ggplot2)
library(dplyr)
library(tidyr)
#library(ComplexHeatmap)
#library(gridExtra)
#library(gtable)

workDir<-getwd()
if(!dir.exists(paste0(workDir,"/plots"))) {
  dir.create(paste0(workDir,"/plots"))
}

hicFeaturePath="/Users/semple/Documents/MeisterLab/otherPeopleProjects/Moushumi/hicFeatures"
RNAseqPath="/Users/semple/Documents/MeisterLab/otherPeopleProjects/Moushumi/2021_RNAseq_MDas"


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


listdf<-list()
grp=pcas$SMC[1]
for (grp in pcas$SMC){
  #pca1<-readRDS(paste0(workDir,"/otherData/",pcas$E1[pcas$SMC==grp]))
  pca2<-import(paste0(hicFeaturePath,"/otherData/",pcas$E2[pcas$SMC==grp]))

  pca2$AB<-ifelse(pca2$score>0,"A",0)
  pca2$AB[pca2$score<0]<-"B"
  pca2A<-pca2[pca2$AB=="A"]
  pca2B<-pca2[pca2$AB=="B"]
  pca20<-pca2[pca2$AB==0]
  pca2A<-reduce(pca2A)
  pca2B<-reduce(pca2B)
  pca20<-reduce(pca20)
  pca2A$AB<-"A"
  pca2B$AB<-"B"
  pca20$AB<-0

  pca2merged<-sort(c(pca2A,pca2B,pca20))
  pca2merged$SMC<-grp
  listdf[[grp]]<-as.data.frame(pca2merged)
}


df<-do.call(rbind,listdf)
df$SMC<-factor(df$SMC,levels=pcas$SMC,labels=pcas$prettyNames)
row.names(df)<-NULL
df$XvA<-ifelse(df$seqnames=="chrX","chrX","Autosomes")
df$AB<-factor(df$AB,levels=c("A","B"))
df<-df[!is.na(df$AB),]

med<-df%>%filter(SMC=="TEVonly") %>% group_by(XvA,AB) %>% summarise(med=median(width))

dfsum<-df %>% group_by(SMC,XvA,AB) %>% summarise(mean=mean(width),median=median(width))
dfsum$percentIncrease<-100*dfsum$mean/dfsum$mean[dfsum$SMC=="TEV only"]

p1<-ggplot(df,aes(x=AB,y=width,fill=AB)) +
  geom_boxplot(outlier.shape=NA,notch=T,varwidth=T) +
  scale_fill_manual(values=c("red","lightblue"))+
  coord_cartesian(ylim = c(0, 65000))+
  facet_grid(rows=vars(XvA),cols=vars(SMC)) + theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        text = element_text(size = 12))+
  #geom_hline(yintercept=8000) +
  stat_summary(fun="mean",geom="point",shape=8,size=2)
p1





###########################-
# compartments - digitized: 366tpm-----
###########################-

### Autosomes

pcas<-data.frame(SMC=c("TEVonly"),
                 strain =c("366"),
                 E1=NA, E2=NA)

tpm366<-import(paste0(RNAseqPath,"/tracks/PMW366_TPM_avr.bedgraph"),format="bedgraph")
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

# p3<-ggplot(df[df$XvA!="X",],aes(x=bin,y=stateFrequency,fill=state)) +
#   geom_col(position="stack") + scale_fill_manual(values=stateClrs) +
#   theme_bw() +
#   theme(legend.position = "none") +
#   facet_grid(rows=vars(compartment),cols=vars(SMC)) + ggtitle(paste0("Frequency of autosomal chromosome states in digitized eigen vector bins"))

p4<-ggplot(df[df$XvA!="X",],aes(x=bin,y=stateWidth,fill=state)) +
  geom_col(position="stack") + scale_fill_manual(values=stateClrs) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(cols=vars(compartment),rows=vars(SMC))+ ggtitle(paste0("Autosomes"))
p4




