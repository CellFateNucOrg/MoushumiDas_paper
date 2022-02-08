library(rtracklayer)
#library(GenomicInteractions)
library(ggplot2)
#library(EnhancedVolcano)
library(BSgenome.Celegans.UCSC.ce11)
library(dplyr)
library(ggpubr)
library(ggpubr)

source("functions.R")
#source("./variableSettings.R")

workDir="."
outputNamePrefix=""



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

ws235toCe11<-function(gr){
  seqlevels(gr)<-c("I","II","III","IV","V","X","MtDNA")
  seqlevels(gr)<-seqlevels(BSgenome.Celegans.UCSC.ce11::Celegans)
  seqinfo(gr)<-seqinfo(BSgenome.Celegans.UCSC.ce11::Celegans)
  return(gr)
}


###########################-
# compartments - digitized ----------------------------------------------------
###########################-

## 366TPM in digitized compartments of different HiCs -----

### Autosomes ------

pcas<-data.frame(SMC=c("wt","TEVonly"),
                 strain =c("N2","366"),
                 E1=NA, E2=NA)

tpm366<-import(paste0(workDir,"/otherData/PMW366_TPM_avr.bedgraph"),
               format="bedgraph")
#tpm366<-import(paste0(workDir,"/otherData/sumFR_366_B_UniqueMultiple.bw"),
#                format="bw")
#tpm366<-ws235toCe11(tpm366)
cov366<-coverage(tpm366,weight="score")

E1files=list.files(paste0(workDir,"/otherData"),
                   pattern="_merge_2000\\.saddle_trans_A_noX_E1\\.digitized\\.rds")
E2files=list.files(paste0(workDir,"/otherData"),
                   pattern="_merge_2000\\.saddle_trans_A_noX_E2\\.digitized\\.rds")
pcas$E1<-E1files[match(pcas$strain,unlist(strsplit(E1files,"_merge_2000\\.saddle_trans_A_noX_E1\\.digitized\\.rds")))]
pcas$E2<-E2files[match(pcas$strain,unlist(strsplit(E2files,"_merge_2000\\.saddle_trans_A_noX_E2\\.digitized\\.rds")))]
listdf<-NULL
grp=pcas$SMC[1]
for (grp in pcas$SMC){

  pca1<-readRDS(paste0(workDir,"/otherData/",pcas$E1[pcas$SMC==grp]))
  pca2<-readRDS(paste0(workDir,"/otherData/",pcas$E2[pcas$SMC==grp]))

  pca1<-GRanges(pca1)
  pca2<-GRanges(pca2)
  start(pca1)<-start(pca1)+1
  start(pca2)<-start(pca2)+1
  seqlevels(pca1)<-seqlevels(BSgenome.Celegans.UCSC.ce11::Celegans)
  seqlevels(pca2)<-seqlevels(BSgenome.Celegans.UCSC.ce11::Celegans)

  pca1<-binnedAverage(pca1,cov366,varname="tpm366",na.rm=T)
  pca2<-binnedAverage(pca2,cov366,varname="tpm366",na.rm=T)

  df1<-data.frame(pca1)
  df2<-data.frame(pca2)

  df1$pca<-"E1"
  df2$pca<-"E2"

  colnames(df1)<-gsub("^E.?\\.d","bin",colnames(df1))
  colnames(df2)<-gsub("^E.?\\.d","bin",colnames(df2))

  df<-rbind(df1,df2)
  df$SMC<-grp

  listdf[[grp]]<-df
}

df<-do.call(rbind,listdf)
df$SMC<-factor(df$SMC,levels=pcas$SMC)
df$bin<-factor(df$bin,levels=1:50)

subdf<-df[df$SMC %in% c("wt","TEVonly"),]
p1<-ggplot(subdf,aes(x=bin,y=log2(tpm366),fill=bin)) +
  geom_boxplot(outlier.shape=NA,size=0.1,fill="lightblue") + facet_grid(SMC~pca)+
  coord_cartesian(ylim=c(-12,12)) + theme_bw()+
  geom_hline(yintercept=0,col="red")+
  ggtitle(paste0("366 TPM in different autosomal bins of digitized pca")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none",axis.text.x=element_blank(),axis.ticks.x=element_blank())

listdf<-NULL
for (grp in pcas$SMC){
  pca1<-readRDS(paste0(workDir,"/otherData/",pcas$E1[pcas$SMC==grp]))
  pca2<-readRDS(paste0(workDir,"/otherData/",pcas$E2[pcas$SMC==grp]))

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

# autosomes
p2<-ggplot(df[df$XvA=="A",],aes(x=bin,y=stateWidth,fill=state)) +
  geom_col(position="stack") + scale_fill_manual(values=stateClrs) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(rows=vars(compartment),cols=vars(SMC))+ ggtitle(paste0("Number of bp per chromosome state in digitized eigen vector bins"))


### chrX -------
pcas<-data.frame(SMC=c("wt","TEVonly"),
                 strain =c("N2","366"),
                 E1=NA, E2=NA)

tpm366<-import(paste0(workDir,"/otherData/PMW366_TPM_avr.bedgraph"),
               format="bedgraph")
#tpm366<-import(paste0(workDir,"/otherData/sumFR_366_B_UniqueMultiple.bw"),
#                format="bw")
#tpm366<-ws235toCe11(tpm366)
cov366<-coverage(tpm366,weight="score")

E1files=list.files(paste0(workDir,"/otherData"),
                   pattern="_merge_2000\\.saddle_cis_X_noA_E1\\.digitized\\.rds")
E2files=list.files(paste0(workDir,"/otherData"),
                   pattern="_merge_2000\\.saddle_cis_X_noA_E2\\.digitized\\.rds")
pcas$E1<-E1files[match(pcas$strain,unlist(strsplit(E1files,"_merge_2000\\.saddle_cis_X_noA_E1\\.digitized\\.rds")))]
pcas$E2<-E2files[match(pcas$strain,unlist(strsplit(E2files,"_merge_2000\\.saddle_cis_X_noA_E2\\.digitized\\.rds")))]
listdf<-NULL
for (grp in pcas$SMC){
  pca1<-readRDS(paste0(workDir,"/otherData/",pcas$E1[pcas$SMC==grp]))
  pca2<-readRDS(paste0(workDir,"/otherData/",pcas$E2[pcas$SMC==grp]))

  pca1<-GRanges(pca1)
  pca2<-GRanges(pca2)
  start(pca1)<-start(pca1)+1
  start(pca2)<-start(pca2)+1
  seqlevels(pca1)<-seqlevels(BSgenome.Celegans.UCSC.ce11::Celegans)
  seqlevels(pca2)<-seqlevels(BSgenome.Celegans.UCSC.ce11::Celegans)

  pca1<-binnedAverage(pca1,cov366,varname="tpm366",na.rm=T)
  pca2<-binnedAverage(pca2,cov366,varname="tpm366",na.rm=T)

  df1<-data.frame(pca1)
  df2<-data.frame(pca2)

  df1$pca<-"E1"
  df2$pca<-"E2"

  colnames(df1)<-gsub("^E.?\\.d","bin",colnames(df1))
  colnames(df2)<-gsub("^E.?\\.d","bin",colnames(df2))

  df<-rbind(df1,df2)
  df$SMC<-grp

  listdf[[grp]]<-df
}

df<-do.call(rbind,listdf)
df$SMC<-factor(df$SMC,levels=pcas$SMC)
df$bin<-factor(df$bin,levels=1:50)

subdf<-df[df$SMC %in% c("wt","TEVonly"),]
p3<-ggplot(subdf,aes(x=bin,y=log2(tpm366),fill=bin)) +
  geom_boxplot(outlier.shape=NA,size=0.1,fill="lightblue") + facet_grid(SMC~pca)+
  coord_cartesian(ylim=c(-12,12)) + theme_bw()+
  geom_hline(yintercept=0,col="red")+
  ggtitle(paste0("366 TPM in different chrX bins of digitized pca")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none",axis.text.x=element_blank(),axis.ticks.x=element_blank())

listdf<-NULL
for (grp in pcas$SMC){
  pca1<-readRDS(paste0(workDir,"/otherData/",pcas$E1[pcas$SMC==grp]))
  pca2<-readRDS(paste0(workDir,"/otherData/",pcas$E2[pcas$SMC==grp]))

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

# chrX
p4<-ggplot(df[df$XvA=="X",],aes(x=bin,y=stateWidth,fill=state)) +
  geom_col(position="stack") + scale_fill_manual(values=stateClrs) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(rows=vars(compartment),cols=vars(SMC))+ ggtitle(paste0("Number of bp per chrX chromosome state in digitized eigen vector bins"))




p<-ggpubr::ggarrange(ggpubr::ggarrange(p1,p3,nrow=1,ncol=2),
                     ggpubr::ggarrange(p2,p4,nrow=1,ncol=2),
                     nrow=2,col=1)

ggpubr::ggexport(p,filename=paste0(workDir, "/plots/",outputNamePrefix,
                         "digitizedCompControls_chrAvX.pdf"),
       device="pdf",width=19,height=19, units="cm")




