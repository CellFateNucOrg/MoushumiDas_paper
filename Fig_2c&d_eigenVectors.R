library(rtracklayer)
library(ggplot2)
library(BSgenome.Celegans.UCSC.ce11)
library(dplyr)
library(ggpubr)

source("functions.R")
#source("./variableSettings.R")


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
          axis.title.x=ggtext::element_markdown(size=9),
          strip.text = element_text(size = 9),
          axis.text=element_text(size=9),
          panel.border = element_rect(colour = "black", fill=NA, size=0.8)
    )
)






states<-import.bed(paste0(publicDataDir,"/chromStates_L3_Evans2016_ce11.bed"))


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
# compartments - correlation ----------------------------------------------------
###########################-
pcas<-data.frame(SMC=c("TEVonly"),
                 strain =c("366"),
                 E1=NA, E2=NA)
tpm366<-import(paste0(projectDir,"/tracks/2021_RNAseq_MDas/PMW366_TPM_avr.bedgraph"),
               format="bedgraph")
cov366<-coverage(tpm366,weight="score")

E1files=list.files(paste0(projectDir,"/otherData"),
                   pattern="_merge_2000\\.oriented_E1\\.vecs\\.bw")
E2files=list.files(paste0(projectDir,"/otherData"),
                   pattern="_merge_2000\\.oriented_E2\\.vecs\\.bw")
pcas$E1<-E1files[match(pcas$strain,unlist(strsplit(E1files,"_merge_2000\\.oriented_E1\\.vecs\\.bw")))]
pcas$E2<-E2files[match(pcas$strain,unlist(strsplit(E2files,"_merge_2000\\.oriented_E2\\.vecs\\.bw")))]
listdf<-NULL
for (grp in pcas$SMC){
  pca1<-import.bw(paste0(projectDir,"/otherData/",pcas$E1[pcas$SMC==grp]))
  pca2<-import.bw(paste0(projectDir,"/otherData/",pcas$E2[pcas$SMC==grp]))

  seqlevels(pca1)<-seqlevels(BSgenome.Celegans.UCSC.ce11::Celegans)
  seqlevels(pca2)<-seqlevels(BSgenome.Celegans.UCSC.ce11::Celegans)

  pca1<-binnedAverage(pca1,cov366,varname="tpm366")
  pca2<-binnedAverage(pca2,cov366,varname="tpm366")

  df1<-data.frame(pca1)
  df2<-data.frame(pca2)

  df1$pca<-"E1"
  df2$pca<-"E2"

  df1$compartment<-NA
  df2$compartment<-NA
  df1$compartment<-ifelse(df1$score>0,"A","B")
  df2$compartment<-ifelse(df2$score>0,"A","B")

  df<-rbind(df1,df2)
  df$compartment<-factor(df$compartment)
  df$SMC<-grp

  listdf[[grp]]<-df
}

df<-do.call(rbind,listdf)
df$SMC<-factor(df$SMC,levels=pcas$SMC)

corMethod="spearman"
tpmThresh=0
allBins<-nrow(df)
df<-df[df$tpm366>=tpmThresh,]
df<-df[df$seqnames!="chrX",]
fracKept<-nrow(df)/allBins # 0.75 of bins
fracKept # 0.75 of bins

with(df[df$SMC=="TEVonly" & df$pca=="E1",], cor.test(score,log2(tpm366),method="spearman"))
with(df[df$SMC=="TEVonly" & df$pca=="E2",], cor.test(score,log2(tpm366),method="spearman"))
nrow(df[df$SMC=="TEVonly" & df$pca=="E1",]) # number of regions

ggplot(df,aes(x=score,y=log2(tpm366))) +
  geom_bin2d(bins=100)+
  facet_grid(cols=vars(pca)) + geom_smooth(method="lm") +
  stat_cor(label.x = -1.5, label.y = 18, size=3,method=corMethod) +
  ggtitle(paste0(corMethod," correlation of PCA eigen value vs PMW366 TPM (for bins > ",
                 formatC(tpmThresh,big.mark=",",format="G"),"tpm, ",
                 round(100*fracKept,1),"% of bins)")) #+

# slight discrepancy between ggplot stat_cor and cor.test(), 0.24 vs 0.25.

###########################-
# compartments - digitized ----------------------------------------------------
###########################-

## 366 TPM in digitized eigen vector bins

### Autosomes ------

pcas<-data.frame(SMC=c("TEVonly"),
                 strain =c("366"),
                 E1=NA, E2=NA)

tpm366<-import(paste0(projectDir,"/tracks/2021_RNAseq_MDas/PMW366_TPM_avr.bedgraph"),
               format="bedgraph")

cov366<-coverage(tpm366,weight="score")

E1files=list.files(paste0(projectDir,"/otherData"),
                   pattern="_merge_2000\\.saddle_trans_A_noX_E1\\.digitized\\.rds")
E2files=list.files(paste0(projectDir,"/otherData"),
                   pattern="_merge_2000\\.saddle_trans_A_noX_E2\\.digitized\\.rds")
pcas$E1<-E1files[match(pcas$strain,unlist(strsplit(E1files,"_merge_2000\\.saddle_trans_A_noX_E1\\.digitized\\.rds")))]
pcas$E2<-E2files[match(pcas$strain,unlist(strsplit(E2files,"_merge_2000\\.saddle_trans_A_noX_E2\\.digitized\\.rds")))]
listdf<-NULL
grp=pcas$SMC[1]
for (grp in pcas$SMC){

  pca1<-readRDS(paste0(projectDir,"/otherData/",pcas$E1[pcas$SMC==grp]))
  pca2<-readRDS(paste0(projectDir,"/otherData/",pcas$E2[pcas$SMC==grp]))

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

subdf<-df[df$SMC %in% c("TEVonly") & df$seqnames!="chrX" & !is.na(df$bin),]
p1<-ggplot(subdf,aes(x=bin,y=log2(tpm366),fill=bin)) +
  geom_boxplot(outlier.shape=NA,size=0.1,fill="lightblue") +
  coord_cartesian(ylim=c(-12,12)) + facet_grid(cols=vars(pca))+
  geom_hline(yintercept=0,col="red")+
  ggtitle(paste0("Autosomes")) + ylab("Log<sub>2</sub>TPM")+
  theme(legend.position="none",axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.title.x=element_blank(),
        plot.title = element_text(hjust = 1))

p1
binCount<-subdf %>% group_by(bin,pca) %>% summarise(count=n())
binCount

# chromatin states-------
listdf<-NULL
for (grp in pcas$SMC){
  pca1<-readRDS(paste0(projectDir,"/otherData/",pcas$E1[pcas$SMC==grp]))
  pca2<-readRDS(paste0(projectDir,"/otherData/",pcas$E2[pcas$SMC==grp]))

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

subdf<-df[df$SMC %in% c("TEVonly") & df$seqnames!="chrX" & !is.na(df$bin),]
# autosomes
p2<-ggplot(df[df$XvA=="A",],aes(x=bin,y=stateWidth/1e6,fill=state)) +
  geom_col(position="stack") + scale_fill_manual(values=stateClrs) +
  theme(legend.position = "none", axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),plot.title = element_text(hjust = 1)) +
  facet_grid(cols=vars(compartment))+ ggtitle(paste0("Autosomes"))+
  ylab("bp per chromatin state (Mb)") + xlab("Eigenvector bin")


p2


p<-ggpubr::ggarrange(p1,p2,nrow=2,col=1,labels=c("c ","d "))

ggpubr::ggexport(p,filename=paste0(finalFigDir, "/Fig_2c&d_eigenVectors.pdf"),
       device="pdf",width=5,height=9, units="cm")




