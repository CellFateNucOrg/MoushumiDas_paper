library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)
library(seqplots)
#library(eulerr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(gridExtra)
library(gtable)
library(alluvial)
library(ggalluvial)
library(ggpubr)

workDir="."
if(!dir.exists(paste0(workDir,"/plots"))){
  dir.create(paste0(workDir,"/plots"))
}


ws235toCe11<-function(gr){
  seqlevels(gr)<-c("I","II","III","IV","V","X","MtDNA")
  seqlevels(gr)<-seqlevels(BSgenome.Celegans.UCSC.ce11::Celegans)
  seqinfo(gr)<-seqinfo(BSgenome.Celegans.UCSC.ce11::Celegans)
  return(gr)
}

####################-
# digitized switching --------
###################-


#' Combine Granges of bins and score switching from High/Medium/Low
#'
#' Two GRanges with the same digitzied eigen vector bins grouped into three quantiles
#' numbered 1-3 with 1 being low and 3 being high
#' These will combined and the switching between high/medium/low between the reference
#' GRanges and the pca GRanges will be scored and stored in a column called "switch"
#' @param pca GRanges of pca for which you want to score the change
#' @param ref GRanges of reference pca from which categories are changing
#' @param namePCA name of pca group. pca GRanges should have a column named "bin.<namePCA>"
#' @param nameREF name of reference group. GRanges should have a column named "bin.<nameREF>"
#' @result Dataframe with switch score as well as original bin and eigen scores.
#' @export
scoreSwitch<-function(pca,ref,namePCA,nameREF){
  pca<-full_join(data.frame(ref),data.frame(pca),
                 by=c("seqnames","start","end","width","strand","pca"),
                 suffix=c(paste0(".",nameREF),paste0(".",namePCA)))
  pca$switch<-factor(paste(pca[,paste0("bin.",nameREF)],pca[,paste0("bin.",namePCA)],sep="."),
                     levels=c("3.3","2.3","3.2","2.2","1.2","2.1","1.1","3.1","1.3"))
  levels(pca$switch)<-c("H->H","M->H","H->M","M->M","L->M","M->L","L->L","H->L","L->H")
  return(pca)
}



#####################-
# alluvial plots -----
#####################-

pcas<-data.frame(SMC=c("TEVonly","dpy26","kle2","scc1","coh1"),
                 strain =c("366","382","775","784","828"),
                 E1=NA, E2=NA)

refPCAs<-c("TEVonly")
otherPCAs<-c("dpy26","kle2","scc1","coh1")
prettyOtherPCAs<-c("dpy-26cs","kle-2cs","scc-1cs","coh-1cs")

#pcaPath<-"/Users/semple/Documents/MeisterLab/otherPeopleProjects/Moushumi/hicFeatures/otherData"
pcaPath=paste0(workDir,"/otherData")
E1files=list.files(paste0(pcaPath),
                   pattern="_10k\\.oriented\\.3quantiles_trans\\.eigs\\.vecs_E1\\.rds")
E2files=list.files(paste0(pcaPath),
                   pattern="_10k\\.oriented\\.3quantiles_trans\\.eigs\\.vecs_E2\\.rds")
pcas$E1<-E1files[match(pcas$strain,unlist(strsplit(E1files,"_10k\\.oriented\\.3quantiles_trans\\.eigs\\.vecs_E1\\.rds")))]
pcas$E2<-E2files[match(pcas$strain,unlist(strsplit(E2files,"_10k\\.oriented\\.3quantiles_trans\\.eigs\\.vecs_E2\\.rds")))]

#listdf<-NULL
refset=refPCAs[1]
ref1<-readRDS(paste0(pcaPath,"/",pcas$E1[pcas$SMC==refset]))
#bedClrEigenBin(ref1,"TEVonly_E1")
ref2<-readRDS(paste0(pcaPath,"/",pcas$E2[pcas$SMC==refset]))
#bedClrEigenBin(ref2,"TEVonly_E2")

bins<-readRDS(paste0(pcaPath,"/",pcas$E1[pcas$SMC==refset]))
mcols(bins)<-NULL


g=1
E1all<-list()
E2all<-list()
statdf<-NULL
for (g in 1:length(otherPCAs)){
  grp<-otherPCAs[g]
  prettyGrp<-prettyOtherPCAs[g]

  pca1<-readRDS(paste0(pcaPath,"/",pcas$E1[pcas$SMC==grp]))
  pca2<-readRDS(paste0(pcaPath,"/",pcas$E2[pcas$SMC==grp]))

  #bedClrEigenBin(pca1,paste0(grp,"_E1"))
  #bedClrEigenBin(pca2,paste0(grp,"_E2"))

  ref1<-readRDS(paste0(pcaPath,"/",pcas$E1[pcas$SMC==refset]))
  ref2<-readRDS(paste0(pcaPath,"/",pcas$E2[pcas$SMC==refset]))

  E1<-scoreSwitch(pca1,ref1,grp,refset)
  E2<-scoreSwitch(pca2,ref2,grp,refset)

  E1[,paste0("bin.",refset)]<-factor(E1[,paste0("bin.",refset)])
  levels(E1[,paste0("bin.",refset)])<-c("L","M","H")
  E1[,paste0("bin.",grp)]<-factor(E1[,paste0("bin.",grp)])
  levels(E1[,paste0("bin.",grp)])<-c("L","M","H")

  E2[,paste0("bin.",refset)]<-factor(E2[,paste0("bin.",refset)])
  levels(E2[,paste0("bin.",refset)])<-c("L","M","H")
  E2[,paste0("bin.",grp)]<-factor(E2[,paste0("bin.",grp)])
  levels(E2[,paste0("bin.",grp)])<-c("L","M","H")

  E1all[[grp]]<-cbind(E1,mcols(bins))
  E2all[[grp]]<-cbind(E2,mcols(bins))

  if(is.null(statdf)){
    statdf<-rbind(data.frame(SMC=prettyGrp,EigenVector="E1",table(E1$switch)),
                  data.frame(SMC=prettyGrp,EigenVector="E2",table(E2$switch)))
  } else {
    statdf<-rbind(statdf,
                  data.frame(SMC=prettyGrp,EigenVector="E1",table(E1$switch)),
                  data.frame(SMC=prettyGrp,EigenVector="E2",table(E2$switch)))
  }
}

statdfwide<-statdf %>% pivot_wider(names_from=Var1,values_from=Freq)

## 1) plot pairwise E1 ref and grp then E2 ref and grp
plotList<-list()
g=1
for (g in 1:length(otherPCAs)){
  grp<-otherPCAs[g]
  prettyGrp<-prettyOtherPCAs[g]

  twoEig<-inner_join(E1all[[grp]],E2all[[grp]],by=c("seqnames","start","end","width","strand"))

  tmp<-twoEig %>% group_by(bin.TEVonly.x,across(paste0("bin.",grp,".x")), bin.TEVonly.y,across(paste0("bin.",grp,".y"))) %>% summarise(count=n(),
                                                                                                                                       na.rm=T)

  p1<-ggplot(tmp,aes(y=count,axis1=bin.TEVonly.x,axis2=get(paste0("bin.",grp,".x"))))+
    geom_alluvium(aes(fill=bin.TEVonly.x),reverse=F)+
    geom_stratum(width = 1/8)+#, fill = "white", color = "black") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum)),reverse=F) +
    scale_x_discrete(limits = paste0("E1.", c("TEVonly",grp)), expand = c(.05, .05)) +
    scale_fill_brewer("E1 tercile",type = "qual", palette = "Set1",na.value = "grey50") +
    theme_bw() +
    ggtitle(paste0("E1 terciles: TEVonly->",prettyGrp))

  p2<-ggplot(tmp,aes(y=count,axis1=bin.TEVonly.y,axis2=get(paste0("bin.",grp,".y"))))+
    geom_alluvium(aes(fill=bin.TEVonly.y),reverse=F)+
    geom_stratum(width = 1/8)+#, fill = "white", color = "black") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum)),reverse=F) +
    scale_x_discrete(limits = paste0("E2.", c("TEVonly",grp)), expand = c(.05, .05)) +
    scale_fill_brewer("E2 tercile",type = "qual", palette = "Set1",na.value = "grey50") +
    theme_bw() +
    ggtitle(paste0("E2 terciles: TEVonly-> ",prettyGrp))
  plotList[[grp]]<-ggarrange(p1,p2,nrow=1)
}


p<-ggarrange(plotlist=plotList,nrow=4,ncol=1)
ggsave(paste0(workDir,"/plots/Figure_S4alluvial1.pdf"),plot=p,height=29,width=19,unit="cm")


## create combined E1-E2 for ref and plot changes to E1 and E2 in different grps
plotList<-list()
g=1
for (g in 1:length(otherPCAs)){
  grp<-otherPCAs[g]
  prettyGrp<-prettyOtherPCAs[g]

  twoEig<-inner_join(E1all[[grp]],E2all[[grp]],by=c("seqnames","start","end","width","strand"))
  twoEig$TEVonly.E1.E2<-factor(paste0(twoEig$bin.TEVonly.x,".",twoEig$bin.TEVonly.y),
                               levels=c("L.L","L.M","L.H","M.L","M.M","M.H","H.L","H.M","H.H"))

  tmp<-twoEig %>% group_by(TEVonly.E1.E2,
                           across(paste0("bin.",grp,".x")),
                           across(paste0("bin.",grp,".y"))) %>%
    summarise(count=n(),na.rm=T)

  p1<-ggplot(tmp,aes(y=count,axis1=get(paste0("bin.",grp,".x")),
                     axis2=TEVonly.E1.E2,
                     axis3=get(paste0("bin.",grp,".y"))))+
    geom_alluvium(aes(fill=TEVonly.E1.E2),reverse=F)+
    geom_stratum(width = 1/8)+#, fill = "white", color = "black") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum)),reverse=F) +
    scale_x_continuous(breaks =1:3, labels= c(paste0(grp,".E1"), "TEVonly.E1.E2",paste0(grp,".E2")))+
    scale_fill_brewer("E1.E2 terciles",type = "qual", palette = "YlOrRd",na.value = "grey50") +
    theme_bw() +
    ggtitle(paste0("E1.E2 terciles: TEVonly->",prettyGrp))

  twoEig$TEVonly.E1.E2<-factor(paste0(twoEig$bin.TEVonly.x,".",twoEig$bin.TEVonly.y),
                               levels=c("L.L","M.L","H.L","L.M","M.M","H.M","L.H","M.H","H.H"))
  tmp<-twoEig %>% group_by(TEVonly.E1.E2,
                           across(paste0("bin.",grp,".x")),
                           across(paste0("bin.",grp,".y"))) %>%
    summarise(count=n(),na.rm=T)
  p2<-ggplot(tmp,aes(y=count,axis1=get(paste0("bin.",grp,".x")),
                     axis2=TEVonly.E1.E2,
                     axis3=get(paste0("bin.",grp,".y"))))+
    geom_alluvium(aes(fill=TEVonly.E1.E2),reverse=F)+
    geom_stratum(width = 1/8)+#, fill = "white", color = "black") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum)),reverse=F) +
    scale_x_continuous(breaks =1:3, labels= c(paste0(grp,".E1"), "TEVonly.E1.E2",paste0(grp,".E2")))+
    scale_fill_brewer("E1.E2 terciles",type = "qual", palette = "YlOrRd",na.value = "grey50") +
    theme_bw() +
    ggtitle(paste0("E1.E2 terciles: TEVonly->",prettyGrp))

  plotList[[grp]]<-ggarrange(p1,p2,nrow=1)
}

p<-marrangeGrob(plotList,nrow=2,ncol=1)
ggsave(paste0(workDir,"/plots/Figure_S4alluvial2.pdf"),plot=p,height=29,width=19,unit="cm")


## 3. plot different combinations of grps for E1 with TEV in the middle. same for E2
pcaPairs<-combn(otherPCAs,2)

plotList<-list()
i=1
for (i in 1:ncol(pcaPairs)){
  grp1<-pcaPairs[1,i]
  grp2<-pcaPairs[2,i]
  prettyGrp1<-prettyOtherPCAs[i]
  prettyGrp2<-prettyOtherPCAs[i]

  twoGrpE1<-inner_join(E1all[[grp1]],E1all[[grp2]],by=c("seqnames","start","end","width","strand",
                                                        "score.TEVonly","bin.TEVonly","pca"))

  tmpE1<-twoGrpE1 %>% group_by(bin.TEVonly, across(paste0("bin.",grp1)),across(paste0("bin.",grp2))) %>% summarise(count=n(),na.rm=T)

  p1<-ggplot(tmpE1,aes(y=count,axis1=get(paste0("bin.",grp1)),axis2=bin.TEVonly,axis3=get(paste0("bin.",grp2))))+
    geom_alluvium(aes(fill=bin.TEVonly),reverse=F)+
    geom_stratum(width = 1/8)+#, fill = "white", color = "black") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum)),reverse=F) +
    scale_x_continuous(breaks = 1:3, labels=paste0("E1.", c(grp1,"TEVonly",grp2))) +
    scale_fill_brewer("E1 tercile",type = "qual", palette = "Set1",na.value = "grey50") +
    theme_bw() +
    ggtitle(paste0("E1 terciles: TEVonly->",grp1,"&",grp2))
  p1

  twoGrpE2<-inner_join(E2all[[grp1]],E2all[[grp2]],by=c("seqnames","start","end","width","strand",
                                                        "score.TEVonly","bin.TEVonly","pca"))

  tmpE2<-twoGrpE2 %>% group_by(bin.TEVonly, across(paste0("bin.",grp1)),across(paste0("bin.",grp2))) %>% summarise(count=n(),na.rm=T)
  p2<-ggplot(tmpE2,aes(y=count,axis1=get(paste0("bin.",grp1)),axis2=bin.TEVonly,axis3=get(paste0("bin.",grp2))))+
    geom_alluvium(aes(fill=bin.TEVonly),reverse=F)+
    geom_stratum(width = 1/8)+#, fill = "white", color = "black") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum)),reverse=F) +
    scale_x_continuous(breaks = 1:3, labels=paste0("E2.", c(grp1,"TEVonly",grp2))) +
    scale_fill_brewer("E2 tercile",type = "qual", palette = "Set1",na.value = "grey50") +
    theme_bw() +
    ggtitle(paste0("E2 terciles: TEVonly->",grp1,"&",grp2))
  p2
  plotList[[paste0(grp1,grp2)]]<-ggarrange(p1,p2,nrow=1)
}

p<-marrangeGrob(plotList,nrow=2,ncol=1)
ggsave(paste0(workDir,"/plots/Figure_S4alluvial3.pdf"),plot=p,height=29,width=19,unit="cm")
