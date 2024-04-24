library(tidyverse)
library(ggplot2)
library("ggpmisc")
library(eulerr)
library(lattice)
library(rtracklayer)
library(GenomicRanges)
library(ggpubr)
library(Cairo)
library(BSgenome.Celegans.UCSC.ce11)

projectDir="."
otherDataDir=paste0(projectDir,"/otherData")
publicDataDir=paste0(projectDir,"/publicData")
finalFigDir=paste0(projectDir,"/finalFigures")
if(!dir.exists(finalFigDir)){
  dir.create(finalFigDir)
}


myTheme<-theme_set(
    theme_classic()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          title=element_text(size=10),
          axis.title.y=ggtext::element_markdown(size=10),
          axis.title.x=ggtext::element_markdown(size=10),
          strip.text = element_text(size = 10),
          axis.text=element_text(size=10)
    )
)




contrastsOI<-c("X.wt.wt.0mM_dpy26cs_vs_wt","X.wt.wt.0mM_kle2cs_vs_wt",
               "X.wt.wt.0mM_scc16cs_vs_wt","X.wt.wt.0mM_coh1cs_vs_wt",
               "X.wt.wt.0mM_scc1coh1cs_vs_wt")[c(3,4,5,1,2)]
useContrasts<-c("dpy26","kle2","scc1","coh1","scc1coh1")[c(3,4,5,1,2)]

prettyNames<-c(substitute(italic(x^cs),list(x="dpy-26")),
               substitute(italic(x^cs),list(x="kle-2")),
               substitute(italic(x^cs),list(x="scc-1")),
               substitute(italic(x^cs),list(x="coh-1")),
               substitute(italic(x^cs*y^cs),list(x="scc-1",y="coh-1")))[c(3,4,5,1,2)]

complexes<-c("condensin~I/I^DC", "condensin~II", "cohesin^{SCC-1}","cohesin^{COH-1}","cohesins")[c(3,4,5,1,2)]
# ?plotmath for parsing details

names(complexes)<-useContrasts

names(contrastsOI)<-useContrasts
# strains<-c("366","382","775","784","828","844")
# strain<-factor(strains,levels=strains)
# SMC<-strain
# levels(SMC)<-c("TEVonly",useContrasts)

#controlGrp<-levels(SMC)[1] # control group
groupsOI<-useContrasts


source(paste0(projectDir,"/functions.R"))



# panel cdf of lfc
lfcVal=0.5
padjVal=0.05
filterPrefix<-"filtCycChrAX"
fileNamePrefix=paste0("p",padjVal,"_lfc",lfcVal,"_",filterPrefix,"/",filterPrefix,"_")
outPath=paste0(projectDir)


########################-
## ECDF of data -----
########################-

# using automatic filtering threshold
sigTables<-list()
localPadj=padjVal
localLFC=0
for (grp in groupsOI){
  print(grp)
  salmon<-readRDS(paste0(outPath,"/",fileNamePrefix,contrastsOI[[grp]], "_DESeq2_fullResults_p",padjVal,".rds"))
  salmon<-salmon[!is.na(salmon$padj),]
  #nrow(filterResults(salmon,padj=0.05,lfc=-0.5,direction="lt",chr="autosomes"))
  print(paste0(nrow(salmon)," genes before filtering"))
  print(paste0(sum(is.na(salmon$padj))," have padj that is NA"))
  #salmon$expressed<-sum(salmon$baseMean>10)
  sigTables[[grp]]<-as.data.frame(salmon) #[salmon$baseMean>10,]
  print(paste0(nrow(sigTables[[grp]])," genes after automatic threshold filter"))
}

SMC<-rep(names(sigTables),lapply(sigTables,nrow))
sig<-do.call(rbind,sigTables)
sig$SMCpretty<-factor(SMC,levels=useContrasts,labels=prettyNames)
sig$SMC<-factor(SMC,levels=useContrasts)
table(sig$SMC)
sig$XvA<-ifelse(sig$chr=="chrX","chrX","Autosomes")
#sig$XvA[sig$chr=="chrX"]<-"chrX"
#sig$XvA<-factor(sig$XvA)
table(sig$XvA,sig$SMC)
sig$upVdown[sig$log2FoldChange<0]<-"down"
sig$upVdown[sig$log2FoldChange>0]<-"up"
sig$upVdown<-factor(sig$upVdown)
table(sig$upVdown)
row.names(sig)<-NULL
SMC<-NULL

lapply(sigTables,function(x){sum(x$padj<0.05)})

options(tibble.width=Inf)
dd1<-sig %>% dplyr::filter(padj<localPadj) %>%
  dplyr::group_by(SMC,upVdown,XvA) %>%
  dplyr::mutate(ecd=ecdf(abs(log2FoldChange))(abs(log2FoldChange)))

dd1$upVdown<-relevel(dd1$upVdown,"up")

p<-ggplot(dd1, aes(x=abs(log2FoldChange),y=ecd,color=SMCpretty,linetype=XvA)) +
  geom_line(size=0.9)+ facet_wrap(vars(upVdown),nrow=2)+
  coord_cartesian(xlim=c(0,1.5)) +
  scale_color_discrete(labels = ggplot2:::parse_safe(levels(dd1$SMCpretty)),
                       name=element_blank())+
  scale_linetype_discrete(name=element_blank()) +
  xlab("Absolute log<sub>2</sub> fold change")+ylab("Fraction of significant genes")


#stat_ecdf(aes(colour=SMC,linetype=XvA),alpha=0.7)
p<-p+geom_vline(aes(xintercept = 0.5), color="grey") +
  annotate("text",label="0.5",size=4, x=0.5, y=1, hjust=-0.05,color="grey") +
  geom_vline(aes(xintercept = 0.15), color="grey") +
  annotate("text",label="0.15",size=4, x=0.15, y=1,hjust=-0.05,color="grey")



# to get fraction between 0.15 and 0.5 lfc
dd2<-sig %>% dplyr::filter(padj<localPadj) %>%
  dplyr::group_by(SMC,upVdown,XvA) %>%
  dplyr::summarise(ecd15=ecdf(abs(log2FoldChange))(c(0.15)),
                   ecd50=ecdf(abs(log2FoldChange))(c(0.5)),
                   q2=100*(ecd50-ecd15),
                   count=n())
dd2[order(dd2$q2),]

sig %>% dplyr::filter(padj<localPadj,abs(log2FoldChange)>0.5) %>%
  dplyr::group_by(SMC,upVdown) %>%
  dplyr::summarise(count=n())

sig %>% dplyr::group_by(SMC) %>% dplyr::mutate(allGenes=n()) %>%
  dplyr::group_by(SMC) %>%
  dplyr::filter(padj<localPadj,abs(log2FoldChange)>0.5) %>%
  dplyr::summarise(count=n(),allGenes=allGenes) %>%
  dplyr::distinct() %>% dplyr::mutate(pc=100*count/allGenes)

sig %>% dplyr::group_by(SMC,XvA) %>% dplyr::mutate(allGenes=n()) %>%
  dplyr::group_by(SMC,XvA,upVdown) %>%
  dplyr::filter(padj<localPadj,abs(log2FoldChange)>0.5) %>%
  dplyr::summarise(count=n(),allGenes=allGenes) %>%
  dplyr::distinct() %>% dplyr::mutate(pc=100*count/allGenes)

# to get total expressed genes
na.omit(sig) %>% group_by(SMC,XvA) %>% summarise(count=n())

# to add tables of number of signifcant genes to plot
ttu<-with(sig[sig$padj<localPadj & sig$upVdown=="up",],table(XvA,SMCpretty))
ttu<-pivot_wider(as.data.frame(ttu),names_from=SMCpretty,values_from=Freq)
names(ttu)<-gsub("XvA","",names(ttu))

ttd<-with(sig[sig$padj<localPadj & sig$upVdown=="down",],table(XvA,SMCpretty))
ttd<-pivot_wider(as.data.frame(ttd),names_from=SMCpretty,values_from=Freq)
names(ttd)<-gsub("XvA","",names(ttd))

tbs<-list("down"=ttd,
          "up"=ttu)

df <- tibble(x = rep(1.5, 2),
             y = rep(0, 2),
             upVdown=factor(c("down","up"),levels=c("up","down")),
             tbl = tbs)

p1<-p + geom_table(data = df, aes(x = x, y = y,label = tbl), parse=T,
                   table.theme = ttheme_gtlight, inherit.aes=F,size=2.5)

p1



##################-
## boxplot of LFC------
##################-
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/",fileNamePrefix,contrastsOI[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
  sigTables[[grp]]<-as.data.frame(salmon)
  sigTables[[grp]]$SMC<-grp
}


allSig<-do.call(rbind,sigTables)
rownames(allSig)<-NULL
allSig$chr<-gsub("chr","",allSig$chr)
allSig$SMC<- factor(allSig$SMC,levels=groupsOI,labels=prettyNames)
allSig$complexes<-complexes[allSig$SMC]
facetLabels<-allSig %>% select(SMC,chr,complexes) %>% distinct()
facetLabels$complexes[facetLabels$chr!="III"]<-""

p2<-ggplot(allSig,aes(x=chr,y=log2FoldChange,fill=chr)) +
  geom_boxplot(outlier.shape=NA,show.legend=F) +
  facet_grid(cols=vars(SMC),labeller=label_parsed) +
  coord_cartesian(ylim=c(-0.7,1.5)) + scale_fill_grey(start=0.8, end=0.4)+
  geom_hline(yintercept=0,linetype="dashed",col="red")+
  xlab("Chromosome") + ylab("Log<sub>2</sub>FC") +
  geom_text(data=facetLabels,aes(label=complexes),parse=T,x=2,y=1.4,size=3.5,
            hjust=0)
p2

# ##################-
# ## boxplot LFC up/down significant------
# ##################-
# localLFCval<-0
# sigTables<-list()
# for (grp in useContrasts){
#   salmon<-readRDS(paste0(outPath,"/",fileNamePrefix,contrastsOI[[grp]],
#                          "_DESeq2_fullResults_p",padjVal,".rds"))
#   salmon<-salmon[!is.na(salmon$chr),]
#   rownames(salmon)<-NULL
#   sigTables[[grp]]<-as.data.frame(getSignificantGenes(salmon,
#                                                       padj=padjVal, lfc=localLFCval,
#                                                       namePadjCol="padj",
#                                                       nameLfcCol="log2FoldChange",
#                                                       direction="both",
#                                                       chr="all", nameChrCol="chr"))
# }
#
# # upregulated
# sigList<-lapply(sigTables, getSignificantGenes,
#                 padj=padjVal,lfc=localLFCval,direction="gt")
#
# sigList<-lapply(sigList, "[", ,c("wormbaseID","chr","log2FoldChange"))
# for(g in names(sigList)){ sigList[[g]]$SMC<-g }
# sigList<-do.call(rbind,sigList)
# sigList$updown<-"up"
# sigTbl<-sigList
#
#
# # downregulated
# sigList<-lapply(sigTables, getSignificantGenes,
#                 padj= padjVal,lfc= -localLFCval,direction="lt")
#
# sigList<-lapply(sigList, "[", ,c("wormbaseID","chr","log2FoldChange"))
# for(g in names(sigList)){ sigList[[g]]$SMC<-g }
# sigList<-do.call(rbind,sigList)
# sigList$updown<-"down"
# sigTbl<-rbind(sigTbl,sigList)
# sigTbl$chr<-as.factor(sigTbl$chr)
# sigTbl$SMC<-factor(sigTbl$SMC, levels=groupsOI,labels=prettyNames)
# sigTbl$updown<-factor(sigTbl$updown, levels=c("up","down"))
# rownames(sigTbl)<-NULL
# sigTbl$XvA<-ifelse(sigTbl$chr=="chrX","chrX","Autosomes")
#
# yminmax=c(0,median(abs(sigTbl$log2FoldChange))+quantile(abs(sigTbl$log2FoldChange))[4]*2)
# p2a<-ggplot(sigTbl,aes(x=updown,y=abs(log2FoldChange),fill=SMC)) +
#   geom_boxplot(notch=T, varwidth=F, position=position_dodge2(padding=0.2),outlier.shape=NA,lwd=0.1,fatten=3) +
#   facet_grid(cols=vars(XvA)) + coord_cartesian(ylim=yminmax) +
#   ggtitle("Absolute LFC of significantly changed genes by chromosome type") +
#   theme_minimal() + xlab("")+ ylab("|log2(FC)|")+
#   scale_fill_brewer(palette="Dark2",labels = ggplot2:::parse_safe(levels(sigTbl$SMC)))




###################-
## volcano plots------
###################-

plotList<-list()
geneTable<-NULL
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/",fileNamePrefix,contrastsOI[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
  p<-plotVolcanoXvA(salmon,addLegend=F) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  plotList[[grp]]<-p+ggtitle(label=prettyNames[[which(groupsOI==grp)]])
}


leg <- get_legend(plotVolcanoXvA(salmon,addLegend=T))
plotList[["legend"]]<-as_ggplot(leg)
p5<-ggarrange(plotlist=plotList,ncol=3,nrow=2)
p5




############### Final assembly #########

p<-ggarrange(p1,p2,p5,nrow=3,heights=c(2.5,1.5,3),labels=c("a ","b ","c "))
p<-annotate_figure(p, top = text_grob("Das et al., Figure S6", size = 14)) #face=bold
ggsave(paste0(finalFigDir,"/Fig_S6_RNAseq_TEV.pdf"),p,device=cairo_pdf,width=21,height=29.7,
       unit="cm")
ggsave(paste0(finalFigDir,"/Fig_S6_RNAseq_TEV.png"),p,device="png",width=21,height=29.7,
       unit="cm",bg="white")


