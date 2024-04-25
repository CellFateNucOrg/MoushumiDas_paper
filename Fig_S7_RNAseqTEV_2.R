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
library(ComplexHeatmap)
library(gridtext)

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

groupsOI<-useContrasts


source(paste0(projectDir,"/functions.R"))

# panel cdf of lfc
lfcVal=0.5
padjVal=0.05
filterPrefix<-"filtCycChrAX"
fileNamePrefix=paste0("p",padjVal,"_lfc",lfcVal,"_",filterPrefix,"/",filterPrefix,"_")



##################-
## venn diagrams------
##################-
eulerLabelsType=c("counts")
numGenes<-list()
## upregulated genes -----
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(projectDir,"/",fileNamePrefix,contrastsOI[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
  sigTables[[grp]]<-as.data.frame(
    getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="gt",
                        chr="all", nameChrCol="chr"))
}

plotTitle<-list(label="Autosomal",fontsize=8)
achr<-lapply(sigTables,function(x) x[x$chr!="chrX",])

subset1<-achr[c("dpy26","kle2","scc1")]
sigGenes<-lapply(subset1, "[[","wormbaseID")
uglyNames<-c("dpy-26cs","kle-2cs","scc-1cs")
names(sigGenes)<-uglyNames
numGenes["Aup"]<-sum(sapply(sigGenes,length))
fit<-euler(sigGenes)
p1a<-plot(fit, quantities=list(type=eulerLabelsType),
          main=plotTitle, labels=list(font=4))
print(p1a)

subset2<-achr[c("scc1","coh1","scc1coh1")]
sigGenes<-lapply(subset2, "[[","wormbaseID")
uglyNames<-c("scc-1cs","coh-1cs","scc-1cscoh-1cs")
names(sigGenes)<-uglyNames
numGenes["Aup"]<-sum(sapply(sigGenes,length))
fit<-euler(sigGenes)
p1b<-plot(fit, quantities=list(type=eulerLabelsType),
          main=plotTitle, labels=list(font=4))
#print(p1b)

plotTitle<-list(label="chrX",fontsize=8)
xchr<-lapply(sigTables,function(x) x[x$chr=="chrX",])
subset1<-xchr[c("dpy26","kle2","scc1")]
sigGenes<-lapply(subset1,"[[","wormbaseID")
uglyNames<-c("dpy-26cs","kle-2cs","scc-1cs")
names(sigGenes)<-uglyNames
numGenes["Xup"]<-sum(sapply(sigGenes,length))
fit<-euler(sigGenes)
p1c<-plot(fit, quantities=list(type=eulerLabelsType),
          main=plotTitle, labels=list(font=4))
#print(p1c)

subset2<-xchr[c("scc1","coh1","scc1coh1")]
sigGenes<-lapply(subset2,"[[","wormbaseID")
uglyNames<-c("scc-1cs","coh-1cs","scc-1cscoh-1cs")
names(sigGenes)<-uglyNames
numGenes["Xup"]<-sum(sapply(sigGenes,length))
fit<-euler(sigGenes)
p1d<-plot(fit, quantities=list(type=eulerLabelsType),
          main=plotTitle, labels=list(font=4))
#print(p1d)


## downregulated genes -----
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(projectDir,"/",fileNamePrefix,contrastsOI[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))

  sigTables[[grp]]<-as.data.frame(
    getSignificantGenes(salmon, padj=padjVal, lfc= -lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="lt",
                        chr="all", nameChrCol="chr"))
}

plotTitle<-list(label="Autosomal",fontsize=8)
achr<-lapply(sigTables,function(x) x[x$chr!="chrX",])
subset1<-achr[c("dpy26","kle2","scc1")]
sigGenes<-lapply(subset1, "[[","wormbaseID")
uglyNames<-c("dpy-26cs","kle-2cs","scc-1cs")
names(sigGenes)<-uglyNames
numGenes["Adown"]<-sum(sapply(sigGenes,length))
fit<-euler(sigGenes)
p1e<-plot(fit, quantities=list(type=eulerLabelsType),
          main=plotTitle, labels=list(font=4))
#print(p1e)

subset2<-achr[c("scc1","coh1","scc1coh1")]
sigGenes<-lapply(subset2, "[[","wormbaseID")
uglyNames<-c("scc-1cs","coh-1cs","scc-1cscoh-1cs")
names(sigGenes)<-uglyNames
numGenes["Adown"]<-sum(sapply(sigGenes,length))
fit<-euler(sigGenes)
p1f<-plot(fit, quantities=list(type=eulerLabelsType),
          main=plotTitle, labels=list(font=4))
#print(p1f)


plotTitle<-list(label="chrX",fontsize=8)
xchr<-lapply(sigTables,function(x) x[x$chr=="chrX",])
subset1<-xchr[c("dpy26","kle2","scc1")]
sigGenes<-lapply(subset1, "[[","wormbaseID")
uglyNames<-c("dpy-26cs","kle-2cs","scc-1cs")
names(sigGenes)<-uglyNames
numGenes["Xdown"]<-sum(sapply(sigGenes,length))
fit<-euler(sigGenes)
p1g<-plot(fit, quantities=list(type=eulerLabelsType),
          main=plotTitle, labels=list(font=4))
#print(p1g)

subset2<-xchr[c("scc1","coh1","scc1coh1")]
sigGenes<-lapply(subset2, "[[","wormbaseID")
uglyNames<-c("scc-1cs","coh-1cs","scc-1cscoh-1cs")
names(sigGenes)<-uglyNames
numGenes["Xdown"]<-sum(sapply(sigGenes,length))
fit<-euler(sigGenes)
p1h<-plot(fit, quantities=list(type=eulerLabelsType),
          main=plotTitle, labels=list(font=4))
#print(p1h)

p1<-ggarrange(ggarrange(p1a,p1c,p1b,p1d,ncol=4, widths=c(0.8,1.3,1,0.5)),
              ggarrange(p1e,p1g,p1f,p1h,ncol=4,widths=c(0.8,0.3,1.2,0.4)),
              nrow=2,labels=c("Up regulated","Down regulated"),
              heights=c(1,0.9))

#p1


##################-
## Correlation------
##################-

geneTable<-NULL
for (grp in groupsOI){
  salmon<-as.data.frame(readRDS(paste0(projectDir,"/",fileNamePrefix,contrastsOI[[grp]],"_DESeq2_fullResults_p",padjVal,".rds")))
  colnames(salmon)[colnames(salmon)=="log2FoldChange"]<-paste0(grp,"_lfc")
  if(is.null(geneTable)){
    geneTable<-as.data.frame(salmon)[,c("wormbaseID","chr",paste0(grp,"_lfc"))]
  } else {
    geneTable<-full_join(geneTable,salmon[,c("wormbaseID","chr",paste0(grp,"_lfc"))], by=c("wormbaseID","chr"))
  }
}

combnTable<-combn(1:length(groupsOI),m=2)

lfcCols<-grep("_lfc$",names(geneTable))
#minScale<-min(geneTable[,lfcCols])*1.05
#maxScale<-max(geneTable[,lfcCols])*1.05
#minScale<- -2
#maxScale<- 2
minScale<-quantile(unlist(geneTable[,lfcCols[c(2,3)]]),0.0001)*1.05
maxScale<-quantile(unlist(geneTable[,lfcCols[c(2,3)]]),0.9999)*1.05

geneTable$XvA<-ifelse(geneTable$chr=="chrX","chrX","Autosomes")
#tmp<-geneTable

geneTable<-na.omit(geneTable)
allContrasts<-NULL
#contrastNm<-NULL
contrastNames<-c()
for (i in 1:ncol(combnTable)){
  grp1<-groupsOI[combnTable[1,i]]
  grp2<-groupsOI[combnTable[2,i]]
  prettyGrp1<-prettyNames[[combnTable[1,i]]]
  prettyGrp2<-prettyNames[[combnTable[2,i]]]
  df<-geneTable[,c(paste0(grp1,"_lfc"),paste0(grp2,"_lfc"),"XvA")]
  names(df)<-c("group1","group2","XvA")
  df$contrast<-deparse(substitute(x~v~y,list(x=prettyGrp1,y=prettyGrp2)))
  contrastNames<-c(contrastNames,df$contrast[1])
  #df$contrast<-paste(grp1,"v",grp2)
  #df$Rval<-Rval
  if(is.null(allContrasts)){
    allContrasts<-df
    #contrastNm<-paste(grp1,"v",grp2)
  } else {
    allContrasts<-rbind(allContrasts,df)
    #contrastNm<-c(contrastNm,paste(grp1,"v",grp2))
  }
}

allContrasts$XvA<-factor(allContrasts$XvA, levels=c("Autosomes","chrX"))
allContrasts$contrast<-factor(allContrasts$contrast,levels=contrastNames,labels=contrastNames)
p2<-ggplot(allContrasts,aes(x=group1,y=group2,col=XvA)) +
  facet_wrap(.~contrast,nrow=2,labeller=label_parsed)+
  geom_point(size=1,alpha=0.4) +
  coord_cartesian(xlim=c(minScale,maxScale),ylim=c(minScale,maxScale)) +
  geom_smooth(method=lm,se=F,fullrange=T, size=0.7,show.legend = T) +
  theme(legend.position="bottom", strip.text.x = element_text(size = 9),
        legend.title=element_blank(),axis.title.x=element_blank(),
        axis.title.y= element_blank()) +
  scale_color_manual(values=c("#111111","#FF1111"))+
  geom_hline(yintercept=0,lty=3,col="grey70",) +
  geom_vline(xintercept=0,lty=3,col="grey70") +
  ggpubr::stat_cor(aes(label = ..r.label..), method="pearson",
                   cor.coef.name = c("R"), output.type = "text",
                   show.legend=F,size=3,
                   label.x=minScale+0.1, label.y=c(maxScale-0.1,maxScale-1))
p2



##########################-
## Manual clicked loops------
##########################-

### new clicked loops
#clickedBatch="366"
#clickedBatch="382"

loopsOrAnchors<-"anchors"
ceTiles<-tileGenome(seqlengths(Celegans),tilewidth=10000,cut.last.tile.in.chrom = T)

if(loopsOrAnchors=="loops"){
  loops<-import(paste0(projectDir,"/otherData/Clicked_loops_",clickedBatch,"_merge.bedpe"),format="bedpe")
  grl<-zipup(loops)
  anchor1<-do.call(c,lapply(grl,"[",1))
  anchor2<-do.call(c,lapply(grl,"[",2))
  mcols(anchor1)<-mcols(loops)
  mcols(anchor2)<-mcols(loops)

  anchor1$loopNum<-paste0("loop",1:length(anchor1))
  anchor2$loopNum<-paste0("loop",1:length(anchor2))
  anchors<-sort(c(anchor1,anchor2))

  #separate anchors from inside tads
  tads_in<-GRanges(seqnames=seqnames(anchor1),IRanges(start=end(anchor1)+1,end=start(anchor2)-1))
  # tads_in<-reduce(tads_in)
  tads_in<-resize(tads_in,width=width(tads_in)-20000,fix="center")

  ol<-findOverlaps(ceTiles,tads_in)
  tenkbInTads<-ceTiles[unique(queryHits(ol))]

  anchors<-resize(anchors,width=10000,fix="center")
  reduce(anchors)

  ol<-findOverlaps(tenkbInTads,anchors)
  tenkbInTads<-tenkbInTads[-queryHits(ol)]
} else {
  anchors<-import(paste0(projectDir,"/otherData/382_X.eigs_cis.vecs_37peaks_p0.65_correct.bed"),format="bed")
  anchors<-resize(anchors,width=10000,fix="center")
  ol<-findOverlaps(ceTiles,anchors)
  tenkbInTads<-ceTiles[-queryHits(ol)]
  clickedBatch="382"
}



width(tenkbInTads)
width(anchors)
dataList<-list()
#grp=useContrasts[3]
statList<-list()
set.seed(34091857)
for (grp in useContrasts){
  salmon<-readRDS(file=paste0(paste0(projectDir,"/",fileNamePrefix,
                                     contrastsOI[[grp]],"_DESeq2_fullResults_p",padjVal,".rds")))

  salmon<-salmon[!is.na(salmon$chr),]
  salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)

  salmongr<-sort(salmongr)

  ol<-findOverlaps(salmongr,tenkbInTads,type="any",minoverlap=100L)
  insideTads<-salmongr[unique(queryHits(ol))]

  ol<-findOverlaps(salmongr,anchors,type="any",minoverlap=100L)
  atAnchors<-salmongr[unique(queryHits(ol))]

  ol<-findOverlaps(insideTads,atAnchors)
  insideTads<-insideTads[-queryHits(ol)]

  insideTads$Loops<-"Not anchor"
  atAnchors$Loops<-"Anchor"

  df<-data.frame(c(insideTads,atAnchors))
  df<-df%>%dplyr::group_by(seqnames,Loops)%>%dplyr::mutate(count=n())
  df$SMC<-grp

  dataList[[grp]]<-df
}



#pvalsDiff<-1-unlist(statList)
#names(pvalsDiff)<-prettyNamesAll

#aux_sdc3bg      dpy26  sdc3dpy26       kle2       scc1       coh1   scc1coh1
#0.8653     0.0001     0.0013     0.3441     0.3545     0.1655     0.2047

## focus on chrX loops
dataTbl<-do.call(rbind,dataList)
xchr<-dataTbl[dataTbl$seqnames=="chrX",]
cntTbl<-xchr %>% dplyr::group_by(SMC,Loops) %>% dplyr::summarise(count=n()) %>%
  filter(SMC=="dpy26")

xchr$SMC<-factor(xchr$SMC,levels=useContrasts,labels=prettyNames)
xchr$complexes<-complexes[xchr$SMC]
facetLabels<-xchr %>% dplyr::select(SMC,Loops,complexes) %>% distinct()

c1<-xchr %>% filter(SMC=="italic(\"dpy-26\"^cs)") %>% group_by(SMC, Loops) %>% summarize(count=n())
c2<-xchr %>% filter(SMC=="italic(\"dpy-26\"^cs)") %>% group_by(SMC, Loops) %>% filter(padj< 0.05, log2FoldChange>0.5) %>% summarize(count=n())

c2$count/c1$count

p3a<-ggplot(xchr,aes(x=Loops,y=log2FoldChange,fill=Loops))+
  geom_boxplot(notch=T,outlier.shape=NA,varwidth=T)+
  #geom_jitter()+
  facet_grid(col=vars(SMC),labeller=label_parsed) +coord_cartesian(ylim=c(-0.5,1.65))+
  ggtitle(paste0("LFC near loop anchors (",
                 cntTbl$count[cntTbl$Loops=="Anchor"]," genes) and not at anchors (",
                 cntTbl$count[cntTbl$Loops=="Not anchor"]," genes) in chrX")) +
  geom_hline(yintercept=0,linetype="dotted",color="grey20") +
  theme(axis.text.x=element_text(angle=45,hjust=1), axis.title.x=element_blank(),
        plot.title=element_text(size=9), legend.position = "none")+
  ylab("Log<sub>2</sub>FC")+
  ggsignif::geom_signif(test=t.test,comparisons = list(c("Anchor", "Not anchor")),
                        map_signif_level = F,tip_length=0.001,y_position=1.4,vjust=-0.1,
                        textsize=3,margin_top=0)+
  geom_text(data=facetLabels,aes(label=complexes),parse=T,x=1.5,y=1.65,size=3,hjust=0.5)

p3a


xchr$measure<-"Expression"
bm<- data.frame(xchr) %>% distinct(wormbaseID,baseMean,measure,Loops)
p3b<-ggplot(bm,aes(x=Loops,y=log2(baseMean),fill=Loops))+
  geom_boxplot(notch=T,outlier.shape=NA,varwidth=T) +
  facet_wrap(.~measure) + ggtitle(" ") +  coord_cartesian(ylim=c(-5,20)) +
  theme(legend.position = "none",axis.text.x=element_text(angle=45,hjust=1),
        plot.title=element_text(size=9),axis.title.x=element_blank()) +
  ylab("Log<sub>2</sub>(base mean counts)")+
  ggsignif::geom_signif(test=t.test,comparisons = list(c("Anchor", "Not anchor")),
                        map_signif_level = F,tip_length=0.001,vjust=-0.1,
                        textsize=3,margin_top=0.1)

p3b

p3<-ggarrange(p3b,p3a,ncol=2,widths=c(1.8,8.7))

p3



p<-ggarrange(p1,p2,p3,nrow=3,heights=c(1,1.2,1),labels=c("a ","b ","c "))
p<-annotate_figure(p, top = text_grob("Das et al., Figure S7", size = 12))
ggsave(paste0(finalFigDir,"/Fig_S7_RNAseqTEV_2.pdf"),p,device="pdf",width=21,height=29.7,
       unit="cm")
ggsave(paste0(finalFigDir,"/Fig_S7_RNAseqTEV_2.png"),p,device="png",width=21,height=29.7,
       unit="cm",bg="white")

