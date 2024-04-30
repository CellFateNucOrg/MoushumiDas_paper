library(tidyverse)
library(ggplot2)
library(DESeq2)
library(ggpubr)
library(eulerr)
library(ComplexHeatmap)
library(circlize)
library(fastcluster)
library(seriation)
library(gridtext)
library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)

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
          panel.border = element_rect(colour = "black", fill=NA, size=0.8)
    )
)


contrastsOI<-c("X.wt.wt.0mM_dpy26cs_vs_wt","X.wt.wt.0mM_kle2cs_vs_wt",
               "X.wt.wt.0mM_scc16cs_vs_wt","X.wt.wt.0mM_coh1cs_vs_wt",
               "X.wt.wt.0mM_scc1coh1cs_vs_wt")
useContrasts<-c("dpy26","kle2","scc1","coh1","scc1coh1")

prettyNames<-c(substitute(italic(x^cs),list(x="dpy-26")),
               substitute(italic(x^cs),list(x="kle-2")),
               substitute(italic(x^cs),list(x="scc-1")),
               substitute(italic(x^cs),list(x="coh-1")),
               substitute(italic(x^cs*y^cs),list(x="scc-1",y="coh-1")))


complexes<-c("condensin~I/I^DC", "condensin~II", "cohesin^{SCC-1}","cohesin^{COH-1}","cohesins")
complexes1<-c("condensin I/I<sup>DC</sup>", "condensin II", "cohesin<sup>SCC-1</sup>",
              "cohesin<sup>COH-1</sup>","cohesin<sup>SCC-1<Br>COH-1</sup>")
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


###################-
# panel A numbers of genes per chr-------
###################-

## up regulated count
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(projectDir,"/",fileNamePrefix,contrastsOI[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))

  sigTables[[grp]]<-as.data.frame(getSignificantGenes(salmon,
                                                                      padj=padjVal, lfc=lfcVal,
                                                                      namePadjCol="padj",
                                                                      nameLfcCol="log2FoldChange",
                                                                      direction="gt",
                                                                      chr="all",
                                                                      nameChrCol="chr"))
}


upPerChr<-lapply(sigTables, "[", ,"chr")
upPerChr<-as.data.frame(do.call(rbind,lapply(upPerChr,table)))
upPerChr$SMC<-rownames(upPerChr)
upPerChr$complexes<-complexes[upPerChr$SMC]
upPerChr<-upPerChr %>% gather(colnames(upPerChr)[1:6],key=chr, value=genes)
upPerChr$chr<-gsub("chr","",upPerChr$chr)
upPerChr$direction<-"up"

## down regulated count
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(projectDir,"/",fileNamePrefix,contrastsOI[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))

  sigTables[[grp]]<-as.data.frame(getSignificantGenes(
    salmon, padj=padjVal,
    lfc= -lfcVal,
    namePadjCol="padj",
    nameLfcCol="log2FoldChange",
    direction="lt",
    chr="all",
    nameChrCol="chr"))
}

downPerChr<-lapply(sigTables, "[", ,"chr")
downPerChr<-as.data.frame(do.call(rbind,lapply(downPerChr,table)))
downPerChr$SMC<-rownames(downPerChr)
downPerChr$complexes<-complexes[downPerChr$SMC]
downPerChr<-downPerChr %>% gather(colnames(downPerChr)[1:6],key=chr, value=genes)
downPerChr$chr<-gsub("chr","",downPerChr$chr)
downPerChr$direction<-"down"

downPerChr$genes<-downPerChr$genes*(-1)

# combine up and down
sigPerChr<-rbind(upPerChr,downPerChr)
sigPerChr$direction<-factor(sigPerChr$direction,levels=c("up","down"))
sigPerChr$SMC<- factor(sigPerChr$SMC,levels=groupsOI[c(3,4,5,1,2)],labels=prettyNames[c(3,4,5,1,2)])


# inspired by: https://www.j4s8.de/post/2018-01-15-broken-axis-with-ggplot2/
maxIdx<-which.max(abs(sigPerChr$genes))
bigMax<-sigPerChr$genes[maxIdx]*1.3
bigMin<-sigPerChr$genes[maxIdx]*0.8
smallMax<-max(abs(sigPerChr$genes[-maxIdx]))*1.1

step<-150
breaks<-seq(floor(-smallMax/step)*step,ceiling(bigMax/step)*step,step)
labels<-abs(seq(floor(-smallMax/step)*step,ceiling(bigMax/step)*step,step))

facetLabels<-sigPerChr %>% select(SMC,chr,complexes) %>% distinct()
facetLabels$complexes[facetLabels$chr!="X"]<-"" # to prevent overplotting of labels

p1  <-  ggplot(sigPerChr, aes(x=chr, y=genes, group=SMC)) +
  facet_grid(direction~SMC,scales="free",space="free",
             labeller=label_parsed) +
  theme(legend.position="bottom",  strip.text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
  axis.text=element_text(size=12), axis.title=element_text(size=12)) +
  geom_bar(stat="identity",position=position_dodge(),aes(fill=chr),
                  show.legend = FALSE) + scale_fill_grey(start=0.8, end=0.4) +
  xlab("Chromosome") + ylab("Number of genes") +
  scale_y_continuous(breaks=breaks,labels=labels,
                            expand = expansion(add = c(70, 70))) +
  geom_text(aes(label=abs(genes)),
                  vjust=ifelse(sign(sigPerChr$genes)>0,-0.2,1.2), color="black",
          position = position_dodge(0.9), size=2.7) +
  geom_text(data=facetLabels,aes(label=complexes),parse=T,x=3.5,y=870,size=3,hjust=0.5)

p1


##########-
# heirarchical clustering of all LFC -------------
##########-

geneTable<-NULL
for (grp in useContrasts){
  salmon<-as.data.frame(readRDS(paste0(projectDir,"/",fileNamePrefix,contrastsOI[[grp]],"_DESeq2_fullResults_p",padjVal,".rds")))
  colnames(salmon)[colnames(salmon)=="log2FoldChange"]<-paste0(grp,"_lfc")
  colnames(salmon)[colnames(salmon)=="padj"]<-paste0(grp,"_padj")
  if(is.null(geneTable)){
    geneTable<-as.data.frame(salmon)[,c("wormbaseID","chr",paste0(grp,"_lfc"),paste0(grp,"_padj"))]
  } else {
    geneTable<-full_join(geneTable,salmon[,c("wormbaseID","chr",paste0(grp,"_lfc"),paste0(grp,"_padj"))], by=c("wormbaseID","chr"))
  }
}

padjCols<-grep("_padj$",colnames(geneTable))
notSigAny<-rowSums(is.na(geneTable[,padjCols]))==length(padjCols)
geneTable<-geneTable[!notSigAny,]
lfcCols<-grep("_lfc$",colnames(geneTable))
colnames(geneTable)<-gsub("_lfc$","",colnames(geneTable))

geneTable$XvA<-ifelse(geneTable$chr=="chrX","chrX","Autosomes")


minQ<- -0.5
maxQ<- 0.5
ht_opt$fast_hclust = TRUE

heatmapCol<-circlize::colorRamp2(c(minQ,0,maxQ),c("blue","white","red"))
heatmapCol<-circlize::colorRamp2(c(minQ,0,maxQ),c("cyan","black","yellow"))



o1 = seriate(as.matrix(geneTable[geneTable$XvA=="Autosomes",lfcCols]), method = "PCA")
hm1<-Heatmap(as.matrix(geneTable[geneTable$XvA=="Autosomes",lfcCols]),
             name="Log2FC", col=heatmapCol,
             row_order = get_order(o1,1), column_order=c(3,4,5,1,2),
             show_row_names=F,row_title="Autosomes",column_names_rot = 90,
             use_raster=T,raster_quality=1,raster_device="CairoPNG",
             heatmap_legend_param = list(
               title = expression(Log[2]~FC)))
o1 = seriate(as.matrix(geneTable[geneTable$XvA=="chrX",lfcCols]), method = "PCA")
hm2<-Heatmap(as.matrix(geneTable[geneTable$XvA=="chrX",lfcCols]), name="NA",
             col=heatmapCol,
             column_labels=gt_render(complexes1,
                                     gp=gpar(fontface="italic", fontsize=12)),
             row_order = get_order(o1,1),  column_order=c(3,4,5,1,2),
             show_row_names=F,row_title="X",column_names_rot = 80,
             show_heatmap_legend=F, use_raster=T,raster_quality=1,
             raster_device="CairoPNG")
htlist=hm1 %v% hm2
ph1<-grid::grid.grabExpr(draw(htlist,padding= unit(c(2, 5, 2, 5), "mm")))
grid.draw(ph1)

table(geneTable$XvA)



##################-
# Correlation------
##################-

geneTable<-NULL
for (grp in groupsOI){
  salmon<-readRDS(paste0(projectDir,"/",fileNamePrefix,contrastsOI[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
  colnames(salmon)[colnames(salmon)=="log2FoldChange"]<-paste0(grp,"_lfc")
  if(is.null(geneTable)){
    geneTable<-as.data.frame(salmon)[,c("wormbaseID","chr",paste0(grp,"_lfc"))]
  } else {
    geneTable<-full_join(geneTable,data.frame(salmon[,c("wormbaseID","chr",paste0(grp,"_lfc"))]), by=c("wormbaseID","chr"))
  }
}

combnTable<-combn(c(3,4,2),m=2)

lfcCols<-grep("_lfc$",names(geneTable))
minScale<-min(geneTable[,lfcCols])*1.05
maxScale<-max(geneTable[,lfcCols])*1.05
minScale<-quantile(unlist(geneTable[,lfcCols[c(2,3)]]),0.0001)*1.05
maxScale<-quantile(unlist(geneTable[,lfcCols[c(2,3)]]),0.9999)*1.05


geneTable$XvA<-ifelse(geneTable$chr=="chrX","chrX","Autosomes")
#tmp<-geneTable


geneTable<-na.omit(geneTable)
allContrasts<-NULL
contrastNames<-c()
#for (i in 3:ncol(combnTable)){
for(i in c(1,2)){
  grp1<-groupsOI[combnTable[1,i]]
  grp2<-groupsOI[combnTable[2,i]]
  prettyGrp1<-prettyNames[[combnTable[1,i]]]
  prettyGrp2<-prettyNames[[combnTable[2,i]]]
  df<-geneTable[,c(paste0(grp1,"_lfc"),paste0(grp2,"_lfc"),"XvA")]
  names(df)<-c("group1","group2","XvA")
  df$contrast<-deparse(substitute(x~v~y,list(x=prettyGrp1,y=prettyGrp2)))
  contrastNames<-c(contrastNames,df$contrast[1])
  if(is.null(allContrasts)){
    allContrasts<-df
  } else {
    allContrasts<-rbind(allContrasts,df)
  }
}

allContrasts$contrast<-factor(allContrasts$contrast,levels=contrastNames,labels=contrastNames)
p3<-ggplot(allContrasts,aes(x=group1,y=group2)) +
  facet_wrap(vars(contrast),nrow=2,labeller=label_parsed)+
  geom_point(col="#11111155",size=1) +
  coord_cartesian(xlim=c(minScale,maxScale),ylim=c(minScale,maxScale)) +
  geom_smooth(method=lm,se=F,fullrange=T, size=0.7) +
  theme(legend.position="right", strip.text.x = element_text(size = 11),
        axis.title.x=element_blank()) +
  geom_hline(yintercept=0,lty=3,col="grey70",) +
  geom_vline(xintercept=0,lty=3,col="grey70") +
  ggpubr::stat_cor(aes(label = ..r.label..), method="pearson",
                   cor.coef.name = c("R"), output.type = "text",
                   label.x=minScale+0.5, label.y=maxScale-0.5)+
  ylab("Log<sub>2</sub>FC")

p3

table(allContrasts$contrast)

##########################-
# Manual clicked loops------
##########################-



ceTiles<-tileGenome(seqlengths(Celegans),tilewidth=10000,cut.last.tile.in.chrom = T)

anchors<-import(paste0(otherDataDir,"/382_X.eigs_cis.vecs_37peaks_p0.65_correct.bed"),format="bed")
anchors<-resize(anchors,width=10000,fix="center")
ol<-findOverlaps(ceTiles,anchors)
tenkbInTads<-ceTiles[-queryHits(ol)]
clickedBatch="382"


width(tenkbInTads)
width(anchors)
dataList<-list()
statList<-list()
set.seed(34091857)
for (grp in useContrasts[1]){
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



## focus on chrX loops
dataTbl<-do.call(rbind,dataList)
xchr<-dataTbl[dataTbl$seqnames=="chrX",]
cntTbl<-xchr %>% dplyr::group_by(SMC,Loops) %>% dplyr::summarise(count=n()) %>%
  filter(SMC=="dpy26")

xchr$SMC<-factor(xchr$SMC,levels=useContrasts,labels=prettyNames)
xchr$complexes<-complexes[xchr$SMC]
facetLabels<-xchr %>% dplyr::select(SMC,Loops,complexes) %>% distinct()
facetLabels$complexes[facetLabels$Loops=="Anchor"]<-""


c1<-xchr %>% filter(SMC=="italic(\"dpy-26\"^cs)") %>% group_by(SMC, Loops) %>% summarize(count=n())
c2<-xchr %>% filter(SMC=="italic(\"dpy-26\"^cs)") %>% group_by(SMC, Loops) %>% filter(padj< 0.05, log2FoldChange>0.5) %>% summarize(count=n())

c2$count/c1$count

p4<-ggplot(xchr,aes(x=Loops,y=log2FoldChange,fill=Loops))+
  geom_boxplot(notch=T,outlier.shape=NA,varwidth=T,show.legend=F)+
  facet_grid(col=vars(SMC),labeller=label_parsed) +coord_cartesian(ylim=c(-0.5,1.7))+
  geom_hline(yintercept=0,linetype="dotted",color="grey20") +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=12),
        axis.text.y=element_text(size=10),
        axis.title.x = element_blank(),
        plot.title=element_text(size=10), legend.title = element_blank(),
        strip.text.x=element_text(size=12))+
  ylab("Log<sub>2</sub>FC")+
  ggsignif::geom_signif(test=wilcox.test,comparisons = list(c("Anchor", "Not anchor")),
                        map_signif_level = F,tip_length=0.001,y_position=1.45,vjust=-0.1,
                        textsize=3,margin_top=0)+
  geom_text(data=facetLabels,aes(label=complexes),parse=T,x=1.5,y=1.65,size=3.5,hjust=0.5)

p4

table(xchr$Loops)


############################### Final arrangement ------
p<-ggarrange(p1,
             ggarrange(ph1, p3,p4,nrow=1,ncol=3,labels=c("b ","c ","d "),
                       widths=c(1.3,1.2,0.9)),nrow=2,heights=c(4,4),
             labels=c("a "))

p
p<-annotate_figure(p, top = text_grob("Das et al., Figure 5", size = 14))
ggsave(paste0(finalFigDir,"/Fig_5_RNAseq_TEV.pdf"),p,device="pdf",width=21,height=21,units="cm")

ggsave(paste0(finalFigDir,"/Fig_5_RNAseq_TEV.png"),p,device="png",width=21,height=21,units="cm",
       bg="white")
  #http://sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page





