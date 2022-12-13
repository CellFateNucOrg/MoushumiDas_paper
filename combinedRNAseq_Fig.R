library(tidyverse)
library(ggplot2)
library(eulerr)
library(ggpubr)
#library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(fastcluster)
library(seriation)
library(gridtext)
library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)
library(showtext)
library(ggtext)
#library(ragg)
#set_null_device("agg")
#par(family="Arial")
font_add(family="Arial","/System/Library/Fonts/Supplemental/Arial.ttf")
showtext_auto()
#theme_set(theme_get() + theme(text = element_text(family = 'Arial')))

workDir<-getwd()
if(!dir.exists(paste0(workDir,"/plots"))) {
  dir.create(paste0(workDir,"/plots"))
}

contrastsOI<-c("X.wt.wt.0mM_scc16cs_vs_wt","X.wt.wt.0mM_coh1cs_vs_wt",
               "X.wt.wt.0mM_scc1coh1cs_vs_wt","X.wt.wt.0mM_dpy26cs_vs_wt",
               "X.wt.wt.0mM_kle2cs_vs_wt")
useContrasts<-c("scc1","coh1","scc1coh1","dpy26","kle2")

prettyNames<-c(substitute(italic(x^cs),list(x="scc-1")),
               substitute(italic(x^cs),list(x="coh-1")),
               substitute(italic(x^cs*y^cs),list(x="scc-1",y="coh-1")),
               substitute(italic(x^cs),list(x="dpy-26")),
               substitute(italic(x^cs),list(x="kle-2")))
#plot(1:100,main=prettyNames[[5]])
# prettyNames<-c("*scc-1*<sup>cs</sup>",
#                "*coh-1*<sup>cs</sup>",
#                "*scc-1*<sup>cs</sup>*coh-1*<sup>cs</sup>",
#                "*dpy-26*<sup>cs</sup>",
#                "*kle-2*<sup>cs</sup>")

prettyNames1<-c("scc-1^cs","coh-1^cs","scc-1^(cs)coh-1^cs","dpy-26^cs","kle-2^cs")
#complexes<-c("cohesin^SCC-1","cohesin^COH-1","cohesins","condensinI/I^DC", "condensinII")
complexes<-c("cohesin<sup>SCC-1</sup>",
             "cohesin<sup>COH-1</sup>",
             "cohesins",
             "condensin I/I<sup>DC</sup>",
             "condensin II")
#complexes1<-c("cohesin^SCC-1","cohesin^COH-1","cohesins^(SCC-1)","condensinI/I^DC", "condensinII")
complexes1<-c("cohesin<sup>SCC-1</sup>",
              "cohesin<sup>COH-1</sup>","cohesins<sup>SCC-1</sup><br><sup>COH-1</sup>","condensin I/I<sup>DC</sup>", "condensin II")

names(complexes)<-useContrasts
names(complexes1)<-useContrasts
names(contrastsOI)<-useContrasts
# strains<-c("366","382","775","784","828","844")
# strain<-factor(strains,levels=strains)
# SMC<-strain
# levels(SMC)<-c("TEVonly",useContrasts)
#
# controlGrp<-levels(SMC)[1] # control group
# groupsOI<-levels(SMC)[-1]
groupsOI<-useContrasts


source(paste0(workDir,"/functions.R"))
# panel cdf of lfc
lfcVal=0.5
padjVal=0.05
filterPrefix<-"filtCycChrAX"
fileNamePrefix=paste0("p",padjVal,"_lfc",lfcVal,"_",filterPrefix,"/",filterPrefix,"_")
outPath=paste0(workDir)


###################-
# panel A numbers of genes per chr-------
###################-

## up regulated count
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/",fileNamePrefix,contrastsOI[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))

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
  salmon<-readRDS(paste0(outPath,"/",fileNamePrefix,contrastsOI[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))

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
sigPerChr$SMC<- factor(sigPerChr$SMC,levels=groupsOI,labels=prettyNames)


# base.plot <- function(data) {
#   p <- ggplot(data, aes(x=chr, y=genes, group=SMC)) + facet_grid(cols=vars(SMC))
#   p <- p + theme_bw()
#   p <- p + theme(legend.position="bottom", panel.grid.major = element_blank(),
#                  panel.grid.minor = element_blank())
#   p <- p + geom_bar(stat="identity",position=position_dodge(),aes(fill=chr),
#                     show.legend = FALSE)
#   p <- p + scale_fill_grey(start=0.8, end=0.4)
#   p <- p + xlab("Chromosome") + ylab("Number of genes")
#   return(p)
# }

# inspired by: https://www.j4s8.de/post/2018-01-15-broken-axis-with-ggplot2/
maxIdx<-which.max(abs(sigPerChr$genes))
bigMax<-sigPerChr$genes[maxIdx]*1.2
bigMin<-sigPerChr$genes[maxIdx]*0.8
smallMax<-max(abs(sigPerChr$genes[-maxIdx]))*1.1

step<-150
breaks<-seq(floor(-smallMax/step)*step,ceiling(bigMax/step)*step,step)
labels<-abs(seq(floor(-smallMax/step)*step,ceiling(bigMax/step)*step,step))

facetLabels<-sigPerChr %>% select(SMC,chr,complexes) %>% distinct()

p1  <-  ggplot(sigPerChr, aes(x=chr, y=genes, group=SMC)) +
  facet_grid(direction~SMC,scales="free",space="free",
             labeller=label_parsed) + theme_bw() +
  theme(legend.position="bottom", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), strip.text = element_text(size = 12),
        axis.text=element_text(size=12), axis.title=element_text(size=12)) +
  geom_bar(stat="identity",position=position_dodge(),aes(fill=chr),
           show.legend = FALSE) + scale_fill_grey(start=0.8, end=0.4) +
  xlab("Chromosome") + ylab("Number of genes") +
  scale_y_continuous(breaks=breaks,labels=labels,
                     expand = expansion(add = c(70, 70))) +
  geom_text(aes(label=abs(genes)),
            vjust=ifelse(sign(sigPerChr$genes)>0,-0.2,1.2), color="black",
            position = position_dodge(0.9), size=2.9) +
  geom_richtext(data=facetLabels,aes(label=complexes),x=3.5,y=870,size=4,hjust=0.5,
                fill = NA, label.color = NA, label.padding = grid::unit(rep(0, 4), "pt"),
                family="Arial")

p1

##########-
# heirarchical clustering of all LFC -------------
##########-

geneTable<-NULL
for (grp in useContrasts){
  salmon<-as.data.frame(readRDS(paste0(outPath,"/",fileNamePrefix,contrastsOI[[grp]],"_DESeq2_fullResults_p",padjVal,".rds")))
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

geneTable$XvA<-ifelse(geneTable$chr=="chrX","X","Autosomes")

#minQ<-quantile(as.matrix(geneTable[,lfcCols]),0.01)
#maxQ<-quantile(as.matrix(geneTable[,lfcCols]),0.99)
minQ<- -0.5
maxQ<- 0.5
ht_opt$fast_hclust = TRUE

heatmapCol<-circlize::colorRamp2(c(minQ,0,maxQ),c("blue","white","red"))
heatmapCol<-circlize::colorRamp2(c(minQ,0,maxQ),c("cyan","black","yellow"))

o1 = seriate(as.matrix(geneTable[geneTable$XvA=="Autosomes",lfcCols]), method = "PCA")
hm1<-Heatmap(as.matrix(geneTable[geneTable$XvA=="Autosomes",lfcCols]),
             col=heatmapCol,
             row_order = get_order(o1,1), column_order=1:length(useContrasts),
             show_row_names=F,row_title="Autosomes",column_names_rot = 90,
             show_heatmap_legend=T,
             heatmap_legend_param = list(title = gt_render("Log<sub>2</sub>FC")),
             use_raster=T,raster_quality=1,
             raster_device="CairoPNG")
o1 = seriate(as.matrix(geneTable[geneTable$XvA=="X",lfcCols]), method = "PCA")
hm2<-Heatmap(as.matrix(geneTable[geneTable$XvA=="X",lfcCols]),
             col=heatmapCol,
             column_labels=gt_render(complexes1,
                                     gp=gpar(fontsize=12)),
             row_order = get_order(o1,1),  column_order=1:5,
             show_row_names=F,row_title="X",column_names_rot = 70,
             show_heatmap_legend=F, use_raster=T,raster_quality=1,
             raster_device="CairoPNG")
htlist=hm1 %v% hm2
ph1<-grid::grid.grabExpr(draw(htlist,padding= unit(c(2, 5, 2, 5), "mm")))
pdf(file=paste0(workDir,"/plots/hclustering_TEV.pdf"),width=5,height=8,
    paper="a4")
draw(htlist,padding= unit(c(2, 0, 2, 0), "mm"))
dev.off()
# clustMethod="pearson"
# hm1<-Heatmap(as.matrix(geneTable[geneTable$XvA=="Autosomes",lfcCols]),name="Log2FC",col=heatmapCol,
#             clustering_distance_rows=clustMethod, column_order=1:5, column_title=clustMethod,
#             show_row_names=F,row_title="Autosomes",column_names_rot = 45)
# hm2<-Heatmap(as.matrix(geneTable[geneTable$XvA=="chrX",lfcCols]), name="NA",col=heatmapCol,
#              clustering_distance_rows=clustMethod, column_order=1:5, column_title=clustMethod,
#              show_row_names=F,row_title="chrX",column_names_rot = 45)
# htlist=hm1 %v% hm2
# draw(htlist)



##################-
# Correlation------
##################-

geneTable<-NULL
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/",fileNamePrefix,contrastsOI[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
  colnames(salmon)[colnames(salmon)=="log2FoldChange"]<-paste0(grp,"_lfc")
  if(is.null(geneTable)){
    geneTable<-as.data.frame(salmon)[,c("wormbaseID","chr",paste0(grp,"_lfc"))]
  } else {
    geneTable<-full_join(geneTable,data.frame(salmon[,c("wormbaseID","chr",paste0(grp,"_lfc"))]), by=c("wormbaseID","chr"))
  }
}

combnTable<-combn(1:length(groupsOI),m=2)

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
for(i in c(1,4)){
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
  geom_smooth(method=lm,se=F,fullrange=T, linewidth=0.7) + theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        legend.position="right", strip.text.x = element_text(size = 12)) +
  geom_hline(yintercept=0,lty=3,col="grey70",) +
  geom_vline(xintercept=0,lty=3,col="grey70") +
  ggpubr::stat_cor(aes(label = after_stat(r.label)), method="pearson",
                   cor.coef.name = c("R"), output.type = "text",
                   label.x=minScale+0.5, label.y=maxScale-0.5)+
  xlab(label=element_blank()) + ylab(label=element_blank())

p3


##########################-
#  loops------
##########################-

### new clicked loops

ceTiles<-tileGenome(seqlengths(Celegans),tilewidth=10000,cut.last.tile.in.chrom = T)

anchors<-import(paste0(outPath,"/otherData/382_X.eigs_cis.vecs_37peaks_p0.65_correct.bed"),format="bed")
anchors<-resize(anchors,width=10000,fix="center")
ol<-findOverlaps(ceTiles,anchors)
tenkbInTads<-ceTiles[-queryHits(ol)]
clickedBatch="382"


width(tenkbInTads)
saveRDS(tenkbInTads[seqnames(tenkbInTads)=="chrX",],paste0(outPath,"/txt/Tenkb_NotAnchor.rds" ))
width(anchors)
anchors<-anchors[-c(9,33),]
saveRDS(anchors,paste0(outPath,"/txt/Tenkb_Anchor.rds" ))
dataList<-list()
#grp=useContrasts[3]
statList<-list()
set.seed(34091857)
for (grp in useContrasts[4]){
  salmon<-readRDS(file=paste0(paste0(outPath,"/",fileNamePrefix,
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
  #  bsAvr<-list()
  #  for(i in 1:10000){
  #    idx<-sample(1:length(insideTads),length(atAnchors))
  #    bsAvr<-c(bsAvr,mean(insideTads$log2FoldChange[idx]))
  #  }
  #  statList[[grp]]<-sum(unlist(bsAvr)>mean(atAnchors$log2FoldChange))/10000
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
  dplyr::filter(SMC=="dpy26")

xchr$SMC<-factor(xchr$SMC,levels=useContrasts,labels=prettyNames)
xchr$complexes<-complexes[xchr$SMC]
facetLabels<-xchr %>% dplyr::select(SMC,Loops,complexes) %>% distinct()

c1<-xchr %>% dplyr::filter(SMC=="italic(\"dpy-26\"^cs)") %>% group_by(SMC, Loops) %>% dplyr::summarize(count=n())
c2<-xchr %>% dplyr::filter(SMC=="italic(\"dpy-26\"^cs)") %>% group_by(SMC, Loops) %>% dplyr::filter(padj< 0.05, log2FoldChange>0.5) %>% summarize(count=n())

c2$count/c1$count

p4<-ggplot(xchr,aes(x=Loops,y=log2FoldChange,fill=Loops))+
  geom_boxplot(notch=T,outlier.shape=NA,varwidth=T,show.legend=F)+
  facet_grid(col=vars(SMC),labeller=label_parsed) +coord_cartesian(ylim=c(-0.5,1.7))+
  geom_hline(yintercept=0,linetype="dotted",color="grey20") + theme_bw()+
  theme(axis.text.x=element_text(angle=70,hjust=1,size=12),
        axis.text.y=element_text(size=12),
        axis.title.y=element_markdown(),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        plot.title=element_text(size=12), legend.title = element_blank(),
        strip.text.x=element_text(size=12))+
  xlab(label=element_blank()) +
  ylab(label="Log<sub>2</sub>FC") +
  #scale_y_continuous(labels=c(expression(Log[2]FC)))+
  ggsignif::geom_signif(test=t.test,comparisons = list(c("Anchor", "Not anchor")),
                        map_signif_level = F,tip_length=0.001,y_position=1.45,vjust=-0.1,
                        textsize=3.5,margin_top=0)+
  geom_richtext(data=facetLabels,aes(label=complexes),x=1.5,y=1.7,size=4,hjust=0.5,
                fill = NA, label.color = NA, label.padding = grid::unit(rep(0, 4), "pt"),
                family="Arial")
  #geom_text(data=facetLabels,aes(label=complexes),parse=T,x=1.5,y=1.65,size=3.5,hjust=0.5)

p4





############################### Final arrangement ------
pnull<-NULL
p<-ggarrange(p1,
             ggarrange(ph1,p3,p4,nrow=1,ncol=3,labels=c("B ","C ","E "),
                       widths=c(1.2,1,0.8),font.label=list(family="Arial")),
             nrow=2,ncol=1,heights=c(1,1),labels=c("A "),font.label=list(family="Arial"))

#ggsave(paste0(workDir,"/plots/RNAseq_TEV.png"),p,device="png",width=10,height=10)
p<-annotate_figure(p, top = text_grob("Das et al., Figure 5", size = 14,family="Arial"))
ggsave(paste0(workDir,"/plots/RNAseq_TEV_Fig5.pdf"),p,device=cairo_pdf,width=21,height=21,
       units="cm",family="Arial")
#embedFonts(file=paste0(workDir,"/plots/RNAseq_TEV_Fig5.pdf"),outfile=paste0(workDir,"/plots/RNAseq_TEV_Fig5a.pdf"),
#           fontpaths="/System/Library/Fonts/Supplemental/Arial.ttf")
ggsave(paste0(workDir,"/plots/RNAseq_TEV_Fig5.png"),p,device=png,width=21,height=21,
       units="cm", bg="white")




