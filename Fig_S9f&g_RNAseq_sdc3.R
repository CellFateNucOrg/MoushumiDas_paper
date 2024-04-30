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
library(seriation)
library(gridtext)
library(rstatix)

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

contrastsOI<-c("wt.TIR1.sdc3deg.X_1mM_vs_0mM",
               "X.wt.wt.0mM_dpy26cs_vs_wt","X.TIR1.X.1mM_dpy26cssdc3deg_vs_wtwt")
useContrasts<-c("aux_sdc3bg","dpy26","sdc3dpy26")

prettyNames<-c(substitute(italic(x^AID),list(x="sdc-3")),
               substitute(italic(x^cs),list(x="dpy-26")),
               substitute(italic(x^AID*y^cs),list(x="sdc-3",y="dpy-26")))

prettyNames1<-c("sdc-3^AID","dpy-26^cs","sdc-3^(AID)dpy-26^(cs)")
prettyNames1<-c("sdc-3<sup>AID</sup>","dpy-26<sup>cs</sup>","sdc-3<sup>AID</sup>dpy-26<sup>cs</sup>")

names(contrastsOI)<-useContrasts

groupsOI<-useContrasts

source(paste0(projectDir,"/functions.R"))

lfcVal=0.5
padjVal=0.05
filterPrefix<-"filtCycChrAX"
fileNamePrefix=paste0("p",padjVal,"_lfc",lfcVal,"_",filterPrefix,"/",filterPrefix,"_")


###################-
# LFC per chr-------
###################-
## all genes
sigTables<-list()
for (grp in groupsOI[1]){
  salmon<-readRDS(paste0(projectDir,"/",fileNamePrefix,contrastsOI[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))

  sigTables[[grp]]<-as.data.frame(salmon)
  sigTables[[grp]]$SMC<-grp
}


allSig<-do.call(rbind,sigTables)
rownames(allSig)<-NULL
allSig$chr<-gsub("chr","",allSig$chr)
allSig$SMC<- factor(allSig$SMC,levels=groupsOI,labels=prettyNames)
allSig$XvA<-factor(ifelse(allSig$chr=="X","X","A"),levels=c("A","X"))

wilcoxt<-wilcox.test(allSig$log2FoldChange[allSig$chr=="X"],allSig$log2FoldChanges[allSig$chr!="X"])


p1<-ggplot(allSig,aes(x=chr,y=log2FoldChange,fill=chr)) +
  geom_boxplot(outlier.shape=NA) +
  facet_grid(cols=vars(SMC),labeller=label_parsed) +
  theme(legend.position="none")+
  coord_cartesian(ylim=c(-0.5,0.5)) + scale_fill_grey(start=0.8, end=0.4)+
  geom_hline(yintercept=0,linetype="dashed",col="red")+
  xlab("Chromosome") +ylab("Log<sub>2</sub>FC") #+
  #geom_signif(y_position=0.49,xmin=NA,xmax=6,annotation="p<2.2%*%10^-16",
  #            parse=T)

p1
table(allSig$chr)




##########-
# heirarchical clustering of all LFC -------------
##########-

geneTable<-NULL
for (grp in useContrasts[1:2]){
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
             name="Log2FC",col=heatmapCol,
             row_order = rev(get_order(o1,1)), column_order=1:length(useContrasts[1:2]),
             show_row_names=F,row_title="Autosomes",column_names_rot = 90,
             heatmap_width = unit(0.7, "npc"), heatmap_legend_param = list(
             title = expression(Log[2]~FC)))
o1 = seriate(as.matrix(geneTable[geneTable$XvA=="chrX",lfcCols]), method = "PCA")
hm2<-Heatmap(as.matrix(geneTable[geneTable$XvA=="chrX",lfcCols]),name="NA",
             col=heatmapCol,
             column_labels=gt_render(prettyNames1[1:2],
                                     gp=gpar(fontface="italic",fontsize=10)),
             row_order = rev(get_order(o1,1)),  column_order=1:length(useContrasts[1:2]),
             show_row_names=F,row_title="chrX",column_names_rot = 90,
             heatmap_width = unit(0.7, "npc"),show_heatmap_legend=F)
htlist=hm1 %v% hm2
ph1<-grid::grid.grabExpr(draw(htlist,padding= unit(c(2, 10, 2, 2), "mm")))

draw(htlist)

table(geneTable$XvA)

############### Final assembly #########

p<-ggarrange(p1,ph1,nrow=1,ncol=2, widths=c(3,3),labels=c("f ","g "))
#p<-annotate_figure(p, top = text_grob("Das et al., Figure S9", size = 14)) #face=bold
ggsave(paste0(finalFigDir,"/Fig_S9_RNAseq_sdc3.pdf"),p,device=cairo_pdf,width=9,height=12,
       unit="cm")



