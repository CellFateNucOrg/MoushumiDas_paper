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

names(contrastsOI)<-useContrasts

groupsOI<-useContrasts


source(paste0(projectDir,"/functions.R"))

##############

contrastsOIall<-c("wt.TIR1.sdc3deg.X_1mM_vs_0mM","X.wt.wt.0mM_dpy26cs_vs_wt",
                  "X.TIR1.X.1mM_dpy26cssdc3deg_vs_wtwt")
useContrastsAll<-c("aux_sdc3bg","dpy26","sdc3dpy26")

prettyNamesAll<-c(substitute(italic(x^AID),list(x="sdc-3")),
                  substitute(italic(x^cs),list(x="dpy-26")),
                  substitute(italic(x^AID*y^cs),list(x="sdc-3",y="dpy-26")))

names(contrastsOIall)<-useContrastsAll




# panel cdf of lfc
lfcVal=0.5
padjVal=0.05
filterPrefix<-"filtCycChrAX"
fileNamePrefix=paste0("p",padjVal,"_lfc",lfcVal,"_",filterPrefix,"/",filterPrefix,"_")




##################-
## barplot of counts per chr------
##################-


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
#   p <- p + theme(legend.position="bottom", panel.grid.major = element_blank(),
#                  panel.grid.minor = element_blank(),)
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


p1  <-  ggplot(sigPerChr, aes(x=chr, y=genes, group=SMC)) +
  facet_grid(direction~SMC,scales="free",space="free",
             labeller=label_parsed) +
  theme(legend.position="bottom", strip.text = element_text(size = 9),
        axis.text=element_text(size=10), axis.title=element_text(size=9)) +
  geom_bar(stat="identity",position=position_dodge(),aes(fill=chr),
           show.legend = FALSE) + scale_fill_grey(start=0.8, end=0.4) +
  xlab("Chromosome") + ylab("Number of genes") +
  scale_y_continuous(breaks=breaks,labels=labels,
                     expand = expansion(add = c(70, 70))) +
  geom_text(aes(label=abs(genes)),
            vjust=ifelse(sign(sigPerChr$genes)>0,-0.2,1.2), color="black",
            position = position_dodge(0.9), size=3)


###################-
## volcano plots------
###################-

plotList<-list()
geneTable<-NULL
for (grp in groupsOI){
  salmon<-readRDS(paste0(projectDir,"/",fileNamePrefix,contrastsOI[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
  p<-plotVolcanoXvA(salmon,addLegend=F)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  plotList[[grp]]<-p+ggtitle(label=prettyNames[[which(groupsOI==grp)]])
}

leg <- get_legend(plotVolcanoXvA(salmon,addLegend=T))
plotList[["legend"]]<-as_ggplot(leg)
p2<-ggarrange(plotlist=plotList,ncol=2,nrow=2)

p2


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
minScale<- -2
maxScale<- 3
#minScale<-quantile(unlist(geneTable[,lfcCols[c(2,3)]]),0.001)*1.05
#maxScale<-quantile(unlist(geneTable[,lfcCols[c(2,3)]]),0.999)*1.05
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
p3<-ggplot(allContrasts,aes(x=group1,y=group2,col=XvA)) +
  facet_wrap(.~contrast,nrow=3,labeller=label_parsed)+
  geom_point(size=1,alpha=0.4) +
  coord_cartesian(xlim=c(minScale,maxScale), ylim=c(minScale,maxScale)) +
  geom_smooth(method=lm,se=F,fullrange=T, size=0.7,show.legend = T) +
  theme(legend.position="bottom", strip.text.x = element_text(size = 9),
        legend.title=element_blank(),axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scale_color_manual(values=c("#111111","#FF1111"))+
  geom_hline(yintercept=0,lty=3,col="grey70",) +
  geom_vline(xintercept=0,lty=3,col="grey70") +
  ggpubr::stat_cor(aes(label = ..r.label..), method="pearson",
                   cor.coef.name = c("R"), output.type = "text",
                   show.legend=F,size=3,
                   label.x=minScale+0.1, label.y=c(maxScale-0.1,maxScale-0.5))
p3

# p4a<-ggplot(allContrasts,aes(x=group1,y=group2)) +
#   facet_grid(cols=vars(XvA),rows=vars(contrast)) +
#   geom_point(col="#11111155",size=1) + xlab(NULL) + ylab(NULL) +
#   coord_cartesian(xlim=c(minScale,maxScale), ylim=c(minScale,maxScale)) +
#   geom_smooth(method=lm,se=F,fullrange=T, size=0.7) + theme_bw() +
#   theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
#         axis.text=element_text(size=12), axis.title=element_text(size=12),
#         strip.text = element_text(size = 12)) +
#   geom_hline(yintercept=0,lty=3,col="grey70",) +
#   geom_vline(xintercept=0,lty=3,col="grey70") +
#   ggpubr::stat_cor(aes(label = ..r.label..), method="pearson",
#                    cor.coef.name = c("R"), output.type = "text")




##################-
## venn diagrams------
##################-

testSigOLgeneSets<-function(DC,chrX=1783){
  set1<-DC[[1]]
  set2<-DC[[2]]
  nums<-c(length(intersect(set1,set2)),length(setdiff(set1,set2)),length(setdiff(set2,set1)))
  mm<-matrix(c(nums, chrX-length(union(set1,set2))),nrow=2)
  tt1<-fisher.test(mm)
  print(names(DC)[c(1,2)])
  print(mm)
  print(tt1)
  set1<-DC[[1]]
  set2<-DC[[3]]
  nums<-c(length(intersect(set1,set2)),length(setdiff(set1,set2)),length(setdiff(set2,set1)))
  mm<-matrix(c(nums, chrX-length(union(set1,set2))),nrow=2)
  tt2<-fisher.test(mm)
  print(names(DC)[c(1,3)])
  print(mm)
  print(tt2)
}


kramer<-as.data.frame(readRDS(file=paste0(projectDir,"/publicData/kramer2015_L3_gr.rds")))

localPadj=0.05
localLFC=0.5


filterList<-list()
filterList[["Cycling_Meeuse"]]<-read.delim(paste0(projectDir,"/publicData/oscillatingGenes.tsv"), header=T, stringsAsFactors=F)$wormbaseID #3739
filterList[["Cycling_Latorre"]]<-read.delim(paste0(projectDir,"/publicData/oscillatingGenes_latorre.tsv"))$wormbaseID #3235
toFilter<-unique(unlist(filterList))

idx<-kramer$wormbaseID %in% toFilter
kramer<-kramer[!idx,]

eulerLabelsType=c("counts")
uglyNames<-c("sdc-3AID","dpy-26cs","sdc-3AIDdpy-26cs")
plotList<-list()
grp=2
for (grp in 1:length(groupsOI)){
  salmon<-data.frame(readRDS(paste0(projectDir,"/",fileNamePrefix,contrastsOI[[groupsOI[grp]]],"_DESeq2_fullResults_p",padjVal,".rds")))

  ###############################-
  ## chrX-----
  ###############################-

  kramerDpy27<-getSignificantGenes(kramer, padj=localPadj, lfc=localLFC,
                                   namePadjCol="dpy27_RNAi_L3_padj",
                                   nameLfcCol="dpy27_RNAi_L3_log2_fold_change",
                                   direction="both",
                                   chr="chrX", nameChrCol="seqnames")

  kramerDpy21<-getSignificantGenes(kramer, padj=localPadj, lfc=localLFC,
                                   namePadjCol="dpy21_mutant_L3_padj",
                                   nameLfcCol="dpy21_mutant_L3_log2_fold_change",
                                   direction="both",
                                   chr="chrX", nameChrCol="seqnames")

  salmondc<-getSignificantGenes(salmon,padjVal,lfcVal,direction="both",chr="chrX")


  DC<-list(salmon=salmondc$wormbaseID, dpy27=kramerDpy27$wormbaseID,
           dpy21=kramerDpy21$wormbaseID)
  names(DC)<-c(uglyNames[grp], "dpy-27 (Kramer 2015)", "dpy-21 (Kramer 2015)")

  # fisher test
  testSigOLgeneSets(DC)

  fit<-euler(DC)

  plotList[[paste(groupsOI[grp],"_","chrX")]]<-plot(fit, quantities=list(type=eulerLabelsType),
                                                    main="chrX",labels=list(font=4))

  ###############################-
  ## Autosomes-----
  ###############################-

  kramerDpy27<-getSignificantGenes(kramer, padj=localPadj, lfc=localLFC,
                                   namePadjCol="dpy27_RNAi_L3_padj",
                                   nameLfcCol="dpy27_RNAi_L3_log2_fold_change",
                                   direction="both",
                                   chr="autosomes", nameChrCol="seqnames")

  kramerDpy21<-getSignificantGenes(kramer, padj=localPadj, lfc=localLFC,
                                   namePadjCol="dpy21_mutant_L3_padj",
                                   nameLfcCol="dpy21_mutant_L3_log2_fold_change",
                                   direction="both",
                                   chr="autosomes", nameChrCol="seqnames")

  salmondc<-getSignificantGenes(salmon,padjVal,lfcVal,direction="both",chr="autosomes")


  DC<-list(salmon=salmondc$wormbaseID, dpy27=kramerDpy27$wormbaseID,
           dpy21=kramerDpy21$wormbaseID)
  names(DC)<-c(uglyNames[grp], "dpy-27 (Kramer 2015)", "dpy-21 (Kramer 2015)")

  fit<-euler(DC)

  plotList[[paste(groupsOI[grp],"_","Autosomes")]]<-plot(fit,
                                                         quantities=list(type=eulerLabelsType),
                                                         main="Autosomes",labels=list(font=4))
}

p4<-ggarrange(plotlist=plotList[c(1,3,5,2,4,6)],ncol=3,nrow=2)
p4




############### Final assembly #########
p<-ggarrange(ggarrange(ggarrange(p1,p2,nrow=2,heights=c(1.2,2),labels=c("a ","b ")),p3,
                       ncol=2,widths=c(2,1),labels=c("","c ")),
             p4,nrow=2,heights=c(2,1),labels=c("","d "))
p<-annotate_figure(p, top = text_grob("Das et al., Figure S10", size = 12))
ggsave(paste0(finalFigDir,"/Fig_S10_RNAseq_sdc3_2.pdf"),p,device=cairo_pdf,width=21,height=29.7,units="cm")
ggsave(paste0(finalFigDir,"/Fig_S10_RNAseq_sdc3_2.png"),p,device=png,width=21,height=29.7,,unit="cm",bg="white")





#------------------------------------------------------------------------------




# contrastsOI<-c("wt.TIR1.X.1mM_sdc3deg_vs_wt","wt.TIR1.sdc3deg.X_1mM_vs_0mM",
#                "X.wt.wt.0mM_dpy26cs_vs_wt","X.TIR1.X.1mM_dpy26cssdc3deg_vs_wtwt",
#                "X.TIR1.sdc3deg.X_dpy26csaux_vs_wt0mM")
# useContrasts<-c("sdc3","aux_sdc3bg","dpy26","sdc3dpy26","aux_sdc3dpy26")
#
# prettyNames<-c(substitute(italic(+-x^AID),list(x="sdc-3")),
#                substitute(italic(x^AID+-IAA),list(x="sdc-3")),
#                substitute(italic(x^cs),list(x="dpy-26")),
#                substitute(italic(+-x^AID~y^cs),list(x="sdc-3",y="dpy-26")),
#                substitute(italic(x^AID+-y^cs~IAA),list(x="sdc-3",y="dpy-26")))

contrastsOI<-c("wt.TIR1.sdc3deg.X_1mM_vs_0mM",
               "X.wt.wt.0mM_dpy26cs_vs_wt","X.TIR1.X.1mM_dpy26cssdc3deg_vs_wtwt")
useContrasts<-c("aux_sdc3bg","dpy26","sdc3dpy26")

prettyNames<-c(substitute(italic(x^AID),list(x="sdc-3")),
               substitute(italic(x^cs),list(x="dpy-26")),
               substitute(italic(x^AID*y^cs),list(x="sdc-3",y="dpy-26")))

prettyNames1<-c("sdc-3^AID","dpy-26^cs","sdc-3^(AID)dpy-26^(cs)")


#plot(1:100,main=prettyNames[[5]])

names(contrastsOI)<-useContrasts
#strains<-c("366","382","775","784","828","844")
#strain<-factor(strains,levels=strains)
#SMC<-strain
#levels(SMC)<-c("TEVonly",useContrasts)

#controlGrp<-levels(SMC)[1] # control group
groupsOI<-useContrasts


source(paste0(projectDir,"/functions.R"))



lfcVal=0.5
padjVal=0.05
filterPrefix<-"filtCycChrAX"
fileNamePrefix=paste0("p",padjVal,"_lfc",lfcVal,"_",filterPrefix,"/",filterPrefix,"_")
projectDir=paste0(projectDir)

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


p11<-ggplot(allSig,aes(x=chr,y=log2FoldChange,fill=chr)) +
  geom_boxplot(outlier.shape=NA) +
  facet_grid(cols=vars(SMC),labeller=label_parsed) +
  theme(legend.position="none")+
  coord_cartesian(ylim=c(-0.5,0.5)) + scale_fill_grey(start=0.8, end=0.4)+
  geom_hline(yintercept=0,linetype="dashed",col="red")+
  xlab("Chromosome") + ylab(bquote(Log[2]~FC))




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

#minQ<-quantile(as.matrix(geneTable[,lfcCols]),0.01)
#maxQ<-quantile(as.matrix(geneTable[,lfcCols]),0.99)
minQ<- -0.5
maxQ<- 0.5
ht_opt$fast_hclust = TRUE

heatmapCol<-circlize::colorRamp2(c(minQ,0,maxQ),c("blue","white","red"))
heatmapCol<-circlize::colorRamp2(c(minQ,0,maxQ),c("cyan","black","yellow"))

o1 = seriate(as.matrix(geneTable[geneTable$XvA=="Autosomes",lfcCols]), method = "PCA")
hm1<-Heatmap(as.matrix(geneTable[geneTable$XvA=="Autosomes",lfcCols]),
             heatmap_legend_param = list(title = gt_render("Log<sub>2</sub>FC")),
             col=heatmapCol,
             row_order = get_order(o1,1), column_order=1:length(useContrasts),
             show_row_names=F,row_title="Autosomes",column_names_rot = 90,
             heatmap_width = unit(0.7, "npc"))
o1 = seriate(as.matrix(geneTable[geneTable$XvA=="chrX",lfcCols]), method = "PCA")
hm2<-Heatmap(as.matrix(geneTable[geneTable$XvA=="chrX",lfcCols]),name="NA",
             col=heatmapCol,
             column_labels=gt_render(prettyNames1,
                                     gp=gpar(fontface="italic",fontsize=9)),
             row_order = get_order(o1,1),  column_order=1:length(useContrasts),
             show_row_names=F,row_title="X",column_names_rot = 75,
             heatmap_width = unit(0.7, "npc"),show_heatmap_legend=F)
htlist=hm1 %v% hm2
ph11<-grid::grid.grabExpr(draw(htlist,padding= unit(c(2, 10, 2, 2), "mm")))

draw(htlist)

pdf(file=paste0(finalFigDir,"/hclustering_deg.pdf"),width=5,height=8,
    paper="a4")
draw(htlist)
dev.off()
# hm1<-Heatmap(as.matrix(geneTable[geneTable$XvA=="Autosomes",lfcCols]),name="Log2FC",col=heatmapCol,
#             clustering_distance_rows=clustMethod, column_order=1:5, column_title=clustMethod,
#             show_row_names=F,row_title="Autosomes",column_names_rot = 45)
# hm2<-Heatmap(as.matrix(geneTable[geneTable$XvA=="chrX",lfcCols]), name="NA",col=heatmapCol,
#              clustering_distance_rows=clustMethod, column_order=1:5, column_title=clustMethod,
#              show_row_names=F,row_title="chrX",column_names_rot = 45)
# htlist=hm1 %v% hm2
# draw(htlist)






############################### Final arrangement ------
pnull<-NULL
p<-ggarrange(p11,ph11,nrow=1,ncol=2,labels=c("F ","G "),
             widths=c(1,1.1))

#ggsave(paste0(workDir,"/plots/RNAseq_TEV.png"),p,device="png",width=10,height=10)
p<-annotate_figure(p, top = text_grob("Das et al., supplementary Figure 9", size = 12))
ggsave(paste0(finalFigDir,"/RNAseqDegwithHiC_SuplFig9.pdf"),p,device="pdf",width=7,height=7,units="cm")

ggsave(paste0(finalFigDir,"/RNAseqDegwithHiC_SuplFig9.png"),p,device="png",width=7,height=7,units="cm",
       bg="white")







