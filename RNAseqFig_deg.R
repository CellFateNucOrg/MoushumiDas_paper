library(tidyverse)
library(ggplot2)
library(eulerr)
library(ggpubr)
#library(gplots)
library(ComplexHeatmap)
library(circlize)
library(fastcluster)
library(seriation)

workDir<-getwd()
if(!dir.exists(paste0(workDir,"/plots"))) {
  dir.create(paste0(workDir,"/plots"))
}

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

#plot(1:100,main=prettyNames[[5]])

names(contrastsOI)<-useContrasts
#strains<-c("366","382","775","784","828","844")
#strain<-factor(strains,levels=strains)
#SMC<-strain
#levels(SMC)<-c("TEVonly",useContrasts)

#controlGrp<-levels(SMC)[1] # control group
groupsOI<-useContrasts


source(paste0(workDir,"/functions.R"))



lfcVal=0.5
padjVal=0.05
filterPrefix<-"filtCycChrAX"
fileNamePrefix=paste0("p",padjVal,"_lfc",lfcVal,"_",filterPrefix,"/",filterPrefix,"_")
outPath=paste0(workDir)

###################-
# LFC per chr-------
###################-
## all genes
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


p1<-ggplot(allSig,aes(x=chr,y=log2FoldChange,fill=chr)) +
  geom_boxplot(outlier.shape=NA) +
  facet_grid(cols=vars(SMC),labeller=label_parsed) + theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylim(c(-0.5,1)) + scale_fill_grey(start=0.8, end=0.4)+
  geom_hline(yintercept=0,linetype="dotted",col="grey40")+
  xlab("Chromosome") +ylab("Log2(Fold Change)")


##################-
## venn diagram------
##################-
eulerLabelsType=c("counts")
## significantly changed genes
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/",fileNamePrefix,contrastsOI[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
  sigTables[[grp]]<-as.data.frame(getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                                                                      namePadjCol="padj",
                                                                      nameLfcCol="log2FoldChange",
                                                                      direction="both",
                                                                      chr="all", nameChrCol="chr"))
}

subset1<-sigTables[c("aux_sdc3bg","dpy26","sdc3dpy26")]
sigGenes<-lapply(subset1,"[[","wormbaseID")
fit<-euler(sigGenes)
p2<-plot(fit, quantities=list(type=eulerLabelsType))#,
#         main=list(label=paste0("All genes: |lfc|>", lfcVal, ", padj<",padjVal,"\n",
#                                paste(lapply(row.names(fit$ellipses), function(x){
#                                  paste(x, sum(fit$original.values[grep(x,names(fit$original.values))]))
#                                }), collapse="  ")), fontsize=8))

p2
print(p2)




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

geneTable$XvA<-ifelse(geneTable$chr=="chrX","chrX","Autosomes")

#minQ<-quantile(as.matrix(geneTable[,lfcCols]),0.01)
#maxQ<-quantile(as.matrix(geneTable[,lfcCols]),0.99)
minQ<- -0.5
maxQ<- 0.5
ht_opt$fast_hclust = TRUE

heatmapCol<-circlize::colorRamp2(c(minQ,0,maxQ),c("blue","white","red"))
heatmapCol<-circlize::colorRamp2(c(minQ,0,maxQ),c("cyan","black","yellow"))

o1 = seriate(as.matrix(geneTable[geneTable$XvA=="Autosomes",lfcCols]), method = "PCA")
hm1<-Heatmap(as.matrix(geneTable[geneTable$XvA=="Autosomes",lfcCols]),name="Log2FC",col=heatmapCol,
             row_order = get_order(o1,1), column_order=1:length(useContrasts),
             show_row_names=F,row_title="Autosomes",column_names_rot = 45)
o1 = seriate(as.matrix(geneTable[geneTable$XvA=="chrX",lfcCols]), method = "PCA")
hm2<-Heatmap(as.matrix(geneTable[geneTable$XvA=="chrX",lfcCols]),name="NA",
             col=heatmapCol, #column_labels=,
             row_order = get_order(o1,1),  column_order=1:length(useContrasts),
             show_row_names=F,row_title="chrX",column_names_rot = 45)
htlist=hm1 %v% hm2
ph1<-grid::grid.grabExpr(draw(htlist))
pdf(file=paste0(workDir,"/plots/hclustering_deg.pdf"),width=5,height=8,
    paper="a4")
draw(htlist)
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
## Correlation------
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
for(i in c(1:3)){
  grp1<-groupsOI[combnTable[1,i]]
  grp2<-groupsOI[combnTable[2,i]]
  df<-geneTable[,c(paste0(grp1,"_lfc"),paste0(grp2,"_lfc"),"XvA")]
  names(df)<-c("group1","group2","XvA")
  Rval<-cor(df[,1],df[,2])
  df$contrast<-paste(prettyGeneName(grp1),"v",prettyGeneName(grp2))
  contrastNames<-c(contrastNames,paste(prettyGeneName(grp1),"v",prettyGeneName(grp2)))
  df$Rval<-Rval
  if(is.null(allContrasts)){
    allContrasts<-df
  } else {
    allContrasts<-rbind(allContrasts,df)
  }
}
#sigPerChr$SMC<- factor(sigPerChr$SMC,levels=useContrasts,labels=prettyNames)

allContrasts$contrast<-factor(allContrasts$contrast,levels=contrastNames)
p3<-ggplot(allContrasts,aes(x=group1,y=group2,col=XvA)) +
  #facet_grid(cols=vars(contrast)) +
  facet_wrap(vars(contrast),nrow=2)+
  geom_point(size=1,alpha=0.4) +
  #xlab(prettyGeneName(grp1)) +
  #ylab(prettyGeneName(grp2)) +
  xlim(c(minScale,maxScale)) + ylim(c(minScale,maxScale)) +
  geom_smooth(method=lm,se=F,fullrange=T, size=0.7) + theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        legend.position="bottom") +
  scale_color_manual(values=c("#111111","#FF1111"))+
  geom_hline(yintercept=0,lty=3,col="grey70",) +
  geom_vline(xintercept=0,lty=3,col="grey70") +
  ggpubr::stat_cor(aes(label = ..r.label..), method="pearson",
                   cor.coef.name = c("R"), output.type = "text")


#p<-ggarrange(p1,ggarrange(ph1,ggarrange(p2,p2a,nrow=2,labels=c("C.","D."),label.x=0.1),ncol=2),
#             p3,nrow=3,ncol=1,labels=c("A.","B.","E."),heights=c(4,4,2.5))

p<-ggarrange(ggarrange(p1,p2,ncol=2,widths=c(3,1),labels=c("A.","B.")),
              ggarrange(ph1,p3,nrow=1,ncol=2,labels=c("C.","D."),widths=c(1,2)),
             nrow=2,heights=c(3,5))


# p<-ggarrange(p1,
#             ggarrange(ggarrange(p2,p2a,nrow=2,labels=c("B.","C."), label.x=0.1),
#                       ph1,p3,ncol=3,labels=c("","D.","E."),widths=c(1.5,1.5,2)),
#              nrow=2,labels=c("A."),heights=c(4,4))

ggsave(paste0(workDir,"/plots/RNAseq_deg.pdf"),p,device="pdf",width=8,height=10)


  #http://sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page





