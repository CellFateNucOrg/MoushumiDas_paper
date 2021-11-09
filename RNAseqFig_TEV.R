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

contrastsOI<-c("X.wt.wt.0mM_dpy26cs_vs_wt","X.wt.wt.0mM_kle2cs_vs_wt",
               "X.wt.wt.0mM_scc16cs_vs_wt","X.wt.wt.0mM_coh1cs_vs_wt",
               "X.wt.wt.0mM_scc1coh1cs_vs_wt")
useContrasts<-c("dpy26","kle2","scc1","coh1","scc1coh1")

prettyNames<-c(substitute(italic(x^cs),list(x="dpy-26")),
               substitute(italic(x^cs),list(x="kle-2")),
               substitute(italic(x^cs),list(x="scc-1")),
               substitute(italic(x^cs),list(x="coh-1")),
               substitute(italic(x^cs*y^cs),list(x="scc-1",y="coh-1")))
#plot(1:100,main=prettyNames[[5]])

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


###################-
# panel A numbers of genes per chr-------
###################-
lfcVal=0.5
padjVal=0.05
filterPrefix<-"filtCycChrAX"
fileNamePrefix=paste0("p",padjVal,"_lfc",lfcVal,"_",filterPrefix,"/",filterPrefix,"_")
outPath=paste0(workDir)

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
downPerChr<-downPerChr %>% gather(colnames(downPerChr)[1:6],key=chr, value=genes)
downPerChr$chr<-gsub("chr","",downPerChr$chr)
downPerChr$direction<-"down"

downPerChr$genes<-downPerChr$genes*(-1)

# combine up and down
sigPerChr<-rbind(upPerChr,downPerChr)
sigPerChr$direction<-factor(sigPerChr$direction,levels=c("up","down"))
sigPerChr$SMC<- factor(sigPerChr$SMC,levels=groupsOI,labels=prettyNames)


base.plot <- function(data) {
  p <- ggplot(data, aes(x=chr, y=genes, group=SMC)) + facet_grid(cols=vars(SMC))
  p <- p + theme_bw()
  p <- p + theme(legend.position="bottom", panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),)
  p <- p + geom_bar(stat="identity",position=position_dodge(),aes(fill=chr),
                    show.legend = FALSE)
  p <- p + scale_fill_grey(start=0.8, end=0.4)
  p <- p + xlab("Chromosome") + ylab("Number of genes")
  return(p)
}

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
             labeller=label_parsed) + theme_bw() +
  theme(legend.position="bottom", panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), strip.text = element_text(size = 12),
  #strip.text.x=element_text(face="italic"),
  axis.text=element_text(size=12), axis.title=element_text(size=12)) +
  geom_bar(stat="identity",position=position_dodge(),aes(fill=chr),
                  show.legend = FALSE) + scale_fill_grey(start=0.8, end=0.4) +
  xlab("Chromosome") + ylab("Number of genes") +
  scale_y_continuous(breaks=breaks,labels=labels,
                            expand = expansion(add = c(70, 70))) +
  geom_text(aes(label=abs(genes)),
                  vjust=ifelse(sign(sigPerChr$genes)>0,-0.2,1.2), color="black",
          position = position_dodge(0.9), size=3)




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

subset1<-sigTables[c("dpy26","kle2","scc1")]
sigGenes<-lapply(subset1,"[[","wormbaseID")
fit<-euler(sigGenes)
p2<-plot(fit, quantities=list(type=eulerLabelsType))#,
#         main=list(label=paste0("All genes: |lfc|>", lfcVal, ", padj<",padjVal,"\n",
#                                paste(lapply(row.names(fit$ellipses), function(x){
#                                  paste(x, sum(fit$original.values[grep(x,names(fit$original.values))]))
#                                }), collapse="  ")), fontsize=8))

p2
print(p2)

subset2<-sigTables[c("scc1","coh1","scc1coh1")]
sigGenes<-lapply(subset2,"[[","wormbaseID")
fit<-euler(sigGenes)
p2a<-plot(fit, quantities=list(type=eulerLabelsType))#,
#         main=list(label=paste0("All genes: |lfc|>", lfcVal, ", padj<",padjVal,"\n",
#                                paste(lapply(row.names(fit$ellipses), function(x){
#                                  paste(x, sum(fit$original.values[grep(x,names(fit$original.values))]))
#                                }), collapse="  ")), fontsize=8))

p2a
print(p2a)

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
             row_order = get_order(o1,1), column_order=1:5,
             show_row_names=F,row_title="Autosomes",column_names_rot = 45)
o1 = seriate(as.matrix(geneTable[geneTable$XvA=="chrX",lfcCols]), method = "PCA")
hm2<-Heatmap(as.matrix(geneTable[geneTable$XvA=="chrX",lfcCols]),name="NA",col=heatmapCol,
             row_order = get_order(o1,1),  column_order=1:5,
             show_row_names=F,row_title="chrX",column_names_rot = 45)
htlist=hm1 %v% hm2
ph1<-grid::grid.grabExpr(draw(htlist))
pdf(file=paste0(workDir,"/plots/hclustering_TEV.pdf"),width=5,height=8,
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
for(i in c(5,8,9,10)){
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
sigPerChr$SMC<- factor(sigPerChr$SMC,labels=prettyNames)

allContrasts$contrast<-factor(allContrasts$contrast,levels=contrastNames)
p3<-ggplot(allContrasts,aes(x=group1,y=group2)) +
  #facet_grid(cols=vars(contrast)) +
  facet_wrap(vars(contrast),nrow=2)+
  geom_point(col="#11111155",size=1) +
  #xlab(prettyGeneName(grp1)) +
  #ylab(prettyGeneName(grp2)) +
  xlim(c(minScale,maxScale)) + ylim(c(minScale,maxScale)) +
  geom_smooth(method=lm,se=F,fullrange=T, size=0.7) + theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  geom_hline(yintercept=0,lty=3,col="grey70",) +
  geom_vline(xintercept=0,lty=3,col="grey70") +
  ggpubr::stat_cor(aes(label = ..r.label..), method="pearson",
                   cor.coef.name = c("R"), output.type = "text")


#p<-ggarrange(p1,ggarrange(ph1,ggarrange(p2,p2a,nrow=2,labels=c("C.","D."),label.x=0.1),ncol=2),
#             p3,nrow=3,ncol=1,labels=c("A.","B.","E."),heights=c(4,4,2.5))

p<-ggarrange(ggarrange(p1,
                 ggarrange(p2,p2a,nrow=2,labels=c("B.","C."),label.x=0.1),ncol=2,widths=c(3,1)),
              ggarrange(ph1, p3,nrow=1,ncol=2,labels=c("D.","E.")),nrow=2,heights=c(4,4), labels=c("A."))


# p<-ggarrange(p1,
#             ggarrange(ggarrange(p2,p2a,nrow=2,labels=c("B.","C."), label.x=0.1),
#                       ph1,p3,ncol=3,labels=c("","D.","E."),widths=c(1.5,1.5,2)),
#              nrow=2,labels=c("A."),heights=c(4,4))

ggsave(paste0(workDir,"/plots/RNAseq_TEV.pdf"),p,device="pdf",width=10,height=10)


  #http://sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page





