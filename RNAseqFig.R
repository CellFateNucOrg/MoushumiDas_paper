library(tidyverse)
library(ggplot2)
library(eulerr)

workDir<-getwd()
if(!dir.exists(paste0(workDir,"/plots"))) {
  dir.create(paste0(workDir,"/plots"))
}

strains<-c("366","382","775","784")
strain<-factor(strains,levels=strains)
SMC<-strain
levels(SMC)<-c("wt","dpy26cs","kle2cs","scc1cs")
controlGrp<-levels(SMC)[1] # control group
groupsOI<-levels(SMC)[-1]

source(paste0(workDir,"/functions.R"))


###################-
# panel A numbers of genes per chr-------
###################-
lfcVal=0.5
padjVal=0.05
filterPrefix<-"prefiltCyc2xChrAX"
fileNamePrefix=paste0("p",padjVal,"_lfc",lfcVal,"/",filterPrefix,"_")
outPath=paste0(workDir,"/",filterPrefix)

## up regulated count
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults_p",padjVal,".rds"))

  sigTables[[prettyGeneName(grp)]]<-as.data.frame(getSignificantGenes(salmon,
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
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults_p",padjVal,".rds"))

  sigTables[[prettyGeneName(grp)]]<-as.data.frame(getSignificantGenes(
    salmon, padj=padjVal,
    lfc= lfcVal,
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
bigMax<-sigPerChr$genes[maxIdx]*1.1
bigMin<-sigPerChr$genes[maxIdx]*0.9
smallMax<-max(abs(sigPerChr$genes[-maxIdx]))*1.1

step<-150
breaks<-seq(floor(-smallMax/step)*step,ceiling(bigMax/step)*step,step)
labels<-abs(seq(floor(-smallMax/step)*step,ceiling(bigMax/step)*step,step))


p1  <-  ggplot(sigPerChr, aes(x=chr, y=genes, group=SMC)) +
  facet_grid(direction~SMC,scales="free",space="free") + theme_bw() +
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
          position = position_dodge(0.9), size=4)




##################-
## venn diagram------
##################-
eulerLabelsType=c("counts")
## significantly changed genes
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults_p",padjVal,".rds"))
  sigTables[[prettyGeneName(grp)]]<-as.data.frame(getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                                                                      namePadjCol="padj",
                                                                      nameLfcCol="log2FoldChange",
                                                                      direction="both",
                                                                      chr="all", nameChrCol="chr"))
}

sigGenes<-lapply(sigTables,"[[","wormbaseID")
fit<-euler(sigGenes)
p2<-plot(fit, quantities=list(type=eulerLabelsType),
         main=list(label=paste0("All genes: |lfc|>", lfcVal, ", padj<",padjVal,"\n",
                                paste(lapply(row.names(fit$ellipses), function(x){
                                  paste(x, sum(fit$original.values[grep(x,names(fit$original.values))]))
                                }), collapse="  ")), fontsize=8))

p2
print(p2)


##################-
## Correlation------
##################-

geneTable<-NULL
for (grp in groupsOI){
  salmon<-as.data.frame(readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults_p",padjVal,".rds")))
  colnames(salmon)[colnames(salmon)=="log2FoldChange"]<-paste0(grp,"_lfc")
  if(is.null(geneTable)){
    geneTable<-as.data.frame(salmon)[,c("wormbaseID","chr",paste0(grp,"_lfc"))]
  } else {
    geneTable<-full_join(geneTable,salmon[,c("wormbaseID","chr",paste0(grp,"_lfc"))], by=c("wormbaseID","chr"))
  }
}

combnTable<-combn(1:length(groupsOI),m=2)

lfcCols<-grep("_lfc$",names(geneTable))
minScale<-min(geneTable[,lfcCols])*1.05
maxScale<-max(geneTable[,lfcCols])*1.05

for (i in 3:ncol(combnTable)){
  grp1<-groupsOI[combnTable[1,i]]
  grp2<-groupsOI[combnTable[2,i]]

  Rval<-round(cor(geneTable[,paste0(grp1,"_lfc")],
                  geneTable[,paste0(grp2,"_lfc")]),2)

  plot(geneTable[,paste0(grp1,"_lfc")],geneTable[,paste0(grp2,"_lfc")],pch=16,
       cex=0.5,col="#11111155",xlab=prettyGeneName(grp1),
       ylab=prettyGeneName(grp2), xlim=c(minScale,maxScale),
       ylim=c(minScale,maxScale))
  abline(v=0,h=0,col="grey60",lty=3)
  bestFitLine<-lm(geneTable[,paste0(grp2,"_lfc")]~geneTable[,paste0(grp1,"_lfc")])
  abline(bestFitLine,col="red")
  title(main=paste0("All genes ",prettyGeneName(grp1)," vs ",
                    prettyGeneName(grp2)," (R=",Rval,")"),
                sub=paste0(nrow(geneTable)," genes"))
  p3<-recordPlot()
}


geneTable<-NULL
for (grp in groupsOI){
  salmon<-as.data.frame(readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults_p",padjVal,".rds")))
  colnames(salmon)[colnames(salmon)=="log2FoldChange"]<-paste0(grp,"_lfc")
  if(is.null(geneTable)){
    geneTable<-as.data.frame(salmon)[,c("wormbaseID","chr",paste0(grp,"_lfc"))]
  } else {
    geneTable<-full_join(geneTable,salmon[,c("wormbaseID","chr",paste0(grp,"_lfc"))], by=c("wormbaseID","chr"))
  }
}

combnTable<-combn(1:length(groupsOI),m=2)

lfcCols<-grep("_lfc$",names(geneTable))
minScale<-min(geneTable[,lfcCols])*1.05
maxScale<-max(geneTable[,lfcCols])*1.05

geneTable$XvA<-ifelse(geneTable$chr=="chrX","chrX","Autosomes")
#tmp<-geneTable


geneTable<-na.omit(geneTable)
allContrasts<-NULL
for (i in 3:ncol(combnTable)){
  grp1<-groupsOI[combnTable[1,i]]
  grp2<-groupsOI[combnTable[2,i]]

  df<-geneTable[,c(paste0(grp1,"_lfc"),paste0(grp2,"_lfc"),"XvA")]
  names(df)<-c("group1","group2","XvA")
  df$contrast<-paste(prettyGeneName(grp1),"vs",prettyGeneName(grp2))
  df$Rval<-Rval
  if(is.null(allContrasts)){
    allContrasts<-df
  } else {
    allContrasts<-rbind(allContrasts,df)
  }
}

p3<-ggplot(allContrasts,aes(x=group1,y=group2)) +
  geom_point(col="#11111155",size=1) + xlab(prettyGeneName(grp1)) +
  ylab(prettyGeneName(grp2)) +
  xlim(c(minScale,maxScale)) + ylim(c(minScale,maxScale)) +
  geom_smooth(method=lm,se=F,fullrange=T, size=0.7) + theme_bw() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  geom_hline(yintercept=0,lty=3,col="grey70",) +
  geom_vline(xintercept=0,lty=3,col="grey70") +
  ggpubr::stat_cor(aes(label = ..r.label..), method="pearson",
                   cor.coef.name = c("R"), output.type = "text")

lay <- rbind(c(1,1,2),
             c(3,NA,NA))
p<-gridExtra::arrangeGrob(p1,p2,p3,layout_matrix=lay)

ggsave(paste0(workDir,"/plots/RNAseq.pdf"),p,device="pdf",width=10,height=7)


