library(tidyverse)
library(ggplot2)
library("ggpmisc")
library(eulerr)
library(lattice)

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



# panel cdf of lfc
lfcVal=0
padjVal=0.05
filterPrefix<-"prefiltCyc2xChrAX"
fileNamePrefix=paste0("p",padjVal,"_lfc",lfcVal,"/",filterPrefix,"_")
outPath=paste0(workDir,"/",filterPrefix)


########################-
## ECDF of data -----
########################-

# using automatic filtering threshold
sigTables<-list()
localPadj=padjVal
localLFC=0
for (grp in groupsOI){
  print(grp)
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults_p",padjVal,".rds"))
  salmon<-salmon[!is.na(salmon$padj),]
  #nrow(filterResults(salmon,padj=0.05,lfc=0.5,direction="lt",chr="autosomes"))
  print(paste0(nrow(salmon)," genes before filtering"))
  print(paste0(sum(is.na(salmon$log2FoldChange))," have log2FoldChange that is NA"))
  #salmon$expressed<-sum(salmon$baseMean>10)
  sigTables[[prettyGeneName(grp)]]<-as.data.frame(salmon) #[salmon$baseMean>10,]
  print(paste0(nrow(sigTables[[prettyGeneName(grp)]])," genes after automatic threshold filter"))
}

SMC<-rep(names(sigTables),lapply(sigTables,nrow))
sig<-do.call(rbind,sigTables)
sig$SMC<-SMC
table(sig$SMC)
sig$XvA<-"Autosomes"
sig$XvA[sig$chr=="chrX"]<-"chrX"
#sig$XvA<-factor(sig$XvA)
table(sig$XvA)
sig$upVdown[sig$log2FoldChange<0]<-"down"
sig$upVdown[sig$log2FoldChange>0]<-"up"
sig$upVdown<-factor(sig$upVdown,levels=c("down","up"))
table(sig$upVdown)
row.names(sig)<-NULL
SMC<-NULL

lapply(sigTables,function(x){sum(x$padj<0.05)})

options(tibble.width=Inf)
dd1<-sig %>% dplyr::filter(padj<localPadj) %>%
  dplyr::group_by(SMC,upVdown,XvA) %>%
  dplyr::mutate(ecd=ecdf(abs(log2FoldChange))(abs(log2FoldChange)))

dd1$upVdown<-relevel(dd1$upVdown,"up")

p<-ggplot(dd1, aes(x=abs(log2FoldChange),y=ecd,color=SMC,linetype=XvA)) +
  geom_line(size=0.9)+ facet_wrap(vars(upVdown),nrow=2)+
  theme_classic() + xlim(c(0,1.5)) +
  xlab("Absolute log2 fold change")+ylab("Fraction of significant genes")
#p

#stat_ecdf(aes(colour=SMC,linetype=XvA),alpha=0.7)
p1<-p+geom_vline(aes(xintercept = 0.5), color="grey") +
  annotate("text",label="0.5",size=4, x=0.5, y=1,hjust=-0.05,color="grey") +
  geom_vline(aes(xintercept = 0.15), color="grey") +
  annotate("text",label="0.15",size=4, x=0.15, y=1,hjust=-0.05,color="grey") +
  theme(strip.text = element_text(size = 12), axis.text=element_text(size=12),
        axis.title=element_text(size=12))
#p1

# to get fraction between 0.15 and 0.5 lfc
dd2<-sig %>% dplyr::filter(padj<localPadj) %>%
  dplyr::group_by(SMC,upVdown,XvA) %>%
  dplyr::summarise(ecd15=ecdf(abs(log2FoldChange))(c(0.15)),
                   ecd50=ecdf(abs(log2FoldChange))(c(0.5)),
                   q2=100*(ecd50-ecd15))
dd2[order(dd2$q2),]

# to get total expressed genes
na.omit(sig) %>% group_by(SMC,XvA) %>% summarise(count=n())

# to add tables of number of signifcant genes to plot
ttu<-with(sig[sig$padj<localPadj & sig$upVdown=="up",],table(XvA,SMC))
ttu<-pivot_wider(as.data.frame(ttu),names_from=SMC,values_from=Freq)
names(ttu)<-gsub("XvA","",names(ttu))

ttd<-with(sig[sig$padj<localPadj & sig$upVdown=="down",],table(XvA,SMC))
ttd<-pivot_wider(as.data.frame(ttd),names_from=SMC,values_from=Freq)
names(ttd)<-gsub("XvA","",names(ttd))

tbs<-list("up"=ttu,
          "down"=ttd)

df <- tibble(x = rep(1.5, 2),
             y = rep(0, 2),
             upVdown= c("up","down"),
             tbl = tbs)

p1<-p1 + geom_table(data = df, aes(x = x, y = y,label = tbl),
                table.theme = ttheme_gtlight)


##################-
## venn diagram------
##################-
lfcVal=0.5
padjVal=0.05
eulerLabelsType=c("counts")
numGenes<-list()
## upregulated genes -----
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults_p",padjVal,".rds"))
  sigTables[[prettyGeneName(grp)]]<-as.data.frame(
    getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="gt",
                        chr="all", nameChrCol="chr"))
}

plotTitle<-list(label=paste("Autosomal up: lfc>", lfcVal, ", padj<",padjVal,"\n",
                            paste(lapply(row.names(fit$ellipses), function(x){
                              paste(x, sum(fit$original.values[grep(x,names(fit$original.values))]))
                            }), collapse="  ")), fontsize=8)
plotTitle<-list(label="Autosomal up",fontsize=8)
achr<-lapply(sigTables,function(x) x[x$chr!="chrX",])
sigGenes<-lapply(achr, "[[","wormbaseID")
numGenes["Aup"]<-sum(sapply(sigGenes,length))
fit<-euler(sigGenes)
p2a<-plot(fit, quantities=list(type=eulerLabelsType),
          main=plotTitle)

#print(p2a)


plotTitle<-list(label=paste("chrX up: lfc>", lfcVal, ", padj<",padjVal,"\n",
                            paste(lapply(row.names(fit$ellipses), function(x){
                              paste(x, sum(fit$original.values[grep(x,names(fit$original.values))]))
                            }), collapse="  ")), fontsize=8)
plotTitle<-list(label="chrX up",fontsize=8)
xchr<-lapply(sigTables,function(x) x[x$chr=="chrX",])
sigGenes<-lapply(xchr, "[[","wormbaseID")
numGenes["Xup"]<-sum(sapply(sigGenes,length))
fit<-euler(sigGenes)
p2b<-plot(fit, quantities=list(type=eulerLabelsType),
         main=plotTitle)
#print(p2b)


## downregulated genes -----
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults_p",padjVal,".rds"))

  sigTables[[prettyGeneName(grp)]]<-as.data.frame(
    getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="lt",
                        chr="all", nameChrCol="chr"))
}

plotTitle<-list(label=paste("Autosomal down: lfc< -", lfcVal, ", padj<",padjVal,"\n",
                            paste(lapply(row.names(fit$ellipses), function(x){
                              paste(x, sum(fit$original.values[grep(x,names(fit$original.values))]))
                            }), collapse="  ")), fontsize=8)
plotTitle<-list(label="Autosomal down",fontsize=8)
achr<-lapply(sigTables,function(x) x[x$chr!="chrX",])
sigGenes<-lapply(achr, "[[","wormbaseID")
numGenes["Adown"]<-sum(sapply(sigGenes,length))
fit<-euler(sigGenes)
p2c<-plot(fit, quantities=list(type=eulerLabelsType),
          main=plotTitle)

#print(p2c)


plotTitle<-list(label=paste("chrX down: lfc< -", lfcVal, ", padj<",padjVal,"\n",
                            paste(lapply(row.names(fit$ellipses), function(x){
                              paste(x, sum(fit$original.values[grep(x,names(fit$original.values))]))
                            }), collapse="  ")), fontsize=8)
plotTitle<-list(label="chrX down",fontsize=8)
xchr<-lapply(sigTables,function(x) x[x$chr=="chrX",])
sigGenes<-lapply(xchr, "[[","wormbaseID")
numGenes["Xdown"]<-sum(sapply(sigGenes,length))
fit<-euler(sigGenes)
p2d<-plot(fit, quantities=list(type=eulerLabelsType),
         main=plotTitle,respect=T)

#print(p2d)

#numGenes<-unlist(numGenes)
# p2<-gridExtra::grid.arrange(gridExtra::arrangeGrob(p2a,p2b,
#                               heights=numGenes[1:2]/sum(numGenes),
#                               nrow=2),
#                             gridExtra::arrangeGrob(p2c,p2d,
#                               heights=(numGenes[3:4]/sum(numGenes)+c(0,0.05)),
#                               nrow=2),ncol=2,respect=T)

lay <- rbind(c(1,1,1,2,2,2),
             c(1,1,1,2,2,2),
             c(1,1,1,2,2,2),
             c(3,3,3,NA,NA,NA),
             c(3,3,3,NA,4,NA),
             c(3,3,3,NA,NA,NA))
p2<-gridExtra::arrangeGrob(p2a,p2b,p2c,p2d,layout_matrix=lay)
#grid::grid.draw(p2)

##cowplot::plot_grid(p2a,p2b,p2c,p2d,align="hv",nrow=2,ncol=2)

##########-
# heirarchical clustering of all LFC -------------
##########-

geneTable<-NULL
for (grp in useContrasts){
  salmon<-as.data.frame(readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults_p",padjVal,".rds")))
  colnames(salmon)[colnames(salmon)=="log2FoldChange"]<-paste0(grp,"_lfc")
  if(is.null(geneTable)){
    geneTable<-as.data.frame(salmon)[,c("wormbaseID","chr",paste0(grp,"_lfc"))]
  } else {
    geneTable<-full_join(geneTable,salmon[,c("wormbaseID","chr",paste0(grp,"_lfc"))], by=c("wormbaseID","chr"))
  }
}


lfcCols<-grep("lfc",colnames(geneTable))
colnames(geneTable)<-gsub("_lfc$","",colnames(geneTable))

#annClrs<-brewer.pal(length(useContrasts), name="Dark2")
#names(annClrs)<-useContrasts
#colClrs<-factor(colnames(geneTable)[lfcCols])
#levels(colClrs)<-annClrs[levels(colClrs)]

heatmapCol<-colorRampPalette(c("cyan","black","pink"))(20)
# Perform the hierarchical clustering with
# A distance based on Pearson-correlation coefficient
# and average linkage clustering as agglomeration criteria
heatmap.2(as.matrix(geneTable[,lfcCols]),
          scale="row",
          hclust=function(x) stats::hclust(x,method="average"),
          distfun=function(x) stats::as.dist((1-cor(t(x)))/2),
          margin=c(6,0),
          trace="none",
          density="none",
          col=heatmapCol,
          labRow="",
          #labCol = names(countTable.kept),
          cexCol=1,
          main=paste0("LFC of all genes (n=",nrow(geneTable),")"),
          #ColSideColors=as.vector(colClrs),
          #colCol=as.vector(colClrs)
)




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

geneTable$XvA<-ifelse(geneTable$chr=="chrX","chrX","Autosomes")
#tmp<-geneTable


geneTable<-na.omit(geneTable)
allContrasts<-NULL
for (i in 1:ncol(combnTable)){
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
  facet_grid(rows=vars(XvA),cols=vars(contrast)) +
  geom_point(col="#11111155",size=1) + xlab(NULL) + ylab(NULL) +
  xlim(c(minScale,maxScale)) + ylim(c(minScale,maxScale)) +
  geom_smooth(method=lm,se=F,fullrange=T, size=0.7) + theme_bw() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=12),
        strip.text = element_text(size = 12)) +
  geom_hline(yintercept=0,lty=3,col="grey70",) +
  geom_vline(xintercept=0,lty=3,col="grey70") +
  ggpubr::stat_cor(aes(label = ..r.label..), method="pearson",
                   cor.coef.name = c("R"), output.type = "text")


lay <- rbind(c(1,1,2),
             c(3,3,3))
p<-gridExtra::arrangeGrob(p1,p2,p3,layout_matrix=lay)

ggsave(paste0(workDir,"/plots/RNAseqSupl.pdf"),p,device="pdf",width=10,height=10)



