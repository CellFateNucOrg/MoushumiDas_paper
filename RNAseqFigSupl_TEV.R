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

complexes<-c("condensinI/I^DC", "condensinII", "cohesin^(SCC-1)","cohesin^(COH-1)","cohesins")
names(complexes)<-useContrasts

names(contrastsOI)<-useContrasts
# strains<-c("366","382","775","784","828","844")
# strain<-factor(strains,levels=strains)
# SMC<-strain
# levels(SMC)<-c("TEVonly",useContrasts)

#controlGrp<-levels(SMC)[1] # control group
groupsOI<-useContrasts


source(paste0(workDir,"/functions.R"))



# panel cdf of lfc
lfcVal=0.5
padjVal=0.05
filterPrefix<-"filtCycChrAX"
fileNamePrefix=paste0("p",padjVal,"_lfc",lfcVal,"_",filterPrefix,"/",filterPrefix,"_")
outPath=paste0(workDir)


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
  theme_classic() + coord_cartesian(xlim=c(0,1.5)) +
  scale_color_discrete(labels = ggplot2:::parse_safe(levels(dd1$SMCpretty)),
                       name=element_blank())+
  scale_linetype_discrete(name=element_blank()) +
  xlab("Absolute log2 fold change")+ylab("Fraction of significant genes")


#stat_ecdf(aes(colour=SMC,linetype=XvA),alpha=0.7)
p<-p+geom_vline(aes(xintercept = 0.5), color="grey") +
  annotate("text",label="0.5",size=4, x=0.5, y=1, hjust=-0.05,color="grey") +
  geom_vline(aes(xintercept = 0.15), color="grey") +
  annotate("text",label="0.15",size=4, x=0.15, y=1,hjust=-0.05,color="grey") +
  theme(strip.text = element_text(size = 12), axis.text=element_text(size=12),
        axis.title=element_text(size=12))


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

#p1



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

p2<-ggplot(allSig,aes(x=chr,y=log2FoldChange,fill=chr)) +
  geom_boxplot(outlier.shape=NA,show.legend=F) +
  facet_grid(cols=vars(SMC),labeller=label_parsed) + theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  coord_cartesian(ylim=c(-0.7,1.5)) + scale_fill_grey(start=0.8, end=0.4)+
  geom_hline(yintercept=0,linetype="dashed",col="red")+
  xlab("Chromosome") + ylab("Log2(Fold Change)") +
  geom_text(data=facetLabels,aes(label=complexes),parse=T,x=2,y=1.4,size=3.5,hjust=0)
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



##################-
## venn diagrams------
##################-
eulerLabelsType=c("counts")
numGenes<-list()
## upregulated genes -----
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/",fileNamePrefix,contrastsOI[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
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
p3a<-plot(fit, quantities=list(type=eulerLabelsType),
          main=plotTitle, labels=list(font=4))
#print(p3a)

subset2<-achr[c("scc1","coh1","scc1coh1")]
sigGenes<-lapply(subset2, "[[","wormbaseID")
uglyNames<-c("scc-1cs","coh-1cs","scc-1cscoh-1cs")
names(sigGenes)<-uglyNames
numGenes["Aup"]<-sum(sapply(sigGenes,length))
fit<-euler(sigGenes)
p3b<-plot(fit, quantities=list(type=eulerLabelsType),
          main=plotTitle, labels=list(font=4))
#print(p3b)

plotTitle<-list(label="chrX",fontsize=8)
xchr<-lapply(sigTables,function(x) x[x$chr=="chrX",])
subset1<-xchr[c("dpy26","kle2","scc1")]
sigGenes<-lapply(subset1,"[[","wormbaseID")
uglyNames<-c("dpy-26cs","kle-2cs","scc-1cs")
names(sigGenes)<-uglyNames
numGenes["Xup"]<-sum(sapply(sigGenes,length))
fit<-euler(sigGenes)
p3c<-plot(fit, quantities=list(type=eulerLabelsType),
         main=plotTitle, labels=list(font=4))
#print(p3c)

subset2<-xchr[c("scc1","coh1","scc1coh1")]
sigGenes<-lapply(subset2,"[[","wormbaseID")
uglyNames<-c("scc-1cs","coh-1cs","scc-1cscoh-1cs")
names(sigGenes)<-uglyNames
numGenes["Xup"]<-sum(sapply(sigGenes,length))
fit<-euler(sigGenes)
p3d<-plot(fit, quantities=list(type=eulerLabelsType),
          main=plotTitle, labels=list(font=4))
#print(p3d)


## downregulated genes -----
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/",fileNamePrefix,contrastsOI[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))

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
p3e<-plot(fit, quantities=list(type=eulerLabelsType),
          main=plotTitle, labels=list(font=4))
#print(p3e)

subset2<-achr[c("scc1","coh1","scc1coh1")]
sigGenes<-lapply(subset2, "[[","wormbaseID")
uglyNames<-c("scc-1cs","coh-1cs","scc-1cscoh-1cs")
names(sigGenes)<-uglyNames
numGenes["Adown"]<-sum(sapply(sigGenes,length))
fit<-euler(sigGenes)
p3f<-plot(fit, quantities=list(type=eulerLabelsType),
          main=plotTitle, labels=list(font=4))
#print(p3f)


plotTitle<-list(label="chrX",fontsize=8)
xchr<-lapply(sigTables,function(x) x[x$chr=="chrX",])
subset1<-xchr[c("dpy26","kle2","scc1")]
sigGenes<-lapply(subset1, "[[","wormbaseID")
uglyNames<-c("dpy-26cs","kle-2cs","scc-1cs")
names(sigGenes)<-uglyNames
numGenes["Xdown"]<-sum(sapply(sigGenes,length))
fit<-euler(sigGenes)
p3g<-plot(fit, quantities=list(type=eulerLabelsType),
          main=plotTitle, labels=list(font=4))
#print(p3g)

subset2<-xchr[c("scc1","coh1","scc1coh1")]
sigGenes<-lapply(subset2, "[[","wormbaseID")
uglyNames<-c("scc-1cs","coh-1cs","scc-1cscoh-1cs")
names(sigGenes)<-uglyNames
numGenes["Xdown"]<-sum(sapply(sigGenes,length))
fit<-euler(sigGenes)
p3h<-plot(fit, quantities=list(type=eulerLabelsType),
          main=plotTitle, labels=list(font=4))
#print(p3h)

p3<-ggarrange(ggarrange(p3a,p3c,p3b,p3d,ncol=4, widths=c(0.8,1.3,1,0.5)),
                        ggarrange(p3e,p3g,p3f,p3h,ncol=4,widths=c(0.8,0.3,1.2,0.4)),
                                  nrow=2,labels=c("Up regulated","Down regulated"),
              heights=c(1,0.9))

#p3


##################-
## Correlation------
##################-

geneTable<-NULL
for (grp in groupsOI){
  salmon<-as.data.frame(readRDS(paste0(outPath,"/",fileNamePrefix,contrastsOI[[grp]],"_DESeq2_fullResults_p",padjVal,".rds")))
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
p4<-ggplot(allContrasts,aes(x=group1,y=group2,col=XvA)) +
  facet_wrap(.~contrast,nrow=2,labeller=label_parsed)+
  geom_point(size=1,alpha=0.4) +
  coord_cartesian(xlim=c(minScale,maxScale),ylim=c(minScale,maxScale)) +
  geom_smooth(method=lm,se=F,fullrange=T, size=0.7,show.legend = T) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        legend.position="bottom", strip.text.x = element_text(size = 9),
        legend.title=element_blank()) +
  scale_color_manual(values=c("#111111","#FF1111"))+
  geom_hline(yintercept=0,lty=3,col="grey70",) +
  geom_vline(xintercept=0,lty=3,col="grey70") +
  ggpubr::stat_cor(aes(label = ..r.label..), method="pearson",
                   cor.coef.name = c("R"), output.type = "text",
                   show.legend=F,size=3,
                   label.x=minScale+0.1, label.y=c(maxScale-0.1,maxScale-1))+
  xlab(label=element_blank()) + ylab(label=element_blank())
p4

###################-
## volcano plots------
###################-

plotList<-list()
geneTable<-NULL
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/",fileNamePrefix,contrastsOI[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
  p<-plotVolcanoXvA(salmon,addLegend=F)
  plotList[[grp]]<-p+ggtitle(label=prettyNames[[which(groupsOI==grp)]])
}

plotList[["legend"]]<-plotVolcanoXvA(salmon,addLegend=T)
p5<-ggarrange(plotlist=plotList,ncol=3,nrow=2)
p5

##########################-
## Manual clicked loops------
##########################-

### new clicked loops
#clickedBatch="366"
#clickedBatch="382"

loopsOrAnchors<-"anchors"
ceTiles<-tileGenome(seqlengths(Celegans),tilewidth=10000,cut.last.tile.in.chrom = T)

if(loopsOrAnchors=="loops"){
  loops<-import(paste0(outPath,"/otherData/Clicked_loops_",clickedBatch,"_merge.bedpe"),format="bedpe")
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
  anchors<-import(paste0(outPath,"/otherData/382_X.eigs_cis.vecs_37peaks_p0.65_correct.bed"),format="bed")
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

p6a<-ggplot(xchr,aes(x=Loops,y=log2FoldChange,fill=Loops))+
  geom_boxplot(notch=T,outlier.shape=NA,varwidth=T)+
  #geom_jitter()+
  facet_grid(col=vars(SMC),labeller=label_parsed) +coord_cartesian(ylim=c(-0.5,1.65))+
  ggtitle(paste0("LFC near ",clickedBatch," loop anchors (",
                 cntTbl$count[cntTbl$Loops=="Anchor"]," genes) and not at anchors (",
                 cntTbl$count[cntTbl$Loops=="Not anchor"]," genes) in chrX")) +
  geom_hline(yintercept=0,linetype="dotted",color="grey20") + theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        plot.title=element_text(size=10), legend.title = element_blank())+
  #scale_fill_discrete(c("darkgreen","darkblue"),labeller=label_parsed)+
  xlab(label=element_blank()) + ylab("Log2(fold change)")+
  ggsignif::geom_signif(test=t.test,comparisons = list(c("Anchor", "Not anchor")),
                        map_signif_level = F,tip_length=0.001,y_position=1.4,vjust=-0.1,
                        textsize=3,margin_top=0)+
  geom_text(data=facetLabels,aes(label=complexes),parse=T,x=1.5,y=1.65,size=3.5,hjust=0.5)

p6a


# #get 366tpm.
# gr<-GRanges(xchr[xchr$SMC=="italic(\"dpy-26\"^cs)",])
# seqlevels(gr)<-seqlevels(Celegans)
# tpm366<-import(paste0(workDir,"/otherData/sumFR_366_B_UniqueMultiple.bw"),
#                format="bw")
# tpm366<-ws235toCe11(tpm366)
# gr$baseMean<-binnedAverage(gr,mcolAsRleList(tpm366,"score"),varname="tpm366",na.rm=T)$tpm366
# gr$measure="Expression"
#bm<- data.frame(gr) #%>% distinct(wormbaseID,baseMean,measure)
bm<- data.frame(xchr) %>% distinct(wormbaseID,baseMean,measure)
p6b<-ggplot(bm,aes(x=Loops,y=log2(baseMean),fill=Loops))+
  geom_boxplot(notch=T,outlier.shape=NA,varwidth=T) +
  facet_wrap(.~measure) + ggtitle(" ") + theme_bw()+ coord_cartesian(ylim=c(-12,15)) +
  theme(legend.position = "none",axis.text.x=element_text(angle=45,hjust=1),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.title=element_text(size=10)) +
  xlab(label=element_blank()) + ylab("Log2(base mean counts)")+
  ggsignif::geom_signif(test=t.test,comparisons = list(c("Anchor", "Not anchor")),
                        map_signif_level = F,tip_length=0.001,vjust=-0.1,
                        textsize=3,margin_top=0.1)

p6<-ggarrange(p6b,p6a,ncol=2,widths=c(1.8,8.7))

p6

#################-
## plot length vs basal expression --------
#################-
chrSubset="autosomes"
localPadj=0.05
localLFC=0
grp=useContrasts[1]
listTbls<-list()
for(grp in useContrasts){
  salmon<-data.frame(readRDS(paste0(outPath,"/",fileNamePrefix,contrastsOI[[grp]],"_DESeq2_fullResults_p",padjVal,".rds")))
  sig<-getSignificantGenes(salmon, padj=localPadj, lfc=localLFC,
                           namePadjCol="padj",
                           nameLfcCol="log2FoldChange",
                           direction="both",
                           chr=chrSubset, nameChrCol="chr")
  sig$geneLength<-sig$end-sig$start
  sig$upVdown<-factor(ifelse(sig$log2FoldChange>0,"up","down"))
  sig$SMC<-grp
  listTbls[[grp]]<-sig
}


allSig<-do.call(rbind,listTbls)
allSig$SMC<-factor(allSig$SMC,levels=useContrasts,labels=prettyNames)
allSig$complexes<-complexes[allSig$SMC]
rownames(allSig)<-NULL
facetLabels<-allSig %>% select(SMC,complexes) %>% distinct()
facetLabels$log2FoldChange<-1

p7a<-ggplot(allSig,aes(x=log2(geneLength),y=log2(baseMean),color=log2FoldChange)) +
  geom_point(size=0.4) +
  scale_color_gradient2(low=scales::muted("#ff000055"),mid="#ffffff22",
                        high=scales::muted("#0000ff55"), na.value="#ffffff22",
                        limits=c(-0.5,0.5),oob=scales::squish,name="Log2FC")+
  facet_grid(rows=vars(SMC),labeller=label_parsed) +theme_bw()+
  ggtitle(paste0("Significantly changed genes on ",chrSubset," p<",localPadj," LFC>",localLFC))+
  theme(legend.position = "bottom", plot.title = element_text(size=12)) +
  xlab("Log2(gene length in bp)") + ylab("Log2(base mean expression)")+
  #geom_text(label=allSig$complexes,parse=T,x=7.5,y=15,size=3.5)
  geom_text(data=facetLabels,aes(label=complexes),parse=T,x=7.5,y=17,size=3.5,color="black",hjust=0)
p7a


uniqGenes<-allSig %>% distinct(wormbaseID,geneLength)

allSig$lengthBin<-cut(allSig$geneLength,quantile(uniqGenes$geneLength,seq(0,1,0.1)),
                      dig.lab=0,ordered_result=T,right=T,include.lowest=T)


labs<-data.frame(lower = factor( as.numeric(sub("(\\(|\\[)(.+),.*", "\\2", levels(allSig$lengthBin) ))),
                 upper = factor( as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", levels(allSig$lengthBin) ))))

levels(allSig$lengthBin)<-paste(levels(labs$lower), levels(labs$upper),sep="-")
#allSig[is.na(allSig$lengthBin),]
p7b<-ggplot(allSig, aes(x=lengthBin,fill=upVdown)) + geom_bar(position="fill") +
  facet_grid(rows=vars(SMC),labeller=label_parsed) + theme_bw() +
  scale_fill_manual(values=c(scales::muted("red"),scales::muted("blue")),name="Log2FC") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "bottom") +
  ylab("Fraction of genes up & down regulated") + xlab("Gene length (bp)") +
  ggtitle(paste0(" "))


# chrSubset="chrX"
# localPadj=0.05
# localLFC=0
# grp=useContrasts[1]
# listTbls<-list()
# for(grp in useContrasts){
#   salmon<-data.frame(readRDS(paste0(outPath,"/",fileNamePrefix,contrastsOI[[grp]],"_DESeq2_fullResults_p",padjVal,".rds")))
#   sig<-getSignificantGenes(salmon, padj=localPadj, lfc=localLFC,
#                            namePadjCol="padj",
#                            nameLfcCol="log2FoldChange",
#                            direction="both",
#                            chr=chrSubset, nameChrCol="chr")
#   sig$geneLength<-sig$end-sig$start
#   sig$upVdown<-factor(ifelse(sig$log2FoldChange>0,"up","down"))
#   sig$SMC<-grp
#   listTbls[[grp]]<-sig
# }
#
#
# allSig<-do.call(rbind,listTbls)
# allSig$SMC<-factor(allSig$SMC,levels=useContrasts,labels=prettyNames)
# p7c<-ggplot(allSig,aes(x=log2(geneLength),y=log2(baseMean),color=log2FoldChange)) +
#   geom_point(size=0.4) +
#   scale_color_gradient2(low=scales::muted("#ff000055"),mid="#ffffff22",
#                         high=scales::muted("#0000ff55"), na.value="#ffffff22",
#                         limits=c(-1,1),oob=scales::squish,name="Log2FC")+
#   facet_grid(rows=vars(SMC),labeller=label_parsed) +theme_bw()+
#   ggtitle(paste0("Significantly changed genes on ",chrSubset," p<",localPadj," LFC>",localLFC))+
#   theme(legend.position = "bottom",plot.title = element_text(size=12)) +
#   xlab("Log2(gene length in bp)") + ylab("Log2(base mean expression)")
# p7c
#
#
# uniqGenes<-allSig %>% distinct(wormbaseID,geneLength)
#
# allSig$lengthBin<-cut(allSig$geneLength,quantile(uniqGenes$geneLength,seq(0,1,0.1)),
#                       dig.lab=0,ordered_result=T,right=T,include.lowest=T)
#
#
# labs<-data.frame(lower = factor( as.numeric(sub("(\\(|\\[)(.+),.*", "\\2", levels(allSig$lengthBin) ))),
#                  upper = factor( as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", levels(allSig$lengthBin) ))))
#
# levels(allSig$lengthBin)<-paste(levels(labs$lower), levels(labs$upper),sep="-")
# #allSig[is.na(allSig$lengthBin),]
# p7d<-ggplot(allSig, aes(x=lengthBin,fill=upVdown)) + geom_bar(position="fill") +
#   facet_grid(rows=vars(SMC),labeller=label_parsed) + theme_bw() +
#   scale_fill_manual(values=c(scales::muted("red"),scales::muted("blue")),name="Log2FC") +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "bottom") +
#   ylab("Fraction of genes up & down regulated") + xlab("Gene length (bp)") +
#   ggtitle(paste0(" "))
#p7d


#################-
## plot length vs 366 TPM --------
#################-
chrSubset="autosomes"
localPadj=0.05
localLFC=0
grp=useContrasts[1]
listTbls<-list()
for(grp in useContrasts){
  salmon<-data.frame(readRDS(paste0(outPath,"/",fileNamePrefix,contrastsOI[[grp]],"_DESeq2_fullResults_p",padjVal,".rds")))
  sig<-getSignificantGenes(salmon, padj=localPadj, lfc=localLFC,
                           namePadjCol="padj",
                           nameLfcCol="log2FoldChange",
                           direction="both",
                           chr=chrSubset, nameChrCol="chr")
  sig$geneLength<-sig$end-sig$start
  sig$upVdown<-factor(ifelse(sig$log2FoldChange>0,"up","down"))
  sig$SMC<-grp
  listTbls[[grp]]<-sig
}


allSig<-do.call(rbind,listTbls)
allSig$SMC<-factor(allSig$SMC,levels=useContrasts,labels=prettyNames)
allSig$complexes<-complexes[allSig$SMC]
rownames(allSig)<-NULL
facetLabels<-allSig %>% select(SMC,complexes) %>% distinct()
facetLabels$log2FoldChange<-1

p7a<-ggplot(allSig,aes(x=log2(geneLength),y=log2(baseMean),color=log2FoldChange)) +
  geom_point(size=0.4) +
  scale_color_gradient2(low=scales::muted("#ff000055"),mid="#ffffff22",
                        high=scales::muted("#0000ff55"), na.value="#ffffff22",
                        limits=c(-0.5,0.5),oob=scales::squish,name="Log2FC")+
  facet_grid(rows=vars(SMC),labeller=label_parsed) +theme_bw()+
  ggtitle(paste0("Significantly changed genes on ",chrSubset," p<",localPadj," LFC>",localLFC))+
  theme(legend.position = "bottom", plot.title = element_text(size=12)) +
  xlab("Log2(gene length in bp)") + ylab("Log2(base mean expression)")+
  #geom_text(label=allSig$complexes,parse=T,x=7.5,y=15,size=3.5)
  geom_text(data=facetLabels,aes(label=complexes),parse=T,x=7.5,y=17,size=3.5,color="black",hjust=0)
p7a


uniqGenes<-allSig %>% distinct(wormbaseID,geneLength)

allSig$lengthBin<-cut(allSig$geneLength,quantile(uniqGenes$geneLength,seq(0,1,0.1)),
                      dig.lab=0,ordered_result=T,right=T,include.lowest=T)


labs<-data.frame(lower = factor( as.numeric(sub("(\\(|\\[)(.+),.*", "\\2", levels(allSig$lengthBin) ))),
                 upper = factor( as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", levels(allSig$lengthBin) ))))

levels(allSig$lengthBin)<-paste(levels(labs$lower), levels(labs$upper),sep="-")
#allSig[is.na(allSig$lengthBin),]
p7b<-ggplot(allSig, aes(x=lengthBin,fill=upVdown)) + geom_bar(position="fill") +
  facet_grid(rows=vars(SMC),labeller=label_parsed) + theme_bw() +
  scale_fill_manual(values=c(scales::muted("red"),scales::muted("blue")),name="Log2FC") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "bottom") +
  ylab("Fraction of genes up & down regulated") + xlab("Gene length (bp)") +
  ggtitle(paste0(" "))





####-
## AB compartment by chromosome 366-----
####-

####### Supplementary functions---------------

#' Collect counts of significantly changed genes per chromosome and per compartment
#'
#' @param listgr List of GRanges for different RNAseq with PCA data in mcols
#' @param namePCAcol Name of PCA column in the listgr object
#' @param padjVal adjusted p value to use as threshold
#' @param lfcVal Log2 fold change value to use as threshold
#' @result Table of significant genes with significant up and down regulated genes
#' counted by chromosome for each sample
#' @export
processCountsPerChr<-function(listgr,namePCAcol,padjVal=0.05,lfcVal=0.5){
  # genes that change significantly
  sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                  padj=padjVal,lfc=lfcVal,direction="both")

  # count genes by category (chr & A/B)
  dfl<-lapply(sigList, function(x){x%>% dplyr::group_by(seqnames, get(namePCAcol)) %>% tally()})
  bgCount<-lapply(lapply(listgr,as.data.frame), function(x){x%>% dplyr::group_by(seqnames,get(namePCAcol)) %>% tally()})
  # add name of SMC protein
  dfl<-do.call(rbind, mapply(cbind,dfl,"SMC"=names(dfl),SIMPLIFY=F))
  bgCount<-do.call(rbind, mapply(cbind,bgCount,"SMC"=names(bgCount),SIMPLIFY=F))
  names(dfl)<-c("seqnames",namePCAcol,"n","SMC")
  names(bgCount)<-c("seqnames",namePCAcol,"n","SMC")
  dfl$seqnames<-gsub("chr","",dfl$seqnames)
  bgCount$seqnames<-gsub("chr","",bgCount$seqnames)

  # do left join to make sure dfl has all the categories required
  dfl<-left_join(bgCount,dfl,by=c("seqnames",namePCAcol,"SMC"),suffix=c("_total",""))
  dfl$n[is.na(dfl$n)]<-0
  dfl$Frac<-dfl$n/dfl$n_total
  return(dfl)
}


#' Plot fractions of significantly changed genes per chromosome and per compartment
#'
#' @param df Data frame with fractions of significant genes per chromosome and compartment
#' @param namePCAcol Name of PCA column in the df object
#' @param padjVal Name of eigen vector to use in plot title
#' @result Plot of significant genes with significant changed genes
#' counted by chromosome for each sample presented as fraction
#' @export
plotFractionPerChrPerCompartment<-function(df,namePCAcol,namePCA){
  p<-ggplot(df,aes(x=seqnames,y=Frac,group=get(namePCAcol))) +
    geom_bar(stat="identity", position=position_dodge(),aes(fill=get(namePCAcol))) +
    facet_grid(cols=vars(SMC),labeller=label_parsed) +
    theme_minimal() + scale_fill_grey(start=0.8, end=0.2) +
    xlab("chr")+ylab("Fraction of genes") +
    #ggtitle(paste0("Fraction changed genes per chromosome by ",namePCA," compartment")) +
    labs(fill=namePCA)
  return(p)
}


####
## 366 compartments -----
####
pca1<-import.bw(paste0(outPath,"/otherData/366_merge_2000.oriented_E1.vecs.bw"))
pca2<-import.bw(paste0(outPath,"/otherData/366_merge_2000.oriented_E2.vecs.bw"))

listgr<-NULL
grp=useContrasts[1]
for (grp in useContrasts){
  salmon<-readRDS(file=paste0(paste0(outPath,"/",fileNamePrefix,
                contrastsOI[[grp]],"_DESeq2_fullResults_p",padjVal,".rds")))

  salmon<-salmon[!is.na(salmon$chr),]
  salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)

  salmongr<-sort(salmongr)

  salmongr<-assignGRtoAB(salmongr,pca1,pcaName="E1")
  salmongr<-assignGRtoAB(salmongr,pca2,pcaName="E2")
  listgr[[grp]]<-salmongr
}

pcaSource="366"

### 366 first Eigenvector ------
dfl<-processCountsPerChr(listgr,namePCAcol="E1_compartment")
dfl$SMC<-factor(dfl$SMC,levels=useContrasts,labels=prettyNames)
p8a<-plotFractionPerChrPerCompartment(dfl, namePCAcol="E1_compartment", namePCA=paste(pcaSource, "E1"))
p8a

### 366 second eigenvector ------
dfl<-processCountsPerChr(listgr,namePCAcol="E2_compartment")
dfl$SMC<-factor(dfl$SMC,levels=useContrasts,labels=prettyNames)
p8b<-plotFractionPerChrPerCompartment(dfl, namePCAcol="E2_compartment", namePCA=paste(pcaSource, "E2"))
p8b

p8<-ggpubr::ggarrange(p8a,p8b,ncol=2,nrow=1)




############### Final assembly #########

p<-ggarrange(p1,p2,p5,nrow=3,heights=c(2.5,1.5,3),labels=c("A ","B ","C "))
p<-annotate_figure(p, top = text_grob("Das et al., Figure S5", size = 14)) #face=bold
ggsave(paste0(workDir,"/plots/RNAseqSupl_TEV1.pdf"),p,device=cairo_pdf,width=21,height=29.7,
       unit="cm")
ggsave(paste0(workDir,"/plots/RNAseqSupl_TEV1.png"),p,device="png",width=21,height=29.7,
       unit="cm",bg="white")

#nullp<-NULL
p<-ggarrange(p3,p4,p6,nrow=3,heights=c(1,1.2,1),labels=c("A ","B ","C "))
p<-annotate_figure(p, top = text_grob("Das et al., Figure S6", size = 14))
ggsave(paste0(workDir,"/plots/RNAseqSupl_TEV2.pdf"),p,device="pdf",width=21,height=29.7,
       unit="cm")
ggsave(paste0(workDir,"/plots/RNAseqSupl_TEV2.png"),p,device="png",width=21,height=29.7,
       unit="cm",bg="white")

# p<-ggarrange(ggarrange(p6a,p6b,ncol=2,widths=c(3,1.5), labels=c("A ","B ")),
#               ggarrange(p6c,p6d,ncol=2,widths=c(3,1.5), labels=c("C ","D ")),nrow=2)

p<-ggarrange(p7a,p7b,ncol=2,widths=c(3,1.5), labels=c("A ","B "))
p<-annotate_figure(p, top = text_grob("Das et al., Figure S8", size = 14))
ggsave(paste0(workDir,"/plots/RNAseqSupl_TEV3.pdf"),p,device="pdf",width=21,height=29.7,units="cm")
ggsave(paste0(workDir,"/plots/RNAseqSupl_TEV3.png"),p,device="png",width=21,height=29.7,units="cm",bg="white")

