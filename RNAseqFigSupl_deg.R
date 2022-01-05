library(tidyverse)
library(ggplot2)
library("ggpmisc")
library(eulerr)
library(lattice)
library(rtracklayer)
library(GenomicRanges)
library(ggpubr)

workDir<-getwd()
if(!dir.exists(paste0(workDir,"/plots"))) {
  dir.create(paste0(workDir,"/plots"))
}


contrastsOI<-c("wt.TIR1.sdc3deg.X_1mM_vs_0mM",
               "X.wt.wt.0mM_dpy26cs_vs_wt","X.TIR1.X.1mM_dpy26cssdc3deg_vs_wtwt")
useContrasts<-c("aux_sdc3bg","dpy26","sdc3dpy26")

prettyNames<-c(substitute(italic(x^AID),list(x="sdc-3")),
               substitute(italic(x^cs),list(x="dpy-26")),
               substitute(italic(x^AID*y^cs),list(x="sdc-3",y="dpy-26")))


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


# ########################-
# ## ECDF of data -----
# ########################-
#
# # using automatic filtering threshold
# sigTables<-list()
# localPadj=padjVal
# localLFC=0
# for (grp in groupsOI){
#   print(grp)
#   salmon<-readRDS(paste0(outPath,"/",fileNamePrefix,contrastsOI[[grp]], "_DESeq2_fullResults_p",padjVal,".rds"))
#   salmon<-salmon[!is.na(salmon$padj),]
#   #nrow(filterResults(salmon,padj=0.05,lfc= -0.5,direction="lt",chr="autosomes"))
#   print(paste0(nrow(salmon)," genes before filtering"))
#   print(paste0(sum(is.na(salmon$padj))," have padj that is NA"))
#   #salmon$expressed<-sum(salmon$baseMean>10)
#   sigTables[[grp]]<-as.data.frame(salmon) #[salmon$baseMean>10,]
#   print(paste0(nrow(sigTables[[grp]])," genes after automatic threshold filter"))
# }
#
# SMC<-rep(names(sigTables),lapply(sigTables,nrow))
# sig<-do.call(rbind,sigTables)
# sig$SMCpretty<-factor(SMC,levels=useContrasts,labels=prettyNames)
# sig$SMC<-factor(SMC,levels=useContrasts)
# table(sig$SMC)
# sig$XvA<-ifelse(sig$chr=="chrX","chrX","Autosomes")
# #sig$XvA[sig$chr=="chrX"]<-"chrX"
# #sig$XvA<-factor(sig$XvA)
# table(sig$XvA)
# sig$upVdown[sig$log2FoldChange<0]<-"down"
# sig$upVdown[sig$log2FoldChange>0]<-"up"
# sig$upVdown<-factor(sig$upVdown,levels=c("down","up"))
# table(sig$upVdown)
# row.names(sig)<-NULL
# SMC<-NULL
#
# lapply(sigTables,function(x){sum(x$padj<0.05)})
#
# options(tibble.width=Inf)
# dd1<-sig %>% dplyr::filter(padj<localPadj) %>%
#   dplyr::group_by(SMC,upVdown,XvA) %>%
#   dplyr::mutate(ecd=ecdf(abs(log2FoldChange))(abs(log2FoldChange)))
#
# dd1$upVdown<-relevel(dd1$upVdown,"up")
#
# p<-ggplot(dd1, aes(x=abs(log2FoldChange),y=ecd,color=SMCpretty,linetype=XvA)) +
#   geom_line(size=0.9)+ facet_wrap(vars(upVdown),nrow=2)+
#   theme_classic() + xlim(c(0,1.5)) +
#   scale_color_discrete(labels = ggplot2:::parse_safe(levels(dd1$SMCpretty)))+
#   xlab("Absolute log2 fold change")+ylab("Fraction of significant genes")
# #p
#
# #stat_ecdf(aes(colour=SMC,linetype=XvA),alpha=0.7)
# p1<-p+geom_vline(aes(xintercept = 0.5), color="grey") +
#   annotate("text",label="0.5",size=4, x=0.5, y=1,hjust=-0.05,color="grey") +
#   geom_vline(aes(xintercept = 0.15), color="grey") +
#   annotate("text",label="0.15",size=4, x=0.15, y=1,hjust=-0.05,color="grey") +
#   theme(strip.text = element_text(size = 12), axis.text=element_text(size=12),
#         axis.title=element_text(size=12))
# #p1
#
# # to get fraction between 0.15 and 0.5 lfc
# dd2<-sig %>% dplyr::filter(padj<localPadj) %>%
#   dplyr::group_by(SMC,upVdown,XvA) %>%
#   dplyr::summarise(ecd15=ecdf(abs(log2FoldChange))(c(0.15)),
#                    ecd50=ecdf(abs(log2FoldChange))(c(0.5)),
#                    q2=100*(ecd50-ecd15))
# dd2[order(dd2$q2),]
#
# # to get total expressed genes
# na.omit(sig) %>% group_by(SMC,XvA) %>% summarise(count=n())
#
# # to add tables of number of signifcant genes to plot
# ttu<-with(sig[sig$padj<localPadj & sig$upVdown=="up",],table(XvA,SMC))
# ttu<-pivot_wider(as.data.frame(ttu),names_from=SMC,values_from=Freq)
# names(ttu)<-gsub("XvA","",names(ttu))
#
# ttd<-with(sig[sig$padj<localPadj & sig$upVdown=="down",],table(XvA,SMC))
# ttd<-pivot_wider(as.data.frame(ttd),names_from=SMC,values_from=Freq)
# names(ttd)<-gsub("XvA","",names(ttd))
#
# tbs<-list("up"=ttu,
#           "down"=ttd)
#
# df <- tibble(x = rep(1.5, 2),
#              y = rep(0, 2),
#              upVdown= c("up","down"),
#              tbl = tbs)
#
# p1<-p1 + geom_table(data = df, aes(x = x, y = y,label = tbl),parse=T,
#                 table.theme = ttheme_gtlight)
#
# #p1



##################-
## boxplot of counts per chr------
##################-


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


p2  <-  ggplot(sigPerChr, aes(x=chr, y=genes, group=SMC)) +
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
#   facet_grid(cols=vars(XvA)) + ylim(yminmax) +
#   ggtitle("Absolute LFC of significantly changed genes by chromosome type") +
#   theme_minimal() + xlab("")+ ylab("|log2(FC)|")+
#   scale_fill_brewer(palette="Dark2",labels = ggplot2:::parse_safe(levels(sigTbl$SMC)))

###################-
## volcano plots------
###################-

plotList<-list()
geneTable<-NULL
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/",fileNamePrefix,contrastsOI[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
  p<-plotVolcanoXvA(salmon)
  plotList[[grp]]<-p+ggtitle(label=prettyNames[[which(groupsOI==grp)]])
}

p5<-ggarrange(plotlist=plotList,ncol=2,nrow=2)




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
minScale<- -2
maxScale<- 2
geneTable$XvA<-ifelse(geneTable$chr=="chrX","chrX","Autosomes")
#tmp<-geneTable

geneTable<-na.omit(geneTable)
allContrasts<-NULL
contrastNm<-NULL
for (i in 1:ncol(combnTable)){
  grp1<-groupsOI[combnTable[1,i]]
  grp2<-groupsOI[combnTable[2,i]]

  df<-geneTable[,c(paste0(grp1,"_lfc"),paste0(grp2,"_lfc"),"XvA")]
  names(df)<-c("group1","group2","XvA")
  df$contrast<-paste(grp1,"v",grp2)
  #df$Rval<-Rval
  if(is.null(allContrasts)){
    allContrasts<-df
    contrastNm<-paste(grp1,"v",grp2)
  } else {
    allContrasts<-rbind(allContrasts,df)
    contrastNm<-c(contrastNm,paste(grp1,"v",grp2))
  }
}


p4a<-ggplot(allContrasts,aes(x=group1,y=group2)) +
  facet_grid(cols=vars(XvA),rows=vars(contrast)) +
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

plotTitle<-list(label="Autosomal up",fontsize=8)
achr<-lapply(sigTables,function(x) x[x$chr!="chrX",])
sigGenes<-lapply(achr, "[[","wormbaseID")
numGenes["Aup"]<-sum(sapply(sigGenes,length))
fit<-euler(sigGenes)
p3a<-plot(fit, quantities=list(type=eulerLabelsType),
          main=plotTitle)
#print(p3a)


plotTitle<-list(label="chrX up",fontsize=8)
xchr<-lapply(sigTables,function(x) x[x$chr=="chrX",])
sigGenes<-lapply(xchr,"[[","wormbaseID")
numGenes["Xup"]<-sum(sapply(sigGenes,length))
fit<-euler(sigGenes)
p3c<-plot(fit, quantities=list(type=eulerLabelsType),
          main=plotTitle)
#print(p3c)


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

plotTitle<-list(label="Autosomal down",fontsize=8)
achr<-lapply(sigTables,function(x) x[x$chr!="chrX",])
sigGenes<-lapply(achr, "[[","wormbaseID")
numGenes["Adown"]<-sum(sapply(sigGenes,length))
fit<-euler(sigGenes)
p3e<-plot(fit, quantities=list(type=eulerLabelsType),
          main=plotTitle)
#print(p3e)


plotTitle<-list(label="chrX down",fontsize=8)
xchr<-lapply(sigTables,function(x) x[x$chr=="chrX",])
sigGenes<-lapply(xchr, "[[","wormbaseID")
numGenes["Xdown"]<-sum(sapply(sigGenes,length))
fit<-euler(sigGenes)
p3g<-plot(fit, quantities=list(type=eulerLabelsType),
          main=plotTitle)
#print(p3g)

p3<-ggarrange(p3a,p3c,p3e,p3g,ncol=2,nrow=2)



# ###################-
# ## AB compartments------
# ###################-
#
# pca2<-import.bw(paste0(outPath,"/otherData/N2_5000b_laminDamID_pca2.bw"))
#
# listgr<-NULL
# for (grp in groupsOI){
#   #grp=groupsOI[1]
#   salmon<-readRDS(file=paste0(outPath,"/",fileNamePrefix,contrastsOI[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
#   salmon<-salmon[!is.na(salmon$chr),]
#   salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)
#   salmongr<-sort(salmongr)
#   salmongr<-assignGRtoAB(salmongr,pca2,pcaName="N2")
#   salmongr$SMC<-grp
#   listgr[[grp]]<-salmongr
# }
#
# # genes that change significantly
# sigList<-lapply(lapply(listgr,data.frame), getSignificantGenes,
#                 padj=padjVal,lfc=lfcVal,direction="both",chr="autosomes",
#                 nameChrCol="seqnames")
#
# compartmentTable<-do.call(rbind,sigList)
# compartmentTable$updown<-ifelse(compartmentTable$log2FoldChange>0,"up","down")
# compartmentTable$updown<-factor(compartmentTable$updown,levels=c("up","down"))
# #compartmentTable$XvA<-ifelse(compartmentTable$seqnames=="chrX","chrX","Autosomes")
# compartmentTable$SMC<-factor(compartmentTable$SMC,levels=useContrasts,labels=prettyNames)
# #table(compartmentTable$updown)
#
# p6a<-ggplot(compartmentTable,aes(x=compartment,fill=updown)) + geom_bar(position=position_dodge()) +
#   facet_wrap(.~SMC, labeller=label_parsed) + theme_bw() +
#   scale_fill_grey() + ggtitle("Number of autosomal genes by compartment that significantly change")
#
#
# allList<-lapply(listgr,data.frame)
# allTable<-do.call(rbind,allList)
# allTable<-allTable[allTable$seqnames!="chrX",]
# allTable$updown<-ifelse(allTable$log2FoldChange>0,"up","down")
# allTable$updown[allTable$log2FoldChange==0]<-NA
# allTable$updown<-factor(allTable$updown,levels=c("up","down"))
# #allTable$XvA<-ifelse(allTable$seqnames=="chrX","chrX","Autosomes")
# allTable$SMC<-factor(allTable$SMC,levels=useContrasts,labels=prettyNames)
#
# pergrp<-countPerGroup<-compartmentTable %>% group_by(compartment,updown,SMC) %>% summarize(count=n())
# all<-allTable %>% group_by(compartment,updown,SMC) %>% summarize(all=n())
#
# pergrp<-left_join(pergrp,all)
# pergrp$fraction<-pergrp$count/pergrp$all
#
# p6b<-ggplot(pergrp,aes(x=compartment,y=fraction,fill=updown)) +
#   geom_bar(position=position_dodge(),stat="identity") +
#   facet_wrap(.~SMC, labeller=label_parsed) + theme_bw() +
#   scale_fill_grey() + ggtitle("Fraction of autosomal genes per compartment that significantly change")



############### Final assembly #########
p<-ggarrange(p1,p2,nrow=2,heights=c(2.5,1),labels=c("A.","B."))
ggsave(paste0(workDir,"/plots/RNAseqSupl_deg1.pdf"),p,device="pdf",width=8,height=11)

p<-ggarrange(p5,p2a,nrow=2,heights=c(3,1),labels=c("A.","B."))
ggsave(paste0(workDir,"/plots/RNAseqSupl_deg2.pdf"),p,device="pdf",width=8,height=11)

p<-ggarrange(p3,p4a,nrow=2,heights=c(1,1.5),labels=c("A.","B."))
ggsave(paste0(workDir,"/plots/RNAseqSupl_deg3.pdf"),p,device="pdf",width=8,height=11)

p<-ggarrange(p6a,p6b,nrow=2,heights=c(1,1),labels=c("A.","B."))
ggsave(paste0(workDir,"/plots/RNAseqSupl_deg4.pdf"),p,device="pdf",width=8,height=6)



