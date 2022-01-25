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

##############

contrastsOIall<-c("wt.TIR1.sdc3deg.X_1mM_vs_0mM","X.wt.wt.0mM_dpy26cs_vs_wt",
               "X.TIR1.X.1mM_dpy26cssdc3deg_vs_wtwt","X.wt.wt.0mM_kle2cs_vs_wt",
               "X.wt.wt.0mM_scc16cs_vs_wt","X.wt.wt.0mM_coh1cs_vs_wt",
               "X.wt.wt.0mM_scc1coh1cs_vs_wt")
useContrastsAll<-c("aux_sdc3bg","dpy26","sdc3dpy26","kle2","scc1","coh1","scc1coh1")

prettyNamesAll<-c(substitute(italic(x^AID),list(x="sdc-3")),
                  substitute(italic(x^cs),list(x="dpy-26")),
                  substitute(italic(x^AID*y^cs),list(x="sdc-3",y="dpy-26")),
               substitute(italic(x^cs),list(x="kle-2")),
               substitute(italic(x^cs),list(x="scc-1")),
               substitute(italic(x^cs),list(x="coh-1")),
               substitute(italic(x^cs*y^cs),list(x="scc-1",y="coh-1")))

names(contrastsOIall)<-useContrastsAll




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
  #nrow(filterResults(salmon,padj=0.05,lfc= -0.5,direction="lt",chr="autosomes"))
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

# to get fraction between 0.15 and 0.5 lfc
dd2<-sig %>% dplyr::filter(padj<localPadj) %>%
  dplyr::group_by(SMC,upVdown,XvA) %>%
  dplyr::summarise(ecd15=ecdf(abs(log2FoldChange))(c(0.15)),
                   ecd50=ecdf(abs(log2FoldChange))(c(0.5)),
                   q2=100*(ecd50-ecd15))
dd2[order(dd2$q2),]

# to get total expressed genes
na.omit(sig) %>% group_by(SMC,XvA) %>% summarise(count=n())

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


forTtest<-sig %>% dplyr::filter(SMC=="aux_sdc3bg")
forTtest %>% dplyr::group_by(XvA) %>% dplyr::summarise(avr=mean(log2FoldChange))

t.test(forTtest$log2FoldChange~forTtest$XvA)


forFisher<-sig %>% dplyr::filter(SMC=="aux_sdc3bg",padj<0.05,abs(log2FoldChange)>0.5)
ff<-forFisher %>% dplyr::group_by(XvA,upVdown) %>% dplyr::summarise(avr=n())
ff

fisher.test(matrix(c(ff$avr),nrow=2,byrow=F))

forTtest<-sig %>% dplyr::filter(SMC=="dpy26")
forTtest %>% dplyr::group_by(XvA) %>% dplyr::summarise(avr=mean(log2FoldChange))

forTtest<-sig %>% dplyr::filter(SMC=="sdc3dpy26")
forTtest %>% dplyr::group_by(XvA) %>% dplyr::summarise(avr=mean(log2FoldChange))

t.test(forTtest$log2FoldChange~forTtest$XvA)

forFisher<-sig %>% dplyr::filter(SMC=="sdc3dpy26",padj<0.05,abs(log2FoldChange)>0.5)
ff<-forFisher %>% dplyr::group_by(XvA,upVdown,.drop=F) %>% dplyr::summarise(avr=n())
ff

fisher.test(matrix(c(ff$avr),nrow=2,byrow=F))


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


p1  <-  ggplot(sigPerChr, aes(x=chr, y=genes, group=SMC)) +
  facet_grid(direction~SMC,scales="free",space="free",
             labeller=label_parsed) + theme_bw() +
  theme(legend.position="bottom", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), strip.text = element_text(size = 10),
        #strip.text.x=element_text(face="italic"),
        axis.text=element_text(size=10), axis.title=element_text(size=10)) +
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
  salmon<-readRDS(paste0(outPath,"/",fileNamePrefix,contrastsOI[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
  p<-plotVolcanoXvA(salmon,addLegend=F)
  plotList[[grp]]<-p+ggtitle(label=prettyNames[[which(groupsOI==grp)]])
}

plotList[["legend"]]<-plotVolcanoXvA(salmon,addLegend=T)
p2<-ggarrange(plotlist=plotList,ncol=2,nrow=2)




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
  xlim(c(minScale,maxScale)) + ylim(c(minScale,maxScale)) +
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
                   show.legend=F,size=3) +
  xlab(label=element_blank()) + ylab(label=element_blank())
#p3

# p4a<-ggplot(allContrasts,aes(x=group1,y=group2)) +
#   facet_grid(cols=vars(XvA),rows=vars(contrast)) +
#   geom_point(col="#11111155",size=1) + xlab(NULL) + ylab(NULL) +
#   xlim(c(minScale,maxScale)) + ylim(c(minScale,maxScale)) +
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

kramer<-as.data.frame(readRDS(file=paste0(outPath,"/publicData/kramer2015_L3_gr.rds")))

localPadj=0.05
localLFC=0.5


filterList<-list()
filterList[["Cycling_Meeuse"]]<-read.delim(paste0(outPath,"/publicData/oscillatingGenes.tsv"), header=T, stringsAsFactors=F)$wormbaseID #3739
filterList[["Cycling_Latorre"]]<-read.delim(paste0(outPath,"/publicData/oscillatingGenes_latorre.tsv"))$wormbaseID #3235
toFilter<-unique(unlist(filterList))

idx<-kramer$wormbaseID %in% toFilter
kramer<-kramer[!idx,]

eulerLabelsType=c("counts")
uglyNames<-c("sdc-3AID","dpy-26cs","sdc-3AIDdpy-26cs")
plotList<-list()
for (grp in 1:length(groupsOI)){
  salmon<-data.frame(readRDS(paste0(outPath,"/",fileNamePrefix,contrastsOI[[groupsOI[grp]]],"_DESeq2_fullResults_p",padjVal,".rds")))

  ###############################-
  ## chrX-----
  ###############################-

  kramerDpy27<-getSignificantGenes(kramer, padj=localPadj, lfc=localLFC,
                                   namePadjCol="dpy27_RNAi_L3_padj",
                                   nameLfcCol="dpy27_RNAi_L3_log2_fold_change",
                                   direction="both",
                                   chr="chrX", nameChrCol="seqnames", outPath=".")

  kramerDpy21<-getSignificantGenes(kramer, padj=localPadj, lfc=localLFC,
                                   namePadjCol="dpy21_mutant_L3_padj",
                                   nameLfcCol="dpy21_mutant_L3_log2_fold_change",
                                   direction="both",
                                   chr="chrX", nameChrCol="seqnames")

  salmondc<-getSignificantGenes(salmon,padjVal,lfcVal,direction="both",chr="chrX")


  DC<-list(salmon=salmondc$wormbaseID, dpy27=kramerDpy27$wormbaseID,
           dpy21=kramerDpy21$wormbaseID)
  names(DC)<-c(uglyNames[grp], "dpy-27 (Kramer 2015)", "dpy-21 (Kramer 2015)")

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
                                   chr="autosomes", nameChrCol="seqnames", outPath=".")

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
#p4


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


# ##########################-
# ## Manual clicked loops------
# ##########################-
#
# ### new clicked loops
# #clickedBatch="366"
# clickedBatch="382"
#
# loopsOrAnchors<-"anchors"
# ceTiles<-tileGenome(seqlengths(Celegans),tilewidth=10000,cut.last.tile.in.chrom = T)
#
# if(loopsOrAnchors=="loops"){
#   loops<-import(paste0(outPath,"/otherData/Clicked_loops_",clickedBatch,"_merge.bedpe"),format="bedpe")
#   grl<-zipup(loops)
#   anchor1<-do.call(c,lapply(grl,"[",1))
#   anchor2<-do.call(c,lapply(grl,"[",2))
#   mcols(anchor1)<-mcols(loops)
#   mcols(anchor2)<-mcols(loops)
#
#   anchor1$loopNum<-paste0("loop",1:length(anchor1))
#   anchor2$loopNum<-paste0("loop",1:length(anchor2))
#   anchors<-sort(c(anchor1,anchor2))
#
#   #separate anchors from inside tads
#   tads_in<-GRanges(seqnames=seqnames(anchor1),IRanges(start=end(anchor1)+1,end=start(anchor2)-1))
#   # tads_in<-reduce(tads_in)
#   tads_in<-resize(tads_in,width=width(tads_in)-20000,fix="center")
#
#   ol<-findOverlaps(ceTiles,tads_in)
#   tenkbInTads<-ceTiles[unique(queryHits(ol))]
#
#   anchors<-resize(anchors,width=10000,fix="center")
#   reduce(anchors)
#
#   ol<-findOverlaps(tenkbInTads,anchors)
#   tenkbInTads<-tenkbInTads[-queryHits(ol)]
# } else {
#   anchors<-import(paste0(outPath,"/otherData/382_X.eigs_cis.vecs_37peaks_p0.65_correct.bed"),format="bed")
#   anchors<-resize(anchors,width=10000,fix="center")
#   ol<-findOverlaps(ceTiles,anchors)
#   tenkbInTads<-ceTiles[-queryHits(ol)]
# }
#
#
#
# width(tenkbInTads)
# width(anchors)
# dataList<-list()
# #grp=useContrasts[3]
# statList<-list()
# set.seed(34091857)
# for (grp in useContrastsAll){
#   salmon<-readRDS(file=paste0(paste0(outPath,"/",fileNamePrefix,
#                                      contrastsOIall[[grp]],"_DESeq2_fullResults_p",padjVal,".rds")))
#
#   salmon<-salmon[!is.na(salmon$chr),]
#   salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)
#
#   salmongr<-sort(salmongr)
#
#   ol<-findOverlaps(salmongr,tenkbInTads,type="any",minoverlap=100L)
#   insideTads<-salmongr[unique(queryHits(ol))]
#
#   ol<-findOverlaps(salmongr,anchors,type="any",minoverlap=100L)
#   atAnchors<-salmongr[unique(queryHits(ol))]
#
#   ol<-findOverlaps(insideTads,atAnchors)
#   insideTads<-insideTads[-queryHits(ol)]
# #  bsAvr<-list()
# #  for(i in 1:10000){
# #    idx<-sample(1:length(insideTads),length(atAnchors))
# #    bsAvr<-c(bsAvr,mean(insideTads$log2FoldChange[idx]))
# #  }
# #  statList[[grp]]<-sum(unlist(bsAvr)>mean(atAnchors$log2FoldChange))/10000
#   insideTads$Loops<-"Non Anchor"
#   atAnchors$Loops<-"Anchor"
#
#
#   df<-data.frame(c(insideTads,atAnchors))
#   df<-df%>%dplyr::group_by(seqnames,Loops)%>%dplyr::mutate(count=n())
#   df$SMC<-grp
#
#   dataList[[grp]]<-df
# }
#
#
#
# #pvalsDiff<-1-unlist(statList)
# #names(pvalsDiff)<-prettyNamesAll
#
# #aux_sdc3bg      dpy26  sdc3dpy26       kle2       scc1       coh1   scc1coh1
# #0.8653     0.0001     0.0013     0.3441     0.3545     0.1655     0.2047
#
# ## focus on chrX loops
# dataTbl<-do.call(rbind,dataList)
# xchr<-dataTbl[dataTbl$seqnames=="chrX",]
# cntTbl<-xchr %>% dplyr::group_by(SMC,Loops) %>% dplyr::summarise(count=n()) %>%
#   filter(SMC=="dpy26")
#
# xchr$SMC<-factor(xchr$SMC,levels=useContrastsAll,labels=prettyNamesAll)
#
# p1<-ggplot(xchr,aes(x=Loops,y=log2FoldChange,fill=Loops))+
#   geom_boxplot(notch=T,outlier.shape=NA,varwidth=T)+
#   #geom_jitter()+
#   facet_grid(col=vars(SMC),labeller=label_parsed) +ylim(c(-0.5,1.5))+
#   ggtitle(paste0("LFC near ",clickedBatch," loop anchors (",
#                  cntTbl$count[cntTbl$Loops=="Anchor"]," genes) and inside loops (",
#                  cntTbl$count[cntTbl$Loops=="Inside loop"]," genes) in chrX")) +
#   geom_hline(yintercept=0,linetype="dotted",color="grey20") + theme_bw()+
#   theme(axis.text.x=element_text(angle=45,hjust=1),
#                                  panel.grid.major=element_blank(),
#                                  panel.grid.minor=element_blank())+
#   #scale_fill_discrete(c("darkgreen","darkblue"),labeller=label_parsed)+
#   xlab(label=element_blank()) +
#   ggsignif::geom_signif(test=wilcox.test,comparisons = list(c("Anchor", "Inside loop")),
#                         map_signif_level = F,tip_length=0.001,y_position=1.4,vjust=-0.1,
#                         textsize=3,margin_top=0)
#
# p1
# #xchr$SMC<-factor(xchr$SMC,levels=useContrastsAll,labels=prettyNamesAll)
# xchr$measure="Expression"
#
# bm<- xchr%>%distinct(wormbaseID,baseMean,measure)
# p2<-ggplot(bm,aes(x=Loops,y=log2(baseMean),fill=Loops))+
#   geom_boxplot(notch=T,outlier.shape=NA,varwidth=T) +
#   facet_wrap(.~measure) + ggtitle("Base mean counts") + theme_bw()+
#   theme(legend.position = "none",axis.text.x=element_text(angle=45,hjust=1),
#         panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank()) +
#   xlab(label=element_blank()) + ylab("Log2(base mean counts)")
#
# p<-ggarrange(p2,p1,ncol=2,widths=c(1.3,8.7))






############### Final assembly #########

p<-ggarrange(ggarrange(ggarrange(p1,p2,nrow=2,heights=c(1.2,2),labels=c("A ","B ")),p3,
          ncol=2,widths=c(2,1),labels=c("","C ")),
          p4,nrow=2,heights=c(2,1),labels=c("","D "))
p<-annotate_figure(p, top = text_grob("Das et al., Figure S10", size = 14))
ggsave(paste0(workDir,"/plots/RNAseqSupl_deg1.pdf"),p,device=cairo_pdf,width=8,height=11)
ggsave(paste0(workDir,"/plots/RNAseqSupl_deg1.png"),p,device=png,width=8,height=11,bg="white")




