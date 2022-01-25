#' Make directories
#'
#' @param path String with path to where the directories should be made
#' @param dirNameList Vector of strings with names of directories to create (can include multilevel directories)
#' @return Creates the directories  listed in dirNameList
#' @examples
#' makeDirs(path=".",dirNameList=c("/txt","/rds/sample1"))
#' @export
makeDirs<-function(path,dirNameList=c()) {
  sub("\\/$","",path) #remove directory slash if present in path string
  for (d in dirNameList) {
    if (!dir.exists(paste0(path,"/",d))){  # for alignments
      dir.create(paste0(path,"/",d), recursive=TRUE, showWarnings=FALSE)
    }
  }
}


#' Filter DESeq2 table results
#'
#' @param resultsTable Table of DESeq results
#' @param padj Adjusted p value threshold
#' @param lfc Log fold change threshold
#' @param direction Whether to find genes that are less than (lt), or greater than (gt) the log fold change threshold, or both extreme tails ("both")
#' @param chr Include all genes in genome ("all") only those on the X chromosome ("chrX"), or only autosomes ("autosomes")
#' @param outPath Path to working directory
#' @param filenamePrefix Text to add to filename
#' @param writeTable Should results table be automatically saved to a file (default: True)
#' @param IDcolumnName Name of column containing feature IDs (default: "ID")
#' @param colsToCopy Vector with names of columns to copy in addition to baseMean,
#' log2FoldChange, padj,
#' @return filtered table of results which is also automatically written to disk
#' @export
filterResults<-function(resultsTable, padj=0.05, lfc=0, direction="both",
                        chr="all", outPath=".", filenamePrefix="",
                        writeTable=T,IDcolName="wormbaseID",
                        colsToCopy=c("chr","start","end","strand")) {
  sigGenes<-getSignificantGenes(resultsTable, padj, lfc, direction=direction,
                                chr=chr)
  idx<-resultsTable[,IDcolName] %in% sigGenes[,IDcolName]
  filtTable<-resultsTable[idx,c("baseMean","log2FoldChange","padj",
                                  IDcolName,colsToCopy)]
  if(writeTable){
    if(!dir.exists(paste0(outPath,"/txt"))){
      dir.create(paste0(outPath,"/txt"))
    }
    write.csv(filtTable,file=paste0(outPath,"/txt/filtResults_p",
                                    padj,"_",direction,"-","lfc",
                                    lfc,"_",chr,".csv"), row.names=F,
              quote=F)
  }
  return(filtTable)
}




#' Get significant genes from  RNAseq results
#'
#' @param resultsTable Table of DESeq results
#' @param padj Adjusted p value threshold
#' @param lfc Log fold change threshold
#' @param namePadjCol Name of column with adjusted P values
#' @param nameLFCcol Name of column with log fold change values
#' @param direction Character value to indicate whether to find genes that are
#' less than (lt), or greater than (gt) the log fold change threshold,
#' or both extreme tails ("both")
#' @param chr Include all genes in genome ("all") only those on the X chromosome
#' ("chrX"), or only autosomes ("autosomes")
#' @param nameChrCol Name of column with chromosome names.
#' @return Filtered table of significant genes at a certain log fold change and adjusted p value.
#' @export
getSignificantGenes<-function(resultsTable, padj=0.05, lfc=0, namePadjCol="padj",
                              nameLfcCol="log2FoldChange", direction="both",
                              chr="all", nameChrCol="chr", outPath="."){
  #remove rows with padj NA value
  idx<-is.na(resultsTable[,namePadjCol])
  resultsTable<-resultsTable[!idx,]
  # do filtering
  if(direction=="both") {
    idx<-!is.na(resultsTable[,namePadjCol]) & resultsTable[,namePadjCol]<padj & abs(resultsTable[,nameLfcCol])>lfc
  } else if(direction=="gt") {
    idx<-!is.na(resultsTable[,namePadjCol]) & resultsTable[,namePadjCol]<padj & resultsTable[,nameLfcCol]>lfc
  } else if(direction=="lt") {
    idx<-!is.na(resultsTable[,namePadjCol]) & resultsTable[,namePadjCol]<padj & resultsTable[,nameLfcCol]<lfc
  } else {
    print("direction must be 'both' to get both tails, \n'gt' to get lfc larger than a specific value, \nor 'lt' to get lfc less than a certain value")
  }
  if(chr=="all"){
    idx<-idx
  } else if(chr=="chrX"){
    idx<-idx & !is.na(resultsTable[,nameChrCol]) & resultsTable[,nameChrCol]=="chrX"
  } else if(chr=="autosomes"){
    idx<-idx & !is.na(resultsTable[,nameChrCol]) & resultsTable[,nameChrCol]!="chrX"
  } else {
    print("chr must be one of 'all', 'chrX' or 'autosomes'")
  }
  filtTable<-resultsTable[idx,]
  return(filtTable)
}


#' Convert short gene name to proper format
#'
#' Takes a gene name like scc1cs and adds hyphen. Function is based
#' on finding one or more digits in the name and putting a hyphen
#' in front of it.
#' @param geneName Gene name without hyphen, e.g. scc1cs
#' @return Gene name with hyphen, e.g. scc-1cs
#' @export
prettyGeneName<-function(geneName){
  gsub("([[:digit:]]+)","-\\1",geneName)
}


#' Assign GRanges to A/B compartment
#'
#' Function takes in a genomic ranges object which you wish to assign,
#' and a genomic ranges for an eigen vector giving the compartments,
#' and then adds pcaScore and compartment (A or B) columns to the
#' metadata of the GRanges object. If strings for the names of the gr
#' and pca are given a bedgraph will be produced with the AB
#' assignments of all the genes.
#' @param gr Genomic ranges object which you wish to assign to A/B compartments
#' @param pcagr Genomic ranges for eigen vector giving the A/B compartments
#' @param grName Name of the data in the gr object. Used to output a bedgraph of AB assignments
#' @param pcaName Name of the data in the pcagr object. Used to output a bedgraph of AB assignments
#' @return The gr genomic ranges that was input with additional metadata columns of pcaScore and compartment
assignGRtoAB<-function(gr, pcagr,grName=NULL,pcaName=NULL,
                       outPath="."){
  ol<-as.data.frame(findOverlaps(gr,pcagr,ignore.strand=T))
  ol$subjectScore<-pcagr$score[ol$subjectHits]
  ol$queryHits<-ol$queryHits
  pcaScore<-ol %>% group_by(queryHits)%>% dplyr::summarise(pcaScore=mean(subjectScore,na.rm=T))

  mcols(gr)[,paste0(pcaName,"_score")]<-NA
  mcols(gr)[,paste0(pcaName,"_compartment")]<-NA
  mcols(gr)[,paste0(pcaName,"_score")][pcaScore$queryHits]<-pcaScore$pcaScore
  mcols(gr)[,paste0(pcaName,"_compartment")]<-as.factor(ifelse(mcols(gr)[,paste0(pcaName,"_score")]>0,"A","B"))
  idx<-is.na(mcols(gr)[,paste0(pcaName,"_score")])
  print(paste0(sum(idx)," genes have no overlapping PCA bin"))
  gr<-gr[! idx ]
  forBG<-gr
  forBG$score<-ifelse(mcols(gr)[,paste0(pcaName,"_compartment")]=="A",1,-1)
  if(!(is.null(grName) | is.null(pcaName))){
    export(forBG,con=paste0(outPath,"/tracks/",grName,"_",
                            "__Compartments_",pcaName, ".bedGraph"),
           format="bedGraph")
  }
  return(gr)
}



#' Calculate number of significant genes at a range of thresholds
#'
#' @param dds DESeq2 object
#' @param contrastOI vector of values for contrast argument for DESeq2 results
#' function in the format of c("condition","treated","untreated")
#' @param padjVals vector of adjusted p value thresholds to explore
#' @param lfcVals vector of log2 fold change thresholds to explore
#' @param direction Vector of strings to indicate whether to find genes that are
#' less than (lt), or greater than (gt) the log fold change threshold,
#' or both extreme tails ("both")
#' @param chr String to indicate whether to analyse all genes in
#' genome ("all") only those on the X chromosome ("chrX"), or only autosomes
#' ("autosomes")
#' @param outPath path to working directory
#' @param filenamePrefix Text to add to filename
#' @return A table of the number of significant genes at different threshold values
#' @export
varyThreshold<-function(dds, contrastOI, padjVals=c(0.05,0.01),
                        lfcVals=c(0,0.25,0.5,1),
                        direction=c("both","lt","gt"), chr="all",
                        outPath=".",fileNamePrefix="",shrink=F){
  # make table of all combinations
  thresholds<-expand.grid(group=grp, lfc=lfcVals, padj=padjVals,
                          direction=c("both","lt","gt"),
                          stringsAsFactors=F)
  thresholds$totalAnalysed<-NA
  thresholds$baseMeanGt10<-NA
  thresholds$numSignificant<-NA
  thresholds$numSigGt10<-NA

  for(pval in padjVals){
    res<-resLFC<-NULL
    res<-results(dds,contrast=contrastOI,alpha=pval)
    if(shrink){
      resLFC<-lfcShrink(dds,coef=paste0(contrastOI[1],"_",contrastOI[2],"_vs_",contrastOI[3]), type="apeglm", res=res)
      res<-resLFC
    }

    pdf(file=paste0(outPath,"/plots/",fileNamePrefix,"_independentFilter_",
                    pval,".pdf"), width=8, height=11, paper="a4")

    plot(metadata(res)$filterNumRej,
         type="b", ylab="number of rejections",
         xlab="quantiles of filter",
         main=paste0("Threshold for independant filtering, alpha=",pval))
    lines(metadata(res)$lo.fit, col="red")
    abline(v=metadata(res)$filterTheta)
    legend("topright",legend=paste0("Read count \nthreshold: ", round(metadata(res)$filterThreshold,2)))
    dev.off()

    ### add metadata
    res$ID<-rownames(res)
    idx<-match(rownames(res),rowData(dds)$gene)
    res$chr<-factor(rowData(dds)$chr,levels=paste0("chr",c("I","II","III","IV","V","X")))[idx]
    if(chr=="chrX"){
      res<-res[na.omit(which(res$chr=="chrX")),]
    }
    if(chr=="autosomes"){
      res<-res[na.omit(which(res$chr!="chrX")),]
    }

    res<-res[!is.na(res$padj),]
    currentPidx<-which(thresholds$padj==pval)
    thresholds$totalAnalysed[currentPidx]<-sum(!is.na(res$padj))
    thresholds$baseMeanGt10[currentPidx]<-sum(res$baseMean>10,na.rm=T)

    for (i in currentPidx){
      thresholds$numSignificant[i]<-nrow(filterResults(res,
                                            padj=thresholds$padj[i],
                                            lfc=thresholds$lfc[i],
                                            direction=thresholds$direction[i],
                                            chr=chr, writeTable=F,
                                            IDcolName="ID",
                                            colsToCopy=c("chr")))
      thresholds$numSigGt10[i]<-nrow(filterResults(res[res$baseMean>10,],
                                              padj=thresholds$padj[i],
                                              lfc=thresholds$lfc[i],
                                              direction=thresholds$direction[i],
                                              chr=chr, writeTable=F,
                                              IDcolName="ID",
                                              colsToCopy=c("chr")))
    }
    thresholds$percentSignificant<-round(100*thresholds$numSignificant/thresholds$totalAnalysed,2)
    thresholds$percentSigGt10<-round(100*thresholds$numSigGt10/thresholds$baseMeanGt10,2)
  }
  return(thresholds)
}


#' Calculate number of significant genes at a range of thresholds
#'
#' @param res DESeq2 results table
#' @param pval vector of adjusted p value thresholds to explore
#' @param lfcVals vector of log2 fold change thresholds to explore
#' @param direction Vector of strings to indicate whether to find genes that are
#' less than (lt), or greater than (gt) the log fold change threshold,
#' or both extreme tails ("both")
#' @param chrs Vector of strings to indicate whether to analyse all genes in
#' genome ("all") only those on the X chromosome ("chrX"), or only autosomes
#' ("autosomes")
#' @return A table of the number of significant genes at different threshold values
#' @export
varyThreshold1<-function(res, pval=0.05,
                        lfcVals=c(0,0.25,0.5,1),
                        direction=c("both","lt","gt"),
                        chrs=c("all","chrX","autosomes")){
  # make table of all combinations
  thresholds<-expand.grid(group=grp, lfc=lfcVals, padj=pval,
                          direction=direction, chr=chrs,
                          stringsAsFactors=F)
  thresholds$totalAnalysed<-NA
  thresholds$baseMeanGt10<-NA
  thresholds$numSignificant<-NA
  thresholds$numSigGt10<-NA
  res<-res[!is.na(res$padj),]

  for(i in 1:nrow(thresholds)){
    if(thresholds$chr[i]=="chrX"){
      res1<-res[na.omit(which(res$chr=="chrX")),]
    } else if(thresholds$chr[i]=="autosomes"){
      res1<-res[na.omit(which(res$chr!="chrX")),]
    } else {
      res1<-res
    }

    thresholds$totalAnalysed[i]<-sum(!is.na(res1$padj))
    thresholds$baseMeanGt10[i]<-sum(res1$baseMean>10,na.rm=T)

    thresholds$numSignificant[i]<-nrow(filterResults(res1,
                                                     padj=thresholds$padj[i],
                                                     lfc=thresholds$lfc[i],
                                                     direction=thresholds$direction[i],
                                                     chr=thresholds$chr[i],
                                                     writeTable=F,
                                                     IDcolName="wormbaseID",
                                                     colsToCopy=c("chr")))
    thresholds$numSigGt10[i]<-nrow(filterResults(res1[res1$baseMean>10,],
                                                 padj=thresholds$padj[i],
                                                 lfc=thresholds$lfc[i],
                                                 direction=thresholds$direction[i],
                                                 chr=thresholds$chr[i],
                                                 writeTable=F,
                                                 IDcolName="wormbaseID",
                                                 colsToCopy=c("chr")))
  }
  thresholds$percentSignificant<-round(100*thresholds$numSignificant/thresholds$totalAnalysed,2)
  thresholds$percentSigGt10<-round(100*thresholds$numSigGt10/thresholds$baseMeanGt10,2)

  return(thresholds)
}




#' Calculate number of significant genes at a range of thresholds
#'
#' @param dds DESeq2 object
#' @param contrastOI Vector of values for contrast argument for DESeq2 results
#' function in the format of c("condition","treated","untreated")
#' @param padjVals Vector of adjusted p value thresholds to explore
#' @param breaks Vector of numbers to serve as breaks for binning the data
#' @param direction String to indicate whether to find genes that are
#' less than (lt), or greater than (gt) the log fold change threshold,
#' or both extreme tails ("both")
#' @param chr String to indicate whether to analyse all genes in
#' genome ("all") only those on the X chromosome ("chrX"), or only autosomes
#' ("autosomes")
#' @param asCounts Logical value to indicate if the function should return counts or proportions
#' @return A table of the number of significant genes at different threshold values
#' @export
getDensity<-function(dds, contrastOI, padjVals=c(0.05,0.01),
                     breaks=c(seq(0,2,0.1),Inf), chr="all",
                     direction="both",asCounts=F,shrink=T){
  groupCounts<-res<-NULL
  breakLabels<-levels(cut(breaks[-1],breaks))
  groupCounts<-data.frame(breaks=breakLabels,
                          pvals=rep(padjVals,each=(length(breaks)-1)), counts=NA)
  for(pval in padjVals){
    res<-NULL
    res<-results(dds,contrast=contrastOI,alpha=pval)
    if(shrink){
      resLFC<-lfcShrink(dds,coef=paste0(contrastOI[1],"_",contrastOI[2],"_vs_",contrastOI[3]),
                        type="apeglm", res=res)
      print(paste0(sum(res$log2FoldChange!=resLFC$log2FoldChange,na.rm=T)," lfc values are different after shrinkage"))
      res<-resLFC
    }
    res$ID<-rowData(dds)$gene
    res$chr<-rowData(dds)$chr
    sig<-filterResults(res, padj=pval, lfc=0, direction=c(direction),
                       chr=chr, writeTable=F, IDcolName="ID", colsToCopy=c("chr"))
    group_tags<-cut(abs(na.omit(sig$log2FoldChange)), breaks=breaks, include.lowest=TRUE, right=TRUE)

    if(asCounts){
      groupCounts[groupCounts$pvals==pval,"counts"]<-summary(group_tags)
    } else {
      groupCounts[groupCounts$pvals==pval,"counts"]<-summary(group_tags)/length(group_tags)
    }
  }
  return(list(groupCounts,sig))
}


#' Calculate number of significant genes at a range of thresholds
#'
#' @param res DESeq2 results table
#' @param pval Adjusted p value threshold to use
#' @param breaks Vector of numbers to serve as breaks for binning the data
#' @param chr String to indicate whether to analyse all genes in
#' genome ("all") only those on the X chromosome ("chrX"), or only autosomes
#' ("autosomes")
#' @param direction String to indicate whether to find genes that are
#' less than (lt), or greater than (gt) the log fold change threshold,
#' or both extreme tails ("both")
#' @param asCounts Logical value to indicate if the function should return counts or proportions
#' @return A table of the number of significant genes at different threshold values
#' @export
getDensity1<-function(res, pval=0.05,
                     breaks=c(seq(0,2,0.1),Inf), chr="all",
                     direction="both",asCounts=F){
  breakLabels<-levels(cut(breaks[-1],breaks))
  groupCounts<-data.frame(breaks=breakLabels,
                          pvals=rep(pval,each=(length(breaks)-1)), counts=NA)

  sig<-filterResults(res, padj=pval, lfc=0, direction=c(direction),
                     chr=chr, writeTable=F, IDcolName="wormbaseID", colsToCopy=c("chr"))
  group_tags<-cut(abs(na.omit(sig$log2FoldChange)), breaks=breaks, include.lowest=TRUE, right=TRUE)

  if(asCounts){
      groupCounts[groupCounts$pvals==pval,"counts"]<-summary(group_tags)
  } else {
      groupCounts[groupCounts$pvals==pval,"counts"]<-summary(group_tags)/length(group_tags)
  }
  return(list(groupCounts,sig))
}




#' Calculate number of significant genes by chormosome
#'
#' @param res DESeq2 results table
#' @param pval Adjusted p value threshold to use
#' @param lfc Log2 fold change threshold to use
#' @return A table of the number of significant genes by chromosome
#' @export
summaryByChr<-function(resLFC,padj,lfc) {
  up<-resLFC[resLFC$padj < padj & resLFC$log2FoldChange > lfc,]
  down<-resLFC[resLFC$padj < padj & resLFC$log2FoldChange < -lfc, ]
  allChr<-as.data.frame(rbind(up=table(up$chr),down=table(down$chr)))
  allChr$autosomes<-rowSums(allChr[,1:5])
  allChr$total<-rowSums(allChr[,1:6])
  rownames(allChr)<-paste0(rownames(allChr),"_p",padj,"_lfc",lfc)
  return(allChr)
}



#' Plot average signal around motifs in particular sized bins
#'
#' @param motif_gr GRanges object for motifs of interest
#' @param bwFiles List of bwFiles whose signal you want to average
#' @winSize Size of the windows on which to avrage (in bp)
#' @numWins Number of windows either side of the motif to look at
#' @return ggplot2 object
#' @export
avrSignalBins<-function(motif_gr, bwFiles, winSize=10000,numWins=10){
  avrbins<-list()
  for (b in 1:length(bwFiles)){
    gr<-GenomicRanges::resize(motif_gr,width=winSize,fix="center")
    bwdata<-rtracklayer::import.bw(bwFiles[[b]])
    cov<-GenomicRanges::coverage(bwdata,weight="score")
    gr<-GenomicRanges::binnedAverage(gr,cov,paste0(names(bwFiles)[b],"__win",0))
    #gr$numGenes<-countOverlaps(gr,bwdata)
    # #extract group name to count genes per window
    # grp<-gsub(paste0(".*\\/",fileNamePrefix),"",bwFiles[[b]])
    # grp<-gsub("_lfc\\.bw","",grp)
    # salmon<-GRanges(readRDS(paste0(outPath,"/rds/",fileNamePrefix,
    #                                      contrastNames[[grp]],
    #                                "_DESeq2_fullResults_p",padjVal,".rds")))
    # gr$expressedGenesInBin<-countOverlaps(gr,salmon[!is.na(salmon$padj)],ignore.strand=T,type="any")
    upstream<-GenomicRanges::resize(motif_gr,width=winSize,fix="center")
    downstream<-GenomicRanges::resize(motif_gr,width=winSize,fix="center")
    for(i in 1:numWins){
      print(i)
      upstream<-GenomicRanges::flank(upstream,width=winSize,start=T)
      upstream<-GenomicRanges::binnedAverage(upstream,cov,paste0(names(bwFiles)[b],"__win-",i))
      #upstream$expressedGenesInBin<-countOverlaps(upstream,salmon[!is.na(salmon$padj)],ignore.strand=T,type="any")
      downstream<-GenomicRanges::flank(downstream,width=winSize,start=F)
      downstream<-GenomicRanges::binnedAverage(downstream,cov,paste0(names(bwFiles)[b],"__win",i))
      #downstream$expressedGenesInBin<-countOverlaps(downstream,salmon[!is.na(salmon$padj)],ignore.strand=T,type="any")
      df<-cbind(data.frame(gr),mcols(upstream),mcols(downstream))
      df<-tidyr::pivot_longer(df,cols=colnames(df)[grep("__win",colnames(df))],names_to="window")
      df$SMC<-do.call(rbind,strsplit(df$window,split="__win"))[,1]
      df$window<-as.numeric(do.call(rbind,strsplit(df$window,split="__win"))[,2])*winSize/1000
      avrbins[[names(bwFiles)[b]]]<-df
    }
  }

  allavrbins<-do.call(rbind,avrbins)
  allavrbins$window<-factor(allavrbins$window,levels=c(-numWins:numWins)*winSize/1000)
  p<-ggplot2::ggplot(allavrbins,ggplot2::aes(x=window,y=value,col=SMC)) +
    ggplot2::facet_grid(SMC~.,space="free_y",shrink=T)+
    ggplot2::ylim(quantile(allavrbins$value,c(0.01,0.99)))+
    ggplot2::geom_boxplot(outlier.shape=NA,col="black",notch=T) +
    ggplot2::geom_jitter(size=0.5,alpha=0.4) + ggplot2::xlab("Window (kb)") +
    ggplot2::theme_bw()+
    ggplot2::ylab("Average score per bin") +
    ggplot2::xlab("Relative distance (kb)")
  #qc plots to be sure number of genes is not strongly different in bins
  # p1<-ggplot2::ggplot(allavrbins[allavrbins$SMC=="dpy26",],
  #                     ggplot2::aes(x=window,y=expressedGenesInBin,
  #                                             col=value)) +
  #   ggplot2::facet_grid(SMC~.,space="free_y",shrink=T)+
  #   ggplot2::geom_jitter()
  # p2<-ggplot2::ggplot(allavrbins[allavrbins$SMC=="dpy26",],
  #                     ggplot2::aes(x=expressedGenesInBin,y=value,
  #                                  col=window)) +
  #   ggplot2::facet_grid(SMC~.,space="free_y",shrink=T)+
  #   ggplot2::geom_jitter()
  return(p)
}



#' Enhanced volcano plot highlighting X chr genes
#'
#' @param resLFC results object from DESeq2
#' @param lfcVal  threshold Log2 fold change for significance (default=0.5)
#' @param padj threshold adjust p value for significance (default=0.05)
#' @param addLegend Should legend be added (default: F)
#' @return Volcano plot object
#' @export
plotVolcanoXvA<-function(resLFC,lfcVal=0.5,padj=0.05, addLegend=F){
  resLFC<-resLFC[!is.na(resLFC$padj),]
  resByChr<-resLFC[order(resLFC$chr),]
  # create custom key-value pairs for 'low', 'chrX', 'autosome' expression by fold-change
  # set the base colour as 'black'
  keyvals <- rep('black', nrow(resByChr))
  # set the base name/label as 'NS'
  names(keyvals) <- rep('NS', nrow(resByChr))

  keyvals[which(resByChr$chr=="chrX")] <- 'red2'
  names(keyvals)[which(resByChr$chr=="chrX")] <- "chrX LFC>0.5"
  nonSigXchr<-which(resByChr$chr=="chrX" & (resByChr$padj>0.05 | abs(resByChr$log2FoldChange)<=0.5))
  keyvals[nonSigXchr]<-"#c3909b"
  names(keyvals)[nonSigXchr]<-"chrX LFC\u22640.5"

  # modify keyvals for variables with fold change < -2.5
  keyvals[which(resByChr$chr!="chrX")] <- 'royalblue'
  names(keyvals)[which(resByChr$chr!="chrX")] <-"Autosomal LFC>0.5"
  nonSigAchr<-which(resByChr$chr!="chrX" & (resByChr$padj>0.05 | abs(resByChr$log2FoldChange)<=0.5))
  keyvals[nonSigAchr]<-"#6b8ba4"
  names(keyvals)[nonSigAchr]<-"Autosomal LFC\u22640.5"

  sigUp<-sum(resByChr$padj<padjVal & resByChr$log2FoldChange>lfcVal,na.rm=T)
  sigDown<-sum(resByChr$padj<padjVal & resByChr$log2FoldChange< -lfcVal,na.rm=T)
  p<-EnhancedVolcano::EnhancedVolcano(resByChr,
                                      lab=rownames(resByChr),
                                      labSize=0.5,
                                      labCol="#11111100",
                                      x="log2FoldChange",y="padj",
                                      selectLab=rownames(resByChr)[12366],
                                      xlim=c(-5.5,5.5),
                                      ylim=c(0,65),
                                      title= NULL,
                                      titleLabSize = 10,
                                      subtitleLabSize = 8,
                                      subtitle=paste0(sum(!is.na(resLFC$padj)), ' tested genes. ',sigUp, " up, ",sigDown," down."),
                                      caption =NULL,
                                      captionLabSize = 0,
                                      pCutoff=padjVal,
                                      FCcutoff=lfcVal,
                                      xlab=bquote(~Log[2]~'FC'),
                                      ylab=bquote(~-Log[10]~adjusted~italic(P)),
                                      legendPosition =ifelse(addLegend,"left","right"),
                                      legendLabSize = ifelse(addLegend,8,0),#8,
                                      legendIconSize = ifelse(addLegend,3,0),#3.0,
                                      axisLabSize=8,
                                      colCustom=keyvals,
                                      colAlpha=0.5,
                                      pointSize = 1.0)
  if(!addLegend){
     p<-p+theme(legend.position="none")
  }
  return(p)
}
