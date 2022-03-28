library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)
library(seqplots)
#library(eulerr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(gridExtra)
library(gtable)

workDir="."
if(!dir.exists(paste0(workDir,"/plots"))){
  dir.create(paste0(workDir,"/plots"))
}


############################-
# overlap of states with digitized compartment---------
############################-

states<-import.bed("./publicData/chromStates_L3_Evans2016_ce11.bed")
#domains<-import.bed("./publicData/chromDomains_L3_Evans2016_ce11.bed")
#seqlevels(domains)<-seqlevels(states)

stateClrs<-c("#fe0003","#f59745","#008100","#74943a",
             "#c4d69c","#05ff7f","#ceff65","#fd0082",
             "#ff70d0","#ffb6c4","#7f00ff","#1900fe",
             "#528dd4","#814006","#b8cde4","#808080",
             "#dcdac3","#c4bc97","#938b54","#141414")

getStateOLtable<-function(gr,states){
  ol<-data.frame(findOverlaps(gr,states,ignore.strand=T,minoverlap=10))
  gr$name<-paste0(seqnames(gr),":",start(gr),"-",end(gr))
  ol$name<-gr$name[ol$queryHits]
  ol$XvA<-ifelse(unlist(strsplit(ol$name,":.*$"))=="chrX","X","A")
  ol$state<-factor(states$score[ol$subjectHits],levels=1:20)
  ol$width<-width(states)[ol$subjectHits]
  df<-ol%>%group_by(state,XvA) %>% summarise(stateFrequency=n(),stateWidth=sum(width))
  allStates<-data.frame(state=factor(1:20, levels=1:20))
  df<-left_join(allStates,df)
  df[is.na(df)]<-0
  return(df)
}

ws235toCe11<-function(gr){
  seqlevels(gr)<-c("I","II","III","IV","V","X","MtDNA")
  seqlevels(gr)<-seqlevels(BSgenome.Celegans.UCSC.ce11::Celegans)
  seqinfo(gr)<-seqinfo(BSgenome.Celegans.UCSC.ce11::Celegans)
  return(gr)
}

####################-
# digitized switching --------
###################-


#' Combine Granges of bins and score switching from High/Medium/Low
#'
#' Two GRanges with the same digitzied eigen vector bins grouped into three quantiles
#' numbered 1-3 with 1 being low and 3 being high
#' These will combined and the switching between high/medium/low between the reference
#' GRanges and the pca GRanges will be scored and stored in a column called "switch"
#' @param pca GRanges of pca for which you want to score the change
#' @param ref GRanges of reference pca from which categories are changing
#' @param namePCA name of pca group. pca GRanges should have a column named "bin.<namePCA>"
#' @param nameREF name of reference group. GRanges should have a column named "bin.<nameREF>"
#' @result Dataframe with switch score as well as original bin and eigen scores.
#' @export
scoreSwitch<-function(pca,ref,namePCA,nameREF){
  pca<-full_join(data.frame(ref),data.frame(pca),
                 by=c("seqnames","start","end","width","strand","pca"),
                 suffix=c(paste0(".",nameREF),paste0(".",namePCA)))
  pca$switch<-factor(paste(pca[,paste0("bin.",nameREF)],pca[,paste0("bin.",namePCA)],sep="."),
                     levels=c("3.3","2.3","3.2","2.2","1.2","2.1","1.1","3.1","1.3"))
  levels(pca$switch)<-c("H->H","M->H","H->M","M->M","L->M","M->L","L->L","H->L","L->H")
  return(pca)
}

#' Get average score of chromatin makrs in bins
#' @param bins A GRanges objec with bins for which to get average score
#' @param bwFiles A data frame with two columns, "filePaths" with complete path
#' to bigwig file, and "shortNames" with short names for each dataset.
#' @return A GRanges object with mcols containing the binned average score for
#' each bigwig
#' @export
getChromMarks<-function(bins,bwFiles){
  for(mark in bwFiles$shortNames){
    bw<-import.bw(bwFiles$filePaths[bwFiles$shortNames==mark])
    seqlevels(bw)<-seqlevels(BSgenome.Celegans.UCSC.ce11::Celegans)
    bwcov<-coverage(bw,weight="score")
    bins<-binnedAverage(bins,numvar=bwcov,varname=mark,na.rm=T)
  }
  return(bins)
}



#' Count chromatin state occurrences in each bin
#'
#' @param bins GRanges object for bins for which you want to counts state occurences
#' @param states GRanges of chromatin states from Evan et al.
#' @return GRanges same as bins but iwth added mcols for states 1-20 countaining
#' counts of their occurrence
#' @export
countChromStates<-function(bins,states){
  ol<-data.frame(findOverlaps(bins,states,ignore.strand=T,minoverlap=10))
  ol$state<-states$score[ol$subjectHits]
  stateCounts<-ol%>% group_by(queryHits,state) %>% summarise(count=n()) %>%
    tidyr::pivot_wider(names_from="state",values_from="count",names_prefix="state")
  stateCounts[is.na(stateCounts)]<-0
  stateCounts<-stateCounts[,c("queryHits",paste0("state",1:20))]
  df<-mcols(bins)
  df[,paste0("state",1:20)]<-0
  df[stateCounts$queryHits,paste0("state",1:20)]<-stateCounts[,paste0("state",1:20)]
  mcols(bins)<-df
  return(bins)
}


#' Count genomic features extracted from gff that overlap each bin
#'
#' @param bins GRanges object for bins for which you want to count features
#' @param featType vector of feature types
#' @param featName text to append to feature types column name
#' @param featureTypePath Path to folder where feature bed files are stored
#' @return GRanges same as bins but with added mcols for counts
#' @export
countGFF<-function(bins,featType,featName,featureTypePath){
  df<-data.frame(matrix(0,nrow=length(bins),ncol=length(featType)))
  colnames(df)<-paste0(featName,"_",featType)
  for(feat in featType){
    featgr<-import.bed(paste0(featureTypePath,"/",feat,".bed"))
    seqlevelsStyle(featgr)<-"ucsc"
    seqlevels(featgr)<-seqlevels(BSgenome.Celegans.UCSC.ce11::Celegans)
    df[,paste0(featName,"_",feat)]<-countOverlaps(bins,featgr,ignore.strand=T,minoverlap=10)
  }
  mcols(bins)<-cbind(mcols(bins),df)
  return(bins)
}


#' Get max or median length of genomic features extracted from gff that overlap each bin
#'
#' @param bins GRanges object for bins for which you want to get max/median length of feature
#' @param featType vector of feature types
#' @param featureTypePath Path to folder where feature bed files are stored
#' @return GRanges same as bins but with added mcols for counts
#' @export
getFeatureLength<-function(bins,featType,featureTypePath){
  df<-data.frame(matrix(0,nrow=length(bins),ncol=length(featType)*2))
  colnames(df)<-c(paste0("maxLength_",featType),paste0("medianLength_",featType))
  for(feat in featType){
    featgr<-import.bed(paste0(featureTypePath,"/",feat,".bed"))
    seqlevelsStyle(featgr)<-"ucsc"
    seqlevels(featgr)<-seqlevels(BSgenome.Celegans.UCSC.ce11::Celegans)
    ol<-data.frame(findOverlaps(bins,featgr,ignore.strand=T,minoverlap=10))
    ol$width<-width(featgr[ol$subjectHits])
    ol<-ol %>% group_by(queryHits) %>% summarise(maxLength=max(width),medianLength=median(width))
    df[ol$queryHits,paste0("maxLength_",feat)]<-ol$maxLength
    df[ol$queryHits,paste0("medianLength_",feat)]<-ol$medianLength
  }
  mcols(bins)<-cbind(mcols(bins),df)
  return(bins)
}



#' Create bed files with colors for bins according to 3 quantiles of the eigen vector
#' @param pca GRanges of bins with a column named "bin" with quantile number
#' @param namePCA name for this data set and eigenvector (e.g. "dpy26_E1") which
#' is used for output file name.
#' @return writes bedfile to disk
#' @export
bedClrEigenBin<-function(pca,namePCA){
  binRGB<-apply(col2rgb(c("#0000CC","#66CC66","#CC6600")),2,paste,collapse=",")
  pca2bed<-pca
  pca2bed$score<-pca$bin
  md<-data.frame(name=c("L","M","H")[pca$bin], score=pca2bed$score,
                 strand=".",
                 thickStart=start(pca2bed), thickEnd=end(pca2bed),
                 itemRGB=binRGB[pca2bed$score], blockCount="1",
                 blockSizes=width(pca2bed), blockStarts=start(pca2bed))
  mcols(pca2bed)<-md
  naidx<-is.na(pca2bed$score)
  pca2bed<-pca2bed[!naidx]
  trackLine<-paste0('track name="bin3q ',namePCA,'" description="3 quantile bins" visibility=1 itemRgb="On"\n')
  export.bed(pca2bed,paste0(workDir,"/otherData/bin3q_",namePCA,"_rgb.bed"))
  pca2bed<-read.delim(paste0(workDir,"/otherData/bin3q_",namePCA,"_rgb.bed"),header=F)
  pca2bed<-cbind(pca2bed,md[!naidx,c("itemRGB")])
  write.table(pca2bed,
              file=paste0(workDir,"/otherData/bin3q_",namePCA,"_rgb.bed"),
              sep="\t",row.names=F,col.names=F,quote=F)
  # add trackline (but not necessary)
  fConn <- file(paste0(workDir,"/otherData/bin3q_",namePCA,"_rgb.bed"), 'r+')
  Lines <- readLines(fConn)
  writeLines(c(trackLine, Lines), con = fConn)
  close(fConn)
}



#' Create bed files with colors for bins that switch up or down
#' @param pca GRanges of bins with a column named "switch" with switch type
#' @param namePCA name for this data set and eigenvector (e.g. "dpy26_E1") which
#' is used for output file name.
#' @return writes bedfile to disk
#' @export
bedClrSwitch<-function(pca,namePCA){
  #switchRGB<-apply(col2rgb(c("#333366","#FFFF66","#FF00FF")),2,paste,collapse=",")
  switchRGB<-apply(col2rgb(c("#333366","#FF9999","#FF0000","#660000","#66FFFF","#00CCCC","#003333")),2,paste,collapse=",")
  pca2bed<-pca
  pca2bed$score<-NA
  pca2bed$score[pca2bed$switch %in% c("H->H","M->M","L->L")]<-1
  pca2bed$score[pca2bed$switch %in% c("L->M")]<-2
  pca2bed$score[pca2bed$switch %in% c("M->H")]<-3
  pca2bed$score[pca2bed$switch %in% c("L->H")]<-4
  pca2bed$score[pca2bed$switch %in% c("H->M")]<-5
  pca2bed$score[pca2bed$switch %in% c("M->L")]<-6
  pca2bed$score[pca2bed$switch %in% c("H->L")]<-7
  #pca2bed$score[pca2bed$switch %in% c("H->M","H->L","M->L")]<-2
  #pca2bed$score[pca2bed$switch %in% c("M->H","L->H","L->M")]<-3
  md<-data.frame(name=pca$switch, score=pca2bed$score,
                 strand=".",
                 thickStart=start(pca2bed), thickEnd=end(pca2bed),
                 itemRGB=switchRGB[pca2bed$score], blockCount="1",
                 blockSizes=width(pca2bed), blockStarts=start(pca2bed))
  mcols(pca2bed)<-md
  naidx<-is.na(pca2bed$score)
  pca2bed<-pca2bed[!naidx]
  trackLine<-paste0('track name="switch ',namePCA,'" description="Switch bin" visibility=1 itemRgb="On"\n')
  export.bed(pca2bed,paste0(workDir,"/otherData/switchBin_",namePCA,"_rgb.bed"))
  pca2bed<-read.delim(paste0(workDir,"/otherData/switchBin_",namePCA,"_rgb.bed"),header=F)
  pca2bed<-cbind(pca2bed,md[!naidx,c("itemRGB")])
  write.table(pca2bed,
              file=paste0(workDir,"/otherData/switchBin_",namePCA,"_rgb.bed"),
              sep="\t",row.names=F,col.names=F,quote=F)
  # add trackline (but not necessary)
  fConn <- file(paste0(workDir,"/otherData/switchBin_",namePCA,"_rgb.bed"), 'r+')
  Lines <- readLines(fConn)
  writeLines(c(trackLine, Lines), con = fConn)
  close(fConn)
}

#' Create complex heatmap of eigen vectors chormatin marks and states
#'
#' @param pca data.frame with various types of data in it
#' @param pcaName Name of sample and eigen vector which is used in output file name
#' @param doRaster whether to rasterize very large matrices (default=T)
#' @return Saves plot to disk
#' @export
doEigenComplexHeatmap<-function(pca, pcaName, doRaster=T){
  unsortedPCA<-pca
  pca<-unsortedPCA %>% group_by(switch) %>% arrange(desc(score.TEVonly))
  scoreCols<-grep("score",colnames(pca))
  tpmCol<-grep("^tpm",colnames(pca))
  lfcCol<-grep("LFC_",colnames(pca))
  chromMarkCols<-grep("H3K",colnames(pca))
  stateCols<-grep("state",colnames(pca))
  naidx<-rowSums(is.na(pca[,c(scoreCols,tpmCol,lfcCol)]))>0
  colnames(pca)<-gsub("^score\\.","",colnames(pca))
  hm1<-Heatmap(t(as.matrix(pca[!naidx,c(scoreCols)])),column_split=pca$switch[!naidx],
               cluster_rows = F,cluster_columns = F, show_row_dend = F, show_column_dend=F,
               col = circlize::colorRamp2(c(-0.6, 0, 0.6), c("darkblue", "yellow", "darkred")),
               column_gap = unit(c(3), "mm"),  border_gp = gpar(col = "black", lty = 1,size=2),
               show_column_names=F, row_names_side="left",use_raster=doRaster,
               raster_quality=5, column_title_rot = 90, #row_labels=list("TEVonly",substitute(italic(x^AID),list(x="sdc-3"))),
               heatmap_legend_param = list(title = "Eigen value"))

  hm2<-Heatmap(t(log2(as.matrix(pca[!naidx,c(tpmCol)])+1)),column_split=pca$switch[!naidx],
               cluster_rows = F, cluster_columns = F, show_row_dend = F, show_column_dend=F,
               col = circlize::colorRamp2(c(0, 5), c("white","darkblue")),
               column_gap = unit(c(3), "mm"), border_gp = gpar(col = "black", lty = 1,size=2),
               show_column_names=F, row_names_side="left",use_raster=doRaster,
               raster_quality=5, column_title_rot = 90,
               show_row_names=T, row_labels=c("PMW366_TPM"),
               heatmap_legend_param = list(title = "log2(TPM+1)"))

  hm3<-Heatmap(t(as.matrix(pca[!naidx,c(lfcCol)])),column_split=pca$switch[!naidx],
               cluster_rows = F, cluster_columns = F, show_row_dend = F, show_column_dend=F,
               col = circlize::colorRamp2(c(-0.5, 0, 0.5), c("darkorange4","white","darkgreen")),
               column_gap = unit(c(3), "mm"), border_gp = gpar(col = "black", lty = 1,size=2),
               show_column_names=F, row_names_side="left",use_raster=doRaster,
               raster_quality=5, column_title_rot = 90,
               show_row_names=T, row_labels=names(pca)[lfcCol],
               heatmap_legend_param=list(title="log2(FC)"))

  hm4<-Heatmap(t(as.matrix(pca[!naidx,c(chromMarkCols)])),column_split=pca$switch[!naidx],
               cluster_rows = F, cluster_columns = F, show_row_dend = F, show_column_dend=F,
               col = circlize::colorRamp2(c(0, 3), c("white","darkred")),
               column_gap = unit(c(3), "mm"), border_gp = gpar(col = "black", lty = 1,size=2),
               show_column_names=F, row_names_side="left",use_raster=doRaster,
               raster_quality=5, column_title_rot = 90,
               heatmap_legend_param = list(title = "Enrichment"))

  hm5<-Heatmap(t(as.matrix(pca[!naidx,c(stateCols)])),column_split=pca$switch[!naidx],
               cluster_rows = F, cluster_columns = F, show_row_dend = F, show_column_dend=F,
               col = circlize::colorRamp2(c(0, 20), c("white","black")),
               column_gap = unit(c(3), "mm"),  border_gp = gpar(col = "black", lty = 1,size=2),
               show_column_names=F, row_names_side="left",use_raster=doRaster,
               raster_quality=5, column_title_rot = 90,
               heatmap_legend_param = list(title = "Count"))

  hmlist=hm1 %v% hm2 %v% hm3 %v% hm4 %v% hm5
  #pdf(file=paste0(workDir,"/plots/eigenHeatmap_",pcaName,".pdf"),width=10,height=10,paper="a4")
  #draw(hmlist,gap = unit(c(3), "mm"))
  #dev.off()
  return(hmlist)
}



#' Create complex heatmap of eigen vectors chormatin marks and states
#'
#' @param pca data.frame with various types of data in it
#' @param pcaName Name of sample and eigen vector which is used in output file name
#' @param doRaster whether to rasterize very large matrices (default=T)
#' @return Saves plot to disk
#' @export
doEigenComplexHeatmap1<-function(pca, pcaName, doRaster=T){
  unsortedPCA<-pca
  pca<-unsortedPCA %>% group_by(switch) %>% arrange(desc(score.TEVonly))
  scoreCols<-grep("score",colnames(pca))
  tpmCol<-grep("^tpm",colnames(pca))
  lfcCol<-grep("LFC_",colnames(pca))
  chromMarkCols<-grep("H3K",colnames(pca))
  lengthCols<-grep("Length_",colnames(pca))
  featureCols1<-grep("feature1_",colnames(pca))
  featureCols2<-grep("feature2_",colnames(pca))
  featureCols3<-grep("feature3_",colnames(pca))
  naidx<-rowSums(is.na(pca[,c(scoreCols,tpmCol,lfcCol)]))>0
  colnames(pca)<-gsub("^score\\.","",colnames(pca))
  hm1<-Heatmap(t(as.matrix(pca[!naidx,c(scoreCols)])),
               column_split=pca$switch[!naidx],
               cluster_rows = F,cluster_columns = F, show_row_dend = F,
               show_column_dend=F,
               col = circlize::colorRamp2(c(-0.6, 0, 0.6),
                                          c("darkblue", "yellow", "darkred")),
               column_gap = unit(c(3), "mm"),
               border_gp = gpar(col = "black", lty = 1,size=2),
               show_column_names=F, row_names_side="left",use_raster=doRaster,
               raster_quality=5, column_title_rot = 90,
               heatmap_legend_param = list(title = "Eigen value"))

  hm2<-Heatmap(t(log2(as.matrix(pca[!naidx,c(tpmCol)])+1)),
               column_split=pca$switch[!naidx],
               cluster_rows = F, cluster_columns = F, show_row_dend = F,
               show_column_dend=F,
               col = circlize::colorRamp2(c(0, 5), c("white","darkblue")),
               column_gap = unit(c(3), "mm"),
               border_gp = gpar(col = "black", lty = 1,size=2),
               show_column_names=F, row_names_side="left",use_raster=doRaster,
               raster_quality=5, column_title_rot = 90,
               show_row_names=T, row_labels=c("PMW366_TPM"),
               heatmap_legend_param = list(title = "log2(TPM+1)"))

  hm3<-Heatmap(t(as.matrix(pca[!naidx,c(lfcCol)])),
               column_split=pca$switch[!naidx],
               cluster_rows = F, cluster_columns = F, show_row_dend = F,
               show_column_dend=F,
               col = circlize::colorRamp2(c(-0.5, 0, 0.5), c("darkorange4","white","darkgreen")),
               column_gap = unit(c(3), "mm"),
               border_gp = gpar(col = "black", lty = 1,size=2),
               show_column_names=F, row_names_side="left",use_raster=doRaster,
               raster_quality=5, column_title_rot = 90,
               show_row_names=T, row_labels=names(pca)[lfcCol],
               heatmap_legend_param=list(title="log2(FC)"))

  hm4<-Heatmap(t(as.matrix(pca[!naidx,c(chromMarkCols)])),
               column_split=pca$switch[!naidx],
               cluster_rows = F, cluster_columns = F, show_row_dend = F,
               show_column_dend=F,
               col = circlize::colorRamp2(c(0, 3), c("white","darkred")),
               column_gap = unit(c(3), "mm"),
               border_gp = gpar(col = "black", lty = 1,size=2),
               show_column_names=F, row_names_side="left",use_raster=doRaster,
               raster_quality=5, column_title_rot = 90,
               heatmap_legend_param = list(title = "Enrichment"))

  hm5<-Heatmap(t(log2(as.matrix(pca[!naidx,c(lengthCols)]))+1),
               column_split=pca$switch[!naidx],
               cluster_rows = F, cluster_columns = F, show_row_dend = F,
               show_column_dend=F,
               col = circlize::colorRamp2(c(0, 25), c("white","black")),
               column_gap = unit(c(3), "mm"),
               border_gp = gpar(col = "black", lty = 1,size=2),
               show_column_names=F, row_names_side="left",
               use_raster=doRaster,
               raster_quality=5, column_title_rot = 90,
               heatmap_legend_param = list(title = "log2(length)"))

  colnames(pca)<-gsub("feature._","",colnames(pca))
  hm6<-Heatmap(t(as.matrix(pca[!naidx,c(featureCols1)])),
               column_split=pca$switch[!naidx],
               cluster_rows = F, cluster_columns = F, show_row_dend = F,
               show_column_dend=F,
               col = circlize::colorRamp2(c(0, 20), c("white","black")),
               column_gap = unit(c(3), "mm"),
               border_gp = gpar(col = "black", lty = 1,size=2),
               show_column_names=F, row_names_side="left",
               use_raster=doRaster,
               raster_quality=5, column_title_rot = 90,
               heatmap_legend_param = list(title = "count1"))

  hm7<-Heatmap(t(as.matrix(pca[!naidx,c(featureCols2)])),
               column_split=pca$switch[!naidx],
               cluster_rows = F, cluster_columns = F, show_row_dend = F,
               show_column_dend=F,
               col = circlize::colorRamp2(c(0, 5), c("white","black")),
               column_gap = unit(c(3), "mm"),
               border_gp = gpar(col = "black", lty = 1,size=2),
               show_column_names=F, row_names_side="left",
               use_raster=doRaster,
               raster_quality=5, column_title_rot = 90,
               heatmap_legend_param = list(title = "count2"))

  hm8<-Heatmap(t(as.matrix(pca[!naidx,c(featureCols3)])),
               column_split=pca$switch[!naidx],
               cluster_rows = F, cluster_columns = F, show_row_dend = F,
               show_column_dend=F,
               col = circlize::colorRamp2(c(0, 0.5), c("white","black")),
               column_gap = unit(c(3), "mm"),
               border_gp = gpar(col = "black", lty = 1,size=2),
               show_column_names=F, row_names_side="left",
               use_raster=doRaster,
               raster_quality=5, column_title_rot = 90,
               heatmap_legend_param = list(title = "count3"))
  hmlist=hm1 %v% hm2 %v% hm3 %v% hm4 %v% hm5 %v% hm6 %v% hm7 %v% hm8
  #pdf(file=paste0(workDir,"/plots/eigenHeatmap1_",pcaName,".pdf"),width=10,height=10,paper="a4")
  #draw(hmlist,gap = unit(c(3), "mm"))
  #dev.off()
  return(hmlist)
}
#######-

####################-
## chromatin feature heatmap------
####################-
evansPath="/Users/semple/Documents/MeisterLab/Datasets/Evans2016_PNAS_L3"
bwFiles<-data.frame(filePaths=paste0(evansPath,"/",
                                     c("GSM624432_WA30634849_H3K27AC_N2_L3_1_ce11.bw",
                                       "GSM624433_WA30634849_H3K27AC_N2_L3_2_ce11.bw",
                                       "GSM562734_HK00013_H3K27ME31E7_N2_L3_1_ce11.bw",
                                       "GSM562735_HK00013_H3K27ME31E7_N2_L3_2_ce11.bw",
                                       "GSM562736_HK00001_H3K36ME313C9_N2_L3_1_ce11.bw",
                                       "GSM562737_HK00001_H3K36ME313C9_N2_L3_2_ce11.bw" )),
                    shortNames=paste0(rep(c("H3K27Ac","H3K27me3","H3K36me3"),each=2),
                                      rep(c("_1","_2"),3)))


pcas<-data.frame(SMC=c("TEVonly","dpy26","kle2","scc1","coh1"),
                 strain =c("366","382","775","784","828"),
                 E1=NA, E2=NA)

refPCAs<-c("TEVonly")
otherPCAs<-c("dpy26","kle2","scc1","coh1")
prettyOtherPCAs<-c("dpy-26cs","kle-2cs","scc-1cs","coh-1cs")

#pcaPath<-"/Users/semple/Documents/MeisterLab/otherPeopleProjects/Moushumi/hicFeatures/otherData"
pcaPath=paste0(workDir,"/otherData")
E1files=list.files(paste0(pcaPath),
                   pattern="_10k\\.oriented\\.3quantiles_trans\\.eigs\\.vecs_E1\\.rds")
E2files=list.files(paste0(pcaPath),
                   pattern="_10k\\.oriented\\.3quantiles_trans\\.eigs\\.vecs_E2\\.rds")
pcas$E1<-E1files[match(pcas$strain,unlist(strsplit(E1files,"_10k\\.oriented\\.3quantiles_trans\\.eigs\\.vecs_E1\\.rds")))]
pcas$E2<-E2files[match(pcas$strain,unlist(strsplit(E2files,"_10k\\.oriented\\.3quantiles_trans\\.eigs\\.vecs_E2\\.rds")))]

#rnaseqTPMpath<-"/Users/semple/Documents/MeisterLab/otherPeopleProjects/Moushumi/2021_RNAseq_MDas/tracks"
#rnaseqLFCpath<-"/Users/semple/Documents/MeisterLab/otherPeopleProjects/Moushumi/2021_RNAseq_MDas/tracks/p0.05_lfc0.5_filtCycChrAX"
rnaseqTPMpath<-paste0(workDir,"/tracks")
rnaseqLFCpath<-paste0(workDir,"/tracks/p0.05_lfc0.5_filtCycChrAX")


rnaTPM<-data.frame(filePaths=paste0(rnaseqTPMpath,"/PMW366_TPM_avr.bw"),
                   shortNames="tpm366")
# rnaTPM<-data.frame(filePaths=paste0(workDir,"/otherData/sumFR_366_B_UniqueMultiple.bw"),
#                    shortNames="tpm366")
rnaLFC<-data.frame(filePaths=paste0(rnaseqLFCpath,"/",
                                    "filtCycChrAX_",otherPCAs,"_lfc.bw"),
                   shortNames=paste0("LFC_",prettyOtherPCAs))


#listdf<-NULL
refset=refPCAs[1]
ref1<-readRDS(paste0(pcaPath,"/",pcas$E1[pcas$SMC==refset]))
#bedClrEigenBin(ref1,"TEVonly_E1")
ref2<-readRDS(paste0(pcaPath,"/",pcas$E2[pcas$SMC==refset]))
#bedClrEigenBin(ref2,"TEVonly_E2")

bins<-readRDS(paste0(pcaPath,"/",pcas$E1[pcas$SMC==refset]))
mcols(bins)<-NULL
bins<-getChromMarks(bins,rnaTPM)
bins<-getChromMarks(bins,bwFiles)
bins<-countChromStates(bins,states)

g=1
E1chrom<-list()
E2chrom<-list()
statdf<-NULL
for (g in 1:length(otherPCAs)){
  grp<-otherPCAs[g]
  prettyGrp<-prettyOtherPCAs[g]

  pca1<-readRDS(paste0(pcaPath,"/",pcas$E1[pcas$SMC==grp]))
  pca2<-readRDS(paste0(pcaPath,"/",pcas$E2[pcas$SMC==grp]))

  #bedClrEigenBin(pca1,paste0(grp,"_E1"))
  #bedClrEigenBin(pca2,paste0(grp,"_E2"))

  ref1<-readRDS(paste0(pcaPath,"/",pcas$E1[pcas$SMC==refset]))
  ref2<-readRDS(paste0(pcaPath,"/",pcas$E2[pcas$SMC==refset]))

  E1<-scoreSwitch(pca1,ref1,prettyGrp,refset)
  E2<-scoreSwitch(pca2,ref2,prettyGrp,refset)

  #pca1<-GRanges(E1)
  #bedClrSwitch(pca1,paste0(grp,"_E1"))
  #pca2<-GRanges(E2)
  #bedClrSwitch(pca2,paste0(grp,"_E2"))


  bins1<-getChromMarks(bins,rnaLFC[rnaLFC$shortNames==paste0("LFC_",prettyGrp),])

  E1<-cbind(E1,mcols(bins1))
  E2<-cbind(E2,mcols(bins1))

  #listdf[[paste0(grp,"_E1")]]<-E1
  #listdf[[paste0(grp,"_E3")]]<-E2
  E1chrom[[paste0(prettyGrp,"_E1")]]<-doEigenComplexHeatmap(E1,paste0(prettyGrp,"_E1"), doRaster=T)
  E2chrom[[paste0(prettyGrp,"_E2")]]<-doEigenComplexHeatmap(E2,paste0(prettyGrp,"_E2"), doRaster=T)

  if(is.null(statdf)){
    statdf<-rbind(data.frame(SMC=prettyGrp,EigenVector="E1",table(E1$switch)),
                  data.frame(SMC=prettyGrp,EigenVector="E2",table(E2$switch)))
  } else {
    statdf<-rbind(statdf,
                  data.frame(SMC=prettyGrp,EigenVector="E1",table(E1$switch)),
                  data.frame(SMC=prettyGrp,EigenVector="E2",table(E2$switch)))
  }
}

statdfwide<-statdf %>% pivot_wider(names_from=Var1,values_from=Freq)

write.table(statdfwide,file=paste0(workDir,"/switchStats.tsv"),sep="\t",
            row.names=F)


####################-
## gene feature heatmap -------
####################-
typesOIwb1<-c("binding_site","TF_binding_site",
              "DNaseI_hypersensitive_site","TSS_region","gene","mRNA",
              "five_prime_UTR","CDS","intron","three_prime_UTR","inverted_repeat",
              "repeat_region", "tandem_repeat", "low_complexity_region")

typesOIwb2<-c("enhancer","promoter","pseudogenic_transcript" , "ncRNA", "piRNA",
              "operon","G_quartet","transposable_element","histone_binding_site")

typesOIwb3<-c( "snoRNA", "nc_primary_transcript","circular_ncRNA" , "pseudogenic_tRNA",
               "tRNA", "antisense_RNA" ,"miRNA" , "snRNA",  "rRNA" ,
               "miRNA_primary_transcript", "pseudogenic_rRNA", "scRNA",
               "regulatory_region")

#listdf<-NULL
refset=refPCAs[1]
#featureTypePath<-"/Users/semple/Documents/MeisterLab/otherPeopleProjects/Moushumi/hicFeatures/tracks"
featureTypePath<-paste0(workDir,"/tracks")
bins<-readRDS(paste0(pcaPath,"/",pcas$E1[pcas$SMC==refset]))
mcols(bins)<-NULL
bins<-getChromMarks(bins,rnaTPM)
bins<-getChromMarks(bins,bwFiles[c(1,3,5),])
bins<-getFeatureLength(bins,c("gene"),featureTypePath)
bins<-countGFF(bins,typesOIwb1,"feature1",featureTypePath)
bins<-countGFF(bins,typesOIwb2,"feature2",featureTypePath)
bins<-countGFF(bins,typesOIwb3,"feature3",featureTypePath)

g=2
E1feat<-list()
E2feat<-list()
for (g in 1:length(otherPCAs)){
  grp<-otherPCAs[g]
  prettyGrp<-prettyOtherPCAs[g]

  pca1<-readRDS(paste0(pcaPath,"/",pcas$E1[pcas$SMC==grp]))
  pca2<-readRDS(paste0(pcaPath,"/",pcas$E2[pcas$SMC==grp]))

  ref1<-readRDS(paste0(pcaPath,"/",pcas$E1[pcas$SMC==refset]))
  ref2<-readRDS(paste0(pcaPath,"/",pcas$E2[pcas$SMC==refset]))

  E1<-scoreSwitch(pca1,ref1,prettyGrp,refset)
  E2<-scoreSwitch(pca2,ref2,prettyGrp,refset)

  bins1<-getChromMarks(bins,rnaLFC[rnaLFC$shortNames==paste0("LFC_",prettyGrp),])

  E1<-cbind(E1,mcols(bins1))
  E2<-cbind(E2,mcols(bins1))

  #listdf[[paste0(grp,"_E1")]]<-E1
  #listdf[[paste0(grp,"_E3")]]<-E2
  E1feat[[paste0(prettyGrp,"_E1")]]<-doEigenComplexHeatmap1(E1,paste0(prettyGrp,"_E1"),doRaster=T)
  E2feat[[paste0(prettyGrp,"_E2")]]<-doEigenComplexHeatmap1(E2,paste0(prettyGrp,"_E2"),doRaster=T)
}


pdf(file=paste0(workDir,"/plots/Figure_S4.pdf"),width=10,height=10,
    paper="a4")#,title="Das et al., Figure S4")
tbltitle <- textGrob("Das et al., Figure S4\nCompartment switching",
                     gp=gpar(fontsize=13))
padding1 <- unit(10,"line")
padding2 <- unit(30,"line")
table <- tableGrob(statdfwide)
table <- gtable_add_rows(table,
                         heights = grobHeight(tbltitle) + padding1,
                         pos = 0)
table <- gtable_add_rows(table,
                         heights = padding2,
                         pos = -1)
table <- gtable_add_grob(table, list(tbltitle),
                         t=1, l=1, r=ncol(table))
grid.newpage()
grid.draw(table)
for(i in 1:length(E1chrom)){
  ComplexHeatmap::draw(E1chrom[[i]],gap = unit(c(3), "mm"),
       column_title=paste0("Das et al., Figure S4\n",names(E1chrom)[i],"\n"),
       column_title_gp = gpar(fontsize = 13))
}
for(i in 1:length(E1feat)){
  ComplexHeatmap::draw(E1feat[[i]],gap = unit(c(3), "mm"),
       column_title=paste0("Das et al., Figure S4\n",names(E1feat)[i],"\n"),
  column_title_gp = gpar(fontsize = 13))
}
for(i in 1:length(E2chrom)){
  ComplexHeatmap::draw(E2chrom[[i]],gap = unit(c(3), "mm"),
       column_title=paste0("Das et al., Figure S4\n",names(E2chrom)[i],"\n"),
       column_title_gp = gpar(fontsize = 13))
  }
for(i in 1:length(E2feat)){
  ComplexHeatmap::draw(E2feat[[i]],gap = unit(c(3), "mm"),
       column_title=paste0("Das et al., Figure S4\n",names(E2chrom)[i],"\n"),
       column_title_gp = gpar(fontsize = 13))
}
dev.off()



