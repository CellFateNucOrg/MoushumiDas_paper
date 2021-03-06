library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)
library(seqplots)
#library(eulerr)
library(ggplot2)
library(dplyr)

workDir="."
if(!dir.exists(paste0(workDir,"/plots"))){
  dir.create(paste0(workDir,"/plots"))
}

hicFeaturePath<-"/Users/semple/Documents/MeisterLab/otherPeopleProjects/Moushumi/hicFeatures"

anchordf<-data.frame(source=c("eigen382"),
                     file=c(paste0(workDir,"/otherData/382_X.eigs_cis.vecs_37peaks_p0.65_correct.bed")))


##### subsitute for getREF function from seqplots thaht has an unfixed bug.
##### Fix comes form  https://github.com/Przemol/seqplots/issues/58
#' Get reference genome
#'
#' @param genome The filename of FASTA file or genome code for BSgenome
#'
#' @return \code{DNAStringSet}
#'
#' @export
#'
getREF <- function(genome) {
  if( file.exists(file.path(Sys.getenv('root'), 'genomes', genome)) ) {
    REF <- Biostrings::readDNAStringSet( file.path(Sys.getenv('root'), 'genomes', genome) )
    names(REF) <- gsub(' .+', '', names(REF))
  } else {

    GENOMES <- BSgenome::installed.genomes(
      splitNameParts=TRUE)$genome
    if( length(GENOMES) )
      names(GENOMES) <- gsub('^BSgenome.', '', BSgenome::installed.genomes())
    if( !length(GENOMES) ) stop('No genomes installed!')

    pkg <- paste0('BSgenome.', names(GENOMES[GENOMES %in% genome]))[[1]]
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE, quietly=TRUE)
    )
    REF <- get(pkg)
  }
  return(REF)
}

assignInNamespace("getREF",getREF,ns="seqplots")
#####################-

# DCC subunits and kleisins -----

bwFiles<-paste0(hicFeaturePath,"/publicData/",c("DPY27_N2_L3_GSE67650_ce11.bw", # Kramer et al (2015)
 "GSE87741_DPY30_N2_MxEmb_MACS141_e5_average_ce11.bw", # Albritton et al (2017)
 "GSE87741_SDC2_N2_MxEmb_MACS141_e5_average_ce11.bw", # Albritton et al (2017)
 "GSE87741_SDC3_N2_MxEmb_MACS141_e5_average_ce11.bw", # Albritton et al (2017)
 "KLE2_N2_L3_GSE45678_ce11.bw", # Kranz et al. (2013)
 "GSM4293382_SCC1_Q0835_mE16_E_N2_L3_rep1_ce11.bw" # Huang et al. (2022)
))

for(anch in 1:nrow(anchordf)){
  anchors<-import(anchordf$file[anch])
  anchorSource<-anchordf$source[anch]

  anchors<-anchors[seqnames(anchors)=="chrX"]
  seqlevels(anchors)<-seqlevels(BSgenome.Celegans.UCSC.ce11::Celegans)
  seqinfo(anchors)<-seqinfo(BSgenome.Celegans.UCSC.ce11::Celegans)

  flankSize<-100000


  p<-getPlotSetArray(tracks=c(bwFiles),
                     features=c(anchors),
                     refgenome="ce11", bin=flankSize/100,
                     xmin=flankSize,
                     xmax=flankSize, type="af",
                     xanchored=max(width(anchors)))

  #bw<-import(bwFiles[6])

  dd<-plotHeatmap(p,plotz=F)
  shortNames<-c("DPY-27",
                "DPY-30","SDC-2",
                "SDC-3",
                "KLE-2","SCC-1")
  heatmapQuantiles<-sapply(dd$HLST,quantile,c(0.01,0.99 ),na.rm=T)
  minVal<-min(heatmapQuantiles[1,])
  maxVal<-max(heatmapQuantiles[2,])
  pdf(file=paste0(workDir,"/plots/Anchors_",anchorSource,"_ChIPSeq_DCC_seqplots_flank",
                  flankSize/1000,"kb.pdf"),
      height=8,width=11,paper="a4r")
  plotAverage(p,main=paste0("ChrX ",anchorSource," anchors"),
              error.estimates=F,labels=shortNames,plotScale="linear")
  x<-plotHeatmap(p, main=paste0("ChrX ",anchorSource," anchors"), plotScale="linear", sortrows=T,
                 clusters=1L,autoscale=F,zmin=-0.7,zmax=9.1, labels=shortNames,
                 indi=F, sort_mids=T,sort_by=c(T,rep(F,length(bwFiles))),
                 ln.v=F,
                 clspace=c("#FFFFFF", "#FF0000"))
  dev.off()
}

# # get number for anchors vs flank
# names(bwFiles)<-c("DPY-27_L3","DPY-30_emb","SDC-2_emb","SDC-3_emb",
#                   "KLE-2_L3","SCC-1_L3")
# anchors<-import(anchordf$file[1])
# seqlevels(anchors)<-seqlevels(Celegans)
# seqinfo(anchors)<-seqinfo(Celegans)
# upflank<-flank(anchors, width=100000,start=T,both=F)
# downflank<-flank(anchors,width=100000,start=F,both=F)
# for(bwfile in 1:length(bwFiles)){
#   bw<-import.bw(bwFiles[bwfile])
#   cov<-coverage(bw,weight="score")
#   anchors<-binnedAverage(bins=anchors,numvar=cov,varname=names(bwFiles)[bwfile])
#   upflank<-binnedAverage(bins=upflank,numvar=cov,varname=names(bwFiles)[bwfile])
#   downflank<-binnedAverage(bins=downflank,numvar=cov,varname=names(bwFiles)[bwfile])
# }
#
# upstream<-colMeans(data.frame(mcols(upflank)[,names(bwFiles)]))
# upstream
# # DPY.27_L3 DPY.30_emb  SDC.2_emb  SDC.3_emb   KLE.2_L3   SCC.1_L3
# # 0.93885903 0.04579297 0.20761486 0.37697418 0.04563449 1.10044808
#
# downstream<-colMeans(data.frame(mcols(downflank)[,names(bwFiles)]))
# downstream
# # DPY.27_L3 DPY.30_emb  SDC.2_emb  SDC.3_emb   KLE.2_L3   SCC.1_L3
# # 0.91662235 0.02374211 0.15455688 0.29675338 0.06287584 1.10569648
#
# atanchors<-colMeans(data.frame(mcols(anchors)[,names(bwFiles)]))
# atanchors
# # DPY.27_L3 DPY.30_emb  SDC.2_emb  SDC.3_emb   KLE.2_L3   SCC.1_L3
# # 9.721283   5.464769  13.075006  14.952511   1.405235   2.923475
#
# atanchors/(0.5*(upstream+downstream))
# # DPY.27_L3 DPY.30_emb  SDC.2_emb  SDC.3_emb   KLE.2_L3   SCC.1_L3
# # 10.478448 157.180196  72.203345  44.387412  25.900476   2.650302

