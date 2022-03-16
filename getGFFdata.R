library(rtracklayer)
library(GenomicRanges)

options(timeout = max(600, getOption("timeout")))
gffPath<-"ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/gff/c_elegans.PRJNA13758.WS280.annotations.gff3.gz"
gffFile<-"c_elegans.PRJNA13758.WS280.annotations.gff3"
if(!file.exists(gffFile)){
  download.file(gffPath,destfile=paste0(gffFile,".gz"))
  system(paste0("gunzip ",gffFile,".gz")) # 7.1Gb

  # remove annotations from BLAT and BLAST
  system(paste0("grep -v BLAT ",gffFile," > tmp.gff3")) # down to 2Gb
  system(paste0("grep -v BLASTX tmp.gff3 > ",gffFile)) # down to 1.9Gb

  #remove SNP and variation data
  system(paste0("grep -v Million_mutation ",gffFile," > tmp.gff3")) # 1.7Gb
  system(paste0("grep -v Variation_project_Polymorphism tmp.gff3 > ",gffFile)) #1.4Gb

  system(paste0("grep -v LASTZ ",gffFile," > tmp.gff3")) # 1.4Gb
  system(paste0("grep -v WormbasePaper tmp.gff3 > ",gffFile)) #1.4Gb
}

gff<-import(gffFile,format="gff3")


typesOI<-c("inverted_repeat","gene","CDS","low_complexity_region","repeat_region",
  "tandem_repeat", "TF_binding_site", "TSS_region","intron", "SL1_acceptor_site" ,
  "SL2_acceptor_site","polyA_site", "snoRNA","nc_primary_transcript",
  "DNaseI_hypersensitive_site", "histone_binding_site", "transcription_end_site",
  "mRNA","three_prime_UTR", "transposable_element",
  "transposable_element_insertion_site", "G_quartet", "five_prime_UTR",
  "binding_site", "pseudogenic_transcript","promoter","ncRNA",
  "protein_coding_primary_transcript", "piRNA", "enhancer","operon",
  "circular_ncRNA" ,"duplication" ,"pseudogenic_tRNA",
  "tRNA", "antisense_RNA","miRNA","snRNA", "regulatory_region" ,"rRNA",
  "miRNA_primary_transcript","pseudogenic_rRNA","scRNA")

typesOIwb<-c( "gene","CDS" ,"intron","snoRNA", "nc_primary_transcript","mRNA","three_prime_UTR", "five_prime_UTR" ,"pseudogenic_transcript" , "ncRNA", "piRNA","circular_ncRNA" , "pseudogenic_tRNA", "tRNA", "antisense_RNA" ,"miRNA" , "snRNA",  "rRNA" ,"miRNA_primary_transcript", "pseudogenic_rRNA", "scRNA")

gff<-gff[gff$type %in% typesOI]

dir.create("./tracks")
for(featType in typesOI){
  curFeat<-gff[gff$type==featType]
  naidx<-is.na(curFeat$score)
  curFeat$score[naidx]<-0
  export(curFeat,paste0("./tracks/",featType,".bed"),format="bed")
}

gff1<-gff[gff$source=="WormBase"]
for(featType in typesOIwb){
  curFeat<-gff1[gff1$type==featType & gff1$source=="WormBase"]
  curFeat$name<-curFeat$ID
  naidx<-is.na(curFeat$name)
  curFeat$name[naidx]<-paste0(curFeat$type[naidx],":",curFeat$Parent[naidx])
  naidx<-is.na(curFeat$score)
  curFeat$score[naidx]<-0
  export(curFeat,paste0("./tracks/",featType,".bed"),format="bed")
}

# curFeat<-gff[gff$type=="TF_binding_site"]
#
# curFeat$name<-curFeat$tf_name
# naidx<-all(is.na(curFeat$name) & !is.na(curFeat$Name))
# curFeat[naidx]$name<-gsub(" [modEncode|(].*$","",gsub("ChIP-Seq TF binding region for ","",unlist(lapply(curFeat[naidx]$Name,"[",1))))
# curFeat<-curFeat[curFeat$name %in% names(table(curFeat$name)[table(curFeat$name)>20])]
# curFeat$Name<-NULL
# grl<-split(curFeat,curFeat$name)
# grl<-lapply(grl,reduce)
# gr<-unlist(as(grl,"GRangesList"))
# gr$name<-names(gr)
# gr$score<-0
#
# export(gr,paste0("./tracks/TF_binding_site.bed"),format="bed")
#
# curFeat<-gff[gff$type=="CDS"]
# table(curFeat$source)
