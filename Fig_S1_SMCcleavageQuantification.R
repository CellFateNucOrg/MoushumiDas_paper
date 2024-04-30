library(ggplot2)
library(cowplot)
library(tidyr)
library(dplyr)
library(readxl)
library(ggtext)
library(grid)
library(ggpubr)
library(ggbeeswarm)

projectDir="."
otherDataDir=paste0(projectDir,"/otherData")
finalFigDir=paste0(projectDir,"/finalFigures")
if(!dir.exists(finalFigDir)){
  dir.create(finalFigDir)
}


theme_set(
  theme_classic()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          title=element_text(size=9),
          axis.title.y=ggtext::element_markdown(size=9),
          axis.title.x=ggtext::element_markdown(size=9)
    )
)


#par(family="ArialMT")
#gpar(fontfamily="ArialMT")


# Western blot data --------
tbl<-read_excel(paste0(otherDataDir,"/wbQuantification.xlsx"),sheet=1)
tbl$strain<-factor(tbl$strain,levels=c("PMW1021","PMW1025","PMW1005","PMW1023"))
tbl$kleisin<-tbl$strain
levels(tbl$kleisin)<-c("SCC-1","COH-1","DPY-26","KLE-2")

# Normalise HS to NHS
tblNormToNHS<-tbl %>% pivot_wider(id_cols=c(wb,kleisin),names_from=temp,values_from=normToTub) %>% mutate(normToNHS=HS/NHS)

tblstatsum<-tblNormToNHS %>% group_by(kleisin) %>%
  summarise(count=n(),avr=mean(normToNHS),sdev=sd(normToNHS),
            se=sdev/sqrt(count),
            label=paste0(round(avr,2),"+-",round(se,2),"\nn=",count))

p1<-ggplot(tblNormToNHS,aes(x=kleisin,y=normToNHS,colour=kleisin)) +
  geom_jitter (width=0.2,size=4)+
  scale_color_manual(values=c("red","orange","darkblue","darkgreen")) +
  ggtitle("Fraction remaining of full length protein after cleavage") +
  ylab("Ratio of protein after heatshock / no heatshock") +
  geom_hline(yintercept=1,linetype="dashed") +
  coord_cartesian(ylim=c(0,1.3))+
  theme(text=element_text(size=10),legend.position="none")+
  geom_text(data=tblstatsum,mapping=aes(x=kleisin,y=1.1,label=label),size=3)

p1
# Normalise to NHS coh-1
tblNormToCOH1<-tbl %>% filter(temp==c("NHS")) %>% group_by(wb) %>% mutate(COH1_NHS=normToTub[kleisin=="COH-1"],normToCOH1=normToTub/COH1_NHS)

tblstatsum<-tblNormToCOH1 %>% group_by(kleisin) %>%
  summarise(count=n(),avr=mean(normToCOH1),sdev=sd(normToCOH1),
            se=sdev/sqrt(count),
            label=paste0(round(avr,2),"+-",round(se,2),"\nn=",count))

p2<-ggplot(tblNormToCOH1,aes(x=kleisin,y=normToCOH1,colour=kleisin)) +
  geom_jitter (width=0.2,size=4)+
  #scale_color_brewer(palette = "Paired")+
  scale_color_manual(values=c("red","orange","darkblue","darkgreen")) +
  ggtitle("Ratio to COH1") +
  ylab("Ratio of protein / COH-1 protein") +
  theme(text=element_text(size=10),legend.position="none")+
  geom_hline(yintercept=1,linetype="dashed") +
  geom_text(data=tblstatsum,mapping=aes(x=kleisin,y=2,label=label),size=3)

p2


# microscopy data ------

microscopyDir<-paste0(otherDataDir,"/microscopy")
dataDirs<-data.frame(dir=c("20230831_941_SMC1GFP_HS","20230907_941-9_SMCGFP1_HS",
                           "20231213_941-9_SMC1GFP_HS"))

dataDirs$date<-unlist(lapply(strsplit(dataDirs$dir,"_"),"[[",1))
min_size=10
max_size=2000
files=NULL
for(i in 1:nrow(dataDirs)){
  tmpdf<-data.frame(file=list.files(paste0(microscopyDir,"/",dataDirs$dir[i],"/csv_",min_size,"_",max_size,"/")))
  tmpdf$file<-paste0(microscopyDir,"/",dataDirs$dir[i],"/csv_",min_size,"_",max_size,"/",tmpdf$file)
  tmpdf$date<-dataDirs$date[i]
  tmpdf$type<-"nuclei"
  if(is.null(files)){
    files<-tmpdf
  } else {
    files<-rbind(files,tmpdf)
  }
}

files$type[grep("_bg",files$file)]<-"notNuclei"
files$type[grep("_head_bg",files$file)]<-"head&background"
files$strain<-unlist(lapply(strsplit(basename(files$file),"_"),"[[",1))
files$HSvNHS<-unlist(lapply(strsplit(basename(files$file),"_"),"[[",2))
files$worm<-unlist(lapply(strsplit(basename(files$file),"[_.]"),"[[",3))
files
nuc=NULL
head=NULL
headbg=NULL
for (f in seq(1,nrow(files),by=3)){
  nuctmp<-read.csv(files$file[f+2])
  bgtmp<-read.csv(files$file[f])
  headtmp<-read.csv(files$file[f+1])
  if(bgtmp$base_name==nuctmp[1,"base_name"]){
    headbgtmp<-headtmp[headtmp$label=="3",]
    headtmp<-headtmp[headtmp$label=="2",]
    nuctmp$intensity_mean_bg<-headbgtmp$intensity_mean
    headtmp$intensity_mean_bg<-headbgtmp$intensity_mean
    headbgtmp$intensity_mean_nonNucBG<-bgtmp$intensity_mean
  }
  headtmp<-cbind(headtmp,files[f+1,c("date","type","strain","HSvNHS","worm")])
  headbgtmp<-cbind(headbgtmp,files[f+1,c("date","type","strain","HSvNHS","worm")])
  nuctmp$date<-files$date[f+2]
  nuctmp$type<-files$type[f+2]
  nuctmp$strain<-files$strain[f+2]
  nuctmp$HSvNHS<-files$HSvNHS[f+2]
  nuctmp$worm<-files$worm[f+2]
  if(is.null(nuc)){
    nuc<-nuctmp
    head<-headtmp
    headbg<-headbgtmp
  } else {
    nuc<-rbind(nuc,nuctmp)
    head<-rbind(head,headtmp)
    headbg<-rbind(headbg,headbgtmp)
  }
}

### calculating per head
head$intensity_meanMinusBg<-head$intensity_mean-head$intensity_mean_bg
head$HSvNHS<-factor(head$HSvNHS,levels=c("NHS","HS"))
counts<-head %>% group_by(strain,HSvNHS) %>% summarise(N=n())
counts$y=10

head<-head[-33,]

counts
stat_box_data <- function(y) {
  df<-data.frame(y = 0.5+1.1*max(y), label =  length(y))
  return(df)
}

# function for number of observations
give.n <- function(x){
  return(c(y = 18, label = length(x)))
}



p3<-ggplot(head,aes(x=HSvNHS,y=intensity_meanMinusBg,color=HSvNHS)) +
  facet_wrap(~strain,nrow=1)+
  geom_beeswarm(size=1.5,alpha=0.8,corral="wrap",corral.width=4)+
  geom_boxplot(notch=F,outlier.shape = NA,fill=NA)+
  scale_color_manual(values=c("blue","red"))+
  theme(axis.text.x = element_text(size=12), legend.position="none") +
  coord_cartesian(ylim=c(0,20)) +
  ggtitle(paste0("Mean intensity per worm head"))+
  xlab("") + ylab("Mean Intensity (A.U.)")+
  stat_summary(fun.data = give.n, geom = "text", fun = 20)
p3




p<-ggpubr::ggarrange(ggpubr::ggarrange(p2,p1,nrow=1,ncol=2, labels=c("a ","b ")),
                     p3,nrow=2,ncol=1,labels=c("","d "))

p
p<-annotate_figure(p, top = text_grob("Das et al., Figure S1",
                                      size = 12,just=0,face="bold"))


p
ggsave(paste0(finalFigDir,"/Fig_S1_SMCcleavageQuantification.pdf"),plot=p,device="pdf",width=19,height=19,
       unit="cm")

