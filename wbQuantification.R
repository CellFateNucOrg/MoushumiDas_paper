library(ggplot2)
library(cowplot)
library(tidyr)
library(dplyr)
library(readxl)
library(ggtext)
library(grid)
library(ggpubr)

par(family="ArialMT")
gpar(fontfamily="ArialMT")

tbl<-read_excel("./wbQuantification.xlsx",sheet=1)
tbl$strain<-factor(tbl$strain,levels=c("PMW1021","PMW1025","PMW1005","PMW1023"))
tbl$kleisin<-tbl$strain
levels(tbl$kleisin)<-c("SCC-1","COH-1","DPY-26","KLE-2")

#Normalise HS to NHS
tblNormToNHS<-tbl %>% pivot_wider(id_cols=c(wb,kleisin),names_from=temp,values_from=normToTub) %>% mutate(normToNHS=HS/NHS)

tblstatsum<-tblNormToNHS %>% group_by(kleisin) %>% summarise(count=n(),avr=mean(normToNHS),sdev=sd(normToNHS),se=sdev/sqrt(count),label=paste0(round(avr,2),"+-",round(se,2),"\nn=",count))

p1<-ggplot(tblNormToNHS,aes(x=kleisin,y=normToNHS,colour=kleisin)) +
  geom_jitter (width=0.2,size=4)+scale_color_brewer(palette = "Paired")+
  ggtitle("Fraction remaining of full length protein after cleavage") +
  ylab("Ratio of protein after heatshock / no heatshock") +
  theme_classic() + geom_hline(yintercept=1,linetype="dashed") +
  coord_cartesian(ylim=c(0,1.3))+
  theme(text=element_text(size=10))+
  geom_text(data=tblstatsum,mapping=aes(x=kleisin,y=1.1,label=label),size=3)


# Normalise to NHS coh-1
tblNormToCOH1<-tbl %>% filter(temp==c("NHS")) %>% group_by(wb) %>% mutate(COH1_NHS=normToTub[kleisin=="COH-1"],normToCOH1=normToTub/COH1_NHS)

tblstatsum<-tblNormToCOH1 %>% group_by(kleisin) %>% summarise(count=n(),avr=mean(normToCOH1),sdev=sd(normToCOH1),se=sdev/sqrt(count),label=paste0(round(avr,2),"+-",round(se,2),"\nn=",count))

p2<-ggplot(tblNormToCOH1,aes(x=kleisin,y=normToCOH1,colour=kleisin)) +
  geom_jitter (width=0.2,size=4)+ scale_color_brewer(palette = "Paired")+
  ggtitle("Ratio to COH1") +
  ylab("Ratio of protein / COH-1 protein") +
  theme_classic() + theme(text=element_text(size=10))+
  geom_hline(yintercept=1,linetype="dashed") +
  geom_text(data=tblstatsum,mapping=aes(x=kleisin,y=2,label=label),size=3)

p2
dflab<-tibble::tibble(kleisin=factor(levels(tbl$kleisin),levels=levels(tbl$kleisin)),
               x=c(0.6,1.6,2.6,3.6),
               y=rep(0.5,4),
               xend=c(1.4,2.4,3.4,4.4),
               yend=rep(0.5,4),
               labels=c("cohesin<sup>SCC-1</sup>",
                        "cohesin<sup>COH-1</sup>",
                        "condensin I/I<sup>DC</sup>",
                        "condensin II"))
p3<-ggplot(dflab,aes(x=kleisin,y=y,colour=kleisin)) +
  coord_cartesian(ylim=c(0,0.5))+ theme_void()+
  scale_color_brewer(palette = "Paired",
                     guide=guide_legend(override.aes = list(color = "white")))+
  theme(legend.text = element_text(color = "white"),
      legend.title = element_text(color = "white"),
      legend.key = element_rect(color = "white"),
      axis.text.y=element_text(color="white")) +
  ylab("kleisin")+
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend),colour="black")+
  aes(x=x+0.4,y=y-0.05,label=labels)+
  geom_richtext(fill = NA, label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt"),size=3 )
p3



p<-ggpubr::ggarrange(NULL,
                     ggarrange(NULL,p1,p2,p3,NULL,nrow=5,ncol=1,
                               heights=c(0.1,1,1,1,0.1),
                               labels=c("","A ","B ")),
                     NULL,ncol=3,widths=c(0.3,1,0.3))

p<-annotate_figure(p, top = text_grob("Das et al., Figure S1",
                                      size = 12,just=0,face="bold"))



ggsave("./Figure_S1_wbQuant.pdf",plot=p,device="pdf",width=19,height=29,
       unit="cm",family="ArialMT")

