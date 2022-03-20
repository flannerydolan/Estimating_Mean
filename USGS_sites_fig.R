library(ggplot2)
library(dplyr)
library(ggallin)
library(tidyr)
library(viridis)

dat<-readRDS("USGS_processed.rds")

dat %>% tidyr::gather(key=Distribution,value=PPCC,LNppcc:PLppcc) %>% mutate(Distribution=substr(Distribution,1,2)) %>%
  group_by(site_no) %>% slice(which.max(PPCC)) ->dat1

dat1 %>% ungroup() %>% mutate(`Record Length`=ifelse(n<5000,"1,000 < n < 10,000",ifelse(n<20000,"10,000 < n < 20,000","n > 20,000")))->dat1

#dat1 %>% mutate(PPCC=ifelse(PPCC>=0.99,">=0.99","<0.99"))->dat2
dat1 %>% filter(PPCC>=0.99)->dat2

dat2 %>% mutate(`Record Length`=factor(`Record Length`,levels=c("1,000 < n < 10,000","10,000 < n < 20,000","n > 20,000")))->dat2

dat2 %>% mutate(Distribution=ifelse(Distribution=="GP","Generalized Pareto",ifelse(Distribution=="LN","Lognormal","Power Law")),
                Distribution=factor(Distribution,levels=c("Lognormal","Generalized Pareto","Power Law")))->dat2


dat2 %>% filter(Distribution=="Lognormal")->ln

dat2 %>% filter(Distribution=="Generalized Pareto")-> gp

dat2 %>% filter(Distribution=="Power Law")->pl

ggplot(dat2,aes(Lskew,perc_err,col=`Record Length`))+geom_point(alpha=0.4,size=2.5)+
  facet_wrap(~Distribution)+ylab("Relative Percent Difference")+theme_bw()+
  scale_color_viridis(discrete=T)+theme(axis.text=element_text(size=11),axis.title = element_text(size=14))+
  scale_y_continuous(trans=ggallin::pseudolog10_trans,breaks=c(-10,0,10,100,1000,10000,100000),labels=c("-10","0","10","100","1,000","10,000","100,000"))



ggsave("usgs_rel_error_n_dist_facet.png",dpi=300)   
