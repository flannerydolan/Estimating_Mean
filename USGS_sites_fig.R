library(ggplot2)
library(dplyr)
library(ggallin)
library(tidyr)
library(viridis)

# read in and process all data

d<-readRDS("USGS_processed_no_date_gaps.rds")

d %>% tidyr::gather(key=Distribution,value=PPCC,LNppcc:PLppcc) %>% mutate(Distribution=substr(Distribution,1,2)) %>%
  group_by(site_no) %>% slice(which.max(PPCC)) ->d1

d1 %>% ungroup() %>%
 mutate(`Record Length`=ifelse(n<5000,"1,000 < n < 10,000",ifelse(n<20000,"10,000 < n < 20,000","n > 20,000")))->d1

d1 %>% filter(PPCC>=0.99)->d2

d2 %>% mutate(`Record Length`=factor(`Record Length`,levels=c("1,000 < n < 10,000","10,000 < n < 20,000","n > 20,000")))->d2

d2 %>% mutate(Distribution=ifelse(Distribution=="GP","Generalized Pareto",ifelse(Distribution=="LN","Lognormal","Power Law")),
                Distribution=factor(Distribution,levels=c("Lognormal","Generalized Pareto","Power Law")))->d2

d2 %>% mutate(Sample_All=Sample) -> all # save sample mean of entire record as 'Sample All'

# read in and process date filtered data

dat<-readRDS("USGS_processed_no_date_gaps_1960-2000.rds")

dat %>% tidyr::gather(key=Distribution,value=PPCC,LNppcc:PLppcc) %>% mutate(Distribution=substr(Distribution,1,2)) %>%
  group_by(site_no) %>% slice(which.max(PPCC)) ->dat1

dat1 %>% ungroup() %>% mutate(`Record Length`=ifelse(n<5000,"1,000 < n < 10,000",ifelse(n<20000,"10,000 < n < 20,000","n > 20,000")))->dat1

dat1 %>% filter(PPCC>=0.99)->dat2

dat2 %>% mutate(`Record Length`=factor(`Record Length`,levels=c("1,000 < n < 10,000","10,000 < n < 20,000","n > 20,000")))->dat2

dat2 %>% mutate(Distribution=ifelse(Distribution=="GP","Generalized Pareto",ifelse(Distribution=="LN","Lognormal","Power Law")),
                Distribution=factor(Distribution,levels=c("Lognormal","Generalized Pareto","Power Law")))->dat2


# calculate percent difference between MOM and sample mean of data in date range over the sample mean of the entire record

dat2 %>% left_join(all,by="site_no") %>% mutate(perc_diff=((MOM.x-Sample.x)/Sample_All)*100) -> j




ggplot(j,aes(Lskew.x,perc_diff))+geom_point(alpha=0.4,size=2.5)+xlab("L-Skew")+
  facet_wrap(~Distribution.x)+ylab("Relative Percent Difference")+theme_bw()+
  scale_color_viridis(discrete=F)+theme(axis.text=element_text(size=11),axis.title = element_text(size=14))+
  scale_y_continuous(trans=ggallin::pseudolog10_trans,breaks=c(-100,-10,-1,0),labels=c("-100","-10","-1","0"))

ggsave("1960-2000_nogaps_overMeanAllpts.png",dpi=300)

  
