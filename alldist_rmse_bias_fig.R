# combine distribution dat

MAE<-readRDS("ln2_mae_n104_m105.rds")

MAE %>% mutate(Type=ifelse(grepl("MLE",Method),Method,"Nonparametric"))->MAE

MAE %>% mutate(Estimator=ifelse(Type=="Nonparametric",Method,"MLE"))->MAE

MAE$Distribution="Lognormal"

MAE %>% filter(!(Estimator=="MLE" & !(grepl("LN2",Type)))) %>%
  filter(Estimator=="MLE" | Estimator=="Median of Means" |Estimator=="Sample Mean") -> MAE


MAE2<-readRDS("GPD_n104_m105_mae.rds")

MAE2$Distribution="Generalized Pareto"


MAE2 %>% filter(!(Estimator=="MLE" & !(grepl("GP",Type)))) %>%
  filter(Estimator=="MLE" | Estimator=="Median of Means" |Estimator=="Sample Mean") -> MAE2


###

MAE3<-readRDS("powerlaw_mae_n104_m105.rds")

MAE3$Distribution="Power Law"


MAE3 %>% filter(!(Estimator=="MLE" & !(grepl("Power Law",Type)))) %>%
  filter(Estimator=="MLE" | Estimator=="Median of Means" |Estimator=="Sample Mean") -> MAE3


rbind(MAE,MAE2,MAE3)->r

r$Metric="MAE"

#########

M<-readRDS("ln2_bias_n104_m105.rds")

M$Distribution="Lognormal"

M %>% filter(!(Estimator=="MLE" & !(grepl("LN2",Type)))) %>%
  filter(Estimator=="MLE" | Estimator=="Median of Means"|Estimator=="Sample Mean") -> M

M2<-readRDS("GPD_n104_m105_bias.rds")

M2$Distribution="Generalized Pareto"

M2 %>% filter(!(Estimator=="MLE" & !(grepl("GP",Type)))) %>%
  filter(Estimator=="MLE" | Estimator=="Median of Means" |Estimator=="Sample Mean") -> M2

M3<-readRDS("powerlaw_bias_n104_m105.rds")

M3$Distribution="Power Law"

M3 %>% filter(!(Estimator=="MLE" & !(grepl("Power Law",Type)))) %>%
  filter(Estimator=="MLE" | Estimator=="Median of Means" |Estimator=="Sample Mean") -> M3

rbind(M,M2,M3)->m

m$Metric="Bias"
###


#Lscale
L<-readRDS("ln2_Lscale_n104_m105.rds")

L %>% mutate(Type=ifelse(grepl("MLE",Method),Method,"Nonparametric"))->L

L %>% mutate(Estimator=ifelse(Type=="Nonparametric",Method,"MLE"))->L

L$Distribution="Lognormal"

L %>% filter(!(Estimator=="MLE" & !(grepl("LN2",Type)))) %>%
  filter(Estimator=="MLE" | Estimator=="Median of Means" |Estimator=="Sample Mean") -> L


L2<-readRDS("GPD_n104_m105_Lscale.rds")

L2$Distribution="Generalized Pareto"


L2 %>% filter(!(Estimator=="MLE" & !(grepl("GP",Type)))) %>%
  filter(Estimator=="MLE" | Estimator=="Median of Means" |Estimator=="Sample Mean") -> L2


###

L3<-readRDS("powerlaw_Lscale_n104_m105.rds")

L3$Distribution="Power Law"


L3 %>% filter(!(Estimator=="MLE" & !(grepl("Power Law",Type)))) %>%
  filter(Estimator=="MLE" | Estimator=="Median of Means" |Estimator=="Sample Mean") -> L3


rbind(L,L2,L3)->l

l$Metric="L-scale"

###
m %>% mutate(value=Mean-1) %>% dplyr::select(LCV,value,Estimator,Distribution,Metric)->m

r %>% rename(value=MAE)%>% dplyr::select(LCV,value,Estimator,Distribution,Metric)->r

l %>% rename(value=Lscale)%>% dplyr::select(LCV,value,Estimator,Distribution,Metric)->l


rbind(m,r,l)->dat





ggplot(dat,aes(LCV,value,col=Estimator))+geom_line(size=1)+facet_grid(Distribution~Metric,scales="free")+theme_bw()+
  scale_color_manual(values=c("black","blue","red"))+ylab("Fraction of True Mean")

ggsave("bias_MAE_lscale_all__n104_m105_free.png",dpi=300)

dat %>% saveRDS("bias_lscale_MAE_all_dists_n104_m105.rds")

#################################
library(ggplot2)
library(dplyr)

dat<-readRDS('bias_lscale_MAE_all_dists.rds')

dat %>% mutate(value=value*100)->dat

dat %>% mutate(Distribution=factor(Distribution,levels=c("Lognormal","Generalized Pareto","Power Law")))->dat

ggplot(dat,aes(LCV,value,col=Estimator))+geom_line(size=1)+facet_grid(Distribution~Metric,scales="free")+theme_bw()+
  scale_color_manual(values=c("black","blue","red"))+ylab("Percent of True Mean")

ggsave("bias_MAE_lscale_all_fixed_percent_free_new.png",dpi=300, height=5,width=7)
