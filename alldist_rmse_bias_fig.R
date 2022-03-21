# combine distribution dat

RMSE<-readRDS("ln2_rmse_n103_m105.rds")

RMSE %>% mutate(Type=ifelse(grepl("MLE",Method),Method,"Nonparametric"))->RMSE

RMSE %>% mutate(Estimator=ifelse(Type=="Nonparametric",Method,"MLE"))->RMSE

RMSE$Distribution="Lognormal"

# Only plot the MLE for that distribution, the MOM and the sample mean

RMSE %>% filter(!(Estimator=="MLE" & !(grepl("LN2",Type)))) %>%
  filter(Estimator=="MLE" | Estimator=="Median of Means" |Estimator=="Sample Mean") -> RMSE


RMSE2<-readRDS("GPD_n104_m105_rmse.rds")

RMSE2$Distribution="Generalized Pareto"


RMSE2 %>% filter(!(Estimator=="MLE" & !(grepl("GP",Type)))) %>%
  filter(Estimator=="MLE" | Estimator=="Median of Means" |Estimator=="Sample Mean") -> RMSE2


###

RMSE3<-readRDS("powerlaw_rmse_n104_m105.rds")

RMSE3$Distribution="Power Law"

# Only plot the MLE for that distribution, the MOM and the sample mean


RMSE3 %>% filter(!(Estimator=="MLE" & !(grepl("Power Law",Type)))) %>%
  filter(Estimator=="MLE" | Estimator=="Median of Means" |Estimator=="Sample Mean") -> RMSE3


rbind(RMSE,RMSE2,RMSE3)->r

r$Metric="RMSE"

#########

M<-readRDS("ln2_bias_n103_m105.rds")

M$Distribution="Lognormal"

# Only plot the MLE for that distribution, the MOM and the sample mean


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

m %>% mutate(value=Mean-1) %>% dplyr::select(LCV,value,Estimator,Distribution,Metric)->m

r %>% rename(value=RMSE)%>% dplyr::select(LCV,value,Estimator,Distribution,Metric)->r


rbind(m,r)->dat





ggplot(dat,aes(LCV,value,col=Estimator))+geom_line(size=1)+facet_grid(Distribution~Metric)+theme_bw()+
  scale_color_manual(values=c("black","blue","red"))+ylab("Fraction of True Mean")

ggsave("bias_rmse_all_fixed.png",dpi=300)

dat %>% saveRDS("bias_and_rmse_all_dists.rds")

#################################
library(ggplot2)
library(dplyr)

dat<-readRDS('bias_and_rmse_all_dists.rds')

dat %>% mutate(value=value*100)->dat

dat %>% mutate(Distribution=factor(Distribution,levels=c("Lognormal","Generalized Pareto","Power Law")))->dat

ggplot(dat,aes(LCV,value,col=Estimator))+geom_line(size=1)+facet_grid(Distribution~Metric,scales="free")+theme_bw()+
  scale_color_manual(values=c("black","blue","red"))+ylab("Percent of True Mean")

ggsave("bias_rmse_all_fixed_percent_free_new.png",dpi=300, height=5,width=6)
