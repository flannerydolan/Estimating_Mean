#################

library(lmom)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

dat=c()
t3=seq(0.1,1,.1)

########################Theoretical Lines#######################################
################################################################################

# Power Law

alpha=(4*t3)/(3*t3-1)
t2pow=1/(2*alpha-3)
t3pow=alpha/(3*alpha-4)
t4pow=(alpha*(2*alpha-1))/((3*alpha-4)*(4*alpha-5))
pow=cbind(t2pow,t3pow,t4pow)
as.data.frame(pow)->pow
names(pow)<-c("t2","t3","t4")
pow$Distribution="Power Law"
pow$Dataset=NA
pow$Type="Distribution"


################################################################################
# Generalized Pareto
# Appendix of Hosking and Wallis

t2gp=0.33299+0.44559*t3+.16641*t3^2+.09111*t3^5-.03625*t3^7


t4gp=0.20196*t3+0.95924*t3^2-0.20096*t3^3+0.04061*t3^4

gp=cbind(t2gp,t3,t4gp)
as.data.frame(gp)->gp

names(gp)<-c("t2","t3","t4")

gp$Distribution="Generalized Pareto"
gp$Dataset=NA
gp$Type="Distribution"

# combine
dat = rbind(pow,gp)

################################################################################
# Lognormal
library(pracma)
library(lamW)

t2ln=1.16008*t3-0.05325*t3^2-0.10501*t3^4-0.00103*t3^6
t4ln=0.12282+0.77518*t3^2+0.12279*t3^4-0.13638*t3^6+0.11368*t3^8


ln=cbind(t2ln,t3,t4ln)
ln=as.data.frame(ln)
names(ln)<-c("t2","t3","t4")
ln$Distribution="Lognormal"
ln$Dataset=NA
ln$Type="Distribution"

# combine
dat = rbind(dat,ln)

################### Datasets ###################################################
################################################################################

dat=as_tibble(dat)

aid<-readRDS('datasets/foreignaid.rds')

Lcur<-as.numeric(lmom::samlmu(aid,nmom=4))
Lcurdat<-c(Lcur[2]/Lcur[1],Lcur[3],Lcur[4],NA,"US Foreign Aid to Low Income Countries (n=664,568)","Dataset")
dat=rbind(Lcurdat,dat)

terr<-read.table('datasets/terrorism_deaths.txt')
terr=as.numeric(unlist(terr))
terr[terr<=0]<-NA
Lterr<-as.numeric(lmom::samlmu(terr,nmom=4))
Lterrdat<-c(Lterr[2]/Lterr[1],Lterr[3],Lterr[4],NA,"Deaths in terrorist attacks 1968-2006 (n=9,101)","Dataset")

dat=rbind(Lterrdat,dat)


crop<-readRDS('datasets/croploss_2017.rds')

croploss=crop$INDEMNITY
croploss[croploss<=0]<-NA

Lcrop<-as.numeric(lmom::samlmu(croploss,nmom=4))
Lcropdat<-c(Lcrop[2]/Lcrop[1],Lcrop[3],Lcrop[4],NA,"US Crop Indemnities per county 2017 (n=115,031)","Dataset")

dat=rbind(Lcropdat,dat)


usgs<-read.table('datasets/USGS/sites/1166.txt')

flow<-usgs$X_00060_00003
flow[flow<=0]<-NA
Lflow<-as.numeric(lmom::samlmu(flow))
Lflowdat<-c(Lflow[2]/Lflow[1],Lflow[3],Lflow[4],NA,"Daily streamflow at USGS site 2357000 (n=27,495)","Dataset")

dat=rbind(Lflowdat,dat)

library(readxl)
honey<-read_excel('datasets/honeycreekdata.xlsx')

honey$`SS, mg/L (suspended solids)`[honey$`SS, mg/L (suspended solids)`==0]<-NA


Lsed<-as.numeric(lmom::samlmu(honey$`SS, mg/L (suspended solids)`))
Lsed<-c(Lsed[2]/Lsed[1],Lsed[3],Lsed[4],NA,"Daily suspended solids in Honey Creek (mg/L) 1976-2019 (n=22,636)","Dataset")
dat=rbind(Lsed,dat)

cities<-read.table('datasets/cities.txt')

cities=as.numeric(unlist(cities))
cities[cities==0]<-NA
Lcities<-as.numeric(lmom::samlmu(cities,nmom=4))
Lcitiesdat<-c(Lcities[2]/Lcities[1],Lcities[3],Lcities[4],NA,"Population of US cities in 2000 (n=19,447)","Dataset")

dat=rbind(Lcitiesdat,dat)

flare<-read.table('datasets/flares.txt')

flare=as.numeric(unlist(flare))
flare[flare==0]<-NA
Lflare<-as.numeric(lmom::samlmu(flare,nmom=4))
Lflaredat<-c(Lflare[2]/Lflare[1],Lflare[3],Lflare[4],NA,"Gamma ray intensity of solar flares 1980-1989 (n=12,773)","Dataset")

dat=rbind(Lflaredat,dat)

fires<-read.table('datasets/fires.txt')

fires=as.numeric(unlist(fires))
fires[fires==0]<-NA
Lfires<-as.numeric(lmom::samlmu(fires,nmom=4))
Lfiresdat<-c(Lfires[2]/Lfires[1],Lfires[3],Lfires[4],NA,"Size of wildfires on US federal land 1986-1996 (acres) (n=203,785)","Dataset")

dat=rbind(Lfiresdat,dat)

dat$t2=as.numeric(dat$t2)
dat$t3=as.numeric(dat$t3)
dat$t4=as.numeric(dat$t4)




################ Plot ##########################################################
################################################################################



p1<-ggplot(dat)+geom_line(data=dat %>%filter(Type=="Distribution"),size=1,aes(t3,t2,linetype=Distribution,group=Distribution))+theme_bw()+
  geom_point(data=dat %>% filter(Type=="Dataset"),aes(t3,t2,col=Dataset),size=4)+scale_color_brewer(palette = "Dark2",labels=function(x) str_wrap(x, width=50))+
  ylab("L-CV")+xlab("L-Skewness")+theme(axis.title = element_text(size=14),legend.text = element_text(margin = margin(t = 10)))+
  guides(color="none",linetype="none")+coord_cartesian(ylim=c(0.5,1),xlim=c(0.5,1))


p2<-ggplot(dat)+geom_line(data=dat %>%filter(Type=="Distribution"),size=1,aes(t3,t4,linetype=Distribution,group=Distribution))+theme_bw()+
  geom_point(data=dat %>% filter(Type=="Dataset"),aes(t3,t4,col=Dataset),size=4)+scale_color_brewer(palette = "Dark2",labels=function(x) str_wrap(x, width=50))+
  ylab("L-kurtosis")+xlab("L-Skewness")+coord_cartesian(xlim=c(0.5,1),ylim=c(0.25,1))+
  theme(axis.title = element_text(size=14),legend.text = element_text(margin = margin(t = 10)))


library(cowplot)

cowplot::plot_grid(p1,p2,labels=c("a)","b)"),rel_widths = c(1,2.25))

ggsave("fig2_Lmomentdiagramscrop_usgs.png",dpi=300,units="in",height=4,width=9.5)
