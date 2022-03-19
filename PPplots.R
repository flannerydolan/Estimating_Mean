### Make PP plots and PPCC coefficients

# libraries

library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(lmom)
library(readxl)
#################

# read in data

honey<-read_excel('datasets/honeycreekdata.xlsx')

honeyss<-honey$`SS, mg/L (suspended solids)`

honeyss[honeyss<=0]<-NA

###############

usgs<-read.table('datasets/USGS/sites/1166.txt')

flow<-usgs$X_00060_00003
flow[flow<=0]<-NA

###############

cities<-read.table('datasets/cities.txt')

cities=as.numeric(unlist(cities))
cities[cities<=0]<-NA

###############

flare<-read.table('datasets/flares.txt')

flare=as.numeric(unlist(flare))
flare[flare<=0]<-NA

###############

fires<-read.table('datasets/fires.txt')

fires=as.numeric(unlist(fires))
fires[fires<=0]<-NA

###############

terr<-read.table('datasets/terrorism_deaths.txt')

terr=as.numeric(unlist(terr))
terr[terr<=0]<-NA

aid<-readRDS('C:/Users/fdolan03/Box/BootstrapMean/datasets/foreignaid.rds')

###############

cropdat<-read.csv('datasets/RMAPayment19892017SouthwestClimateHub.csv')
cropdat %>% filter(YEAR==2017)->crop17
crop17 %>% saveRDS('datasets/croploss_2017.rds')

crop17$INDEMNITY->croploss
###############################

# For each needed plot, set x = dataset

flare->x
cities->x
fires->x
terr->x
honeyss->x
aid->x
flow->x
croploss->x

##########################################
x[x<=0]<-NA

na.omit(x)->x
x=as.numeric(x)
Lx<-as.numeric(lmom::samlmu(x,nmom=4))



# Generalized Pareto

k=(Lx[1]-2*Lx[2])/Lx[2] ## Hosking and Wallis Appendix

a=Lx[1]*(1+k)

y=seq(1:length(x))/(length(x)+1)

Fgp = 1 - (1-(k*x)/a)^(1/k) # Hosking and Wallis 1987


# Power Law

xmin=min(x,na.rm=T)

if (xmin==0){
  xmin=sort(unique(x))[2]
}

LCV=Lx[2]/Lx[1]

alpha=(1+3*LCV)/(2*LCV)


C = (alpha-1)*xmin^(alpha-1)

Fpl = 1-(C/(alpha-1))*x^-(alpha-1)



# Lognormal


vy=(sqrt(2)*qnorm((1+LCV)/2))^2

uy=log(Lx[1])-vy/2


Fln=plnorm(x, meanlog=uy, sdlog=sqrt(vy))


#

df=as.data.frame(cbind(y,sort(Fgp),sort(Fpl),sort(Fln)))

names(df)<-c("phat","Generalized Pareto","Power Law","Lognormal")

df %>% as_tibble() %>%
  gather(Distribution, value, `Generalized Pareto`:Lognormal) -> df

ggplot(df,aes(value,phat,col=Distribution))+geom_point()+theme_bw()+
  scale_color_viridis(discrete=TRUE)+ylab("Empirical cumulative distribution")+xlab("Theoretical cumulative distribution")


ggsave("figures/SI/croploss_PPplot.png",dpi=300)

# ppcc
cor(sort(y),sort(Fln))
cor(sort(y),sort(Fgp))
cor(sort(y),sort(Fpl))

cor(y,Fgp)

y[is.na(Fgp)]<-NA

na.omit(y)->y
na.omit(Fgp)->Fgp
