library(revdbayes)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
# libraries

#####
# read in dat

honey<-read_excel('datasets/honeycreekdata.xlsx')

honeyss<-honey$`SS, mg/L (suspended solids)`
honeyf<-honey$`Flow, CFS`

honeyss[honeyss<=0]<-NA
honeyf[honeyf<=0]<-NA
#
cities<-read.table('datasets/cities.txt')

cities=as.numeric(unlist(cities))
cities[cities<=0]<-NA
#
flare<-read.table('datasets/flares.txt')

flare=as.numeric(unlist(flare))
flare[flare<=0]<-NA
#
fires<-read.table('datasets/fires.txt')

fires=as.numeric(unlist(fires))
fires[fires<=0]<-NA
#
terr<-read.table('datasets/terrorism_deaths.txt')

terr=as.numeric(unlist(terr))
terr[terr<=0]<-NA

aid<-readRDS('C:/Users/fdolan03/Box/BootstrapMean/datasets/foreignaid.rds')

#waves<-read.csv('C:/Users/fdolan03/Box/BootstrapMean/datasets/indonesianwave_heights_GPA.txt',skip=9)
#wh<-waves$HS[waves$HS>0]
#######################################
names=c("Cities","Wildfire Area","Solar Flare Intensity","Deaths from Terrorism","Sediment Load","Foreign Aid","Normal","Streamflow")




for (j in 2:9){
  ##
  if (j==1){cities->x
            x=sample(x)} #1
  
  else if (j==2){fires->x
                  x=sample(x)} #2
  
  else if(j==3){flare->x} #3
  
  else if(j==4){terr->x
                x=sample(x)} #4
  
  else if(j==5){honeyss->x} #5
  
  else if(j==6){aid->x}
  
  else if(j==7){a=0.11
    k=-0.88
    p=runif(10000,0,1)
    (a*(1-(1-p)^k)/k)->x}
  
  else if(j==8){flow->x}
  
  else{croploss->x}
  
  na.omit(x)->x
  
  mandel=matrix(NA,ncol=2,nrow=length(x))
  
  for (i in 1:length(x)){
    
    #sample=sample(x,size=i,replace=TRUE)
    
    mandel[i,1]=i
    
    #x=rev(x)
    
    mandel[i,2]=sd(x[1:i],na.rm=T)
    

  }
  message(j)
  
  mandel %>% saveRDS(paste0("datasets/MandelPlot/","gpa",".rds"))
  
  }


Mnfracs=cbind(seq(1,100),Mnfracs)

Mnfracs=as.data.frame(Mnfracs)

names(Mnfracs)=c("p","US City Population","Wildfire Area","Solar Flare Intensity",
                 "Deaths from Terrorism","Sediment Load","Foreign Aid","Normal",
                 "Streamflow","Wave Heights")


Mnfracs %>% gather(key=Dataset,value=value,-p) ->dat

lmom::samlmu(wh,nmom=3)[3]

dat_text <- data.frame(
  label = c("L-Skew=0.90", "L-Skew=0.81", "L-Skew=0.97","L-Skew=0.92","L-Skew=0.78",
            "L-Skew=0.61","L-Skew=0","L-Skew=0.60","L-Skew=0.03"),
  Dataset   = c("Foreign Aid","US City Population","Wildfire Area","Solar Flare Intensity",
                "Deaths from Terrorism","Sediment Load","Normal","Streamflow","Wave Heights") ,
  Col = c("0","0","0","0","0","0","1","0","0"))

dat %>% mutate(Col=ifelse(Dataset=="Normal","1","0"))->dat

p<-ggplot(dat,aes(p,value,col=Col))+geom_line(size=1)+theme_bw()+
  facet_wrap(~factor(Dataset,levels=c("Normal","Wave Heights","Streamflow","Sediment Load","Deaths from Terrorism",
                                      "US City Population","Foreign Aid","Solar Flare Intensity","Wildfire Area")))+ylab("Fraction of Empirical Mean")+
  xlab("Percentile Bin Dropped")+scale_color_manual(values=c("black","red"))+
  guides(color="none")

P<- p + geom_text(
  data    = dat_text,
  mapping = aes(x = -Inf, y = -Inf, label = label),
  hjust   = -0.1,
  vjust   = -1
)


######

sites<-readRDS("newusgssites.rds")
sites %>% ungroup() %>% arrange(desc(perc_err))

sites %>% ungroup() %>% arrange(desc(PLppcc))

high_gpppcc <- sites %>% filter(GPppcc>0.999 & Lskew>0.75)

high_PLppcc <- sites %>% filter(PLppcc>=0.999)


higherr<-sites %>% filter(perc_err>1000)

sites %>% filter(site_no=="2357000")

files=list.files('datasets/USGS/sites/',full.names=TRUE)


for (i in 1:length(files)){
  
  dat=read.table(files[i])
  
  if (nrow(dat)>1000){
    
    message(i)
    
  
  if ( dat$site_no %in% high_gpppcc$site_no){
    
    mandel=matrix(NA,nrow=2,ncol=nrow(dat))
    
    
    for (j in 1:nrow(dat)){   
      
    #sample=sample(dat$X_00060_00003,size=j,replace=TRUE)
    
    mandel[1,j]=j
    
    mandel[2,j]=sd((dat$X_00060_00003)[1:j],na.rm=T)
    
    
    }
    
    mandel %>% saveRDS(paste0("datasets/MandelPlot/croploss.rds"))}
}
}

mandel<-readRDS('datasets/MandelPlot/2357000.rds')

plot(mandel[1,],mandel[2,],type="l",ylab="Standard Deviation",xlab="Sample Size",main="USGS site 2357000")


mandel<-readRDS('datasets/MandelPlot/Wildfire Area.rds')

plot(mandel[1,],mandel[2,],type="l",ylab="Standard Deviation",xlab="Sample Size",main="Deaths from Terrorism")



######



##

terr<-readRDS('datasets/MandelPlot/Deaths from Terrorism.rds')
terr=t(terr)
terr=cbind(terr,"Deaths from Terrorism")
terr=as_tibble(terr)
names(terr)=c("Sample Size","Standard Deviation","Dataset")

##

cities<-readRDS('datasets/MandelPlot/Cities.rds')
cities=t(cities)
cities=cbind(cities,"US City Population")
cities=as_tibble(cities)
names(cities)=c("Sample Size","Standard Deviation","Dataset")

##
fire=readRDS('datasets/MandelPlot/Wildfire Area.rds')
fire=t(fire)
fire=cbind(fire,"Wildfire Area")
fire=as_tibble(fire)
names(fire)=c("Sample Size","Standard Deviation","Dataset")

##

flare=readRDS('datasets/MandelPlot/Solar Flare Intensity.rds')
flare=t(flare)
flare=cbind(flare,"Solar Flare Intensity")
flare=as_tibble(flare)
names(flare)=c("Sample Size","Standard Deviation","Dataset")

##

flow=readRDS('datasets/MandelPlot/2357000.rds')
flow=t(flow)
flow=cbind(flow,"Streamflow")
flow=as_tibble(flow)
names(flow)=c("Sample Size","Standard Deviation","Dataset")

##

aid=readRDS('datasets/MandelPlot/Foreign Aid.rds')
aid=t(aid)
aid=cbind(aid,"Foreign Aid")
aid=as_tibble(aid)
names(aid)=c("Sample Size","Standard Deviation","Dataset")

##

sed=readRDS('datasets/MandelPlot/Sediment Load.rds')
sed=t(sed)
sed=cbind(sed,"Sediment Load")
sed=as_tibble(sed)
names(sed)=c("Sample Size","Standard Deviation","Dataset")

##


crop=readRDS('datasets/MandelPlot/croploss.rds')
crop=cbind(crop,"US Crop Losses")
crop=as_tibble(crop)
names(crop)=c("Sample Size","Standard Deviation","Dataset")


##


gpa=readRDS('datasets/MandelPlot/gpa.rds')
gpa=cbind(gpa,"Generalized Pareto LCV=0.9")
gpa=as_tibble(gpa)
names(gpa)=c("Sample Size","Standard Deviation","Dataset")

###

m=rbind(gpa,crop,sed,terr,cities,flare,fire,flow,aid)

m<-m %>% mutate(`Sample Size`=as.numeric(`Sample Size`, `Standard Deviation`=as.numeric(`Standard Deviation`)))


ggplot(m,aes(`Sample Size`,`Standard Deviation`,group=Dataset))+geom_line(size=1)+theme_bw()+
  facet_grid(~factor(Dataset,levels=c("Generalized Pareto LCV=0.9","Streamflow","Sediment Load","US Crop Losses","Deaths from Terrorism",
                                      "US City Population","Foreign Aid","Solar Flare Intensity","Wildfire Area")))+ylab("Fraction of Empirical Mean")+
  xlab("Percentile Bin Dropped")#scale_color_manual(values=c("black","red"))+
 # guides(color="none")

gpa$`Sample Size`=as.numeric(gpa$`Sample Size`)
gpa$`Standard Deviation`=as.numeric(gpa$`Standard Deviation`)

p1<-ggplot(gpa,aes(`Sample Size`,`Standard Deviation`,group=Dataset))+geom_line(size=1,col="red")+theme_bw()+
  ggtitle("GPA LCV=0.9")

crop$`Sample Size`=as.numeric(crop$`Sample Size`)
crop$`Standard Deviation`=as.numeric(crop$`Standard Deviation`)

p2<-ggplot(crop,aes(`Sample Size`,`Standard Deviation`,group=Dataset))+
  geom_line(size=1)+theme_bw()+ggtitle("US Crop Losses")

sed$`Sample Size`=as.numeric(sed$`Sample Size`)
sed$`Standard Deviation`=as.numeric(sed$`Standard Deviation`)

p3<-ggplot(sed,aes(`Sample Size`,`Standard Deviation`,group=Dataset))+
  geom_line(size=1)+theme_bw()+ggtitle("Sediment Load")


terr$`Sample Size`=as.numeric(terr$`Sample Size`)
terr$`Standard Deviation`=as.numeric(terr$`Standard Deviation`)

p4<-ggplot(terr,aes(`Sample Size`,`Standard Deviation`,group=Dataset))+
  geom_line(size=1)+theme_bw()+ggtitle("Terrorism Deaths")


cities$`Sample Size`=as.numeric(cities$`Sample Size`)
cities$`Standard Deviation`=as.numeric(cities$`Standard Deviation`)

p5<-ggplot(cities,aes(`Sample Size`,`Standard Deviation`,group=Dataset))+
  geom_line(size=1)+theme_bw()+ggtitle("US City Population")

flare$`Sample Size`=as.numeric(flare$`Sample Size`)
flare$`Standard Deviation`=as.numeric(flare$`Standard Deviation`)

p6<-ggplot(flare,aes(`Sample Size`,`Standard Deviation`,group=Dataset))+
  geom_line(size=1)+theme_bw()+ggtitle("Solar Flare Intensity")


fire$`Sample Size`=as.numeric(fire$`Sample Size`)
fire$`Standard Deviation`=as.numeric(fire$`Standard Deviation`)

p7<-ggplot(fire,aes(`Sample Size`,`Standard Deviation`,group=Dataset))+
  geom_line(size=1)+theme_bw()+ggtitle("Wildfire Area")

flow$`Sample Size`=as.numeric(flow$`Sample Size`)
flow$`Standard Deviation`=as.numeric(flow$`Standard Deviation`)

p8<-ggplot(flow,aes(`Sample Size`,`Standard Deviation`,group=Dataset))+
  geom_line(size=1)+theme_bw()+ggtitle("Streamflow")

aid$`Sample Size`=as.numeric(aid$`Sample Size`)
aid$`Standard Deviation`=as.numeric(aid$`Standard Deviation`)

p9<-ggplot(aid,aes(`Sample Size`,`Standard Deviation`,group=Dataset))+
  geom_line(size=1)+theme_bw()+ggtitle("Foreign Aid")


cowplot::plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,labels="AUTO",ncol=3)

ggsave("mandelplot_all.png",dpi=300)
