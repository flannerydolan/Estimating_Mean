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
  


#######################################


fmean=c()
Mnfracs=matrix(NA,nrow=100,ncol=9)


for (j in 1:9){
  ##
  if (j==1){cities->x} 
  
  else if (j==2){fires->x} 
  
  else if(j==3){flare->x} 
  
  else if(j==4){terr->x} 
  
  else if(j==5){honeyss->x} 
  
  else if(j==6){aid->x}
  
  else if(j==7){rnorm(10000,1,1)->x}
  
  else if(j==8){flow->x}
  
  else{croploss->x}
  
  na.omit(x)->x
  
for (i in 1:100){
  
  
  q2=quantile(x,i/100,na.rm=T)
  

  x1=x[x<=q2]
  
  fmean[i]=mean(x1,na.rm=T)
  
  frac[i]=fmean[i]/mean(x,na.rm=T)

}
Mnfracs[,j]=frac
}


Mnfracs=cbind(seq(1,100),Mnfracs)

Mnfracs=as.data.frame(Mnfracs)

names(Mnfracs)=c("p","US City Population","Wildfire Area","Solar Flare Intensity",
                 "Deaths from Terrorism","Sediment Load","Foreign Aid","Normal",
                 "Streamflow","US Crop Losses")


Mnfracs %>% gather(key=Dataset,value=value,-p) ->dat

lmom::samlmu(croploss,nmom=3)[3]
lmom::samlmu(x,nmom=3)[3]
  
dat_text <- data.frame(
    label = c("L-Skew=0.90", "L-Skew=0.81", "L-Skew=0.97","L-Skew=0.92","L-Skew=0.78",
              "L-Skew=0.61","L-Skew=0","L-Skew=0.51","L-Skew=0.76"),
    Dataset   = c("Foreign Aid","US City Population","Wildfire Area","Solar Flare Intensity",
              "Deaths from Terrorism","Sediment Load","Normal","Streamflow","US Crop Losses") ,
    Col = c("0","0","0","0","0","0","1","0","0"))

dat %>% mutate(Col=ifelse(Dataset=="Normal","1","0"))->dat

p<-ggplot(dat,aes(p,value,col=Col))+geom_line(size=1)+theme_bw()+
  facet_wrap(~factor(Dataset,levels=c("Normal","Streamflow","Sediment Load","US Crop Losses","Deaths from Terrorism",
                               "US City Population","Foreign Aid","Solar Flare Intensity","Wildfire Area")))+ylab("Fraction of Empirical Mean")+
  xlab("Percentile Bin Dropped")+scale_color_manual(values=c("black","red"))+
  guides(color="none")

P<- p + geom_text(
  data    = dat_text,
  mapping = aes(x = -Inf, y = -Inf, label = label),
  hjust   = -0.25,
  vjust   = -1
  )

P


ggsave("Fig1_percentiledropped_crop_usgs.png",P,dpi=300)
