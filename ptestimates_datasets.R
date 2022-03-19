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
  
  crop17 %>% saveRDS('datasets/croploss_2017.rds')
  
  crop17$INDEMNITY->croploss
  #######################################
  MnEst=matrix(NA,nrow=8,ncol=7)
  
  for (i in 1:8){
  ##
  if (i==1){cities->x} #1
  
  else if (i==2){fires->x} #2
  
  else if(i==3){flare->x} #3
  
  else if(i==4){terr->x} #4
  
  else if(i==5){honeyss->x} #5
  
  else if (i==6){croploss->x}
  
  else if (i==7){aid->x}
    
  else {flow->x}
  
  ###
  n=length(x)
  Means=c()
  
  x[x<=0]<-NA
  x<-na.omit(x)
  # LN2 MLE
  
  meanreal=function(dist){
            u=log(dist)
            sum(u,na.rm=T)/n->my
            sqrt(sum((u-my)^2,na.rm = T)/n)->sdu
            exp(my+(sdu^2)/2)}
  
  Means[1]=meanreal(x)
  
  # GP MLE
  gpmle<-function(x){
           gp=grimshaw_gp_mle(x)
            gp$a/(gp$k+1) }
  
  Means[2]=gpmle(x)
  
  ##????????? negative k?
  
  # Power Law MLE
  
  plmle<-function(x){alphahat=1+n*(1/sum(log(x/min(x)),na.rm = T))
                    (min(x,na.rm=T)/(alphahat-2))+min(x,na.rm=T)}
  
  Means[3]=plmle(x)
  
  # Median of Means
  k=4
  
  mom<-function(x){index<-sample(seq(1:n),n,replace=FALSE)
                    f=rep(1:k,n/k)
                    groups<-split(x[index],f,drop=FALSE)
                    groups<-sapply(groups,na.omit)
                    means<-sapply(groups,mean)
                    median(means,na.rm=T) }
  
  Means[4]=mom(x)
  
  # Trimmed Mean
  # 
  # delta=0.01
  # 
  # trimmed<-function(x){#x2<-sample(x,n/2,replace =FALSE)
  #                     #epsilon= (16*log(8/delta))/(3*(n/2))
  #                     #A=sort(x2)[epsilon*(n/2)]
  #                     #B=sort(x2)[(1-epsilon)*(n/2)]
  #                     #x2[x2<A]<-NA
  #                     #x2[x2>B]<-NA
  #                     #mean(x2,na.rm=T)
  #   q=quantile(x,c(.05,.95))
  #   x[x<q[1]]<-NA
  #   x[x>q[2]]<-NA
  #   mean(x,na.rm=T)
  #   }
  # 
  # Means[5]=trimmed(x)
  # 
  # # winsorized mean
  # 
  # winsor<-function(x){q=quantile(x,c(.05,.95))
  #                     x[x<q[1]]<-q[1]
  #                     x[x>q[2]]<-q[2]
  #                     mean(x,na.rm=T)}
  # 
  # Means[6]=winsor(x)
  
  # Empirical Mean
  
  Means[7]=mean(x,na.rm=T)
  
  ##################################
  
  
  Means-> MnEst[i,]
  
  }
  
  MnEst
  
  
  MnEst=as.data.frame(MnEst)
  
  names(MnEst)<-c("LN2 MLE","GP MLE","Power Law MLE","Median of Means","Trimmed Mean","Winsorized Mean","Empirical Mean")
  rownames(MnEst)<-c("City Population","Wildfire Area","Solar Flare Intensity","Terrorism Deaths","Honey Creek SS","Honey Creek Flow","Foreign Aid","Wave Heights")
  
  MnEst$Dataset=rownames(MnEst)
  
  MnEst %>% gather(key=Estimator,value=value,`LN2 MLE`:`Empirical Mean`) ->MnEst2

MnEst2 %>% mutate(fill=as.factor(ifelse(Estimator=="Empirical Mean",1,0)))->MnEst2

ggplot(MnEst2,aes(Estimator,value,fill=fill))+geom_col()+facet_wrap(~Dataset,scales="free")+scale_fill_manual(values=c("darkgrey","red"))+
  theme(axis.text=element_text(angle=90),legend.position = "none")+ylab("")

#######################################

# 
# # bootstrap sample
# 
# meanboot=c()
# R=c()
# 
# # LN2 MLE
# 
# boot=c()
# 
# for(b in 1:1000){
#   
#   boot[b]=meanreal(sample(x,n,replace=T))
#   
# }
# 
# meanboot[1]=mean(boot)
# R[1]=var(boot)
# 
# # GP MLE
# 
# boot=c()
# 
# for(b in 1:1000){
#   
#   
#   boot[b]=gpmle(sample(x,n,replace=T))
#   
# }
# 
# meanboot[2]=mean(boot)
# R[2]=var(boot)
# 
# # Power Law
# 
# boot=c()
# 
# for(b in 1:1000){
#   
#   
#   boot[b]=plmle(sample(x,n,replace=T))
#   
# }
# 
# meanboot[3]=mean(boot)
# R[3]=var(boot)
# 
# # Median of Means
# 
# boot=c()
# 
# for(b in 1:1000){
#   
#   
#   boot[b]=mom(sample(x,n,replace=T))
#   
# }
# 
# meanboot[4]=mean(boot)
# R[4]=var(boot)
# 
# # Trimmed Mean
# 
# boot=c()
# 
# for(b in 1:1000){
#   
#   
#   boot[b]=trimmed(sample(x,n,replace=T))
#   
# }
# 
# meanboot[5]=mean(boot)
# R[5]=var(boot)
# 
# # winsorized
# 
# boot=c()
# 
# for(b in 1:1000){
#   
#   
#   boot[b]=winsor(sample(x,n,replace=T))
#   
# }
# 
# meanboot[6]=mean(boot)
# R[6]=var(boot)
# 
# 
# # Empirical Mean
# boot=c()
# 
# for(b in 1:1000){
#       boot[b]=mean(sample(x,n,replace=T))
# }
# 
# meanboot[7]=mean(boot)
# R[7]=var(boot)
# ############
# 
# names(meanboot)<-c("LN2 MLE","GP MLE","Power Law MLE","Median of Means","Trimmed Mean","Winsorized Mean","Empirical Mean")
# 
# as.data.frame(meanboot)->meanboot
# 
# meanboot$Estimator=rownames(meanboot)
# 
# meanboot %>% rename(Mean=meanboot)->meanboot
# 
# 
# ## var
# 
# names(R)<-c("LN2 MLE","GP MLE","Power Law MLE","Median of Means","Trimmed Mean","Winsorized Mean","Empirical Mean")
# 
# as.data.frame(R)->R
# 
# R$Estimator=rownames(R)
# 
# R %>% rename(Variance=R)->R
# 
# 
# R %>% left_join(meanboot,by="Estimator") %>% gather(key=Statistic,value=Value,Variance,Mean) -> bootdat
# 
# 
# bootdat %>% mutate(Fill=as.factor(ifelse(Estimator=="Empirical Mean",1,0)))->bootdat
# 
# ggplot(bootdat,aes(Estimator,Value,fill=Fill))+geom_col()+facet_wrap(~Statistic,scales="free")+scale_fill_manual(values=c("darkgrey","red"))+
#   theme(plot.title=element_text(size=16),axis.text=element_text(angle=90),legend.position = "none")+ylab("")+ylab("Honey Creek flow (cfs)")+
#   coord_cartesian(ylim=c(-200,500))


fmean=c()
Mnfracs=matrix(NA,nrow=100,ncol=9)


for (j in 1:9){
  ##
  if (j==1){cities->x} #1
  
  else if (j==2){fires->x} #2
  
  else if(j==3){flare->x} #3
  
  else if(j==4){terr->x} #4
  
  else if(j==5){honeyss->x} #5
  
  else if(j==6){aid->x}
  
  else if(j==7){rnorm(10000,mean=1,sd=1)->x}
  
  else if(j==8){honeyf->x}
  
  else{wh->x}
  
  na.omit(x)->x
  
for (i in 1:100){
  
  
  q1=quantile(x,(i-1)/100,na.rm = T)
  q2=quantile(x,(i)/100,na.rm=T)
  
  x1=x[x<=q1 | x>=q2]
  
  fmean[i]=mean(x1,na.rm=T)
  
  frac=fmean/mean(x,na.rm=T)

}
Mnfracs[,j]=frac
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


ggsave("Fig1_percentiledropped_new.png",P,dpi=300)
