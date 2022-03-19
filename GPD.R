library(dplyr)
library(tidyr)
library(ggplot2)
library(MASS)
library(revdbayes)
library(RColorBrewer)
# declare constants
n=10000
m=100000
#b=100
ux=1

LCV=seq(0.1,.9,.1)

#kappa=(((1/LCV)^2)-1)/2

#alpha=ux*(kappa+1)

kappa=(1-2*LCV)/LCV

alpha=ux*(kappa+1)

set.seed(123)


M=matrix(data=NA,nrow=length(LCV),ncol=7)
RMSE=matrix(data=NA,nrow=length(LCV),ncol=7)


rmse<-function(x){sqrt(sum((x-ux)^2)/m)}


####################
# loop over standard deviations
for (s in 1:length(LCV)){
  
  # make distribution with quantile function
  # set.seed(123)
  
  gpd=sapply(1:m, function(x) {p=runif(n,0,1)
  if (kappa[s]==0){
    -alpha[s]*log(1-p)}
  else{alpha[s]*(1-(1-p)^kappa[s])/kappa[s]} })
  
  
  # define expected value
  
  meanreal=function(y){dist=gpd[,y]
  u=log(dist)
  sum(u)/n->my
  sqrt(sum((u-my)^2)/n)->sdu
  exp(my+(sdu^2)/2)}
  
  
  varreal=function(y){
    dist=gpd[,y]
    u=log(dist)
    vu=var(u)
    (exp(vu)-1)*exp(2*mean(u)+vu)}
  
  # sample mean via Monte Carlo
  
  mc=colMeans(gpd)
  
  RMSE[s,1]=rmse(mc)
  
  M[s,1]=mean(mc) 
  
  # MVUE
  
  lnmle = sapply(1:ncol(gpd), meanreal)
  
  RMSE[s,2] = rmse(lnmle)
  M[s,2] = mean(lnmle)
  
  gpmle<-sapply(1:m, function(y) {gp=grimshaw_gp_mle(gpd[,y])
  gp$a/(gp$k+1) })
  
  #  pwm=sapply(1:m, function(y) {pwMoment(gpd[,y],k=0)})
  
  RMSE[s,3] = rmse(gpmle)
  M[s,3] = mean(gpmle)
  
  
  
  # 
  ####
  # median of Means
  
  # larger groups do better
  k=4
  
  # sample numbers between 1 and # of runs without replacement k times. these will be the indices of the samples that go into each group
  
  # shuffle the data
  
  mom=sapply(1:ncol(gpd), function(y){
    index<-sample(seq(1:n),n,replace=FALSE)
    f=rep(1:k,n/k)
    groups<-split(gpd[index,y],f,drop=FALSE)
    means<-sapply(groups,mean)
    median(means)})
  
  M[s,4]=mean(mom) 
  RMSE[s,4]=rmse(mom)
  
  # trimmed mean
  trimmedmean=c()
  trimmedmean=sapply(1:m, function(y){q=quantile(gpd[,y],c(.05,.95),na.rm = T)
  x=gpd[,y]
  x[x<q[1]]<-NA
  x[x>q[2]]<-NA
  mean(x,na.rm=T)})
  
  
  
  M[s,5]=mean(trimmedmean)
  RMSE[s,5]=rmse(trimmedmean)
  
  winsor=c()
  
  winsor=sapply(1:m, function(y){q=quantile(gpd[,y],c(.05,.95),na.rm = T)
  x=gpd[,y]
  x[x<q[1]]<-q[1]
  x[x>q[2]]<-q[2]
  mean(x,na.rm=T)})
  
  M[s,6]=mean(winsor)
  RMSE[s,6]=rmse(winsor)
  
  
  plmle=c()
  
  plmle=sapply(1:m, function(y){alphahat=1+n*(1/sum(log(gpd[,y]/min(gpd[,y]))))
  (min(gpd[,y])/(alphahat-2))+min(gpd[,y])
  })
  
  M[s,7]=mean(plmle)
  RMSE[s,7]=rmse(plmle)
  
  
  message(s)
}


as.data.frame(M)->M  
names(M)=c("Sample Mean","LN2 MLE","GP MLE","Median of Means","Trimmed Mean","Winsorized Mean","Power Law MLE")

M$LCV=LCV

as.data.frame(RMSE)->RMSE
names(RMSE)=c("Sample Mean","LN2 MLE","GP MLE","Median of Means","Trimmed Mean","Winsorized Mean","Power Law MLE")

RMSE$LCV=LCV

######

RMSE %>% #select(`Sample Mean`,`LN2 MLE`,`GP MLE`,`MOM k=23`,`Trimmed Mean`,LCV) %>%
  gather(Method,RMSE,-LCV)->RMSE

RMSE %>% saveRDS("/cluster/tufts/lamontagnelab/fdolan03/EstimatingMean/GPD_n104_m105_rmse.rds")


RMSE %>% mutate(Type=ifelse(grepl("MLE",Method),Method,"Nonparametric"))->RMSE

RMSE %>% mutate(Estimator=ifelse(Type=="Nonparametric",Method,"MLE"))->RMSE

p1<-ggplot(RMSE,aes(LCV,RMSE,col=Estimator))+geom_point(aes(shape=Type,group=Method),size=3)+scale_shape_manual(values=c(8,15,6,17))+
  geom_line(aes(group=Method),size=1)+ggtitle("GPA n=10,000, mc=100,000")+guides(shape=guide_legend(""))+
  xlab("LCV distribution")+ylab("RMSE of Sample Mean")+theme_bw()+scale_color_manual(values=c("Black","Green","Red","Blue","darkgrey"))+
  #scale_y_continuous(trans = 'log10',breaks = trans_breaks('log10', function(x) 10^x),labels = trans_format('log10', math_format(10^.x)))+#++
  coord_cartesian(ylim=c(0,2 ))

ggsave("/cluster/tufts/lamontagnelab/fdolan03/EstimatingMean/GPD_n104_m105_rmse.png",plot=p1,units="in",dpi=300,height=5,width=7)


M %>% gather(Method,Mean,-LCV)->M

M %>% mutate(Type=ifelse(grepl("MLE",Method),Method,"Nonparametric"))->M

M %>% mutate(Estimator=ifelse(Type=="Nonparametric",Method,"MLE"))->M

p2<-ggplot(M,aes(LCV,Mean,col=Estimator))+geom_point(aes(shape=Type,group=Method),size=3)+scale_shape_manual(values=c(8,15,6,17))+
  geom_line(aes(group=Method),size=1)+ggtitle("GPA n=1,000, mc=100,000")+guides(shape=guide_legend(""))+
  xlab("LCV distribution")+ylab("RMSE of Sample Mean")+theme_bw()+scale_color_manual(values=c("Black","Green","Red","Blue","darkgrey"))+
  #scale_y_continuous(trans = 'log10',breaks = trans_breaks('log10', function(x) 10^x),labels = trans_format('log10', math_format(10^.x)))+#++
  coord_cartesian(ylim=c(0,2 ))

ggsave("/cluster/tufts/lamontagnelab/fdolan03/EstimatingMean/GPD_n104_m105_bias.png",plot=p2,units="in",dpi=300,height=5,width=7)

M %>% saveRDS("/cluster/tufts/lamontagnelab/fdolan03/EstimatingMean/GPD_n104_m105_bias.rds")
RMSE %>% saveRDS("/cluster/tufts/lamontagnelab/fdolan03/EstimatingMean/GPD_n104_m105_rmse.rds")
