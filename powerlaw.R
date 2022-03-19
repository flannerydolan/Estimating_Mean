library(dplyr)
library(tidyr)
library(ggplot2)
library(MASS)
library(revdbayes)



n=10000
m=100000
ux=1


LCV=seq(.1,.9,.1)

alpha=(1+3*LCV)/(2*LCV)


xmin=ux*(alpha-2)/(alpha-1)

rmse<-function(x){sqrt(sum((x-ux)^2)/m)}


RMSE=matrix(data=NA,nrow=length(LCV),ncol=7)
M=matrix(data=NA,nrow=length(LCV),ncol=7)

####################
# loop over standard deviations
ptm <- proc.time()


for (s in 1:length(LCV)){
  
  # make distribution with quantile function
  pow=c()
  
  pow=sapply(1:m, function(x) {p=runif(n,0,1)
  #xmin*p^(1/(1-alpha[s]))
  #(p/(xmin^(alpha-1)))^(1/(-(alpha-1)))
  exp((-log((-p+1)/exp((alpha[s]-1)*log(xmin[s]))))/(alpha[s]-1))
  })
  
  meanreal=function(y){dist=pow[,y]
  u=log(dist)
  sum(u)/n->my
  sqrt(sum((u-my)^2)/n)->sdu
  exp(my+(sdu^2)/2)}
  
  
  varreal=function(y){
    dist=pow[,y]
    u=log(dist)
    vu=var(u)
    
    (exp(vu)-1)*exp(2*mean(u)+vu)}
  
  # sample mean via Monte Carlo
  
  mc=colMeans(pow)
  
  RMSE[s,1]=rmse(mc)
  
  M[s,1]=mean(mc) 
  
  # LN2 MLE
  
  lnmle = sapply(1:ncol(pow), meanreal)
  
  RMSE[s,2] = rmse(lnmle)
  M[s,2] = mean(lnmle)
  
  gpmle<-sapply(1:m, function(y) {gp=grimshaw_gp_mle(pow[,y])
  gp$a/(gp$k+1) })
  
  RMSE[s,5] = rmse(gpmle)
  M[s,5] = mean(gpmle)
  
  
  ####
  # median of Means
  
  # larger groups do better
  k=4
  
  # sample numbers between 1 and # of runs without replacement k times. these will be the indices of the samples that go into each group
  
  # shuffle the data
  
  mom23=sapply(1:ncol(pow), function(y){
    index<-sample(seq(1:n),n,replace=FALSE)
    f=rep(1:k,n/k)
    groups<-split(pow[index,y],f,drop=FALSE)
    means<-sapply(groups,mean)
    median(means)})
  
  M[s,3]=mean(mom23) 
  RMSE[s,3]=rmse(mom23)
  
  mle=c()
  
  mle=sapply(1:m, function(y){alphahat=1+n*(1/sum(log(pow[,y]/min(pow[,y]))))
  (min(pow[,y])/(alphahat-2))+min(pow[,y])
  })
  
  M[s,4]=mean(mle)
  RMSE[s,4]=rmse(mle)
  
  
  # trimmed mean
  trimmedmean=c()
  trimmedmean=sapply(1:m, function(y){q=quantile(pow[,y],c(.05,.95),na.rm = T)
  x=pow[,y]
  x[x<q[1]]<-NA
  x[x>q[2]]<-NA
  mean(x,na.rm=T)})
  
  
  
  M[s,7]=mean(trimmedmean)
  RMSE[s,7]=rmse(trimmedmean)
  
  winsor=c()
  
  winsor=sapply(1:m, function(y){q=quantile(pow[,y],c(.05,.95),na.rm = T)
  x=pow[,y]
  x[x<q[1]]<-q[1]
  x[x>q[2]]<-q[2]
  mean(x,na.rm=T)})
  
  M[s,6]=mean(winsor)
  RMSE[s,6]=rmse(winsor)
  
  message(s)
}


as.data.frame(M)->M  
names(M)=c("Sample Mean","LN2 MLE","Median of Means","Power Law MLE","GP MLE","Winsorized Mean","Trimmed Mean")

M$LCV=LCV

M %>% #select(`Sample Mean`,Bagged,Bragged,`MOM k=23`,`LN2 MLE`,CV,`Power Law MLE`) %>%
  gather(Method,Mean,-LCV)->M


M %>% mutate(Type=ifelse(grepl("MLE",Method),Method,"Nonparametric"))->M

M %>% mutate(Estimator=ifelse(Type=="Nonparametric",Method,"MLE"))->M

p1<-ggplot(M,aes(LCV,Mean-1,col=Estimator))+geom_point(aes(shape=Type,group=Method),size=3)+scale_shape_manual(values=c(8,15,6,17))+
  geom_line(aes(group=Method),size=1)+ggtitle("Power Law n=10,000, mc=1,000")+guides(shape=guide_legend(""))+
  xlab("LCV distribution")+ylab("Bias")+theme_bw()+scale_color_manual(values=c("Green","Black","Red","Blue","darkgrey"))+
  #scale_y_continuous(trans = 'log10',breaks = trans_breaks('log10', function(x) 10^x),labels = trans_format('log10', math_format(10^.x)))+#++
  coord_cartesian(ylim=c(-1,1 ))

ggsave("/cluster/tufts/lamontagnelab/fdolan03/EstimatingMean/pl_n104_m105_bias.png",plot=p1,units="in",dpi=300,height=5,width=7)

###################################################

as.data.frame(RMSE)->RMSE
names(RMSE)=c("Sample Mean","LN2 MLE","Median of Means","Power Law MLE","GP MLE","Winsorized Mean","Trimmed Mean")

RMSE$LCV=LCV

######

RMSE %>% #select(`Sample Mean`,Bagged,Bragged,`MOM k=23`,`LN2 MLE`,CV,`Power Law MLE`) %>%
  gather(Method,RMSE,-LCV)->RMSE


RMSE %>% mutate(Type=ifelse(grepl("MLE",Method),Method,"Nonparametric"))->RMSE

RMSE %>% mutate(Estimator=ifelse(Type=="Nonparametric",Method,"MLE"))->RMSE

p2<-ggplot(RMSE,aes(LCV,RMSE,col=Estimator))+geom_point(aes(shape=Type,group=Method),size=3)+scale_shape_manual(values=c(8,15,6,17))+
  geom_line(aes(group=Method),size=1)+ggtitle("Power Law n=10,000, mc=1,000")+guides(shape=guide_legend(""))+
  xlab("LCV distribution")+ylab("RMSE of Sample Mean")+theme_bw()+scale_color_manual(values=c("Green","Black","Red","Blue","darkgrey"))+
  #scale_y_continuous(trans = 'log10',breaks = trans_breaks('log10', function(x) 10^x),labels = trans_format('log10', math_format(10^.x)))+#++
  coord_cartesian(ylim=c(0,.7 ))

ggsave("/cluster/tufts/lamontagnelab/fdolan03/EstimatingMean/pl_n104_m105_rmse.png",plot=p2,units="in",dpi=300,height=5,width=7)

M %>% saveRDS("/cluster/tufts/lamontagnelab/fdolan03/EstimatingMean/powerlaw_bias_n104_m105.rds")
RMSE %>% saveRDS("/cluster/tufts/lamontagnelab/fdolan03/EstimatingMean/powerlaw_rmse_n104_m105.rds")
