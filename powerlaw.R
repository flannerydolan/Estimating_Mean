library(dplyr)
library(tidyr)
library(ggplot2)
library(MASS)
library(revdbayes)

# declare constants

n=100
m=100000
ux=1

# Initialize sequence of LCVs

LCV=seq(.1,.9,.1)

# Distribution parameters in terms of LCV

alpha=(1+3*LCV)/(2*LCV)

xmin=ux*(alpha-2)/(alpha-1)

###
set.seed(123)

# functions

lscale<-function(x){lmom::samlmu(x,nmom=2)[2]}

mae<-function(x){mean(abs(x-ux),na.rm=T)}

bias<-function(x){mean(x,na.rm=T)-ux}


#####
# Initialize matrices


L=matrix(data=NA,nrow=length(LCV),ncol=7)
MAE=matrix(data=NA,nrow=length(LCV),ncol=7)
B=matrix(data=NA,nrow=length(LCV),ncol=7)

####################
# loop over LCVs


for (s in 1:length(LCV)){
  
  # make distribution with quantile function
  pow=c()
  
  pow=sapply(1:m, function(x) {p=runif(n,0,1)
		exp((-log((-p+1)/exp((alpha[s]-1)*log(xmin[s]))))/(alpha[s]-1))
 		 })
  
  meanreal=function(y){dist=pow[,y]
  		u=log(dist)
  		sum(u)/n->my
  		sqrt(sum((u-my)^2)/n)->sdu
  		exp(my+(sdu^2)/2)}
  
  
  # sample mean via Monte Carlo
  
  mc=colMeans(pow)
  
  
  L[s,1]=lscale(mc)
  MAE[s,1]=mae(mc) 
  B[s,1]=bias(mc)

  # LN2 MLE
  
  lnmle = sapply(1:ncol(pow), meanreal)
  

  L[s,2]=lscale(lnmle)
  MAE[s,2]=mae(lnmle) 
  B[s,2]=bias(lnmle)
    
  gpmle<-sapply(1:m, function(y) {gp=grimshaw_gp_mle(pow[,y])
  		gp$a/(gp$k+1) })
  
  
  L[s,5]=lscale(gpmle)
  MAE[s,5]=mae(gpmle)
  
  ####
  # median of Means
  
  k=4
  
  # sample numbers between 1 and # of runs without replacement k times. these will be the indices of the samples that go into each group
  
  # shuffle the data
  
  mom=sapply(1:ncol(pow), function(y){
    		index<-sample(seq(1:n),n,replace=FALSE)
    		f=rep(1:k,n/k)
   		groups<-split(pow[index,y],f,drop=FALSE)
    		means<-sapply(groups,mean)
    		median(means)})
  

  L[s,3]=lscale(mom)
  MAE[s,3]=mae(mom)
  B[s,3]=bias(mom)
  
  ##
   
  mle=c()
  
  mle=sapply(1:m, function(y){alphahat=1+n*(1/sum(log(pow[,y]/min(pow[,y]))))
  		(min(pow[,y])/(alphahat-2))+min(pow[,y])
  		})
  

  L[s,4]=lscale(mle)
  MAE[s,4]=mae(mle)   
  B[s,4]=bias(mle)

  # trimmed mean

  trimmedmean=c()
  trimmedmean=sapply(1:m, function(y){q=quantile(pow[,y],c(.05,.95),na.rm = T)
  		     x=pow[,y]
  		     x[x<q[1]]<-NA
  		     x[x>q[2]]<-NA
  		     mean(x,na.rm=T)})
  
    
  L[s,7]=lscale(trimmedmean)
  MAE[s,7]=mae(trimmedmean) 
  B[s,7]=bias(trimmedmean)

  # winsorized mean

  winsor=c()
  
  winsor=sapply(1:m, function(y){q=quantile(pow[,y],c(.05,.95),na.rm = T)
  		x=pow[,y]
  		x[x<q[1]]<-q[1]
  		x[x>q[2]]<-q[2]
  		mean(x,na.rm=T)})

 
  L[s,6]=lscale(winsor)
  MAE[s,6]=mae(winsor) 
  B[s,6]=bias(winsor)
   
  message(s)
}


as.data.frame(MAE)->MAE  
names(MAE)=c("Sample Mean","LN2 MLE","Median of Means","Power Law MLE","GP MLE","Winsorized Mean","Trimmed Mean")

MAE$LCV=LCV

MAE %>% gather(Method,MAE,-LCV)->MAE


MAE %>% mutate(Type=ifelse(grepl("MLE",Method),Method,"Nonparametric"))->MAE

MAE %>% mutate(Estimator=ifelse(Type=="Nonparametric",Method,"MLE"))->MAE

p1<-ggplot(MAE,aes(LCV,MAE,col=Estimator))+geom_point(aes(shape=Type,group=Method),size=3)+scale_shape_manual(values=c(8,15,6,17))+
  geom_line(aes(group=Method),size=1)+ggtitle("Power Law n=100, mc=100,000")+guides(shape=guide_legend(""))+
  xlab("LCV distribution")+ylab("Mean Absolute Error")+theme_bw()+scale_color_manual(values=c("Green","Black","Red","Blue","darkgrey"))+
  coord_cartesian(ylim=c(-1,1 ))

ggsave("/cluster/tufts/lamontagnelab/fdolan03/EstimatingMean/pl_n100_m105_mae.png",plot=p1,units="in",dpi=300,height=5,width=7)

###################################################

as.data.frame(L)->L
names(L)=c("Sample Mean","LN2 MLE","Median of Means","Power Law MLE","GP MLE","Winsorized Mean","Trimmed Mean")

L$LCV=LCV

######

L %>% #select(`Sample Mean`,Bagged,Bragged,`MOM k=23`,`LN2 MLE`,CV,`Power Law MLE`) %>%
  gather(Method,Lscale,-LCV)->L


L %>% mutate(Type=ifelse(grepl("MLE",Method),Method,"Nonparametric"))->L

L %>% mutate(Estimator=ifelse(Type=="Nonparametric",Method,"MLE"))->L

p2<-ggplot(L,aes(LCV,Lscale,col=Estimator))+geom_point(aes(shape=Type,group=Method),size=3)+scale_shape_manual(values=c(8,15,6,17))+
  geom_line(aes(group=Method),size=1)+ggtitle("Power Law n=100, mc=100,000")+guides(shape=guide_legend(""))+
  xlab("LCV distribution")+ylab("L-Scale")+theme_bw()+scale_color_manual(values=c("Green","Black","Red","Blue","darkgrey"))+
  coord_cartesian(ylim=c(0,.7 ))

ggsave("/cluster/tufts/lamontagnelab/fdolan03/EstimatingMean/pl_n100_m105_Lscale.png",plot=p2,units="in",dpi=300,height=5,width=7)

MAE %>% saveRDS("/cluster/tufts/lamontagnelab/fdolan03/EstimatingMean/powerlaw_mae_n100_m105.rds")
L %>% saveRDS("/cluster/tufts/lamontagnelab/fdolan03/EstimatingMean/powerlaw_Lscale_n100_m105.rds")


##
as.data.frame(B)->B  
names(B)=c("Sample Mean","LN2 MLE","Median of Means","Power Law MLE","GP MLE","Winsorized Mean","Trimmed Mean")

B$LCV=LCV

B %>% gather(Method,Bias,-LCV)->B


B %>% mutate(Type=ifelse(grepl("MLE",Method),Method,"Nonparametric"))->B

B %>% mutate(Estimator=ifelse(Type=="Nonparametric",Method,"MLE"))->B

p3<-ggplot(B,aes(LCV,Bias,col=Estimator))+geom_point(aes(shape=Type,group=Method),size=3)+scale_shape_manual(values=c(8,15,6,17))+
  geom_line(aes(group=Method),size=1)+ggtitle("Power Law n=100, mc=100,000")+guides(shape=guide_legend(""))+
  xlab("LCV distribution")+ylab("Bias")+theme_bw()+scale_color_manual(values=c("Green","Black","Red","Blue","darkgrey"))+
  coord_cartesian(ylim=c(-1,1 ))

ggsave("/cluster/tufts/lamontagnelab/fdolan03/EstimatingMean/pl_n100_m105_bias.png",plot=p3,units="in",dpi=300,height=5,width=7)

B %>% saveRDS("/cluster/tufts/lamontagnelab/fdolan03/EstimatingMean/powerlaw_bias_n100_m105.rds")

