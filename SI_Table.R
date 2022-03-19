
#######################################
MnEst=matrix(NA,nrow=8,ncol=7)

for (i in 1:8){
  
  if (i==1){cities->x} #1
  
  else if (i==2){fires->x} #2
  
  else if(i==3){flare->x} #3
  
  else if(i==4){terr->x} #4
  
  else if(i==5){honeyss->x} #5
  
  else if (i==6){honeyf->x}
  
  else if (i==7){aid->x}
  
  else {wh->x}
  
  ###
  n=length(x)
  Means=c()
  
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
  
  delta=0.01
  
  trimmed<-function(x){
    q=quantile(x,c(.05,.95))
    x[x<q[1]]<-NA
    x[x>q[2]]<-NA
    mean(x,na.rm=T)
  }
  
  Means[5]=trimmed(x)
  
  # winsorized mean
  
  winsor<-function(x){q=quantile(x,c(.05,.95))
  x[x<q[1]]<-q[1]
  x[x>q[2]]<-q[2]
  mean(x,na.rm=T)}
  
  Means[6]=winsor(x)
  
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