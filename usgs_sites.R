library(dplyr)
library(lmom)
library(ggplot2)
library(purrr)
# Median of Means
k=4

mom<-function(x){
n=length(x)
index<-sample(seq(1:n),n,replace=FALSE)
k=4
f=rep(1:k,n/k)
groups<-split(x[index],f,drop=FALSE)
groups<-sapply(groups,na.omit)
means<-sapply(groups,mean)
median(means,na.rm=T) }

SF=c()

files=list.files(path="datasets/USGS",pattern="streamflow",full.names=T)

files %>% map(read.csv) %>% reduce(rbind)->SF

SF->sf
#sf<-read.csv('datasets/USGS/result_01_streamflow_cfs.csv')

SF=na.omit(SF)

SF %>% group_by(staid) %>% summarize(MOM=mom(obsQ),Sample=mean(obsQ,na.rm=T),
                                     Lskew=lmom::samlmu(obsQ,nmom=3)[3],n=length(obsQ),
                                     perc_err=((Sample-MOM)/MOM)*100) ->sf

sf=distinct(sf)

sf %>% tidyr::gather(key=Estimator, value=value, MOM:Sample)->sf

ggplot(sf,aes(Lskew,(value)))+geom_point(aes(col=Estimator),alpha=0.3,size=1)+geom_line(aes(group=staid))+theme_bw()+
  ylab("Estimate of Mean Streamflow (cfs)")+scale_color_manual(values=c("red","blue"))+scale_y_continuous(trans="log10")

ggsave("usgs_MOM_sample_logged_scale.png",dpi=300)


ggplot(sf,aes(n,(value)))+geom_point(aes(col=Estimator),alpha=0.3,size=1)+geom_line(aes(group=staid))+theme_bw()+
  ylab("Estimate of Mean Streamflow (cfs)")+scale_color_manual(values=c("red","blue"))+scale_y_continuous(trans="log10")+xlab("Sample Size (n)")

ggsave("usgs_MOM_sample_vs_n.png",dpi=300)

ggplot(sf,aes(P0,value))+geom_point(aes(col=Estimator),alpha=0.3,size=1)+geom_line(aes(group=staid))+theme_bw()+
  ylab("Estimate of Mean Streamflow (cfs)")+scale_color_manual(values=c("red","blue"))+scale_y_continuous(trans="log10")+xlab("Probability of Zeros")

ggsave("usgs_MOM_sample_vs_n.png",dpi=300)


ggplot(sf,aes(Lskew,(value)))+geom_point(aes(col=Estimator),alpha=0.2)+geom_line(aes(group=staid))+theme_bw()+
  ylab("Estimate of Mean Streamflow (cfs)")+scale_color_manual(values=c("red","blue"))

ggsave("usgs_MOM_sample_realspace.png",dpi=300)



#######################################
SF->sf
sf=na.omit(sf)

SF %>% group_by(staid) %>% dplyr::summarize(MOM=mom(obsQ),Sample=mean(obsQ,na.rm=T),Lskew=lmom::samlmu(obsQ,nmom=3)[3],
                                    LCV=lmom::samlmu(obsQ,nmom=2)[2]/mean(obsQ), n=length(obsQ),zero=ifelse(obsQ==0,1,0),
                                    P0=sum(zero)/n, LNppcc=CDFln(obsQ), GPppcc=CDFgp(obsQ), PLppcc=CDFpl(obsQ),
                                    perc_err=((Sample-MOM)/MOM)*100,Lkurtosis=lmom::samlmu(obsQ,nmom=4)[4]) ->sf1


#sf1 %>% tidyr::gather(key=Estimator, value=value, MOM:Sample)->sf

sf1 %>% distinct() ->sf1

sf1 %>% write.csv("usgs_site_dat_newer.csv")

library(ggExtra)
library(scales)

ggplot(sf,aes(Lskew,perc_err))+geom_point(size=2,alpha=.6,aes(col=n))+theme_bw()+
  ylab("Percent Error")+scale_color_viridis(option="D")+scale_y_continuous(trans="log10",labels=trans_format("log10",math_format()))+xlab("Lskew")+
  guides(color=guide_legend(title="Record Length"))

ggsave("errorY_recordlengthX_lskew.png",dpi=300)


sf %>% select(-Estimator,-value) %>% distinct() %>% ungroup() %>% mutate(bad=ifelse(perc_err>=100,1,0)) %>% summarize(badperc=sum(bad)/nrow(.))

# 14% of sites exhibit errors of 10% or over, 8% had 100% error or over, 0.7% had 1000% error or over

ggplot(sf,aes(Lskew,value))+geom_point(aes(col=Estimator),alpha=0.3,size=1)+geom_line(aes(group=staid))+theme_bw()+
  ylab("Estimate of Mean Streamflow (cfs)")+scale_color_manual(values=c("red","blue"))+scale_y_continuous(trans="log10")+xlab("L-Skew")


p1<-ggMarginal(p, type="density",margins='x',size=4,color="purple")

sf1 %>% mutate(good=ifelse(LNppcc>=.99,1,0)) %>% ungroup() %>% summarize(perc_good=sum(good)/n())

distinct(sf1)->sf1

sf %>% tidyr::gather(key=Distribution,value=PPCC,LNppcc:PLppcc)->sf2


ggplot(sf1,aes(Lskew,perc_err,col=Lkurtosis))+geom_point()+scale_color_viridis()+theme_bw()+
  scale_y_continuous(trans="log10")+ylab("Relative Percent Error")

ggsave("usgs_lskew_percerr_lkurt.png",dpi=300)

sf2 %>% group_by(staid) %>% slice(which.max(PPCC)) %>% 
  mutate(Distribution=substr(Distribution,1,2)) ->sf3

sf3 %>% mutate(ppcc = ifelse(PPCC>=0.99,">=0.99",ifelse(PPCC>=0.9,">=0.9","<0.9")))->sf3

sf3 %>% mutate(good=ifelse(PPCC>=.99,1,0)) %>% ungroup() %>% summarize(percgood=sum(good)/1168)

ggplot(sf3, aes(Lskew, perc_err,col=log10(MOM)))+geom_point(size=2,alpha=.6)+theme_bw()+scale_y_continuous(trans="log10")+
  scale_color_viridis(discrete=F)+ylab("Relative Percent Error")#+guides(color=guide_legend(reverse=TRUE))

ggsave("usgs_lskew_err_lcv.png",dpi=300)



CDFln<-function(x){

x=na.omit(x)
uy=mean(log(x),na.rm=T)
vy=var(log(x),na.rm=T)
ux=exp(uy+vy/2)
vx=(exp(vy)-1)*(exp(2*uy+vy))
Fln=plnorm(x, meanlog=uy, sdlog=sqrt(vy))

cor(seq(1:length(x))/(length(x)+1), sort(Fln))
}

CDFpl<-function(x){
  
  x=na.omit(x)
  xmin=min(x,na.rm=T)
  
  if (xmin==0){
    xmin=sort(unique(x))[2]
  }
  
  alpha=1+length(x)*(1/(sum(log(x/xmin),na.rm=T))) # Newman eq 5
  
  C = (alpha-1)*xmin^(alpha-1)
  
  Fpl = (C/(alpha-1))*x^-(alpha-1)
  
  cor(seq(1:length(x))/(length(x)+1), sort(Fpl))
}


CDFgp<-function(x){
  
  x=na.omit(x)
  x=as.numeric(x)
  Lx<-as.numeric(lmom::samlmu(x,nmom=4))
  k=(Lx[1]-2*Lx[2])/Lx[2] ## Hosking and Wallis Appendix
  
  a=Lx[1]*(1+k)
  
  y=seq(1:length(x))/(length(x)+1)
  
  Fgp = 1 - (1-(k*x)/a)^(1/k) # Hosking and Wallis 1987
  
  y[is.na(Fgp)]<-NA
  
  na.omit(y)->y
  na.omit(Fgp)->Fgp
  
  as.numeric(y)->y
  as.numeric(Fgp)->Fgp
  
  cor(y, sort(Fgp),use = "complete.obs")
}


ggplot(sf1,aes(Lskew,perc_err,col=n))+geom_point()




#####

files=list.files('datasets/USGS/sites',full.names = T)

#sites=c()

for (i in 1:length(files)){
  
  
  file=read.table(files[i])
  
  if (nrow(file)>1){
 # sites[i]=unique(file$site_no)
    
    file %>% write.table(paste0('datasets/USGS/newnames/',unique(file$site_no),'.txt'))
  
  message(i)
  }
}

which(sites=="2357000")
