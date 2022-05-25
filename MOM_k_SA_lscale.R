# set parameters 
ux=1
m=1000
# change parameters manually
n=2000
LCV=0.1

library(dplyr)
library(tidyr)
library(ggplot2)
# parameters varying with LCV

#ln2
vy=(sqrt(2)*qnorm((1+LCV)/2))^2
uy=log(ux)-vy/2

# gpa
kappa=(1-2*LCV)/LCV
alpha=ux*(kappa+1)

# power law
a=(1+3*LCV)/(2*LCV)
xmin=ux*(a-2)/(a-1)

# functions
library(lmom)

#mae<-function(x){mean(abs(x-ux),na.rm=T)}

lscale<-function(x){lmom::samlmu(x,nmom=2)[2]}
#rmse<-function(x){sqrt(sum((x-ux)^2)/m)}  

# initialize array

E=matrix(NA, nrow=(n/2), ncol=3)

# Run Mean absolute error first

# loop through distributions and k groups

for (i in 1:3){
  
  if (i==1){dist=sapply(1:m, function(x) {p=runif(n,0,1)
  exp(uy+sqrt(vy)*qnorm(p))    })
  message("Lognormal")}
  
  else if (i==2){dist=sapply(1:m, function(x) {p=runif(n,0,1)
  if (kappa==0){
    -alpha*log(1-p)}
  else{alpha*(1-(1-p)^kappa)/kappa} })
  message("GPA")}
  
  else{dist=sapply(1:m, function(x) {p=runif(n,0,1)
  exp((-log((-p+1)/exp((a-1)*log(xmin))))/(a-1))
  })
  message("Power") }
  
  for (k in 1:(n/2)){
    
    
    mom=sapply(1:ncol(dist), function(y){
      index<-sample(seq(1:n),n,replace=FALSE)
      f=rep(1:k,n/k)
      groups<-split(dist[index,y],f,drop=FALSE)
      means<-sapply(groups,mean)
      median(means)})
    
    message(k)
    
    E[k,i]=lscale(mom)
  }
}

E=as.data.frame(E)

names(E)=c("Lognormal","Generalized Pareto","Power Law")
E$Groups=seq(1,(n/2))
E$LCV="0.1"

##################################
LCV=0.5

# parameters varying with LCV

#ln2
vy=(sqrt(2)*qnorm((1+LCV)/2))^2
uy=log(ux)-vy/2

# gpa
kappa=(1-2*LCV)/LCV
alpha=ux*(kappa+1)

# power law
a=(1+3*LCV)/(2*LCV)
xmin=ux*(a-2)/(a-1)



# initialize array

E2=matrix(NA, nrow=(n/2), ncol=3)

# loop through distributions and k groups

for (i in 1:3){
  
  if (i==1){dist=sapply(1:m, function(x) {p=runif(n,0,1)
  exp(uy+sqrt(vy)*qnorm(p))    })
  message("Lognormal")}
  
  else if (i==2){dist=sapply(1:m, function(x) {p=runif(n,0,1)
  if (kappa==0){
    -alpha*log(1-p)}
  else{alpha*(1-(1-p)^kappa)/kappa} })
  message("GPA")}
  
  else{dist=sapply(1:m, function(x) {p=runif(n,0,1)
  exp((-log((-p+1)/exp((a-1)*log(xmin))))/(a-1))
  })
  message("Power") }
  
  for (k in 1:(n/2)){
    
    
    mom=sapply(1:ncol(dist), function(y){
      index<-sample(seq(1:n),n,replace=FALSE)
      f=rep(1:k,n/k)
      groups<-split(dist[index,y],f,drop=FALSE)
      means<-sapply(groups,mean)
      median(means)})
    
    message(k)
    
    E2[k,i]=lscale(mom)
  }
}

E2=as.data.frame(E2)

names(E2)=c("Lognormal","Generalized Pareto","Power Law")
E2$Groups=seq(1,(n/2))
E2$LCV="0.5"


#############################################
LCV=0.9

# parameters varying with LCV

#ln2
vy=(sqrt(2)*qnorm((1+LCV)/2))^2
uy=log(ux)-vy/2

# gpa
kappa=(1-2*LCV)/LCV
alpha=ux*(kappa+1)

# power law
a=(1+3*LCV)/(2*LCV)
xmin=ux*(a-2)/(a-1)

# functions


# initialize array

E3=matrix(NA, nrow=(n/2), ncol=3)

# loop through distributions and k groups

for (i in 1:3){
  
  if (i==1){dist=sapply(1:m, function(x) {p=runif(n,0,1)
  exp(uy+sqrt(vy)*qnorm(p))    })
  message("Lognormal")}
  
  else if (i==2){dist=sapply(1:m, function(x) {p=runif(n,0,1)
  if (kappa==0){
    -alpha*log(1-p)}
  else{alpha*(1-(1-p)^kappa)/kappa} })
  message("GPA")}
  
  else{dist=sapply(1:m, function(x) {p=runif(n,0,1)
  exp((-log((-p+1)/exp((a-1)*log(xmin))))/(a-1))
  })
  message("Power") }
  
  for (k in 1:(n/2)){
    
    
    mom=sapply(1:ncol(dist), function(y){
      index<-sample(seq(1:n),n,replace=FALSE)
      f=rep(1:k,n/k)
      groups<-split(dist[index,y],f,drop=FALSE)
      means<-sapply(groups,mean)
      median(means)})
    
    message(k)
    
    E3[k,i]=lscale(mom)
  }
}

E3=as.data.frame(E3)

names(E3)=c("Lognormal","Generalized Pareto","Power Law")
E3$Groups=seq(1,(n/2))
E3$LCV="0.9"

rbind(E,E2,E3)->lscale


lscale %>% gather(key=Distribution,value=lscale,Lognormal:`Power Law`)->lscale

lscale %>% saveRDS("MOMkSA_n2000_lscale.rds")

ggplot(lscale,aes(Groups,lscale,col=LCV,linetype=Distribution))+theme_bw()+coord_cartesian(ylim=c(0,1))+ylab("L-Scale")+
  geom_line(size=1)+ggtitle("n=2000, m=1,000")+scale_color_manual(values=c("black","grey","red"))+scale_x_continuous(trans='log10')

ggsave("MOMkSA_n2000_logged_lscale.png",dpi=300)


