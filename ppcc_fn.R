
CDFln<-function(x){
  
  x=na.omit(x)
  x=as.numeric(x)
  Lx<-as.numeric(lmom::samlmu(x,nmom=4))
  
  LCV=Lx[2]/Lx[1]
  
  vy=(sqrt(2)*qnorm((1+LCV)/2))^2
  
  uy=log(Lx[1])-vy/2
  
  
  Fln=plnorm(x, meanlog=uy, sdlog=sqrt(vy))
  
  cor(seq(1:length(x))/(length(x)+1), sort(Fln))
}



CDFpl<-function(x){
  
  x=na.omit(x)
  xmin=min(x,na.rm=T)
  x=as.numeric(x)
  Lx<-as.numeric(lmom::samlmu(x,nmom=4))
  
  if (xmin==0){
    xmin=sort(unique(x))[2]
  }
  
  LCV=Lx[2]/Lx[1]
  
  alpha=(1+3*LCV)/(2*LCV)
  
  C = (alpha-1)*xmin^(alpha-1)
  
  Fpl = 1-(C/(alpha-1))*x^-(alpha-1)
  
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
