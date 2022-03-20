mom<-function(x){
  n=length(x)
  index<-sample(seq(1:n),n,replace=FALSE)
  k=4
  f=rep(1:k,n/k)
  groups<-split(x[index],f,drop=FALSE)
  groups<-sapply(groups,na.omit)
  means<-sapply(groups,mean)
  median(means,na.rm=T) }