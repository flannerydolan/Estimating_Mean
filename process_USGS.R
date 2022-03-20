library(dataRetrieval)



sites <- read.table('USGSsites.txt',sep='\t')

sites<-sites$V2

sites=as.character(sites)
# site information

for (i in 1:length(sites)){
  
  siteNumber<-sites[i]
  
  if (nchar(siteNumber<8)){siteNumber=paste0("0",siteNumber)}
  
  # parameter codes from "https://help.waterdata.usgs.gov/code/parameter_cd_query?fmt=rdb&group_cd=PHY&inline=true"
  #"00060" means daily mean discharge and the unit is cubic feet per second
  # Raw daily data:
  rawDailyData <- readNWISdv(siteNumber,'00060')
  
  # write out site data
  
  write.table(rawDailyData, file = paste0("datasets/USGS/sites/",i,".txt"))
  
}


#####################################################################################################################

source("ppcc_fn.R")
source("medianofmeans.R")
library(lmom)
library(dplyr)

files=list.files('datasets/USGS/sites/',full.names=TRUE)

sf=c()

for (i in 1:length(files)){
  
  dat=read.table(files[i])
  
  # don't use sites with less than 1000 measurements
  
  if (nrow(dat)>=1000){
    
    if (names(dat)[4]=="X_00060_00003"){
      dat %>% rename(Q=X_00060_00003) %>% select(site_no,Q)->dat }
    
    else{
      q1=c()
      q=c()
      
      for (name in names(dat)){
        
        if (grepl("_00060",name) & !(grepl("_cd",name))){
          
          dat %>% rename(Q=name) %>% select(site_no,Date,Q) -> q
          
          q1=rbind(q1,q)}}
      
      dat<-q1
      
    }
    
    dat %>% na.omit() ->dat
    
    
    dat %>% group_by(site_no) %>% filter(Q>0) %>% summarize(MOM=mom(Q),Sample=mean(Q,na.rm=T),Lskew=lmom::samlmu(Q,nmom=3)[3],
                                                            LCV=lmom::samlmu(Q,nmom=2)[2]/mean(Q), n=length(Q), LNppcc=CDFln(Q), GPppcc=CDFgp(Q), PLppcc=CDFpl(Q),
                                                            perc_err=((Sample-MOM)/MOM)*100,Lkurtosis=lmom::samlmu(Q,nmom=4)[4]) %>% 
      distinct()->sf1
    

    sf=rbind(sf1,sf)
    
    message(i)
  }
  
}

sf %>% saveRDS("USGS_processed.rds")
