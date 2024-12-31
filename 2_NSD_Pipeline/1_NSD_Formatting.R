#set home dir of pipeline
home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline"

# Purpose ----------------------------------------------------------------------

#Format datasets to run NSD calculations for analysis

# In/Out -----------------------------------------------------------------------

#geo.tox, geo.aer, geo.trap inputs
#geo.aerd.wk, geo.toxd.wk, geo.trapd.wk outputs

# Setup ------------------------------------------------------------------------

#load libraries
library(amt)
library(dplyr)
library(stringr)
library(sf)
library(forcats)

#Set dirs
input=file.path(home,"1_Data","Input",fsep=.Platform$file.sep)
objdir=file.path(home,"1_Data","Objects",fsep=.Platform$file.sep)

#Read geoloc data
geo.tox<-readRDS(file.path(objdir,"geotox.rds"))
geo.aer<-readRDS(file.path(objdir,"geoaer.rds"))
geo.trap<-readRDS(file.path(objdir,"geotrap.rds"))

#source needed functions
func.list=list.files(file.path(home,"2_Scripts","Functions",fsep=.Platform$file.sep),full.names=TRUE)
for(f in 1:length(func.list)){
  source(func.list[f])
}

# Trim incomplete weeks ------------------------------------------------------

#Make function to trim incomplete weeks
trimwks=function(geo.rem,mindays){
  remchk=geo.rem %>% 
    group_by(animalid,removal.period.akdecalc,week) %>% 
    dplyr::summarise(minday=n_distinct(jDate)) %>%
    as.data.frame()
  
  remchk2=remchk[remchk$minday<mindays,]
  
  for(i in 1:nrow(remchk2)){
    geo.rem=geo.rem[!(geo.rem$animalid==remchk2$animalid[i]&
              geo.rem$removal.period.akdecalc==remchk2$removal.period.akdecalc[i]&
              geo.rem$week==remchk2$week[i]),]
    
  }
  
  return(geo.rem)
  
}

#Make another function to adjust week intervals
adjust_intvals<-function(geo.rem){
  ids=unique(geo.rem$animalid)
  for(i in 1:length(ids)){
    geo.rem_i=geo.rem[geo.rem$animalid==ids[i],]
    key1=unique(geo.rem_i$week)
    key2=1:length(key1)
    
    geo.rem_
    
  }
  
  
}

#Trim incomplete weeks
geo.aer=trimwks(geo.aer,6)
geo.trap=trimwks(geo.trap,6)
geo.tox=trimwks(geo.tox,6)

#View summary
nrow(geo.tox) #315490
nrow(geo.trap) #588735
nrow(geo.aer) #225757

# Calculate NSD ----------------------------------------------------------------

#Make functions to calculate NSD
#pigex is pig/period subset from geo.rem dfs
do.NSDcalcs<-function(pigex){
  pigex=pigex[order(pigex$datetime),]
  
  #lagged loc
  pigex$Xstart=pigex$X[1]
  pigex$Ystart=pigex$Y[1]
  #pigex$datetime2=dplyr::lead(pigex$datetime)
  
  pigex$displacement=sqrt((pigex$X-pigex$Xstart)^2+(pigex$Y-pigex$Ystart)^2)
  pigex$NSD=(pigex$displacement)^2
  return(pigex)
}

#loop through each pig, get NSD, rbind together, output geo.rem with NSD
do.NSDcalcs.georem=function(geo.rem){
  pigID=unique(geo.rem$animalid)
  
  for(i in 1:length(pigID)){
    
    print(paste0("calculating NSD for ",pigID[i]," (pig ",i," of ",length(pigID),")"))
    geo.pig=geo.rem[geo.rem$animalid==pigID[i],]
    geo.pig[order(geo.pig$datetime),]
    
    geo.pig2=do.NSDcalcs(geo.pig)
    
    if(i==1){
      geo.rem.out=geo.pig2
    } else{
      geo.rem.out=rbind(geo.rem.out,geo.pig2)
    }
  }
  
  return(geo.rem.out)
  
}

geo.aerd=do.NSDcalcs.georem(geo.aer)
geo.trapd=do.NSDcalcs.georem(geo.trap)
geo.toxd=do.NSDcalcs.georem(geo.tox)

# Summarize NSD by weekly median -----------------------------------------------

#summarize by week
geo.aerd.wk=geo.aerd %>% group_by(animalid, Removal.Type, removal.period.akdecalc,sex, week) %>% dplyr::summarise(mNSD=median(NSD),mX=mean(X),mY=mean(Y)) %>% as.data.frame()
geo.trapd.wk=geo.trapd %>% group_by(animalid, Removal.Type, removal.period.akdecalc,sex, week) %>% dplyr::summarise(mNSD=median(NSD),mX=mean(X),mY=mean(Y)) %>% as.data.frame()
geo.toxd.wk=geo.toxd %>% group_by(animalid, Removal.Type, removal.period.akdecalc,sex, week) %>% dplyr::summarise(mNSD=median(NSD),mX=mean(X),mY=mean(Y)) %>% as.data.frame()

#Relevel periods
geo.aerd.wk$removal.period.akdecalc<-fct_relevel(geo.aerd.wk$removal.period.akdecalc,c("before","during"))
geo.trapd.wk$removal.period.akdecalc<-fct_relevel(geo.trapd.wk$removal.period.akdecalc,c("before","during"))
geo.toxd.wk$removal.period.akdecalc<-fct_relevel(geo.toxd.wk$removal.period.akdecalc,c("before","during"))

#Relevel removal types
geo.aerd.wk$Removal.Type<-fct_relevel(geo.aerd.wk$Removal.Type,"ctrl")
geo.trapd.wk$Removal.Type<-fct_relevel(geo.trapd.wk$Removal.Type,"ctrl")
geo.toxd.wk$Removal.Type<-fct_relevel(geo.toxd.wk$Removal.Type,"ctrl")

# Write out data ---------------------------------------------------------------

saveRDS(geo.aerd.wk,file.path(objdir,"NSDgeoaer.rds",fsep=.Platform$file.sep))
saveRDS(geo.trapd.wk,file.path(objdir,"NSDgeotrap.rds",fsep=.Platform$file.sep))
saveRDS(geo.toxd.wk,file.path(objdir,"NSDgeotox.rds",fsep=.Platform$file.sep))
#replace geo_toxd_wk.rds
