#set home dir of pipeline
home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline"

# Purpose ----------------------------------------------------------------------
#The purpose of this script is to format outputs to use for analysis in removals
#responses pipeline

#1. AKDE output formatting
  #relevel, remove missing periods, etc.
#2. Geolocation data formatting
  #General tidying
  #set periods
#3. Calculate NSD
#4. Summarize weekly median NSD

# In/Out -----------------------------------------------------------------------

#for NSD calcs:
#geo_remtype.csv --> geo_aerd_wk.rds, geo_toxd_wk.rds, geo_trapd_wk.rds

#for AKDE calcs
#outdf_akde_aer_corrected.rds --> outdf_akde_aer_corrected_f.rds
#outdf_akde_tox_corrected.rds --> outdf_akde_tox_corrected_f.rds
#outdf_akde_trap_corrected.rds --> outdf_akde_trap_corrected_f.rds

# Setup ------------------------------------------------------------------------

#load libraries

#Set dirs
input=file.path(home,"1_Data","Input",fsep=.Platform$file.sep)
objdir=file.path(home,"1_Data","Objects",fsep=.Platform$file.sep)

#read activities file, for getting cutoff dates during field work
activities.raw=readRDS(file.path(input,"activities.rds",fsep=.Platform$file.sep))

#read home range area outputs
akaer=readRDS(file.path(objdir,"outdf_akde_aer_corrected.rds",fsep=.Platform$file.sep))
aktrap=readRDS(file.path(objdir,"outdf_akde_trap_corrected.rds",fsep=.Platform$file.sep))
aktox=readRDS(file.path(objdir,"outdf_akde_tox_corrected.rds",fsep=.Platform$file.sep))

#read geolocations with removal type designations:
geo=read.csv(file.path(objdir,"geo_remtyp.csv",fsep=.Platform$file.sep))

#source needed functions
func.list=list.files(file.path(home,"2_Scripts","Functions",fsep=.Platform$file.sep),full.names=TRUE)
for(f in 1:length(func.list)){
  source(func.list[f])
}

# Area df object data formatting --------------------------------------------------

#fix units-- sq meters/hectares-- want all to be sq km
Fix.Area.Units<-function(akrem){
  akrem[grep("square meters",akrem$akde.hrarea.units),][3:5]=akrem[grep("square meters",akrem$akde.hrarea.units),][3:5]/1e6
  akrem[grep("hectares",akrem$akde.hrarea.units),][3:5]=akrem[grep("hectares",akrem$akde.hrarea.units),][3:5]/100
  return(akrem)
}

aktox=Fix.Area.Units(aktox)
aktrap=Fix.Area.Units(aktrap)
akaer=Fix.Area.Units(akaer)

#Remove pigs that were killed during aerial
killed.aer=akaer[is.na(akaer$area.CI.low),]$animalid
akaer=akaer[!(akaer$animalid%in%killed.aer),]

#Remove NA 'after' periods for tox-killed pigs
aktox=aktox[!is.na(aktox$area.CI.low),]

#Relevel periods
akaer$period<-forcats::fct_relevel(akaer$period,c("before","after"))
aktrap$period<-forcats::fct_relevel(aktrap$period,c("before","during","after"))
aktox$period<-forcats::fct_relevel(aktox$period,c("before","after"))

#Relevel removal types
akaer$Removal.Type<-forcats::fct_relevel(akaer$Removal.Type,c("ctrl","aer"))
aktrap$Removal.Type<-forcats::fct_relevel(aktrap$Removal.Type,c("ctrl","trap"))
aktox$Removal.Type<-forcats::fct_relevel(aktox$Removal.Type,c("ctrl","tox"))

#write out formatted versions
saveRDS(akaer,file.path(objdir,"outdf_akde_aer_corrected_f.rds",fsep=.Platform$file.sep))
saveRDS(aktrap,file.path(objdir,"outdf_akde_trap_corrected_f.rds",fsep=.Platform$file.sep))
saveRDS(aktox,file.path(objdir,"outdf_akde_aktox_corrected_f.rds",fsep=.Platform$file.sep))

# Geolocation data formatting --------------------------------------------------

#Basic tidying
geo=geo[,-1] #First X sequence col remove
geo=geo[!is.na(geo$Removal.Type),] #remove undesignated pigs
geo$date_only<-Neat.Dates.POSIXct(geo$date_only,tz="UTC")
  
#Set dates for before, during, after periods for each removal type
#trap dates
trap.act=activities.raw[activities.raw$activity=="trap",]
trap.start.date=as.Date(min(trap.act$datetime))
trap.end.date=as.Date(max(trap.act$datetime))

#tox dates
tox.act=activities.raw[activities.raw$activity=="toxic",]
tox.start.date=as.Date(min(tox.act$datetime))
tox.end.date=as.Date("2023-03-09")

#aerial dates
aer.start.date=as.Date("2023-03-27")
aer.end.date=as.Date("2023-03-29")

#Get cutoffs for each removal type
#Aerial 
#determined 4 weeks from sensitivity analysis for how much time needed to get accurate akde est
aer.cutoff1=aer.start.date-37
aer.cutoff2=aer.end.date+37

#Trap
#use length of time in trap after period (smaller than during period)
trap.len=as.Date(max(geo[geo$Removal.Type=="trap",]$date_only))-as.Date(trap.end.date)-1
#use as trap start date, 56 days before trap end date
#think this is when trapping really started more intensively, anyways
trap.start.date.akde=trap.end.date-trap.len-1
trap.cutoff1=trap.start.date-trap.len-1
trap.cutoff2=trap.end.date+trap.len+1
  
#Tox
#use length of time in tox period
#tox.len=tox.end.date-tox.start.date
tox.len=36
tox.cutoff1=tox.start.date-tox.len-1
tox.cutoff2=tox.end.date+tox.len+1
  
Do.Date.Cutoffs<-function(geo,removal,cutoff1,cutoff2,startdate,enddate){
    geo.rem=geo[geo$Removal.Type==removal|geo$Removal.Type=="ctrl",]
    geo.rem.before=geo.rem[geo.rem$date_only>=cutoff1&geo.rem$date_only<startdate,]
    geo.rem.after=geo.rem[geo.rem$date_only>enddate&geo.rem$date_only<=cutoff2,]
    geo.rem.during=geo.rem[geo.rem$date_only>=startdate&geo.rem$date_only<=enddate,]
    geo.rem.before$removal.period.akdecalc="before"
    geo.rem.after$removal.period.akdecalc="after"
    geo.rem.during$removal.period.akdecalc="during"
    geo.rem=rbind(geo.rem.before,geo.rem.during,geo.rem.after)
    return(geo.rem)
  }
geo.aer=Do.Date.Cutoffs(geo,"aer",aer.cutoff1,aer.cutoff2,aer.start.date,aer.end.date)
geo.tox=Do.Date.Cutoffs(geo,"tox",tox.cutoff1,tox.cutoff2,tox.start.date,tox.end.date)
  
#do trap ones manually
geo.rem=geo[geo$Removal.Type=="trap"|geo$Removal.Type=="ctrl",]
geo.rem.before=geo.rem[geo.rem$date_only>=trap.cutoff1&geo.rem$date_only<trap.start.date,]
geo.rem.after=geo.rem[geo.rem$date_only>trap.end.date&geo.rem$date_only<=trap.cutoff2,]
geo.rem.during=geo.rem[geo.rem$date_only>trap.start.date.akde&geo.rem$date_only<=trap.end.date,]
geo.rem.before$removal.period.akdecalc="before"
geo.rem.after$removal.period.akdecalc="after"
geo.rem.during$removal.period.akdecalc="during"
geo.rem=rbind(geo.rem.before,geo.rem.during,geo.rem.after)
geo.trap=geo.rem
  
#set datetime format
geo.aer$datetime<-Neat.Dates.POSIXct(geo.aer$datetime,tz="UTC")
geo.trap$datetime<-Neat.Dates.POSIXct(geo.trap$datetime,tz="UTC")
geo.tox$datetime<-Neat.Dates.POSIXct(geo.tox$datetime,tz="UTC")
  
#make join key
geo.aer$joinkey=paste(geo.aer$animalid,geo.aer$removal.period.akdecalc,sep="_")
geo.trap$joinkey=paste(geo.trap$animalid,geo.trap$removal.period.akdecalc,sep="_")
geo.tox$joinkey=paste(geo.tox$animalid,geo.tox$removal.period.akdecalc,sep="_")
  
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

#loop through each pig, get disps, rbind together, output geo.rem with disp
#geo.rem=geo.aer
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

# Summarize NSD by weekly median ----------------------------------------------------------------

#redo week num summary
#confirm all start on same date
geo.aerd %>% group_by(animalid) %>% 
  dplyr::summarise(strtdt=min(date_only)) %>%
  as.data.frame()
geo.trapd %>% group_by(animalid) %>% 
  dplyr::summarise(strtdt=min(date_only)) %>%
  as.data.frame()
geo.toxd %>% group_by(animalid) %>% 
  dplyr::summarise(strtdt=min(date_only)) %>%
  as.data.frame()

#1-get start date of removal period
#date 2-5 weeks/35 days prior to that

#set origins
#same periods as akde area calcs
#trim to be multiple of 7
trap.origin=trap.start.date-56
trap.origin.after=trap.end.date+1
trap.origin.during=trap.end.date-56

tox.origin=tox.start.date-35
tox.origin.after=tox.end.date+1
tox.origin.during=tox.start.date

aer.origin=aer.start.date-35
aer.origin.after=aer.end.date+1

#do week summaries for each, starting at day 1
#separate each period
Do.Week.Split<-function(geo.aerd,removal.str,origin.vector){
  #c(aer.origin,aer.origin.after)
  
  geo.aerd.before=geo.aerd[geo.aerd$removal.period.akdecalc=="before",]
  geo.aerd.after=geo.aerd[geo.aerd$removal.period.akdecalc=="after",]
  if(removal.str!="aer"){
    geo.aerd.during=geo.aerd[geo.aerd$removal.period.akdecalc=="during",]
  }
  
  #set jDate according to origins
  geo.aerd.before$jDate <- julian(as.Date(geo.aerd.before$datetime), origin = origin.vector[1])
  geo.aerd.before$jDate=geo.aerd.before$jDate+1
  geo.aerd.after$jDate <- julian(as.Date(geo.aerd.after$datetime), origin = origin.vector[2])
  geo.aerd.after$jDate=geo.aerd.after$jDate+1
  if(removal.str!="aer"){
    geo.aerd.during$jDate=julian(as.Date(geo.aerd.during$datetime), origin = origin.vector[3])
  }
  #set weeks
  geo.aerd.before=transform(geo.aerd.before, week = cut(jDate, c(0, seq(7, max(jDate) + 7, 7)), labels = FALSE))
  geo.aerd.before=geo.aerd.before[!is.na(geo.aerd.before$week),]
  geo.aerd.after=transform(geo.aerd.after, week = cut(jDate, c(0, seq(7, max(jDate) + 7, 7)), labels = FALSE))
  
  if(removal.str!="aer"){
    geo.aerd.during=transform(geo.aerd.during, week = cut(jDate, c(0, seq(7, max(jDate) + 7, 7)), labels = FALSE))
    geo.aerd.during$week=geo.aerd.during$week+max(geo.aerd.before$week)
    geo.aerd.during=geo.aerd.during[!is.na(geo.aerd.during$week),]
    geo.aerd.after$week=geo.aerd.after$week+max(geo.aerd.during$week)
    
  } else{
    geo.aerd.after$week=geo.aerd.after$week+max(geo.aerd.before$week)
  }
  
  #bind back together
  if(removal.str!="aer"){
    geo.aerd2=rbind(geo.aerd.before,geo.aerd.after,geo.aerd.during)
  } else{
    geo.aerd2=rbind(geo.aerd.before,geo.aerd.after)
  }
  
  return(geo.aerd2)
}

geo.aerd2=Do.Week.Split(geo.aerd,"aer",c(aer.origin,aer.origin.after))
geo.trapd2=Do.Week.Split(geo.trapd,"trap",c(trap.origin,trap.origin.after,trap.origin.during))
geo.toxd2=Do.Week.Split(geo.toxd,"tox",c(tox.origin,tox.origin.after,tox.origin.during))

#verify that no weeks are very short

#Aerial
aerchk=geo.aerd2 %>% 
  group_by(animalid,removal.period.akdecalc,week) %>% 
  dplyr::summarise(n_distinct(jDate)) %>%
  as.data.frame()
#remove week 11, only two days
geo.aerd2=geo.aerd2[geo.aerd2$week!=11,]

#Trap
trapchk=geo.trapd2 %>% 
  group_by(animalid,removal.period.akdecalc,week) %>% 
  dplyr::summarise(n_distinct(jDate)) %>%
  as.data.frame()
#remove week 25, only 1 days
geo.trapd2=geo.trapd2[geo.trapd2$week!=25,]

#Tox
toxchk=geo.toxd2 %>% 
  group_by(animalid,removal.period.akdecalc,week) %>% 
  dplyr::summarise(n_distinct(jDate)) %>%
  as.data.frame()
#remove week 14
geo.toxd2=geo.toxd2[geo.toxd2$week!=14,]

#remove pig/weeks with very short weeks that died during toxicant
#say any with <5 days of geolocations on week 8
rem.tox.wk8.id=toxchk[toxchk$week==8&toxchk$`n_distinct(jDate)`<5,]$animalid
geo.toxd2=geo.toxd2[!(geo.toxd2$animalid%in%rem.tox.wk8.id&geo.toxd2$week==8),]

#recheck
toxchk=geo.toxd2 %>% 
  group_by(animalid,removal.period.akdecalc,week) %>% 
  dplyr::summarise(n_distinct(jDate)) %>%
  as.data.frame()

#summarize by week
geo.aerd.wk=geo.aerd2 %>% group_by(animalid, Removal.Type, removal.period.akdecalc,sex, week) %>% dplyr::summarise(mNSD=median(NSD),mX=mean(X),mY=mean(Y)) %>% as.data.frame()
geo.trapd.wk=geo.trapd2 %>% group_by(animalid, Removal.Type, removal.period.akdecalc,sex, week) %>% dplyr::summarise(mNSD=median(NSD),mX=mean(X),mY=mean(Y)) %>% as.data.frame()
geo.toxd.wk=geo.toxd2 %>% group_by(animalid, Removal.Type, removal.period.akdecalc,sex, week) %>% dplyr::summarise(mNSD=median(NSD),mX=mean(X),mY=mean(Y)) %>% as.data.frame()

#Remove pigs from aer that were accidentally culled
aer.cull.sum=geo.aerd.wk %>% 
  group_by(animalid) %>% 
  dplyr::summarise(np=n_distinct(removal.period.akdecalc))
aer.culled=aer.cull.sum[aer.cull.sum$np<2,]$animalid
geo.aerd.wk=geo.aerd.wk[!(geo.aerd.wk$animalid%in%aer.culled),]

#Write out data to use for NSD analysis:
saveRDS(geo.aerd.wk,file.path(objdir,"geo_aerd_wk.rds"))
saveRDS(geo.trapd.wk,file.path(objdir,"geo_trapd_wk.rds"))
saveRDS(geo.toxd.wk,file.path(objdir,"geo_toxd_wk.rds"))


