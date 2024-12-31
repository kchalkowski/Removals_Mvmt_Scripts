
########Sensitivity analyses for Removals/Mvmt ms
#1. subsampling analysis
#2. leave one out analysis

########Subsampling
#Steps:
#1. get sample of 10 pigs with 20 min fix rates
#2. resample with amt- 20 mins and 4 hrs
#3. run akde and NSD calcs
#4. output df with full and shortened versions

library(amt)
library(ctmm)
library(ggplot2)
library(segclust2d)
library(NIFApackagev1.1)
library(stringr)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggeffects)
library(hrbrthemes)
library(glmmTMB)

#source needed functions
funcdir<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/2_Scripts/Functions"
func.list=list.files(funcdir,full.names=TRUE)
for(f in 1:length(func.list)){
  source(func.list[f])
}

#read objects from input
input="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/1_Data/Input/"
sitelocs=readRDS(paste0(input,"trapntox_sites.raw.rds")) #trap/tox site locs
fp.chulls=readRDS(paste0(input,"fp.chulls.rds")) #flight path convex hulls
activities.raw=readRDS(paste0(input,"activities.rds"))

#read in geo
home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/"
geo=read.csv(paste0(home,"geo_remtyp.csv"))
geo=geo[,-1]
geo=geo[!is.na(geo$Removal.Type),]
geo$datetime<-Neat.Dates.POSIXct(geo$datetime,tz="UTC")
geo$date_only<-Neat.Dates.POSIXct(geo$date_only,tz="UTC")

tox.act=activities.raw[activities.raw$activity=="toxic",]
tox.act$tox_datetime=Neat.Dates.POSIXct(tox.act$tox_datetime,tz="UTC")
tox.start.date=as.Date(min(tox.act$tox_datetime))
tox.end.date=as.Date("2023-03-09")

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

geo.tox=Do.Date.Cutoffs(geo,"tox",tox.cutoff1,tox.cutoff2,tox.start.date,tox.end.date)

#geo[geo$type=="tox",]
#geot=geo.tox[geo.tox$Removal.Type=="tox",]

IDsf20=unique(geo.tox[geo.tox$Removal.Type=="tox",]$animalid)

#IDsf20=unique(id)
#for i in 1:length(IDsf20)
Convert.Telemetry.fromtk<-function(geolocs){
  geolocs<-as.data.frame(geolocs)
  id.n=which(colnames(geolocs)=="animalid")
  id.t=which(colnames(geolocs)=="t_")
  id.lon=which(colnames(geolocs)=="latitude")
  id.lat=which(colnames(geolocs)=="longitude")
  tk=geolocs[,c(id.n,id.t,id.lon,id.lat)]
  colnames(tk)<-c("ID","timestamp","longitude","latitude")
  crs_str="+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  te1 <- ctmm::as.telemetry(tk,projection=crs_str)
  return(te1)
}

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


#for(i in 1:2){
for(i in 1:length(IDsf20)){
  print(paste0("starting pig ",IDsf20[i],": ",i," of ",length(IDsf20)))
  pig=geo.tox[geo.tox$animalid==IDsf20[i],]
  pigtk=mk_track(pig,.x=X,.y=Y,.t=datetime,all_cols=TRUE) #
  pigtk_20=track_resample(pigtk,rate=minutes(20),tolerance=minutes(2))
  pigtk_240=track_resample(pigtk,rate=minutes(240),tolerance=minutes(24))
  
  #Get akde hr area of pigtk full, 20, 240

  tk.list=list(pigtk,pigtk_20,pigtk_240)
  tk.names=c("full","20min","240min")
  #do akde stuff
  for(t in 1:3){
  print(paste0("starting akde ests for subsample ",tk.names[t]))
  #convert to telemetry format
  te1=Convert.Telemetry.fromtk(tk.list[[t]])
  GUESS1 <- ctmm.guess(te1, interactive = FALSE)
  print("fitting ctmm model for full path")
  FIT1_pHREML <- ctmm.select(te1, GUESS1, method = 'pHREML')
  Neff.t=summary(FIT1_pHREML)$DOF["area"]
  UD1_pHREML <- akde(te1, FIT1_pHREML)
  area.est.full.t=summary(UD1_pHREML)$CI

  dfa.t=as.data.frame(t(c(IDsf20[i],tk.names[t],area.est.full.t,Neff.t)))
  colnames(dfa.t)=c("animalid","sstype","low","est","high","Neff")
  
  #Do NSD calcs
  pigex=tk.list[[t]]
  pigex=pigex[order(pigex$t_),]
  #lagged loc
  pigex$Xstart=pigex$x_[1]
  pigex$Ystart=pigex$y_[1]
  #pigex$datetime2=dplyr::lead(pigex$datetime)
  
  pigex$displacement=sqrt((pigex$x_-pigex$Xstart)^2+(pigex$y_-pigex$Ystart)^2)
  pigex$NSD=(pigex$displacement)^2
  pwNSD.t=pigex %>% group_by(week) %>% dplyr::summarise(ndays=n_distinct(jDate),mNSD=median(NSD))

  #keep only full weeks
  pwNSD.t=pwNSD.t[pwNSD.t$ndays==7,]
  
  #remove days col
  pwNSD.t=pwNSD.t[,c(1,3)]
  
  pwNSD.t$animalid=IDsf20[i]
  pwNSD.t$sstype=tk.names[t]
  
  if(t==1&i==1){
    dfa=dfa.t
    pwNSD=pwNSD.t
  } else{
    dfa=rbind(dfa,dfa.t)
    pwNSD=rbind(pwNSD,pwNSD.t)
  }
  
  }
  
  }

saveRDS(dfa,"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Output/SA_Results/SA_areas.rds")
saveRDS(pwNSD,"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Output/SA_Results/SA_disp.rds")

ggplot(dfa,aes(y=animalid,x=as.numeric(est),color=sstype))+
  geom_point()+
  geom_segment(aes(x=as.numeric(low),xend=as.numeric(high)))

#remove 86101_j3_j3 and 86075_y5_y5 are atypical (still not issue with respect to diff fix rates)
dfa2=dfa[dfa$animalid!="86101_J3_J3"&dfa$animalid!="86075_Y5_Y5",]
#change levels to make full plot on top
dfa2$sstype<-forcats::fct_relevel(dfa2$sstype,c("full","20min","240min"))

ggplot(dfa2,aes(y=animalid,x=as.numeric(est),color=sstype))+
  geom_point()+
  geom_segment(aes(x=as.numeric(low),xend=as.numeric(high)))

#do sep plots for 240/20min-- too crowded to viz
ggplot(dfa2[dfa2$sstype!="20min",],aes(y=animalid,x=as.numeric(est),color=as.factor(sstype)))+
  geom_point(size=5,alpha=0.7)+
  #scale_alpha_discrete(values=c(0.2,1))+
  geom_segment(aes(x=as.numeric(low),xend=as.numeric(high)),linewidth=1,alpha=0.7)+
  xlab("Area est. (m)")
filename="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Output/SA_Results/"
ggsave(paste0(filename,"fix240SA_areas.png"),bg="white")


#Confirmed that area calcs within 240mins were within confidence interval of original
ggplot(dfa2[dfa2$sstype!="240min",],aes(y=animalid,x=as.numeric(est),color=as.factor(sstype)))+
  geom_point(size=5,alpha=0.7)+
  geom_segment(aes(x=as.numeric(low),xend=as.numeric(high)),
               linewidth=1,alpha=0.7)+
  xlab("Area est (m)")
filename="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Output/SA_Results/"
ggsave(paste0(filename,"fix20SA_areas.png"),bg="white")


dfav20=dfa2 %>% 
  group_by(animalid) %>% 
  dplyr::summarise("full_est"=first(est[sstype=="full"]),
                   "full_high"=first(high[sstype=="full"]),
                   "full_low"=first(low[sstype=="full"]),
                   "min20_est"=first(est[sstype=="20min"]),
                   "min240_est"=first(est[sstype=="240min"]))
dfav20$full_est<-as.numeric(dfav20$full_est)
dfav20$full_low<-as.numeric(dfav20$full_low)
dfav20$full_high<-as.numeric(dfav20$full_high)
dfav20$min20_est<-as.numeric(dfav20$min20_est)
dfav20$min240_est<-as.numeric(dfav20$min240_est)

dfav20$bool_20min=(dfav20$full_low<=dfav20$min20_est)&(dfav20$min20_est<=dfav20$full_high)
dfav20$bool_240min=(dfav20$full_low<=dfav20$min240_est)&(dfav20$min240_est<=dfav20$full_high)

dfav20$diff_20=dfav20$full_est-dfav20$min20_est
dfav20$diff_240=dfav20$full_est-dfav20$min240_est

dfav20$prop_20=dfav20$min20_est/dfav20$full_est
dfav20$prop_240=dfav20$min240_est/dfav20$full_est

#####Now for displacement

pwNSD2=pwNSD[pwNSD$animalid!="86101_J3_J3",]
ggplot(pwNSD2,aes(x=week,y=mNSD,color=sstype))+
  geom_point(position=position_jitter())+
  facet_wrap(~animalid)+
  theme_ipsum()
filename="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Output/SA_Results/"
ggsave(paste0(filename,"fix_SA_disps.png"),bg="white",width=15,height=15,units="in")

View(pwNSD)

pwNSDs=pwNSD %>% 
  group_by(animalid,week) %>% 
  dplyr::summarise("mNSD"=first(mNSD[sstype=="full"]),
                   "mNSD20"=first(mNSD[sstype=="20min"]),
                   "mNSD240"=first(mNSD[sstype=="240min"]))


pwNSDw=pwNSD %>% tidyr::pivot_wider(names_from=sstype,values_from=mNSD)
View(pwNSDw)
colnames(pwNSDw)[4:5]<-c("min20","min240")

pwNSDw$diff20=sqrt(pwNSDw$full)-sqrt(pwNSDw$min20)
pwNSDw$diff240=sqrt(pwNSDw$full)-sqrt(pwNSDw$min240)

ggplot(pwNSDw,aes(x=diff20))+
  geom_histogram()+
  theme_ipsum()+
  xlab("Displacement diff (m), 20 min subsample")
ggsave(paste0(filename,"fix_SA_hist_20.png"),bg="white")


ggplot(pwNSDw,aes(x=diff240))+
  geom_histogram()+
  theme_ipsum()+
  xlab("Displacement diff (m), 240 min subsample")
ggsave(paste0(filename,"fix_SA_hist_240.png"),bg="white")

max(pwNSDw$diff240)
min(pwNSDw$diff240)

max(pwNSDw$diff20)
min(pwNSDw$diff20)

#Get akde hr areas of pigtk full, 20, and 240
#place output in df

#Get displacement vals by week for each
#place output in sep df


#################################################################################
########Plan leave one out analysis
#Area first

library(glmmTMB)
library(DHARMa)
library(ctmm)
library(fitdistrplus)
library(ggplot2)
library(plyr)
library(dplyr)
library(hrbrthemes)

aktox=readRDS("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/outdf_akde_tox_corrected.rds")

#Remove NA after periods for tox-killed pigs
aktox=aktox[!is.na(aktox$area.CI.low),]

#Relevel period
aktox$period<-forcats::fct_relevel(aktox$period,c("before","after"))

#Relevel removal type
aktox$Removal.Type<-forcats::fct_relevel(aktox$Removal.Type,c("ctrl","tox"))

#first, does it change which is top model
#after assess that, see what top model parms ests look like

###for each animalid in akaer, remove one at a time and rerun models
tox.IDs=unique(aktox$animalid)

res.rps=glmmTMB(area.est ~ Removal.Type*period+sex + (1|animalid), data=aktox,family=Gamma(link=log))
parms=broom.mixed::tidy(res.rps)
parms$minusID="full"

for(i in 1:length(tox.IDs)){
print(paste0(tox.IDs[i],": ",i," of ",length(tox.IDs)))
  aktox.i=aktox[aktox$animalid!=tox.IDs[i],]
#####Run model
  res.rps=glmmTMB(area.est ~ Removal.Type*period+sex + (1|animalid), data=aktox.i,family=Gamma(link=log))
  
#combine with others

parms.i=broom.mixed::tidy(res.rps)
parms.i$minusID=tox.IDs[i]

  parms=rbind(parms,parms.i)
  
}

#Plot showing whether interpretation changes-- 
intxn=parms[parms$term=="Removal.Typetox:periodafter",]
ggplot(intxn,aes(y=minusID,x=estimate))+
  geom_point()+
  geom_segment(aes(x=(estimate-2*std.error),xend=(estimate+2*std.error)))+
  geom_vline(xintercept=0,linetype="dashed",color="red")+
  theme_ipsum()+
  xlab("Parameter estimate, Removal type: tox * Period: after")
ggsave("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Output/SA_Results/L1O_Area_Tox.png",width=5,height=8,units="in")


#################################################################################
#####Now displacement SA
source("~/Library/CloudStorage/OneDrive-USDA/Projects/Automation/NeatDates.R")
library(stringr)
library(NIFApackagev1.1)
#Setup, copied from Disp_Model_wDiagnostics.R
{
  home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/"
  geo=read.csv(paste0(home,"geo_remtyp.csv"))
  geo=geo[,-1]
  geo=geo[!is.na(geo$Removal.Type),]
  geo$date_only<-Neat.Dates.POSIXct(geo$date_only,tz="UTC")
  
  #Load df outputs from Get_Period_AKDEs
  akaer=readRDS("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/outdf_akde_aer_corrected.rds")
  aktrap=readRDS("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/outdf_akde_trap_corrected.rds")
  aktox=readRDS("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/outdf_akde_tox_corrected.rds")
  
  #Fix area units-- sometimes area rec'd as hectares, meters
  Fix.Area.Units<-function(akrem){
    akrem[grep("square meters",akrem$akde.hrarea.units),][3:5]=akrem[grep("square meters",akrem$akde.hrarea.units),][3:5]/1e6
    akrem[grep("hectares",akrem$akde.hrarea.units),][3:5]=akrem[grep("hectares",akrem$akde.hrarea.units),][3:5]/100
    return(akrem)
  }
  
  aktox=Fix.Area.Units(aktox)
  aktrap=Fix.Area.Units(aktrap)
  akaer=Fix.Area.Units(akaer)
  
  
  #Formatting
  akaer$period<-forcats::fct_relevel(akaer$period,c("before","after"))
  akaer$Removal.Type<-forcats::fct_relevel(akaer$Removal.Type,c("ctrl","aer"))
  
  #Remove pigs that were killed during aerial
  killed.aer=akaer[is.na(akaer$area.CI.low),]$animalid
  akaer=akaer[!(akaer$animalid%in%killed.aer),]
  
  #Remove NA after periods for tox-killed pigs
  aktox=aktox[!is.na(aktox$area.CI.low),]
  
  #Relevel period
  akaer$period<-forcats::fct_relevel(akaer$period,c("before","after"))
  aktrap$period<-forcats::fct_relevel(aktrap$period,c("before","during","after"))
  aktox$period<-forcats::fct_relevel(aktox$period,c("before","after"))
  
  #Relevel removal type
  akaer$Removal.Type<-forcats::fct_relevel(akaer$Removal.Type,c("ctrl","aer"))
  aktrap$Removal.Type<-forcats::fct_relevel(aktrap$Removal.Type,c("ctrl","trap"))
  aktox$Removal.Type<-forcats::fct_relevel(aktox$Removal.Type,c("ctrl","tox"))
  
  aktox.sumperiods=aktox %>% group_by(animalid,period) %>% dplyr::summarise(n()) %>% tidyr::pivot_wider(names_from="period",values_from=`n()`) %>% as.data.frame()
  died.tox=aktox.sumperiods[which(is.na(aktox.sumperiods$after)),1]
  aktox$died.tox=0
  aktox[aktox$animalid%in%died.tox,]$died.tox=1
  
  #make sure died tox is factor
  aktox$died.tox<-as.factor(aktox$died.tox)
  aktox$died.tox<-forcats::fct_relevel(aktox$died.tox, c("0","1"))
  aktox2=aktox[aktox$died.tox==0,]
  
  #View(akaer)
  {
    #trap dates
    trap.act=activities.raw[activities.raw$activity=="trap",]
    trap.act$trap_datetime=Neat.Dates.POSIXct(trap.act$trap_datetime,tz="UTC")
    trap.start.date=as.Date(min(trap.act$trap_datetime))
    trap.end.date=as.Date(max(trap.act$trap_datetime))
    
    #tox dates
    tox.act=activities.raw[activities.raw$activity=="toxic",]
    tox.act$tox_datetime=Neat.Dates.POSIXct(tox.act$tox_datetime,tz="UTC")
    tox.start.date=as.Date(min(tox.act$tox_datetime))
    tox.end.date=as.Date("2023-03-09")
    
    #aerial dates
    aer.start.date=as.Date("2023-03-27")
    aer.end.date=as.Date("2023-03-29")
    
  }
  
  #Get cutoffs for each removal type
  #Aerial 
  #determined 4 weeks from sensitivity analysis
  #geo.aer=geo2[geo2$Removal.Type=="aer"|geo2$Removal.Type=="ctrl",]
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
  
  #removal="tox"
  #cutoff1=tox.cutoff1
  #cutoff2=tox.cutoff2
  #startdate=tox.start.date
  #enddate=tox.end.date
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
  #geo.trap=Do.Date.Cutoffs(geo,"trap",trap.cutoff1,trap.cutoff2,trap.start.date.akde,trap.end.date)
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
  
  
  
}

################## Calc NSD and summarise by week
{
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
  
  #ggplot(geo.aerd.wk,aes(y=mNSD,x=week,group=animalid,color=removal.period.akdecalc))+
  #  geom_line()+facet_wrap(~Removal.Type)
  #ggplot(geo.trapd.wk,aes(y=mNSD,x=week,group=animalid,color=removal.period.akdecalc))+
  #  geom_line()+facet_wrap(~Removal.Type)
  #ggplot(geo.toxd.wk,aes(y=mNSD,x=week,group=animalid,color=removal.period.akdecalc))+
  #  geom_line()+facet_wrap(~Removal.Type)
  
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
  #origin.vector=c(trap.origin,trap.origin.after,trap.origin.during)
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
  
  #test=Do.Week.Split(geo.aerd,"aer",c(aer.origin,aer.origin.after))
  #test=test[test$week!=11,]
  #test==geo.aerd2
  
  geo.aerd2=Do.Week.Split(geo.aerd,"aer",c(aer.origin,aer.origin.after))
  geo.trapd2=Do.Week.Split(geo.trapd,"trap",c(trap.origin,trap.origin.after,trap.origin.during))
  geo.toxd2=Do.Week.Split(geo.toxd,"tox",c(tox.origin,tox.origin.after,tox.origin.during))
  
  
  #verify that no weeks are very short
  aerchk=geo.aerd2 %>% 
    group_by(animalid,removal.period.akdecalc,week) %>% 
    dplyr::summarise(n_distinct(jDate)) %>%
    as.data.frame()
  #View(aerchk)
  #remove week 11, only two days
  geo.aerd2=geo.aerd2[geo.aerd2$week!=11,]
  
  trapchk=geo.trapd2 %>% 
    group_by(animalid,removal.period.akdecalc,week) %>% 
    dplyr::summarise(n_distinct(jDate)) %>%
    as.data.frame()
  #View(trapchk)
  #remove week 25, only 1 days
  geo.trapd2=geo.trapd2[geo.trapd2$week!=25,]
  
  toxchk=geo.toxd2 %>% 
    group_by(animalid,removal.period.akdecalc,week) %>% 
    dplyr::summarise(n_distinct(jDate)) %>%
    as.data.frame()
  #View(toxchk)
  
  #remove week 11, only 2 days
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
  
  #visualize to check
  #ggplot(geo.aerd.wk,aes(y=mNSD,x=week,group=animalid,color=removal.period.akdecalc))+
  #  geom_line()+facet_wrap(~Removal.Type)+theme_ipsum()
  #ggplot(geo.trapd.wk,aes(y=mNSD,x=week,group=animalid,color=removal.period.akdecalc))+
  #  geom_line()+facet_wrap(~Removal.Type)+theme_ipsum()
  #ggplot(geo.toxd.wk,aes(y=mNSD,x=week,group=animalid,color=removal.period.akdecalc))+
  #  geom_line()+facet_wrap(~Removal.Type)+theme_ipsum()
}

#tox- res.rpdto
#Separate data by who died during tox treatment
geo.toxd.wk$died_tox=NA
geo.toxd.wk[geo.toxd.wk$animalid%in%died.tox,]$died_tox=1
geo.toxd.wk[is.na(geo.toxd.wk$died_tox),]$died_tox<-0
geo.toxd.wk$died_tox<-as.factor(geo.toxd.wk$died_tox)

geo.wkdf=geo.toxd.wk
geo.wkdf$Removal.Type<-forcats::fct_relevel(geo.wkdf$Removal.Type,c("ctrl","tox"))
geo.wkdf$removal.period.akdecalc<-forcats::fct_relevel(geo.wkdf$removal.period.akdecalc,c("before","during","after"))
geo.wkdf$sex<-forcats::fct_relevel(geo.wkdf$sex,c("Female","Male"))

res.rp.dto=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + Removal.Type*removal.period.akdecalc + died_tox, data=geo.wkdf,family=Gamma(link=log))
parmsxd=broom.mixed::tidy(res.rp.dto)
parmsxd$minusID="full"
tox.IDs=unique(geo.toxd.wk$animalid)

for(i in 1:length(tox.IDs)){
  print(paste0(tox.IDs[i],": ",i," of ",length(tox.IDs)))
  geo.wkdf.i=geo.wkdf[geo.wkdf$animalid!=tox.IDs[i],]
  #####Run model
  res.rp.dto=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + Removal.Type*removal.period.akdecalc + died_tox, data=geo.wkdf.i,family=Gamma(link=log))
  
  #combine with others
  parmsxd.i=broom.mixed::tidy(res.rp.dto)
  parmsxd.i$minusID=tox.IDs[i]
  
  parmsxd=rbind(parmsxd,parmsxd.i)
  
}

intxn_xd_d=parmsxd[parmsxd$term=="Removal.Typetox:removal.period.akdecalcduring",]
intxn_xd_a=parmsxd[parmsxd$term=="Removal.Typetox:removal.period.akdecalcafter",]

ggplot(intxn_xd_d,aes(y=minusID,x=estimate))+
  geom_point()+
  geom_segment(aes(x=(estimate-2*std.error),xend=(estimate+2*std.error)))+
  geom_vline(xintercept=0,linetype="dashed",color="red")+
  theme_ipsum()+
  xlab("Parameter estimate, Removal type: tox * Period: during")
ggsave("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Output/SA_Results/L1O_Disp_Tox_1.png",width=5,height=8,units="in")

ggplot(intxn_xd_a,aes(y=minusID,x=estimate))+
  geom_point()+
  geom_segment(aes(x=(estimate-2*std.error),xend=(estimate+2*std.error)))+
  geom_vline(xintercept=0,linetype="dashed",color="red")+
  theme_ipsum()+
  xlab("Parameter estimate, Removal type: tox * Period: after")
ggsave("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Output/SA_Results/L1O_Disp_Tox2.png",width=5,height=8,units="in")

#trap- res.rps###################################
geo.wkdf=geo.trapd.wk
geo.wkdf$Removal.Type<-forcats::fct_relevel(geo.wkdf$Removal.Type,c("ctrl","trap"))
geo.wkdf$removal.period.akdecalc<-forcats::fct_relevel(geo.wkdf$removal.period.akdecalc,c("before","during","after"))
geo.wkdf$sex<-forcats::fct_relevel(geo.wkdf$sex,c("Female","Male"))

res.rps2=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + Removal.Type*removal.period.akdecalc*sex, data=geo.wkdf,family=Gamma(link=log))
parmstd=broom.mixed::tidy(res.rps2)
parmstd$minusID="full"
trap.IDs=unique(geo.trapd.wk$animalid)

for(i in 1:length(trap.IDs)){
  print(paste0(trap.IDs[i],": ",i," of ",length(trap.IDs)))
  geo.wkdf.i=geo.wkdf[geo.wkdf$animalid!=trap.IDs[i],]
  #####Run model
  res.rps2=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + Removal.Type*removal.period.akdecalc*sex, data=geo.wkdf.i,family=Gamma(link=log))
  
  #combine with others
  parmstd.i=broom.mixed::tidy(res.rps2)
  parmstd.i$minusID=trap.IDs[i]
  
  parmstd=rbind(parmstd,parmstd.i)
  
}


intxn_td_d=parmstd[parmstd$term=="Removal.Typetrap:removal.period.akdecalcduring",]
intxn_td_a=parmstd[parmstd$term=="Removal.Typetrap:removal.period.akdecalcafter",]
intxn_td_ds=parmstd[parmstd$term=="Removal.Typetrap:removal.period.akdecalcduring:sexMale",]
intxn_td_as=parmstd[parmstd$term=="Removal.Typetrap:removal.period.akdecalcafter:sexMale",]

ggplot(intxn_td_d,aes(y=minusID,x=estimate))+
  geom_point()+
  geom_segment(aes(x=(estimate-2*std.error),xend=(estimate+2*std.error)))+
  geom_vline(xintercept=0,linetype="dashed",color="red")+
  theme_ipsum()+
  xlab("Parameter estimate, Removal type: trap * Period: during")
ggsave("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Output/SA_Results/L1O_Disp_Trap1.png",width=5,height=8,units="in")

ggplot(intxn_td_a,aes(y=minusID,x=estimate))+
  geom_point()+
  geom_segment(aes(x=(estimate-2*std.error),xend=(estimate+2*std.error)))+
  geom_vline(xintercept=0,linetype="dashed",color="red")+
  theme_ipsum()+
  xlab("Parameter estimate, Removal type: trap * Period: after")
ggsave("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Output/SA_Results/L1O_Disp_Trap2.png",width=5,height=8,units="in")


ggplot(intxn_td_ds,aes(y=minusID,x=estimate))+
  geom_point()+
  geom_segment(aes(x=(estimate-2*std.error),xend=(estimate+2*std.error)))+
  geom_vline(xintercept=0,linetype="dashed",color="red")+
  theme_ipsum()+
  xlab("Parameter estimate, Removal type: trap * Period: during * Sex: male")
ggsave("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Output/SA_Results/L1O_Disp_Trap3.png",width=5,height=8,units="in")

ggplot(intxn_td_as,aes(y=minusID,x=estimate))+
  geom_point()+
  geom_segment(aes(x=(estimate-2*std.error),xend=(estimate+2*std.error)))+
  geom_vline(xintercept=0,linetype="dashed",color="red")+
  theme_ipsum()+
  xlab("Parameter estimate, Removal type: trap * Period: after * Sex: male")
ggsave("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Output/SA_Results/L1O_Disp_Trap4.png",width=5,height=8,units="in")




