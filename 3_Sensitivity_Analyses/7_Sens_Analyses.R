
# Purpose --------------------------
#The purpose of this script is to run several sensitivity analyses for analyses on effects of removal methods on movement/contact responses

# Overview -------------------------
#1. subsampling analysis
#2. toxicant fate analysis
#3. leave one out analysis

# Setup -------------------------
home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline"

#Set dirs
input=file.path(home,"1_Data","Input",fsep=.Platform$file.sep)
objdir=file.path(home,"1_Data","Objects",fsep=.Platform$file.sep)
outdir=file.path(home,"3_Output",fsep=.Platform$file.sep)

#load libraries
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
library(DHARMa)
library(fitdistrplus)
library(stringr)

#read objects from input
input="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/1_Data/Input/"
sitelocs=readRDS(paste0(input,"trapntox_sites.raw.rds")) #trap/tox site locs
fp.chulls=readRDS(paste0(input,"fp.chulls.rds")) #flight path convex hulls
activities.raw=readRDS(paste0(input,"activities.rds"))

## Read in data -----------
#for subsampling analysis:
home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/"
geo=read.csv(paste0(home,"geo_remtyp.csv"))
geo=geo[,-1]
geo=geo[!is.na(geo$Removal.Type),]
geo$datetime<-Neat.Dates.POSIXct(geo$datetime,tz="UTC")
geo$date_only<-Neat.Dates.POSIXct(geo$date_only,tz="UTC")

#For leave one out analysis:
#load NSD tox data
geo.toxd.wk=readRDS(file.path(objdir,"NSDgeotox.rds",fsep=.Platform$file.sep))
#load distance data
dist<- readRDS(paste0(objdir,"/pig_weekly_distance_ctmm.rds"))

## Load relevant dates ---------------------
tox.act=activities.raw[activities.raw$activity=="toxic",]
tox.act$tox_datetime=Neat.Dates.POSIXct(tox.act$tox_datetime,tz="UTC")
tox.start.date=as.Date(min(tox.act$tox_datetime))
tox.end.date=as.Date("2023-03-09")

tox.len=36
tox.cutoff1=tox.start.date-tox.len-1
tox.cutoff2=tox.end.date+tox.len+1

## Set needed functions ---------------------
#source needed functions
funcdir<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/2_Scripts/Functions"
func.list=list.files(funcdir,full.names=TRUE)
for(f in 1:length(func.list)){
  source(func.list[f])
}

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

source("~/Library/CloudStorage/OneDrive-USDA/Projects/Automation/NeatDates.R")

#colors for plots
tox.hex="#FF985A"
trap.hex="#FF57A9"
aer.hex="#822AFF"
ctrl.hex="#ffc400"

# Subsampling analysis -------------------------

#Steps:
#1. get sample of 10 pigs with 20 min fix rates
#2. resample with amt- 20 mins and 4 hrs
#3. run akde and NSD calcs
#4. output df with full and shortened versions

geo.tox=Do.Date.Cutoffs(geo,"tox",tox.cutoff1,tox.cutoff2,tox.start.date,tox.end.date)

IDsf20=unique(geo.tox[geo.tox$Removal.Type=="tox",]$animalid)

## AKDE Areas -------------------------------------------

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

#remove 86101_j3_j3 and 86075_y5_y5 are atypical
dfa2=dfa[dfa$animalid!="86101_J3_J3"&dfa$animalid!="86075_Y5_Y5",]

#change levels to make full plot on top
dfa2$sstype<-forcats::fct_relevel(dfa2$sstype,c("full","20min","240min"))

ggplot(dfa2,aes(y=animalid,x=as.numeric(est),color=sstype))+
  geom_point()+
  geom_segment(aes(x=as.numeric(low),xend=as.numeric(high)))

#do sep plots for 240/20min-- too crowded to viz
ggplot(dfa2[dfa2$sstype!="20min",],aes(y=animalid,x=as.numeric(est),color=as.factor(sstype)))+
  geom_point(size=5,alpha=0.7)+
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

## NSD -------------------------------------------

pwNSD2=pwNSD[pwNSD$animalid!="86101_J3_J3",]
ggplot(pwNSD2,aes(x=week,y=mNSD,color=sstype))+
  geom_point(position=position_jitter())+
  facet_wrap(~animalid)+
  theme_ipsum()
filename="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Output/SA_Results/"
ggsave(paste0(filename,"fix_SA_disps.png"),bg="white",width=15,height=15,units="in")

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

# Leave one out analysis -------------------------------------------------------

#Doing only for interaction effect sizes >20 where significant in overall model, but not when testing male/female interaction
#Justification: 
  #is likely that loss of power results in loss of significance in m/f
  #suggests that could be small number of individuals causing the sig effect, because is otherwise very strong
  #informative to know if effect caused by high inter-individual variability, or broader population-level effect but low sample size of one sex

#Effects of interest:
  #Tox, NSD, after
  #Trap, distance, after

## Make function to do the leave one out thing ---------------------------------
#inputs:
  #dataset, subset to period, response and removal type of interest
  #model call

#data=geo.toxd.wk
#call=res.spat$call
#NSD_l1o_parms=LeaveOneOut(geo.toxd.wk,call)
LeaveOneOut<-function(data,call){
  
  #get vector of IDs
  IDs=unique(data$animalid)
  
  #evaluate model with full dataset
  print("running full model")
  res=eval(call)
  
  #tidy output
  parms=broom.mixed::tidy(res)
  
  #add minusID col
  parms$minusID="full"
  
  #rename data to reuse call for removals
  dat=data
  
  #loop through each id, remove from data, rerun model, output results
  for(i in 1:length(IDs)){
    
    print(paste0(IDs[i],": ",i," of ",length(IDs)))
    
    data=dat[dat$animalid!=IDs[i],]
    res=eval(call)
    
    parms.i=broom.mixed::tidy(res)
    parms.i$minusID=IDs[i]
    
    parms=rbind(parms,parms.i)
    
  }
  
  return(parms)
  
}


## Tox, NSD, after ---------------------------------
#geo.toxd.wk=readRDS(file.path(objdir,"NSDgeotox.rds",fsep=.Platform$file.sep))

#Formatting needed to run model
colnames(geo.toxd.wk)[c(2,3)]<-c("trt_ctrl","period")
geo.toxd.wk$trt_ctrl<-as.character(geo.toxd.wk$trt_ctrl)
geo.toxd.wk$trt_ctrl[geo.toxd.wk$trt_ctrl!="ctrl"]<-"trt"
geo.toxd.wk$trt_ctrl<-as.factor(geo.toxd.wk$trt_ctrl)
geo.toxd.wk$trt_ctrl<-forcats::fct_relevel(geo.toxd.wk$trt_ctrl,"ctrl","trt")
geo.toxd.wk<-geo.toxd.wk[geo.toxd.wk$period!="during",]

#Run model to get call for function
data=geo.toxd.wk
res=glmmTMB(mNSD ~ trt_ctrl*period, data=data,family=Gamma(link=log))
call=res$call

#Run function for leave one out process
NSD_l1o_parms=LeaveOneOut(data,call)

#subset to interaction effect of interest
intxn=NSD_l1o_parms[NSD_l1o_parms$term=="trt_ctrltrt:periodafter",]

#join minusID to info about pigs
  #sex, trt vs control
NSDids=unique(geo.toxd.wk[,c("animalid","trt_ctrl","sex")])
colnames(NSDids)[1]<-"minusID"
intxn=left_join(intxn,NSDids,by="minusID")
intxn$trt_ctrl<-as.character(intxn$trt_ctrl)
intxn[intxn$minusID=="full",]$sex="Ref"
intxn[intxn$minusID=="full",]$trt_ctrl="Ref"

#Relevel for plot
intxn$sex<-factor(intxn$sex,levels=c("Ref","Female","Male"))
intxn$trt_ctrl<-factor(intxn$trt_ctrl,levels=c("Ref","ctrl","trt"))

#Plot showing whether interpretation changes-- 
ggplot(intxn,aes(y=minusID,x=estimate))+
  geom_point(aes(color=trt_ctrl,
                 shape=sex),
             size=3)+
  geom_segment(aes(x=(estimate-2*std.error),
                   xend=(estimate+2*std.error),
                   color=trt_ctrl))+
  scale_color_manual(values=c("black","#d6a400",tox.hex))+
  geom_vline(xintercept=0,linetype="dashed",color="red")+
  theme_ipsum()+
  theme(axis.text.y=element_blank())+
  ylab("Individual removed")
  xlab("Parameter estimate, Removal type: tox * Period: after")
#ggsave("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Output/SA_Results/L1O_Area_Tox.png",width=5,height=8,units="in")

## Trap, distance, after ---------------------------------
  dist<- readRDS(paste0(objdir,"/pig_weekly_distance_ctmm.rds"))
  
  #Formatting needed to run model
  colnames(dist)[4]<-"period"
  dist<-dist[dist$period!="during",]
  dist<-dist[dist$Removal.Type=="trap",]
  
  #relevel factors
  dist$period<-factor(dist$period,levels=c("before","after"))
  dist$trt_ctrl<-factor(dist$trt_ctrl,levels=c("ctrl","trt"))
  
  #Run model to get call for function
  data=dist
  res=glmmTMB(weekly_dist_km~(1|animalid)+trt_ctrl*period,
                       data=data,family=Gamma(link="log"))
  call=res$call
  
  #Run function for leave one out process
  dist_l1o_parms=LeaveOneOut(data,call)
  
  #subset to interaction effect of interest
  intxn=dist_l1o_parms[dist_l1o_parms$term=="trt_ctrltrt:periodafter",]
  
  #join minusID to info about pigs
  #sex, trt vs control
  distids=unique(dist[,c("animalid","trt_ctrl","sex")])
  colnames(distids)[1]<-"minusID"
  intxn=left_join(intxn,distids,by="minusID")
  intxn$sex<-as.character(intxn$sex)
  intxn$trt_ctrl<-as.character(intxn$trt_ctrl)
  intxn[intxn$minusID=="full",]$sex="Ref"
  intxn[intxn$minusID=="full",]$trt_ctrl="Ref"
  
  #Relevel for plot
  intxn$sex<-factor(intxn$sex,levels=c("Ref","Female","Male"))
  intxn$trt_ctrl<-factor(intxn$trt_ctrl,levels=c("Ref","ctrl","trt"))
  
  #Plot showing whether interpretation changes-- 
  intxn=intxn[order(intxn$trt_ctrl),]
  ggplot(intxn,aes(y=minusID,x=estimate))+
    geom_point(aes(color=trt_ctrl,
                   shape=sex),
               size=3)+
    geom_segment(aes(x=(estimate-2*std.error),
                     xend=(estimate+2*std.error),
                     color=trt_ctrl))+
    scale_color_manual(values=c("black","#d6a400",trap.hex))+
    geom_vline(xintercept=0,linetype="dashed",color="red")+
    theme_ipsum()+
    theme(axis.text.y=element_blank())+
    ylab("Individual removed")
  xlab("Parameter estimate, Removal type: tox * Period: after")
  #ggsave("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Output/SA_Results/L1O_Area_Tox.png",width=5,height=8,units="in")
  


  
  
