#set home dir of pipeline
home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline"

# Purpose ----------------------------------------------------------------------

#The purpose of this script is to make figures showing change in each response variable over time

# In/Out -----------------------------------------------------------------------

#inputs:
#outdf_akde_aer_corrected_f.rds, 
#outdf_akde_tox_corrected_f.rds, 
#outdf_akde_trap_corrected_f.rds

#outputs: 
#.png figures in X folder

# Setup ------------------------------------------------------------------------

#load libraries
library(glmmTMB)
library(DHARMa)
library(ctmm)
library(fitdistrplus)
library(ggplot2)
library(plyr)
library(dplyr)
library(hrbrthemes)
library(sf)
library(mapview)
library(cowplot)

#Set dirs
input=file.path(home,"1_Data","Input",fsep=.Platform$file.sep)
objdir=file.path(home,"1_Data","Objects",fsep=.Platform$file.sep)
if(!dir.exists(file.path(home,"3_Output","Area_GLM_Results",fsep=.Platform$file.sep))){
  dir.create(file.path(home,"3_Output","Area_GLM_Results",fsep=.Platform$file.sep))
}
outdir=file.path(home,"3_Output","Area_GLM_Results",fsep=.Platform$file.sep)

#Read in area data
akaer=readRDS(file.path(objdir,"outdf_akde_aer_corrected_f.rds",fsep=.Platform$file.sep))
aktrap=readRDS(file.path(objdir,"outdf_akde_trap_corrected_f.rds",fsep=.Platform$file.sep))
aktox=readRDS(file.path(objdir,"outdf_akde_aktox_corrected_f.rds",fsep=.Platform$file.sep))

#Read in NSD data:
geo.aerd.wk=readRDS(file.path(objdir,"NSDgeoaer.rds"))
geo.trapd.wk=readRDS(file.path(objdir,"NSDgeotrap.rds"))
geo.toxd.wk=readRDS(file.path(objdir,"NSDgeotox.rds"))

#get IDs of satellite pigs-- these will have produced bad CTMMS
geo=readRDS(file.path(input,"geo.rds",fsep=.Platform$file.sep)) #geolocations, tidied and filtered
satpigs=unique(geo[geo$data_from=="satellite",]$animalid)

#Read in distance data:
#distaer=readRDS(file.path(objdir,"distaer.rds"))
#disttrap=readRDS(file.path(objdir,"distrap.rds"))
#disttox=readRDS(file.path(objdir,"disttox.rds"))
dist<- readRDS(paste0(objdir,"/pig_weekly_distance_ctmm.rds"))

#subset dfs
distaer <- dist %>% filter(Removal.Type%in%c('aer')) %>% filter(removal.period.akdecalc!='during') 
disttox <- dist %>% filter(Removal.Type%in%c('tox')) 
disttrap <- dist %>% filter(Removal.Type%in%c('trap')) 

#relevel when during isn't recorded
distaer$removal.period.akdecalc <- factor(distaer$removal.period.akdecalc,
                                          levels=c('before','after'))
disttox$removal.period.akdecalc <- factor(disttox$removal.period.akdecalc,
                                          levels=c('before','during','after'))
disttrap$removal.period.akdecalc <- factor(disttrap$removal.period.akdecalc,
                                           levels=c('before','during','after'))


#Read in speed data:
speed<- readRDS(paste0(objdir,"/pig_weekly_speed_ctmm.rds"))

speedaer <- speed %>% filter(Removal.Type%in%c('aer')) %>% 
  filter(removal.period.akdecalc!='during') 
speedtox <- speed %>% filter(Removal.Type%in%c('tox')) 
speedtrap <- speed %>% filter(Removal.Type%in%c('trap')) 

#relevel when during isn't recorded
speedaer$removal.period.akdecalc <- factor(speedaer$removal.period.akdecalc,
                                           levels=c('before','after'))
speedtox$removal.period.akdecalc <- factor(speedtox$removal.period.akdecalc,
                                           levels=c('before','during','after'))
speedtrap$removal.period.akdecalc <- factor(speedtrap$removal.period.akdecalc,
                                            levels=c('before','during','after'))


# Format and join data ---------------------------------------------------------

#extract columns from each area set
akaer2=akaer[,c("animalid","Removal.Type","period", "area.est")]
aktrap2=aktrap[,c("animalid","Removal.Type","period", "area.est")]
aktox2=aktox[,c("animalid","Removal.Type","period", "area.est")]

#extract cols from NSD sets
geo.aerd.wk2=geo.aerd.wk[,c("animalid","Removal.Type","removal.period.akdecalc","week","mNSD")]
geo.trapd.wk2=geo.trapd.wk[,c("animalid","Removal.Type","removal.period.akdecalc","week","mNSD")]
geo.toxd.wk2=geo.toxd.wk[,c("animalid","Removal.Type","removal.period.akdecalc","week","mNSD")]

#extract cols from distance sets
distaer2=distaer[,c("animalid","trt_ctrl","removal.period.akdecalc","week","weekly_dist_km")]
disttrap2=disttrap[,c("animalid","trt_ctrl","removal.period.akdecalc","week","weekly_dist_km")]
disttox2=disttox[,c("animalid","trt_ctrl","removal.period.akdecalc","week","weekly_dist_km")]

#extract cols from speed sets
speedaer2=speedaer[,c("animalid","trt_ctrl","removal.period.akdecalc","week","weekly_md_km_hr")]
speedtrap2=speedtrap[,c("animalid","trt_ctrl","removal.period.akdecalc","week","weekly_md_km_hr")]
speedtox2=speedtox[,c("animalid","trt_ctrl","removal.period.akdecalc","week","weekly_md_km_hr")]

#make data frames for formatting
distaer2<-as.data.frame(distaer2)
disttrap2<-as.data.frame(disttrap2)
disttox2<-as.data.frame(disttox2)
speedaer2<-as.data.frame(speedaer2)
speedtrap2<-as.data.frame(speedtrap2)
speedtox2<-as.data.frame(speedtox2)

#rename cols for ease
colnames(akaer2)<-c("animalid","trt","per","area")
colnames(aktrap2)<-c("animalid","trt","per","area")
colnames(aktox2)<-c("animalid","trt","per","area")
colnames(geo.aerd.wk2)<-c("animalid","trt","per","wk","nsd")
colnames(geo.trapd.wk2)<-c("animalid","trt","per","wk","nsd")
colnames(geo.toxd.wk2)<-c("animalid","trt","per","wk","nsd")
colnames(distaer2)<-c("animalid","trt","per","wk","dist")
colnames(disttrap2)<-c("animalid","trt","per","wk","dist")
colnames(disttox2)<-c("animalid","trt","per","wk","dist")
colnames(speedaer2)<-c("animalid","trt","per","wk","speed")
colnames(speedtrap2)<-c("animalid","trt","per","wk","speed")
colnames(speedtox2)<-c("animalid","trt","per","wk","speed")

#make akaer ones as character to change the levels
akaer2$trt=as.character(akaer2$trt)
aktrap2$trt=as.character(aktrap2$trt)
aktox2$trt=as.character(aktox2$trt)

#character
distaer2$trt=as.character(distaer2$trt)
disttrap2$trt=as.character(disttrap2$trt)
disttox2$trt=as.character(disttox2$trt)

#character
speedaer2$trt=as.character(speedaer2$trt)
speedtrap2$trt=as.character(speedtrap2$trt)
speedtox2$trt=as.character(speedtox2$trt)
geo.aerd.wk2$trt=as.character(geo.aerd.wk2$trt)
geo.toxd.wk2$trt=as.character(geo.toxd.wk2$trt)
geo.trapd.wk2$trt=as.character(geo.trapd.wk2$trt)

#change removal type cols to be treatment or control
akaer2[,c("trt")][akaer2[,c("trt")]=="aer"]<-"trt"
aktrap2[,c("trt")][aktrap2[,c("trt")]=="trap"]<-"trt"
aktox2[,c("trt")][aktox2[,c("trt")]=="tox"]<-"trt"
geo.aerd.wk2[,c("trt")][geo.aerd.wk2[,c("trt")]=="aer"]<-"trt"
geo.trapd.wk2[,c("trt")][geo.trapd.wk2[,c("trt")]=="trap"]<-"trt"
geo.toxd.wk2[,c("trt")][geo.toxd.wk2[,c("trt")]=="tox"]<-"trt"

#left join to get one df for each removal type
aer=left_join(geo.aerd.wk2,akaer2,by=c("animalid","trt","per"))
trap=left_join(geo.trapd.wk2,aktrap2,by=c("animalid","trt","per"))
tox=left_join(geo.toxd.wk2,aktox2,by=c("animalid","trt","per"))

#left join to get one df for each removal type
aer2=left_join(aer,distaer2,by=c("animalid","trt","per","wk"))
trap2=left_join(trap,disttrap2,by=c("animalid","trt","per","wk"))
tox2=left_join(tox,disttox2,by=c("animalid","trt","per","wk"))

#left join to get one df for each removal type
aer3=left_join(aer2,speedaer2,by=c("animalid","trt","per","wk"))
trap3=left_join(trap2,speedtrap2,by=c("animalid","trt","per","wk"))
tox3=left_join(tox2,speedtox2,by=c("animalid","trt","per","wk"))

#Check incomplete cases
trap_gone=c("48479_T2_T2","85411_C3_C3","85454_1W_1W")
tox_gone=c("86070_H2_H2")
View(trap3[!complete.cases(trap3)&
             !(trap3$animalid%in%satpigs)&
             !(trap3$animalid%in%trap_gone),])
View(tox3[!complete.cases(tox3)&
            !(tox3$animalid%in%satpigs)&
            !(tox3$animalid%in%tox_gone),])

#Add removal type columns
aer3$rem<-"aer"
trap3$rem<-"trap"
tox3$rem<-"tox"

#rbind into one df
dat=rbind(aer3,trap3,tox3)

#reorder cols, want removal type in front
dat=dat[,c(9,1:8)]

# Summarize data ---------------------------------------------------------------

#want data summarized for each removal type, treatment, period, week
#get median and 25/75th percentiles for each response
quant_25<-function(x){
  x=quantile(x,probs=0.25,na.rm=TRUE)
}

quant_75<-function(x){
  quantile(x,probs=0.75,na.rm=TRUE)
}

med<-function(x){
  median(x,na.rm=TRUE)
}

dat2w=dat %>% dplyr::group_by(rem,trt,per,wk) %>% 
  dplyr::summarise_at(vars("nsd","dist","speed"),list(med=med,q25=quant_25,q75=quant_75))

dat2p=dat %>% dplyr::group_by(rem,trt,per) %>% 
  dplyr::summarise_at(vars("area"),list(med=med,q25=quant_25,q75=quant_75))

#remove NA row with tox during
dat2p=dat2p[complete.cases(dat2p),] 

#format cols for ggplot

# Set variables for plotting ---------------------------------------------------
tox.hex="#FF985A"
trap.hex="#FF57A9"
aer.hex="#822AFF"
ctrl.hex="#FFD065"

dat2w$hex=NA
dat2w[,c("hex")][dat2w[,c("rem")]=="aer"]<-aer.hex
dat2w[,c("hex")][dat2w[,c("rem")]=="trap"]<-trap.hex
dat2w[,c("hex")][dat2w[,c("rem")]=="tox"]<-tox.hex
dat2w[,c("hex")][dat2w[,c("trt")]=="ctrl"]<-ctrl.hex

dat2p$hex=NA
dat2p[,c("hex")][dat2p[,c("rem")]=="aer"]<-aer.hex
dat2p[,c("hex")][dat2p[,c("rem")]=="trap"]<-trap.hex
dat2p[,c("hex")][dat2p[,c("rem")]=="tox"]<-tox.hex
dat2p[,c("hex")][dat2p[,c("trt")]=="ctrl"]<-ctrl.hex

#Get start dates for each removal type in week terms for plot
#(below code from ./NRM_Descriptions.R)
trap.data.start=as.Date("2022-12-23")
trap.start.date=as.Date("2023-02-17")
trap.end.date=as.Date("2023-05-04")
trap.data.end=as.Date("2023-06-30")

aer.start.date=as.Date("2023-03-27")
aer.data.start=aer.start.date-37
aer.end.date=as.Date("2023-03-29")
aer.data.end=aer.end.date+37

tox.start.date=as.Date("2023-02-17")
tox.data.start=as.Date("2023-01-11")
tox.end.date=as.Date("2023-03-09")
tox.data.end=as.Date("2023-04-15")

#week 1 starts at start date
#make df with start date and day sequence through end date+1
#make col for week, numeric
#then find week corresponding with start and end dates for vline week numbers

#Match weeks with dates
#wk_dates=function(dstart,start,end,dend){
wk_dates=function(dstart,dend){
  dates=seq(from=dstart,to=dend+1,by=1)
  days=difftime(time1 = dates, time2 = min(dates), units = "days")+1
  days=difftime(time1 = dates, time2 = min(dates), units = "days")+1
  wks_num=days/7
  wks_int=as.integer(floor(days/7))
  #df=data.frame("dt"=dates,"nwks"=as.numeric(wks),"iwks"=wks_int)
  df=data.frame("wk"=unique(wks_int)+1,"dt"=seq(from=dstart,to=dend+1,by=7))
  return(df)
}

aer_wk_df=wk_dates(aer.data.start,aer.data.end)
trap_wk_df=wk_dates(trap.data.start,trap.data.end)
tox_wk_df=wk_dates(tox.data.start,tox.data.end)

#join with weeks to get actual start dates of week for timeline figs as x label
aer_w=left_join(dat2w[dat2w$rem=="aer",],aer_wk_df,by="wk")
trap_w=left_join(dat2w[dat2w$rem=="trap",],trap_wk_df,by="wk")
tox_w=left_join(dat2w[dat2w$rem=="tox",],tox_wk_df,by="wk")

#bind back
dat3w=rbind(aer_w,trap_w,tox_w)

#Add start/end date columns by removal type
dat3w$strt=NA
dat3w[,c("strt")][dat3w[,c("rem")]=="aer"]<-aer.start.date
dat3w[,c("strt")][dat3w[,c("rem")]=="trap"]<-trap.start.date
dat3w[,c("strt")][dat3w[,c("rem")]=="tox"]<-tox.start.date

dat3w$end=NA
dat3w[,c("end")][dat3w[,c("rem")]=="aer"]<-aer.end.date
dat3w[,c("end")][dat3w[,c("rem")]=="trap"]<-trap.end.date
dat3w[,c("end")][dat3w[,c("rem")]=="tox"]<-tox.end.date

#Relevel factors for dat2p
dat2p$per=forcats::fct_relevel(dat2p$per,c("before","during","after"))

# Make timeline figures --------------------------------------------------------

#* Make NSD figures ------------------------------------------------------------

nsd_aer=dat3w %>% filter(rem=="aer") %>%
ggplot(., aes(x=dt,y=nsd_med,ymin=nsd_q25,ymax=nsd_q75,color=hex,fill=hex))+
  scale_color_identity()+
  scale_fill_identity()+
  geom_line()+
  geom_ribbon(alpha=0.2,color=NA)+
  theme_ipsum()+
  geom_vline(aes(xintercept=strt),linetype="dashed")+
  geom_vline(aes(xintercept=end),linetype="dashed")+
  xlab("date")+
  scale_x_date(date_breaks = "1 week", date_labels =  "%b-%d")+ 
  theme(axis.text.x=element_text(angle=60, hjust=1))+
  ylab("NSD (m^2)")

nsd_trap=dat3w %>% filter(rem=="trap") %>%
  ggplot(., aes(x=dt,y=nsd_med,ymin=nsd_q25,ymax=nsd_q75,color=hex,fill=hex))+
  scale_color_identity()+
  scale_fill_identity()+
  geom_line()+
  geom_ribbon(alpha=0.2,color=NA)+
  theme_ipsum()+
  geom_vline(aes(xintercept=strt),linetype="dashed")+
  geom_vline(aes(xintercept=end),linetype="dashed")+
  xlab("date")+
  scale_x_date(date_breaks = "2 week", date_labels =  "%b-%d")+ 
  theme(axis.text.x=element_text(angle=60, hjust=1))+
  ylab("NSD (m^2)")

nsd_tox=dat3w %>% filter(rem=="tox") %>%
  ggplot(., aes(x=dt,y=nsd_med,ymin=nsd_q25,ymax=nsd_q75,color=hex,fill=hex))+
  scale_color_identity()+
  scale_fill_identity()+
  geom_line()+
  geom_ribbon(alpha=0.2,color=NA)+
  theme_ipsum()+
  geom_vline(aes(xintercept=strt),linetype="dashed")+
  geom_vline(aes(xintercept=end),linetype="dashed")+
  xlab("date")+
  scale_x_date(date_breaks = "1 week", date_labels =  "%b-%d")+ 
  theme(axis.text.x=element_text(angle=60, hjust=1))+
  ylab("NSD (m^2)")

#* Make Area figures ------------------------------------------------------------

hr_aer=dat2p %>% filter(rem=="aer") %>%
  ggplot(., aes(x=per,ymin=q25,yend=q75,color=hex,fill=hex))+
  geom_point(aes(y=med),size=3,position = position_dodge(width = 0.1))+
  geom_segment(aes(y=q25),linewidth=1,position = position_dodge(width = 0.1))+
  scale_color_identity()+
  theme_ipsum()+
  theme(legend.position="none")+
  xlab("period")+
  ylab("HR area (km^2)")

hr_trap=dat2p %>% filter(rem=="trap") %>%
  ggplot(., aes(x=per,ymin=q25,yend=q75,color=hex,fill=hex))+
  geom_point(aes(y=med),size=3,position = position_dodge(width = 0.1))+
  geom_segment(aes(y=q25),linewidth=1,position = position_dodge(width = 0.1))+
  scale_color_identity()+
  theme_ipsum()+
  theme(legend.position="none")+
  xlab("period")+
  ylab("HR area (km^2)")

hr_tox=dat2p %>% filter(rem=="tox") %>%
  ggplot(., aes(x=per,ymin=q25,yend=q75,color=hex,fill=hex))+
  geom_point(aes(y=med),size=3,position = position_dodge(width = 0.1))+
  geom_segment(aes(y=q25),linewidth=1,position = position_dodge(width = 0.1))+
  scale_color_identity()+
  theme_ipsum()+
  theme(legend.position="none")+
  xlab("period")+
  ylab("HR area (km^2)")


#* Make distance figures ------------------------------------------------------------

dist_aer=dat3w %>% filter(rem=="aer") %>%
  ggplot(., aes(x=dt,y=dist_med,ymin=dist_q25,ymax=dist_q75,color=hex,fill=hex))+
  scale_color_identity()+
  scale_fill_identity()+
  geom_line()+
  geom_ribbon(alpha=0.2,color=NA)+
  theme_ipsum()+
  geom_vline(aes(xintercept=strt),linetype="dashed")+
  geom_vline(aes(xintercept=end),linetype="dashed")+
  xlab("date")+
  scale_x_date(date_breaks = "1 week", date_labels =  "%b-%d")+ 
  theme(axis.text.x=element_text(angle=60, hjust=1))+
  ylab("dist (m)")

dist_trap=dat3w %>% filter(rem=="trap") %>%
  ggplot(., aes(x=dt,y=dist_med,ymin=dist_q25,ymax=dist_q75,color=hex,fill=hex))+
  scale_color_identity()+
  scale_fill_identity()+
  geom_line()+
  geom_ribbon(alpha=0.2,color=NA)+
  theme_ipsum()+
  geom_vline(aes(xintercept=strt),linetype="dashed")+
  geom_vline(aes(xintercept=end),linetype="dashed")+
  xlab("date")+
  scale_x_date(date_breaks = "2 week", date_labels =  "%b-%d")+ 
  theme(axis.text.x=element_text(angle=60, hjust=1))+
  ylab("dist (m)")

dist_tox=dat3w %>% filter(rem=="tox") %>%
  ggplot(., aes(x=dt,y=dist_med,ymin=dist_q25,ymax=dist_q75,color=hex,fill=hex))+
  scale_color_identity()+
  scale_fill_identity()+
  geom_line()+
  geom_ribbon(alpha=0.2,color=NA)+
  theme_ipsum()+
  geom_vline(aes(xintercept=strt),linetype="dashed")+
  geom_vline(aes(xintercept=end),linetype="dashed")+
  xlab("date")+
  scale_x_date(date_breaks = "1 week", date_labels =  "%b-%d")+ 
  theme(axis.text.x=element_text(angle=60, hjust=1))+
  ylab("dist (m)")

#* Make speed figures ------------------------------------------------------------

speed_aer=dat3w %>% filter(rem=="aer") %>%
  ggplot(., aes(x=dt,y=speed_med,ymin=speed_q25,ymax=speed_q75,color=hex,fill=hex))+
  scale_color_identity()+
  scale_fill_identity()+
  geom_line()+
  geom_ribbon(alpha=0.2,color=NA)+
  theme_ipsum()+
  geom_vline(aes(xintercept=strt),linetype="dashed")+
  geom_vline(aes(xintercept=end),linetype="dashed")+
  xlab("date")+
  scale_x_date(date_breaks = "1 week", date_labels =  "%b-%d")+ 
  theme(axis.text.x=element_text(angle=60, hjust=1))+
  ylab("speed (km/hr)")

speed_trap=dat3w %>% filter(rem=="trap") %>%
  ggplot(., aes(x=dt,y=speed_med,ymin=speed_q25,ymax=speed_q75,color=hex,fill=hex))+
  scale_color_identity()+
  scale_fill_identity()+
  geom_line()+
  geom_ribbon(alpha=0.2,color=NA)+
  theme_ipsum()+
  geom_vline(aes(xintercept=strt),linetype="dashed")+
  geom_vline(aes(xintercept=end),linetype="dashed")+
  xlab("date")+
  scale_x_date(date_breaks = "1 week", date_labels =  "%b-%d")+ 
  theme(axis.text.x=element_text(angle=60, hjust=1))+
  ylab("speed (km/hr)")

speed_tox=dat3w %>% filter(rem=="tox") %>%
  ggplot(., aes(x=dt,y=speed_med,ymin=speed_q25,ymax=speed_q75,color=hex,fill=hex))+
  scale_color_identity()+
  scale_fill_identity()+
  geom_line()+
  geom_ribbon(alpha=0.2,color=NA)+
  theme_ipsum()+
  geom_vline(aes(xintercept=strt),linetype="dashed")+
  geom_vline(aes(xintercept=end),linetype="dashed")+
  xlab("date")+
  scale_x_date(date_breaks = "1 week", date_labels =  "%b-%d")+ 
  theme(axis.text.x=element_text(angle=60, hjust=1))+
  ylab("speed (km/hr)")

dist_aer
dist_trap
dist_tox

speed_aer
speed_trap
speed_tox

nsd_aer
nsd_trap
nsd_tox

# Combine plots --------------------------------------------------------

cowplot::plot_grid(hr_aer,nsd_aer,ncol=2)
cowplot::plot_grid(hr_trap,nsd_trap,ncol=2)
cowplot::plot_grid(hr_tox,nsd_tox,ncol=2)


