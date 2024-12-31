
#Need revise this a bit more
#Pipeline in proc data folder:
  #1_Prep_Geolocation_Data
    #-does removal designations, periods, week splits
    #-remove week trimming from this script!
  #2_Get_Period_AKDEs
    #-trims pig period data as needed
    #-fits ctmm for each pig/period/removal
    #-outputs area df

#Next steps for tomorrow:
#3-finish cleaning up 4_Preanalysisformatting (should just be NSD calcs, move to NSD folder)

#set home dir of pipeline
home<-"/Users/kayleigh.chalkowski/OneDrive - USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline"

# Purpose ----------------------------------------------------------------------

#The purpose of this script is to format geolocation data with all designations:
  #removal type
  #period
  #week

# Process ----------------------------------------------------------------------

#In: geo.csv, geolocation data
#Out: geotox.rds, geoaer.rds, geotrap.rds 

# Removal designations
#1. pull geolocs, needed objects 
#2. get trap/tox chulls and flight paths
#3. get mcps for all pigs, respective to each treatment area and summarise overlap
    #a-50% mcps, to determine treatment pigs
    #b-95% mcps, to determine control pigs
#4. summarize pigs not included in any treatment
#5. write out geolocation data with treatment/control designations

# Period

# Week

# Setup ----------------------------------------------------------------------

#load libraries
library(amt)
library(dplyr)
library(stringr)
library(sf)

#set directories
input=file.path(home,"1_Data","Input",fsep=.Platform$file.sep)
objdir=file.path(home,"1_Data","Objects",fsep=.Platform$file.sep)
funcdir<-file.path(home,"2_Scripts","Functions",fsep=.Platform$file.sep)

#read objects from input
trapntox_sites.raw=readRDS(file.path(input,"trapntox_sites.raw.rds",fsep=.Platform$file.sep)) #trap/tox site locs
geo=readRDS(file.path(input,"geo.rds",fsep=.Platform$file.sep)) #geolocations, tidied and filtered
sd=readRDS(file.path(input,"sd.rds",fsep=.Platform$file.sep)) #study site polygons list
fp.chulls=readRDS(file.path(input,"fp.chulls.rds",fsep=.Platform$file.sep)) #flight path convex hulls

#source needed functions
func.list=list.files(funcdir,full.names=TRUE)
for(f in 1:length(func.list)){
  source(func.list[f])
}

#Set dates for before, during, after periods for each removal type
#trap dates
activities.raw=readRDS(file.path(input,"activities.rds",fsep=.Platform$file.sep))
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

#Format sites data
tox.sites=trapntox_sites.raw[trapntox_sites.raw$activity=="toxic"&!is.na(trapntox_sites.raw$trap_status),]
trap.sites=trapntox_sites.raw[trapntox_sites.raw$activity=="trap"&!is.na(trapntox_sites.raw$trap_status),]

tox.sites=tox.sites[!is.na(tox.sites$x),]
trap.sites=trap.sites[!is.na(trap.sites$x),]

tox.sites=st_as_sf(tox.sites,coords=c(6,7),crs=st_crs(32614))
trap.sites=st_as_sf(trap.sites,coords=c(6,7),crs=st_crs(32614))

# Determine removal designations -----------------------------------------------

# * Get removal area convex hulls ----------------------------------------------

trap.chull=st_convex_hull(st_union(trap.sites))
tox.chull=st_convex_hull(st_union(tox.sites))

trap.chull = trap.chull %>%
  st_sf %>%
  st_cast

tox.chull = tox.chull %>%
  st_sf %>%
  st_cast

#aerial convex hull

# * Get pig MCPs respective to treatments ----------------------------------------

#d-get week grouping
geo$jDate <- julian(as.Date(geo$datetime), origin = min(as.Date(geo$datetime)))
geo$jDate=geo$jDate+1
geo=transform(geo, week = cut(jDate, c(0, seq(7, max(jDate) + 7, 7)), labels = FALSE))
geo$wk_id=paste(geo$animalid,geo$week,sep="_")
geo$datetime<-Neat.Dates.POSIXct(geo$datetime,tz="UTC")

geo.sf<-st_as_sf(geo,coords=c("latitude","longitude"),crs=st_crs(4326))
geo.sf<-st_transform(geo.sf,st_crs(32614))
geo$X<-st_coordinates(geo.sf)[,1]
geo$Y<-st_coordinates(geo.sf)[,2]

#get into track format
tkmt=mk_track(geo,.x=X,.y=Y,.t=datetime,all_cols=TRUE,crs=32614)

#Get dateonly format
geo$date_only=as.Date(geo$datetime)

#make summary to get week number within start/end dates for treatments
wkdates=geo %>% group_by(week) %>% summarise(mindt=min(date_only),maxdt=max(date_only)) %>% as.data.frame()
jdates=geo %>% group_by(jDate) %>% summarise(dates=min(date_only)) %>% as.data.frame()

#nest tracks by animalid
tkmt.n=tkmt |> nest(data=-animalid)

#filter tkmt tracks by each cutoff interval
tkmt.n=tkmt.n |>
  mutate(trap.trks = map(data, function(x)
    x |>
      filter(jDate>=which(jdates$dates==trap.start.date)&jDate<=which(jdates$dates==trap.end.date))
  ))

tkmt.n=tkmt.n |>
  mutate(aer.trks = map(data, function(x)
    x |>
      filter(jDate>=which(jdates$dates==aer.start.date)&
               jDate<=which(jdates$dates==aer.end.date))
  ))

tkmt.n=tkmt.n |>
  mutate(tox.trks = map(data, function(x)
    x |>
      filter(jDate>=which(jdates$dates==tox.start.date)&
               jDate<=which(jdates$dates==tox.end.date))
  ))

#Get MCPs for each set of tracks, cutoff by trtmnt area dates
mcp.iso=0.5
trap.mcps=GetAreaMetrics(tkmt.n$trap.trks,"mcp",mcp.iso,tkmt.n$animalid)
aer.mcps=GetAreaMetrics(tkmt.n$aer.trks,"mcp",mcp.iso,tkmt.n$animalid)
tox.mcps=GetAreaMetrics(tkmt.n$tox.trks,"mcp",mcp.iso,tkmt.n$animalid)

#make data frame with ID, proportion overlap, area overlap
trap.overlaps=GetRemovalOverlaps(trap.mcps$mcp,trap.mcps$mcp$animalid,trap.chull)
aer.overlaps=GetRemovalOverlaps(aer.mcps$mcp,aer.mcps$mcp$animalid,st_union(fp.chulls)) #error
tox.overlaps=GetRemovalOverlaps(tox.mcps$mcp,tox.mcps$mcp$animalid,tox.chull)

trap.trt=trap.overlaps[which(trap.overlaps$prop_overlap!=0),]
aer.trt=aer.overlaps[which(aer.overlaps$prop_overlap!=0),]
tox.trt=tox.overlaps[which(tox.overlaps$prop_overlap!=0),]

trt.sums=data.frame(
  removal=c("trap","aerial","tox"),
  pigs=c(nrow(trap.trt),nrow(aer.trt),nrow(tox.trt)),
  avg_prop_overlap=c(mean(trap.trt$prop_overlap),mean(aer.trt$prop_overlap),mean(tox.trt$prop_overlap)),
  avg_area_overlap_km=c(mean(trap.trt$area_overlap)/1e6,mean(aer.trt$area_overlap)/1e6,mean(tox.trt$area_overlap)/1e6)
)

# * Control pig designations ---------------------------------------------------
#Control site designations
#Use same tracks from before subset to each treatment area
#Get 95% MCPs
#Determine which to call controls
#a-are there enough pigs that didn't overlap with any treatment, at any period?
#b-what pigs can be considered control for each treatment period


#Get MCPs for each set of tracks, cutoff by trtmnt area dates
mcp.iso=0.95
trap.mcps.95=GetAreaMetrics(tkmt.n$trap.trks,"mcp",mcp.iso,tkmt.n$animalid)
aer.mcps.95=GetAreaMetrics(tkmt.n$aer.trks,"mcp",mcp.iso,tkmt.n$animalid)
tox.mcps.95=GetAreaMetrics(tkmt.n$tox.trks,"mcp",mcp.iso,tkmt.n$animalid)

#make data frame with ID, proportion overlap, area overlap
trap.overlaps.95=GetRemovalOverlaps(trap.mcps.95$mcp,tox.mcps.95$mcp$animalid,trap.chull)
aer.overlaps.95=GetRemovalOverlaps(aer.mcps.95$mcp,aer.mcps.95$mcp$animalid,st_union(fp.chulls)) #error
tox.overlaps.95=GetRemovalOverlaps(tox.mcps.95$mcp,tox.mcps.95$mcp$animalid,tox.chull)

trap.ctrl=trap.overlaps.95[which(trap.overlaps.95$prop_overlap==0),]
aer.ctrl=aer.overlaps.95[which(aer.overlaps.95$prop_overlap==0),]
tox.ctrl=tox.overlaps.95[which(tox.overlaps.95$prop_overlap==0),]

trap.ctrl$trt="trap"
aer.ctrl$trt="aer"
tox.ctrl$trt="tox"

ctrls=rbind(trap.ctrl[,c(1,4)],aer.ctrl[,c(1,4)],tox.ctrl[,c(1,4)])

ctrlsw=xtabs( ~ animalid + trt, ctrls)
ctrlsw=as.data.frame.matrix(ctrlsw)

ctrlsw$all=rowSums(ctrlsw)

ct=ctrlsw

nrow(ct[ct$all==3,])
nrow(ct[ct$all==2,])
nrow(ct[ct$all==1,])

######Show MCPs of pigs which do not overlap with any treatment areas
ctrl.all.1=trap.mcps.95$mcp[which(trap.mcps.95$mcp$animalid%in%rownames(ct[ct$all==3,])),]
ctrl.all.2=aer.mcps.95$mcp[which(aer.mcps.95$mcp$animalid%in%rownames(ct[ct$all==3,])),]
ctrl.all.3=tox.mcps.95$mcp[which(tox.mcps.95$mcp$animalid%in%rownames(ct[ct$all==3,])),]

# * Summarize pigs not incl in any treatment/ctrl ------------------------------

#IDs
ctrl.pig.IDs=rownames(ct[ct$all==3,])
trap.pig.IDs=trap.trt$animalid
aer.pig.IDs=aer.trt$animalid
tox.pig.IDs=tox.trt$animalid

#sanity check, verify that control pigs not also in treatments
any(ctrl.pig.IDs%in%trap.pig.IDs) #FALSE
any(ctrl.pig.IDs%in%aer.pig.IDs) #FALSE
any(ctrl.pig.IDs%in%tox.pig.IDs) #FALSE

#any overlapping treatments?
any(trap.pig.IDs%in%aer.pig.IDs) #FALSE
any(trap.pig.IDs%in%tox.pig.IDs) #FALSE
any(tox.pig.IDs%in%trap.pig.IDs) #FALSE
any(tox.pig.IDs%in%aer.pig.IDs) #FALSE
any(aer.pig.IDs%in%trap.pig.IDs) #FALSE
any(aer.pig.IDs%in%tox.pig.IDs) #FALSE

all.selected.IDs=c(ctrl.pig.IDs,trap.pig.IDs,tox.pig.IDs,aer.pig.IDs)

not.selected=unique(geo$animalid)[which(!unique(geo$animalid)%in%all.selected.IDs)]

#cutoff by date:
#48462_E6_E6
#48467_H6_H6
#48476_S1_S1
#85396_W2_W2
#85401_A1_A1
#85412_X1_X1
#85429_U3_U4
#85434_2H_2H
#85736_6X_6X

#pigs not selected:
length(not.selected) #14 total
no.dt.overlap=c("48462_E6_E6", #9 were not collared during any treatment periods
                "48467_H6_H6",
                "48476_S1_S1",
                "85396_W2_W2",
                "85401_A1_A1",
                "85412_X1_X1",
                "85429_U3_U4",
                "85434_2H_2H",
                "85736_6X_6X")

not.sel.spat=not.selected[-which(not.selected%in%no.dt.overlap)]
gns=geo[geo$animalid%in%not.sel.spat,]
gnsf=st_as_sf(gns,coords=c("X","Y"),crs=st_crs(32614))

#base+mapview(gnsf,zcol="animalid")

gns %>% group_by(animalid) %>% summarise(mindt=min(date_only),maxdt=max(date_only))
#85403_S7_S7 near aerial, but data ends 2/28
#85396_2_L2_L2 near aerial, data ends 2/28

#48464_E5_E5 overlaps, but must not have been overlap with 50% MCP core
#base+chulls+mapview(gnsf[gnsf$animalid=="48464_E5_E5",] |>
#                      filter(jDate>=which(jdates$dates==trap.start.date)&jDate<=which(jdates$dates==trap.end.date)))

#85434_2_3K_3K didn't have enough points during aerial gunning to properly estimate MCP
#no points during that time period available
#base+chulls+mapview(gnsf[gnsf$animalid=="85434_2_3K_3K",] |>
#                      filter(jDate>=which(jdates$dates==aer.start.date)&jDate<=which(jdates$dates==aer.end.date)))

#because was the satellite one, poor resolution, points every 4 hours
#chulls+mapview(gnsf[gnsf$animalid=="85434_2_3K_3K",])

#48461_B4_B4 overlaps, but must have been no overlap with 50% MCP core
#base+chulls+mapview(gnsf[gnsf$animalid=="48461_B4_B4",] |>
                     #filter(jDate>=which(jdates$dates==trap.start.date)&jDate<=which(jdates$dates==trap.end.date)))

#*48447_3Y_3Y overlaps a bit, but must have been no overlap with 50% MCP core
#base+chulls+mapview(gnsf[gnsf$animalid=="48447_3Y_3Y",] |>
#                      filter(jDate>=which(jdates$dates==trap.start.date)&jDate<=which(jdates$dates==trap.end.date)))
#base+chulls+mapview(gnsf[gnsf$animalid=="48447_3Y_3Y",] |>
#                      filter(jDate>=which(jdates$dates==tox.start.date)&jDate<=which(jdates$dates==tox.end.date)))
#* Note: ^^was added to treatment area after geolocation tidying fix

#Summary:
#15 pigs not used for any treatment or ctrls
#Of these 15, 
#9 of these pigs were not collared during any of the treatment dates
#2 pigs overlapped with aerial but collar data ends before aerial treatment dates
#2 pigs did not have 50% core MCP home ranges that overlapped with any treatment area
#1 pig did not have enough points to calculate 50% MCP core range overlap in treatment area of interest
#1 pig removed earlier during tidying/filtering process-- malfunctioning mortality signal

#Note: actually 14, 48447_3Y_3Y moved to trap after geolocation tidying fix

#Removing three more pigs from aer: 
#aerial pigs that were accidentally culled
aer.pig.IDs=aer.pig.IDs[aer.pig.IDs!="48476_2_4Y_4Y"&
              aer.pig.IDs!="85401_2_U_U"&
              aer.pig.IDs!="85440_E2_E2"]
not.selected=c(not.selected,"48476_2_4Y_4Y","85401_2_U_U","85440_E2_E2")

#tox pig not trapped in before tox period:
tox.pig.IDs=tox.pig.IDs[tox.pig.IDs!="86070_H2_H2"]
not.selected=c(not.selected,"86070_H2_H2")

# * Link designations to geolocs -----------------------------------------------

geo$Removal.Type=NA

geo[geo$animalid%in%ctrl.pig.IDs,]$Removal.Type="ctrl"
geo[geo$animalid%in%trap.pig.IDs,]$Removal.Type="trap"
geo[geo$animalid%in%aer.pig.IDs,]$Removal.Type="aer"
geo[geo$animalid%in%tox.pig.IDs,]$Removal.Type="tox"

#Remove pigs with NA removal types-- not designated
geo=geo[!is.na(geo$Removal.Type),]
#nrow(geo) #1416606

#Clean up unneeded columns used for removal designations
#will redo week/day designations later relative to each removal area
geo=geo[,-c(which(colnames(geo)=="jDate"),
        which(colnames(geo)=="week"),
        which(colnames(geo)=="wk_id"))]

# Determine period designations ------------------------------------------------

# * Get cutoff dates for each removal type--------------------------------------
#Trim data before/after certain point to keep period timespans consistent

#Aerial 
#determined 4 weeks from sensitivity analysis for how much time needed to get accurate akde est
aer.cutoff1=aer.start.date-37
aer.cutoff2=aer.end.date+37

#Trap
#use length of time in trap 'after' period (smaller than during period)
trap.len=as.Date(max(geo[geo$Removal.Type=="trap",]$date_only))-as.Date(trap.end.date)-1
#use as trap start date: 56 days before trap end date
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

# * Set periods according to date cutoffs --------------------------------------

#Note: this step separates data for each removal type into sep df's
  #Needed bc each removal type has different dates for each period, and need 
  #associated ctrls to be same respective timespans

#Make function to do it
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

#Use function for aer/tox
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

# Determine week designations --------------------------------------------------

# * Confirm all start on same date ---------------------------------------------

geo.aer %>% group_by(animalid) %>% 
  dplyr::summarise(strtdt=min(date_only)) %>%
  as.data.frame()
geo.trap %>% group_by(animalid) %>% 
  dplyr::summarise(strtdt=min(date_only)) %>%
  as.data.frame()
geo.tox %>% group_by(animalid) %>% 
  dplyr::summarise(strtdt=min(date_only)) %>%
  as.data.frame()

# * Get dates for day 1 for each removal period/type ---------------------------
#Trap
trap.origin=trap.cutoff1 #before period start date
trap.origin.after=trap.end.date+1 #after period start date
trap.origin.during=trap.end.date-56 #during period start date

#Tox
tox.origin=tox.cutoff1 #before period start date
tox.origin.after=tox.end.date+1 #after period start date
tox.origin.during=tox.start.date #during period start date

#Aerial
aer.origin=aer.cutoff1 #before period start date
aer.origin.after=aer.end.date+1 #after period start date

# * Run week splits using period origins ---------------------------------------

#Make function to do splits by week
Do.Week.Split<-function(geo.aerd,removal.str,origin.vector){
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

#Split by week
geo.aer=Do.Week.Split(geo.aer,"aer",c(aer.origin,aer.origin.after))
geo.trap=Do.Week.Split(geo.trap,"trap",c(trap.origin,trap.origin.after,trap.origin.during))
geo.tox=Do.Week.Split(geo.tox,"tox",c(tox.origin,tox.origin.after,tox.origin.during))

# Write out data and summaries -------------------------------------------------

#write out csv
out.dir<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/"
#write.csv(geo,paste0(out.dir,"geo_remtyp.csv"))
saveRDS(geo.tox,file.path(objdir,"geotox.rds"))
saveRDS(geo.trap,file.path(objdir,"geotox.rds"))
saveRDS(geo.aer,file.path(objdir,"geotox.rds"))


