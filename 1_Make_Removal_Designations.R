
##############################################
################## Overview ################## 
##############################################

#1. pull geolocs, needed objects 
#2. get trap/tox chulls and flight paths
#3. get mcps for all pigs, respective to each treatment area and summarise overlap
    #a-50% mcps, to determine treatment pigs
    #b-95% mcps, to determine control pigs
#4. summarize pigs not included in any treatment
#5. write out geolocation data with treatment/control designations

###########################################
################## Setup ################## 
###########################################

#load libraries
library(amt)
library(dplyr)
library(stringr)
library(sf)

#set directories
#sdir<-"/Volumes/Projects/MUDD/ASF_NIFA/Datasets/Movement data/Tidied_WP_Geolocations"
#filtered<-"/Volumes/Projects/MUDD/ASF_NIFA/Datasets/Movement data/Tidied_WP_Geolocations/4_filtered/"
input="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/1_Data/Input/"

#read objects from input
trapntox_sites.raw=readRDS(paste0(input,"trapntox_sites.raw.rds")) #trap/tox site locs
geo=readRDS(paste0(input,"geo.rds")) #geolocations, tidied and filtered
sd=readRDS(paste0(input,"sd.rds")) #study site polygons list
fp.chulls=readRDS(paste0(input,"fp.chulls.rds")) #flight path convex hulls

#source needed functions
funcdir<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/2_Scripts/Functions"
func.list=list.files(funcdir,full.names=TRUE)
for(f in 1:length(func.list)){
  source(func.list[f])
}

#a-set removal start/end dates
#trap dates
trap.act=activities.raw[activities.raw$activity=="trap",]
trap.act$trap_datetime=Neat.Dates.POSIXct(trap.act$trap_datetime,tz="UTC")
trap.act[trap.act$trap_datetime==min(trap.act$trap_datetime),]
trap.start.date=as.Date(min(trap.act$trap_datetime))
trap.end.date=as.Date(max(trap.act$trap_datetime))

#tox dates
tox.act=activities.raw[activities.raw$activity=="toxic",]
tox.act$tox_datetime=Neat.Dates.POSIXct(tox.act$tox_datetime,tz="UTC")
tox.act[tox.act$tox_datetime==min(tox.act$tox_datetime),]
tox.start.date=as.Date(min(tox.act$tox_datetime))
tox.end.date=as.Date("2023-03-09")

#aerial dates
aer.start.date=as.Date("2023-03-27")
aer.end.date=as.Date("2023-03-29")

#load site data
#trapntox_sites.raw=readRDS("/Volumes/Projects/MUDD/ASF_NIFA/Datasets/Daily_Activities_Trap_Locations/Activities_Tidying_Pipeline_KK/Data/Final/site.locations_tidied.RDS")

#################################
#### Read and format objects ####
#################################
#tidied<-"/Volumes/Projects/MUDD/ASF_NIFA/Datasets/Movement data/Tidied_WP_Geolocations/3_tidied/"

#load data
#geolocations
#filtered<-"/Volumes/Projects/MUDD/ASF_NIFA/Datasets/Movement data/Tidied_WP_Geolocations/4_filtered/"
#geo=read.csv(paste0(filtered,"geolocs.csv"))

#study site polygons
#sd=loadNIFAspatial(c("cd_outline","aerial_pastures","ctrl_pastures","trap_pastures","eup_pastures"))
ap=sd$aerial_pastures

#flight path convex hulls
#fp.chulls=readRDS("/Volumes/Projects/MUDD/ASF_NIFA/Datasets/Aerial_Flight_Paths/tidy_flightpaths/aerial_flightpaths_chull.rds")


###################################################################
################## Get removal area convex hulls ################## 
###################################################################

#fp.chull.adj=st_convex_hull(st_union(fpts2)) %>% st_sf %>% st_cast

#assign to old variable name to match code below
#fp.chulls=readRDS("smb://aapcoftc3fp13/Projects/MUDD/ASF_NIFA/Datasets/Aerial_Flight_Paths/tidy_flightpaths/aerial_flightpaths_chull.rds")

tox.sites=trapntox_sites.raw[trapntox_sites.raw$activity=="toxic"&!is.na(trapntox_sites.raw$trap_status),]
trap.sites=trapntox_sites.raw[trapntox_sites.raw$activity=="trap"&!is.na(trapntox_sites.raw$trap_status),]

tox.sites=tox.sites[!is.na(tox.sites$x),]
trap.sites=trap.sites[!is.na(trap.sites$x),]

tox.sites=st_as_sf(tox.sites,coords=c(6,7),crs=st_crs(32614))
trap.sites=st_as_sf(trap.sites,coords=c(6,7),crs=st_crs(32614))

trap.chull=st_convex_hull(st_union(trap.sites))
tox.chull=st_convex_hull(st_union(tox.sites))

trap.chull = trap.chull %>%
  st_sf %>%
  st_cast

tox.chull = tox.chull %>%
  st_sf %>%
  st_cast

############################################################################### 
################## Get Mcps for pigs, respective to treatments ################ 
############################################################################### 

#d-get week grouping
geo$jDate <- julian(as.Date(geo$datetime), origin = min(as.Date(geo$datetime)))
geo$jDate=geo$jDate+1
geo=transform(geo, week = cut(jDate, c(0, seq(7, max(jDate) + 7, 7)), labels = FALSE))
geo$wk_id=paste(geo$animalid,geo$week,sep="_")

#geowk=geo %>% group_by(wk_id) %>% summarise(wkmed_latitude=median(latitude),wkmed_longitude=median(longitude)) %>% as.data.frame()
geo$datetime<-Neat.Dates.POSIXct(geo$datetime,tz="UTC")

geo.sf<-st_as_sf(geo,coords=c("latitude","longitude"),crs=st_crs(4326))
geo.sf<-st_transform(geo.sf,st_crs(32614))
geo$X<-st_coordinates(geo.sf)[,1]
geo$Y<-st_coordinates(geo.sf)[,2]

#I can do mcps but for now going to use hr_area, have code for it alreaday, is based on a kde (95%)
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


#########################################################################
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


###############################################################################
################## Summarise pigs not incl in any trtment/ctrl ################ 
###############################################################################

#IDs
ctrl.pig.IDs=rownames(ct[ct$all==3,])
trap.pig.IDs=trap.trt$animalid
aer.pig.IDs=aer.trt$animalid
tox.pig.IDs=tox.trt$animalid

#sanity check, verify that control pigs not also in treatments
any(ctrl.pig.IDs%in%trap.pig.IDs)
any(ctrl.pig.IDs%in%aer.pig.IDs)
any(ctrl.pig.IDs%in%tox.pig.IDs)

#any overlapping treatments?
any(trap.pig.IDs%in%aer.pig.IDs)
any(trap.pig.IDs%in%tox.pig.IDs)
any(tox.pig.IDs%in%trap.pig.IDs)
any(tox.pig.IDs%in%aer.pig.IDs)
any(aer.pig.IDs%in%trap.pig.IDs)
any(aer.pig.IDs%in%tox.pig.IDs)

all.selected.IDs=c(ctrl.pig.IDs,trap.pig.IDs,tox.pig.IDs,aer.pig.IDs)

not.selected=unique(geo$animalid)[which(!unique(geo$animalid)%in%all.selected.IDs)]

geo[geo$animalid%in%not.selected,] %>% group_by(animalid) %>% summarise(mindt=min(date_only),maxdt=max(date_only))

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
length(not.selected) #15 total
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
gnsf=st_as_sf(gns,coords=c(1,15),crs=st_crs(32614))

base+mapview(gnsf,zcol="animalid")

gns %>% group_by(animalid) %>% summarise(mindt=min(date_only),maxdt=max(date_only))
#85403_S7_S7 near aerial, but data ends 2/28
#85396_2_L2_L2 near aerial, data ends 2/28

#48464_E5_E5 overlaps, but must not have been overlap with 50% MCP core
base+chulls+mapview(gnsf[gnsf$animalid=="48464_E5_E5",] |>
                      filter(jDate>=which(jdates$dates==trap.start.date)&jDate<=which(jdates$dates==trap.end.date)))

#85434_2_3K_3K didn't have enough points during aerial gunning to properly estimate MCP
#no points during that time period available
base+chulls+mapview(gnsf[gnsf$animalid=="85434_2_3K_3K",] |>
                      filter(jDate>=which(jdates$dates==aer.start.date)&jDate<=which(jdates$dates==aer.end.date)))

#because was the satellite one, poor resolution, points every 4 hours
chulls+mapview(gnsf[gnsf$animalid=="85434_2_3K_3K",])

#48461_B4_B4 overlaps, but must have been no overlap with 50% MCP core
base+chulls+mapview(gnsf[gnsf$animalid=="48461_B4_B4",] |>
                      filter(jDate>=which(jdates$dates==trap.start.date)&jDate<=which(jdates$dates==trap.end.date)))

#48447_3Y_3Y overlaps a bit, but must have been no overlap with 50% MCP core
#base+chulls+mapview(gnsf[gnsf$animalid=="48447_3Y_3Y",] |>
#                      filter(jDate>=which(jdates$dates==trap.start.date)&jDate<=which(jdates$dates==trap.end.date)))
#base+chulls+mapview(gnsf[gnsf$animalid=="48447_3Y_3Y",] |>
#                      filter(jDate>=which(jdates$dates==tox.start.date)&jDate<=which(jdates$dates==tox.end.date)))
#^^was added to treatment area after geolocation tidying fix


#Summary:
#15 pigs not used for any treatment or ctrls
#Of these 15, 
#9 of these pigs were not collared during any of the treatment dates
#2 pigs overlapped with aerial but collar data ends before aerial treatment dates
#2 pigs did not have 50% core MCP home ranges that overlapped with any treatment area
#1 pig did not have enough points to calculate 50% MCP core range overlap in treatment area of interest
#1 pig removed earlier during tidying/filtering process-- malfunctioning mortality signal

#################################################################################################
#link designations to geolocations

geo$Removal.Type=NA

geo[geo$animalid%in%ctrl.pig.IDs,]$Removal.Type="ctrl"
geo[geo$animalid%in%trap.pig.IDs,]$Removal.Type="trap"
geo[geo$animalid%in%aer.pig.IDs,]$Removal.Type="aer"
geo[geo$animalid%in%tox.pig.IDs,]$Removal.Type="tox"

#write out csv
out.dir<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/"
write.csv(geo,paste0(out.dir,"geo_remtyp.csv"))
write.csv(trap.overlaps.95,paste0(out.dir,"trap_mcp_overlaps_95_summary.csv"))
write.csv(tox.overlaps.95,paste0(out.dir,"tox_mcp_overlaps_95_summary.csv"))
write.csv(aer.overlaps.95,paste0(out.dir,"aer_mcp_overlaps_95_summary.csv"))
write.csv(trap.overlaps,paste0(out.dir,"trap_mcp_overlaps_50_summary.csv"))
write.csv(tox.overlaps,paste0(out.dir,"tox_mcp_overlaps_50_summary.csv"))
write.csv(aer.overlaps,paste0(out.dir,"aer_mcp_overlaps_50_summary.csv"))





