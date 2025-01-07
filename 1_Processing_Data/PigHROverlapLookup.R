#set home dir of pipeline
home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline"

# Purpose ----------------------------------------------------------------------
#The purpose of this script is to create a df of pig pairs with hr overlap status

# In/Out -----------------------------------------------------------------------

#inputs:geo_remtype.rds

#outputs: pairs.rds

# Setup ------------------------------------------------------------------------

#load libraries
library(glmmTMB)
library(DHARMa)
library(ctmm)
library(fitdistrplus)
library(ggplot2)
library(tidyr)
library(plyr)
library(dplyr)
library(hrbrthemes)
library(sf)
library(purrr)
library(mapview)
library(amt)

#Set dirs
input=file.path(home,"1_Data","Input",fsep=.Platform$file.sep)
objdir=file.path(home,"1_Data","Objects",fsep=.Platform$file.sep)
if(!dir.exists(file.path(home,"3_Output","Area_GLM_Results",fsep=.Platform$file.sep))){
  dir.create(file.path(home,"3_Output","Area_GLM_Results",fsep=.Platform$file.sep))
}
outdir=file.path(home,"3_Output","Area_GLM_Results",fsep=.Platform$file.sep)

#Read input objects
geo=readRDS(file.path(objdir,"geo_remtype.rds"))

#source needed functions
funcdir<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/2_Scripts/Functions"
func.list=list.files(funcdir,full.names=TRUE)
func.list=func.list[-grep("cpp",func.list)]
for(f in 1:length(func.list)){
  source(func.list[f])
}


# Pull IDs by type --------------------------------------------------------------------
trapctrl=unique(geo[geo$Removal.Type=="trap"|geo$Removal.Type=="ctrl",c("animalid","Removal.Type")])
aerctrl=unique(geo[geo$Removal.Type=="aer"|geo$Removal.Type=="ctrl",c("animalid","Removal.Type")])
toxctrl=unique(geo[geo$Removal.Type=="tox"|geo$Removal.Type=="ctrl",c("animalid","Removal.Type")])

# Set removal dates --------------------------------------------------------------------

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

# Format geolocations --------------------------------------------------------------------

geo$jDate <- julian(as.Date(geo$datetime), origin = min(as.Date(geo$datetime)))
geo$jDate=geo$jDate+1

# Keep all IDs, but trim and divide datasets by cutoffs
tkmt=amt::mk_track(geo,.x=X,.y=Y,.t=datetime,all_cols=TRUE,crs=32614)

jdates=geo %>% group_by(jDate) %>% summarise(dates=min(date_only)) %>% as.data.frame()


#Next by animalid
tkmt.n=tkmt |> nest(data=-animalid)

#filter tkmt tracks by each cutoff interval
tkmt.n=tkmt.n |>
  mutate(trap.trks = map(data, function(x)
    x |>
      filter(jDate>=which(jdates$dates==trap.cutoff1)&jDate<=which(jdates$dates==trap.start.date))
  ))

tkmt.n=tkmt.n |>
  mutate(aer.trks = map(data, function(x)
    x |>
      filter(jDate>=which(jdates$dates==aer.cutoff1)&
               jDate<=which(jdates$dates==aer.start.date))
  ))

tkmt.n=tkmt.n |>
  mutate(tox.trks = map(data, function(x)
    x |>
      filter(jDate>=which(jdates$dates==tox.cutoff1)&
               jDate<=which(jdates$dates==tox.start.date))
  ))

# Make MCPs --------------------------------------------------------------------

#Get MCPs for each set of tracks, cutoff by trtmnt area dates
mcp.iso=0.95
trap.mcps=GetAreaMetrics(tkmt.n$trap.trks,"mcp",mcp.iso,tkmt.n$animalid)
aer.mcps=GetAreaMetrics(tkmt.n$aer.trks,"mcp",mcp.iso,tkmt.n$animalid)
tox.mcps=GetAreaMetrics(tkmt.n$tox.trks,"mcp",mcp.iso,tkmt.n$animalid)

# Identify MCP overlaps --------------------------------------------------------------------
#mcp_sf=trap.mcps$mcp
#idrem=trapctrl
#prop=0.1 #minimum proportion of hr overlap allowable, 0.1=10% overlap
Pig_Overlaps=function(mcp_sf,idrem,prop){
ids=idrem$animalid
  
intxn=sf::st_intersection(mcp_sf,mcp_sf) %>% 
  dplyr::mutate(
    area       = sf::st_area(.),
    proportion = area / area.1
  ) %>%
  tibble::as_tibble() %>%
  dplyr::select(
    id_1 = animalid,
    id_2 = animalid.1,
    proportion,
    area
  ) %>% as.data.frame()


intxn=intxn[intxn$id_1%in%ids,]
intxn=intxn[intxn$id_2%in%ids,]

intxn$proportion=as.numeric(intxn$proportion)
intxn=intxn[intxn$proportion>=prop,]
intxn=intxn[intxn$id_1!=intxn$id_2,]

intxn$trt_ctrl_1=NA
intxn$trt_ctrl_2=NA
for(i in 1:nrow(intxn)){
  intxn$trt_ctrl_1[i]=idrem[grep(intxn[i,"id_1"],ids),"Removal.Type"]
  intxn$trt_ctrl_2[i]=idrem[grep(intxn[i,"id_2"],ids),"Removal.Type"]
}

intxn$trt_ctrl_1[intxn$trt_ctrl_1!="ctrl"]<-"trt"
intxn$trt_ctrl_2[intxn$trt_ctrl_2!="ctrl"]<-"trt"

#remove intxn crossovers
intxn=intxn[intxn$trt_ctrl_1==intxn$trt_ctrl_2,]

intxn$rem=unique(idrem$Removal.Type[idrem$Removal.Type!="ctrl"])

intxn=intxn[,c(7,1:6)]

return(intxn)
}

#Run pig overlap functions, get pair data frames
trap_pairs=Pig_Overlaps(trap.mcps$mcp,trapctrl,0.01)
aer_pairs=Pig_Overlaps(aer.mcps$mcp,aerctrl,0.01)
tox_pairs=Pig_Overlaps(tox.mcps$mcp,toxctrl,0.01)

#bind together
pairs=dplyr::bind_rows(trap_pairs,aer_pairs,tox_pairs)

#save output
saveRDS(pairs,file.path(objdir,"pairs.rds"))


