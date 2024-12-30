#The purpose of this script is to draft descriptive 
#figures/data frames for removals/mvmt manuscript. 

#Sections:
#I. Script setup
#II. Methods descriptive outputs 1
#III. Methods descriptive outputs 2 (this fig removed)
#IV. Data descriptive output figure

# Script Setup -----------------------------------------------------

#Load libraries
{
library("ggplot2")
library("gridExtra")
library("carData")
library("car")
library(factoextra)
library(ggpattern)
library(viridis)
library(ggfortify)
library(amt)
#install.packages("GGally")
library(GGally)
library(NIFApackagev1.1)
library(NIFADataPackage)
library(ggthemes)
library(sf)
library(terrainr)
#install.packages("ggspatial")  
library(ggspatial)
library(sp)
library(hrbrthemes)
library(vistime)
library(cowplot)
library(lubridate)
library(tidyr)
}

#Load and format spatial overview objects
{
ls=readRDS("/Volumes/Projects/MUDD/ASF_NIFA/Shapefiles/R_Spatial_Files/landsat_agg_10_df.rds")
cdo=loadNIFAspatial("cd_outline")
ap=loadNIFAspatial("aerial_pastures")
ctp=loadNIFAspatial("ctrl_pastures")
tp=loadNIFAspatial("trap_pastures")
xp=loadNIFAspatial("eup_pastures")
t1=loadNIFAspatial("trap.z1")
t2=loadNIFAspatial("trap.z2")
t3=loadNIFAspatial("trap.z3")
a1=loadNIFAspatial("aerial.z1")
a2=loadNIFAspatial("aerial.z2")
a3=loadNIFAspatial("aerial.z3")
afp.path="/Volumes/Projects/MUDD/ASF_NIFA/Datasets/Aerial_Flight_Paths/tidy_flightpaths/aerial_flightpaths.rds"
afp=readRDS(afp.path)

#need add this to main flight path processing....
#binds together and turns into line
for(i in 1:length(afp)){
  afp.i=afp[[i]] %>%
    dplyr::summarize(do_union=FALSE) %>%  # do_union=FALSE doesn't work as well
    st_cast("LINESTRING") 
  if(i==1){
    afp.ls=afp.i
  } else{
    afp.ls=rbind(afp.ls,afp.i)
  }
}

#flight path convex hulls
fp.chulls=readRDS("/Volumes/Projects/MUDD/ASF_NIFA/Datasets/Aerial_Flight_Paths/tidy_flightpaths/aerial_flightpaths_chull.rds")
st_area(fp.chulls)/1000000 #area of fp chulls

#load site data
trapntox_sites.raw=readRDS("/Volumes/Projects/MUDD/ASF_NIFA/Datasets/Daily_Activities_Trap_Locations/Activities_Tidying_Pipeline_KK/Data/Final/site.locations_tidied.RDS")

tox.sites=trapntox_sites.raw[trapntox_sites.raw$activity=="toxic"&!is.na(trapntox_sites.raw$trap_status),]
trap.sites=trapntox_sites.raw[trapntox_sites.raw$activity=="trap"&!is.na(trapntox_sites.raw$trap_status),]

tox.sites=tox.sites[!is.na(tox.sites$x),]
trap.sites=trap.sites[!is.na(trap.sites$x),]
nrow(tox.sites) #27 tox sites total
nrow(trap.sites) #78 trap sites total

nrow(tox.sites[tox.sites$trap_status=="never_set",]) #16 prebait only
nrow(tox.sites[tox.sites$trap_status=="activated",]) #11 toxicant bait sites

nrow(trap.sites[trap.sites$trap_status=="never_set",]) #26 prebait only
nrow(trap.sites[trap.sites$trap_status=="activated",]) #52 toxicant bait sites

trapsf=st_as_sf(trap.sites,coords=c("x","y"),crs=st_crs(32614))
toxsf=st_as_sf(tox.sites,coords=c("x","y"),crs=st_crs(32614))

trap.chull = st_convex_hull(st_union(trapsf))
tox.chull = st_convex_hull(st_union(toxsf))

trap.chull = trap.chull %>% st_sf %>% st_cast
tox.chull = tox.chull %>% st_sf %>% st_cast

st_area(trap.chull)/1000000
st_area(tox.chull)/1000000

}

#Load other needed objects
{
  dir<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/"
  geor=read.csv(paste0(dir,"geo_remtyp.csv"))
  trap.overlaps.95=read.csv(paste0(dir,"trap_mcp_overlaps_summary.csv"))
  tox.overlaps.95=read.csv(paste0(dir,"tox_mcp_overlaps_summary.csv"))
  aer.overlaps.95=read.csv(paste0(dir,"aer_mcp_overlaps_summary.csv"))
  
  #capture data for pigs (want capture location, X/Y, for each animalid)
  caps=read.csv("/Volumes/Projects/MUDD/ASF_NIFA/Datasets/Capture_Data/NIFA_EUP_capture_combined.csv")
  caps=caps[,c(4,5,8,13,14,33,34)]
  caps.sf=st_as_sf(caps,coords=c(7,6),crs=st_crs(4326))
  caps.sf=st_transform(caps.sf,crs=st_crs(32614))
  
}

#Get start/end dates for removal activities
#This code grabbed from Make_Removal_Designations.R, 
#Removal start/end dates here also verified with N. Snow
{
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
}


# Methods Descriptive Outputs ----------------------------------

#Needed outputs for methods descriptive figures/data frames:
#1. Pig geolocation data summary
    #a. Total number
    #b. X pigs not included in analysis because X
    #c. Average total number of fixes and fix rates (collar download/sat)
#2. Maps
    #a. Collared pig capture locs with removal chulls
    #b. Close-up of each zone with removal sites
#3. Timeline
#4. Combined figure

# Pig geolocation data summary
{
#Remove pigs not included in study
geor=geor[!is.na(geor$Removal.Type),]
geor$date_only<-Neat.Dates.POSIXct(geor$date_only,"UTC")
geor$datetime<-Neat.Dates.POSIXct(geor$datetime,"UTC")
geor.sums=geor %>% group_by(animalid) %>% dplyr::summarise(sex=first(sex),remtype=first(Removal.Type),startdt=min(date_only),enddt=max(date_only),numfix=n()) %>% as.data.frame()
#write.csv(geor.sums,"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Descriptive_Outputs/Supplementary_Data/pig_geoloc_fixes_summary.csv")
#Get average fix rate

IDs=unique(geor$animalid)
for(i in 1:length(IDs)){
geor.id=mk_track(geor[geor$animalid==IDs[i],],.x="X",.y="Y",.t="datetime",crs=st_crs(32614))
ssr.i=amt::summarize_sampling_rate(geor.id,time_unit="min") %>% as.data.frame()
ssr.i$animalid=IDs[i]
if(i==1){
  ssr=ssr.i
} else{
  ssr=rbind(ssr,ssr.i)
}
}

geo.fixsums=left_join(geor.sums,ssr,by="animalid")
geo.fixsums=geo.fixsums[,c(1:6,9)]
geo.fixsums$numdays=as.numeric(difftime(geo.fixsums$enddt,geo.fixsums$startdt))

colnames(geo.fixsums)[7]<-"medianfixrate"

write.csv(geo.fixsums,paste0(dir,"geo.fixsummaries.csv"))

length(IDs)
min(as.numeric(geo.fixsums$numdays)) #54 days
max(as.numeric(geo.fixsums$numdays)) #243 days

colnames(geor)

geor.datafrom=unique(geor[,c(5,6)])
geor.sums.datafrom=left_join(geor.sums,unique(geor[,c(5,6)]),by="animalid")

#Fix rates were 15-20 for collars that were retrieved (Supp Table X). 
#For unretrieved collars, median fix rates for available data was 240 minutes (Supp. Table X)
#Duration of geolocation data (after filtering), for each pig, ranges from 54-243 days. (Supp Table X)

#Notes from Make_Removal_Designations:
#Summary:
#122 pigs collared, removed one because lost collar very early
#121 pigs with geolocation data
#106 collared pigs included in analysis
#15 pigs not used for any treatment or ctrls
#Of these 15, 
#9 of these pigs were not collared during any of the treatment dates
#2 pigs overlapped with aerial but collar data ends before aerial treatment dates
#4 pigs did not have 50% core MCP home ranges that overlapped with any treatment area
}

#### Maps maps maps: overview, each zone
{
#Study design map figures
  #Map with flight paths, chull of tox/trap sites, pts for tox/trap sites
  #Locations of where pigs were collared
  #Timeline, before/during/after for each treatment zone

base=ggplot() +
  terrainr::geom_spatial_rgb(
    data = ls,
    mapping = aes(
      x = x,
      y = y,
      r = Red,
      g = Green,
      b = Blue
    ),alpha=0.7)+coord_fixed()#+#theme_map()+
  #ggspatial::annotation_north_arrow(which_north = "true",location = "tr",
  #                       style = north_arrow_fancy_orienteering(text_size=30), 
  #                       height=unit(3, "cm"),
  #                       width=unit(3,"cm"),
  #                       pad_x = unit(3, "cm"),
  #                       pad_y = unit(3, "cm"))

#Add study area polygons
ctr.locs.df=readRDS("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Output/AllCtrPts/ctrpts.rds")
ctrs=st_as_sf(ctr.locs.df,coords=c("X","Y"),crs=st_crs(32614))

base2=base+
  geom_sf(data=cdo$cd_outline,alpha=0.3,fill="#FFFFFF",linewidth=0)+
  geom_sf(data=ap$aerial_pastures,alpha=0,color="#822AFF",linewidth=2,show.legend = TRUE)+
  geom_sf(data=tp$trap_pastures,alpha=0,color="#FF57A9",linewidth=2,show.legend = TRUE)+
  geom_sf(data=xp$eup_pastures,alpha=0,color="#FF985A",linewidth=2,show.legend = TRUE)+
  geom_sf(data=ctp$ctrl_pastures,alpha=0,color="#FFD065",linewidth=2,show.legend = TRUE)
  
map.overview=base2+
  geom_sf(data=tox.chull,fill="#FF985A",color="#FF985A",linewidth=1,alpha=0.3,show.legend = TRUE)+
  geom_sf(data=trap.chull,fill="#FF57A9",color="#FF57A9",linewidth=1,alpha=0.3,show.legend = TRUE)+
  geom_sf(data=fp.chulls,fill="#822AFF",color="#822AFF",linewidth=1,alpha=0.3,show.legend = TRUE)+
  geom_sf(data=ctrs,color="#FFFFFF",fill="#000000",shape=22,size=6,show.legend = TRUE)+ggthemes::theme_map()
#map.overview

#make ls raster
lsr=raster::rasterFromXYZ(ls,crs="EPSG:32614")

ap.b500=st_buffer(ap$aerial_pastures,500)
tp.b500=st_buffer(tp$trap_pastures,500)
xp.b500=st_buffer(xp$eup_pastures,500)

lsr.ap=raster::mask(lsr,ap.b500)
lsr.tp=raster::mask(lsr,tp.b500)
lsr.xp=raster::mask(lsr,xp.b500)

rasptap <- raster::rasterToPoints(lsr.ap)
lsr2ap=as.data.frame(rasptap)

raspttp <- raster::rasterToPoints(lsr.tp)
lsr2tp=as.data.frame(raspttp)

rasptxp <- raster::rasterToPoints(lsr.xp)
lsr2xp=as.data.frame(rasptxp)

base.ap=ggplot() +
  terrainr::geom_spatial_rgb(
    data = lsr2ap,
    mapping = aes(
      x = x,
      y = y,
      r = Red,
      g = Green,
      b = Blue
    ),alpha=0.7)+coord_fixed()

base.tp=ggplot() +
  terrainr::geom_spatial_rgb(
    data = lsr2tp,
    mapping = aes(
      x = x,
      y = y,
      r = Red,
      g = Green,
      b = Blue
    ),alpha=0.7)+coord_fixed()

base.xp=ggplot() +
  terrainr::geom_spatial_rgb(
    data = lsr2xp,
    mapping = aes(
      x = x,
      y = y,
      r = Red,
      g = Green,
      b = Blue
    ),alpha=0.7)+coord_fixed()

#Close-up of tox zone
tox.map=base.xp+
  geom_sf(data=xp$eup_pastures,alpha=0.2,fill="#FFFFFF",color="#FF985A",linewidth=2,show.legend = TRUE)+
  geom_sf(data=toxsf[toxsf$trap_status=="never_set",],shape=21,fill="#FFFFFF",color="#000000",size=7)+
  geom_sf(data=toxsf[toxsf$trap_status=="activated",],color="#000000",size=7)+theme_map()

#Close-up of trap zone
trap.map=base.tp+
  geom_sf(data=tp$trap_pastures,alpha=0.2,fill="#FFFFFF",color="#FF57A9",linewidth=2,show.legend = TRUE)+
  geom_sf(data=trapsf[trapsf$trap_status=="never_set",],shape=21,fill="#FFFFFF",color="#000000",size=7)+
  geom_sf(data=trapsf[trapsf$trap_status=="activated",],color="#000000",size=7)+theme_map()

#Close-up of aerial zone
aer.map=base.ap+
  geom_sf(data=ap$aerial_pastures,alpha=0.2,fill="#FFFFFF",color="#822AFF",linewidth=2,show.legend = TRUE)+
  #geom_sf(data=fp,shape=1,color="#000000",size=4)
  geom_sf(data=afp.ls,alpha=1,linewidth=1,color="black",size=0.1)+theme_map()
  
}

# Timeline of removal activities
{
#source NRM_Descriptions, same folder as this script
#object name: timeline
source("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/2_Scripts/SummarizingScripts/NRM_Descriptions.R")  
}

# Combined method description figure
{
#objects:
#map.overview, tox.map, trap.map, aer.map, timeline
#tox.spaces=plot_grid(NULL,tox.map,NULL,ncol=1,rel_heights=c(0.0,0.7,0.3))
#tox.overview=plot_grid(NULL,map.overview,NULL,tox.spaces,nrow=1,rel_widths=c(-0.4,1.2,-0.5,0.3))
#trap.spaces=plot_grid(NULL,trap.map,ncol=1,rel_heights=c(0.01,0.99))
#aer.trap=plot_grid(NULL,aer.map,NULL,trap.spaces,NULL,nrow=1,rel_widths=c(-0.28,0.7,-0.30,0.13,0.1))
#all_maps=plot_grid(tox.overview,NULL,aer.trap,ncol=1,rel_heights=c(0.5,-0.08,0.5))
#test=plot_grid(all_maps,timeline,ncol=1,rel_heights=c(0.85,0.15))
#ggsave(plot=test,"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Descriptive_Outputs/Methods_Descriptive_Fig/Methods_Description.png",width=16.67,height=19.44,units="in")

aer.spaces=plot_grid(NULL,aer.map,
                     ncol=1,
                     rel_heights=c(-0.15,1.15))
aer.trap.tox=plot_grid(NULL,aer.spaces,NULL,trap.map,NULL,tox.map,
                       nrow=1,
                       rel_widths=c(-0.1,0.5,-0.18,0.5,-0.12,0.5))
allmap.spaces=plot_grid(map.overview,NULL,NULL,
                        nrow=1,
                        rel_widths=c(0.6,0.2,0.2))
all_maps=plot_grid(allmap.spaces,
                   NULL,
                   aer.trap.tox,
                   ncol=1,
                   rel_heights=c(0.7,-0.15,0.3))

maps.tl=plot_grid(all_maps,
                  timeline,
                  ncol=1,
                  rel_heights=c(0.8,0.2))

ggsave("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Descriptive_Outputs/Methods_Descriptive_Fig/Methods_Description.png",
       plot=maps.tl,
       bg="white",
       height=14,
       width=12,
       units="in")

}


# Results Descriptive Outputs ------------------------------------
#1. Removal designations, bar chart
    #a. Number in each treatment, number in control
    #b. Average overlaps, prop and km2
#2. Collared pigs by treatment, M/F
#3. Spatial overlaps of removals on collared pig 50% MCPs, M/F
    #proportion of range overlap
    #area (km2) of range overlap
#4. Temporal overlap of treatment
    #ridgeline plot, one histo for each treatment area
    #y axis is how many weeks of overlap for each pig
#5. Combined figures

# Removal numbers barchart
{
  #removed per site
  #trapping: daily activities estimate: most recent update from Kelly, 296
  #aerial: estimate from operators/nate in heli kills, 256
  #toxicant: estimate from Nate, 58
  trap.removals=296
  aer.removals=256
  tox.removals=58
  
  rems=data.frame(name=c("aer","tox","trap"),removals=c(aer.removals,tox.removals,trap.removals))
  
  set.text.size=10
  set.title.size=12
  
  rembars=ggplot() +
    geom_bar(rems, mapping=aes(x=name, y=removals, fill=name),stat="identity")+
    scale_fill_manual(values=c("#822AFF","#FF985A","#FF57A9"))+
    geom_bar() + 
    theme_ipsum()+
    theme(axis.text.x = element_text(size = set.text.size,angle=90),
          axis.text.y = element_text(size = set.text.size),
          axis.title.x = element_text(size = set.title.size),
          axis.title.y = element_text(size = set.title.size),
          legend.text = element_text(size = set.text.size),
          legend.title = element_text(size = set.title.size),
          legend.position = "none")+
    labs(fill="Method",x="Method",y="Num. Removals")
}

# Collared pig Treatment/control barchart by designation, M/F
{
  geor.sums$sex[geor.sums$sex=="Male"]<-"male"
  geor.sums$sex[geor.sums$sex=="Female"]<-"female"
  
  tc.sums=geor.sums %>% group_by(remtype,sex) %>% dplyr::summarise(num=n()) %>% as.data.frame() 
  set.text.size=10
  set.title.size=12
  tc.bars=ggplot() +
    geom_bar(tc.sums, mapping=aes(x=remtype, y=num, color=remtype,fill=sex),stat="identity",linewidth=2,position=position_dodge())+
    scale_color_manual(values=c("#822AFF","#FF985A","#FF57A9","#FFD065"))+
    scale_fill_manual(values=c("#c9c9c9","#4a4a4a"))+
    #scale_pattern_manual(values=c("stripe","circle"))+
    #scale_color_manual(values=c("#822AFF","#FF985A","#FF57A9","#FFD065"))+
    theme_ipsum()+
    theme(axis.text.x = element_text(size = set.text.size,angle=90),
          axis.text.y = element_text(size = set.text.size),
          axis.title.x = element_text(size = set.title.size),
          axis.title.y = element_text(size = set.title.size),
          legend.position = "none"
    )+
    labs(x="Method",y="Number collared pigs")
}

# Spatial overlaps of removals on collared pig 50% MCPs, M/F
{
#ridgeline plot, one histo for each treatment area
#one for proportion of range overlap and one for km2
#trap.overlaps,tox.overlaps,aer.overlaps
trap.overlaps$trt="trap"
tox.overlaps$trt="tox"
aer.overlaps$trt="aer"

geor.sums.trap=geor.sums[geor.sums$remtype=="trap",]
geor.sums.tox=geor.sums[geor.sums$remtype=="tox",]
geor.sums.aer=geor.sums[geor.sums$remtype=="aer",]

geor.sums.trap=left_join(geor.sums.trap,trap.overlaps,by="animalid")
geor.sums.tox=left_join(geor.sums.tox,tox.overlaps,by="animalid")
geor.sums.aer=left_join(geor.sums.aer,aer.overlaps,by="animalid")

geor.overlaps=rbind(geor.sums.trap,geor.sums.tox,geor.sums.aer)

ggplot()+
  geom_density_ridges2(geor.overlaps,
                      mapping=aes(x=prop_overlap,y=trt,fill=sex),
                      alpha=0.5)+
  theme_ipsum()

geor.overlaps$km_overlap=geor.overlaps$area_overlap/1e6

prop.over.ridge=ggplot()+
  geom_density_ridges(geor.overlaps,
                       mapping=aes(x=prop_overlap,y=trt,fill=sex,color=trt),
                       alpha=0.5,linewidth=1)+
  scale_color_manual(values=c("#822AFF","#FF985A","#FF57A9","#FFD065"))+
  scale_fill_manual(values=c("#c9c9c9","#4a4a4a"))+
  theme_ipsum()+
  theme(legend.position = "none")

km.over.ridge=ggplot()+
  geom_density_ridges(geor.overlaps,
                      mapping=aes(x=km_overlap,y=trt,fill=sex,color=trt),
                      alpha=0.5,linewidth=1)+
  scale_color_manual(values=c("#822AFF","#FF985A","#FF57A9","#FFD065"))+
  scale_fill_manual(values=c("#c9c9c9","#4a4a4a"))+
  theme_ipsum()+
  theme(legend.position = "none")
}

# Temporal overlap of treatment
{
geor.dumbell=geor.sums[,1:5]
geor.dumbell=pivot_longer(geor.dumbell,cols=c(4,5)) %>% as.data.frame()

#Need get week form of trap start/end dates
weekdtkey=unique(geor[,c(15,17)])
colnames(geor.dumbell)[5]<-"date_only"
geor.dumbell=left_join(geor.dumbell,weekdtkey,by="date_only")

#Need get removal start/end dates converted to week numbers
trap.start.wk=weekdtkey[weekdtkey$date_only==trap.start.date,]$week
trap.end.wk=weekdtkey[weekdtkey$date_only==trap.end.date,]$week
aer.start.wk=weekdtkey[weekdtkey$date_only==aer.start.date,]$week
aer.end.wk=weekdtkey[weekdtkey$date_only==aer.end.date,]$week
tox.start.wk=weekdtkey[weekdtkey$date_only==tox.start.date,]$week
tox.end.wk=weekdtkey[weekdtkey$date_only==tox.end.date,]$week

Starts <- geor.dumbell %>%
  filter(name == "startdt")
Ends <- geor.dumbell %>%
  filter(name == "enddt")
geor.dumbell$animalid<-as.factor(geor.dumbell$animalid)


trap.db=geor.dumbell[geor.dumbell$remtype=="trap",]
tox.db=geor.dumbell[geor.dumbell$remtype=="tox",]
aer.db=geor.dumbell[geor.dumbell$remtype=="aer",]
ctrl.db=geor.dumbell[geor.dumbell$remtype=="ctrl",]


p.trap.dumbell <- ggplot(trap.db)+
  geom_segment(data = Starts[Starts$remtype=="trap",],
               aes(x = week, y = animalid,
                   yend = Ends[Ends$remtype=="trap",]$animalid, xend = Ends[Ends$remtype=="trap",]$week), #use the $ operator to fetch data from our "Females" tibble
               color = "#FF57A9",
               size = 4.5, #Note that I sized the segment to fit the points
               alpha = .25) +
  geom_point(aes(x = week, y = animalid, color = name), size = 4, show.legend = TRUE)+
  scale_color_manual(values=c("#000000","#000000"))+
  theme_ipsum()+
  theme(axis.text.y=element_blank(),
        legend.position = "none")+
  geom_vline(xintercept=trap.start.wk,linetype="dashed",color="#FF57A9", linewidth=2)+
  geom_vline(xintercept=trap.end.wk,linetype="dashed",color="#FF57A9", linewidth=2)

p.tox.dumbell <- ggplot(tox.db)+
  geom_segment(data = Starts[Starts$remtype=="tox",],
               aes(x = week, y = animalid,
                   yend = Ends[Ends$remtype=="tox",]$animalid, xend = Ends[Ends$remtype=="tox",]$week), #use the $ operator to fetch data from our "Females" tibble
               color = "#FF985A",
               size = 4.5, #Note that I sized the segment to fit the points
               alpha = .5) +
  geom_point(aes(x = week, y = animalid, color = name), size = 4, show.legend = TRUE)+
  scale_color_manual(values=c("#000000","#000000"))+
  theme_ipsum()+
  #theme(axis.text.y=element_blank(),
  #      legend.position = "none")+
  theme(
        legend.position = "none")+
  geom_vline(xintercept=tox.start.wk,linetype="dashed",color="#FF985A", linewidth=2)+
  geom_vline(xintercept=tox.end.wk,linetype="dashed",color="#FF985A", linewidth=2)

p.aer.dumbell <- ggplot(aer.db)+
  geom_segment(data = Starts[Starts$remtype=="aer",],
               aes(x = week, y = animalid,
                   yend = Ends[Ends$remtype=="aer",]$animalid, xend = Ends[Ends$remtype=="aer",]$week), #use the $ operator to fetch data from our "Females" tibble
               color = "#822AFF",
               size = 4.5, #Note that I sized the segment to fit the points
               alpha = .5) +
  geom_point(aes(x = week, y = animalid, color = name), size = 4, show.legend = TRUE)+
  scale_color_manual(values=c("#000000","#000000"))+
  theme_ipsum()+
  theme(axis.text.y=element_blank(),
        legend.position = "none")+
  geom_vline(xintercept=aer.start.wk,linetype="dashed",color="#822AFF", linewidth=2)+
  geom_vline(xintercept=aer.end.wk,linetype="dashed",color="#822AFF", linewidth=2)

p.ctrl.dumbell <- ggplot(ctrl.db)+
  geom_segment(data = Starts[Starts$remtype=="ctrl",],
               aes(x = week, y = animalid,
                   yend = Ends[Ends$remtype=="ctrl",]$animalid, xend = Ends[Ends$remtype=="ctrl",]$week), #use the $ operator to fetch data from our "Females" tibble
               color = "#FFD065",
               size = 4.5, #Note that I sized the segment to fit the points
               alpha = 1) +
  geom_point(aes(x = week, y = animalid, color = name), size = 4, show.legend = TRUE)+
  scale_color_manual(values=c("#000000","#000000"))+
  theme_ipsum()+
  theme(axis.text.y=element_blank(),
        legend.position = "none")+
  geom_vline(xintercept=aer.start.wk,linetype="dashed",color="#822AFF", linewidth=2)+
  geom_vline(xintercept=aer.end.wk,linetype="dashed",color="#822AFF", linewidth=2)+
  geom_vline(xintercept=trap.start.wk,linetype="dashed",color="#FF57A9", linewidth=2)+
  geom_vline(xintercept=trap.end.wk,linetype="dashed",color="#FF57A9", linewidth=2)+
  geom_vline(xintercept=tox.start.wk,linetype="dashed",color="#FF985A", linewidth=2)+
  geom_vline(xintercept=tox.end.wk,linetype="dashed",color="#FF985A", linewidth=2)
}

# Combine results figures
{
#Stitch results figure.
#Leaving out timelines for separate figure, maybe, since already a lot here
#object names:
#rembars, tc.bars, prop.over.ridge, km.over.ridge
#p.trap.dumbell, p.tox.dumbell, p.aer.dumbell, p.ctrl.dumbell

bars=plot_grid(rembars,tc.bars,nrow=1,rel_widths=c(0.4,0.6),labels=c("A","B"))
overlaps=plot_grid(prop.over.ridge,km.over.ridge,ncol=1,labels=c("C","D"))
dumbells=plot_grid(p.trap.dumbell,p.tox.dumbell,p.aer.dumbell,p.ctrl.dumbell,labels=c("A","B","C","D"))
bar_overlaps=plot_grid(bars,overlaps,ncol=1,rel_heights=c(0.3,0.7))

#save bar_overlaps and dumbells
}

# Data descriptive output -----------------------------------------------------

#Timeline figures for each response
  #pull in outdf corrected for each treatment
  #find scripts with timeline figures....
outaer=readRDS("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/1_Data/Objects/outdf_akde_aer_corrected.rds")



#Scrap code
{
p=ggplot(data=geo.summaries[geo.summaries$date>=x.start&geo.summaries$date<=x.end,], aes(x=date, y=medx))+
  geom_line(colour="grey",alpha=0.8)+
  geom_ribbon(aes(ymin=q.25x, ymax=q.75x), linetype=2, alpha=0.1)+
  scale_x_date(labels = date_format("%Y-%m-%d"),
               breaks = seq(x.start, x.end+intervaln, intervaln),
               limits = c(x.start, x.end+intervaln))+
  geom_vline(xintercept = startdate,linetype=2)+
  geom_vline(xintercept = enddate, linetype=2)+
  geom_segment(data=geo.meds.e,
               mapping=aes(x=xs,xend=xe,y=median.entire,yend=median.entire),colour="red")+
theme_ipsum()+
  theme(axis.text.x = element_text(angle = 90))+
  labs(y=y.lab)+ylim(0,y.max)+
  facet_wrap(~Removal.Type,ncol=1)+
  theme(strip.text.x = element_text(size = 25, colour = "grey", face="bold"),
        axis.title.y = element_text(size = 15))
}









