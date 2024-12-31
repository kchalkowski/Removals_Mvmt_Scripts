
#set home dir of pipeline
home<-"/Users/kayleigh.chalkowski/OneDrive - USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline"

# Purpose ----------------------------------------------------------------------

#The purpose of this script is to output AKDE objects by period/remoavl type

# Process ----------------------------------------------------------------------

#In: geotox.rds, geoaer.rds, geotrap.rds 

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
library(ctmm)
library(beepr)

#set directories
input=file.path(home,"1_Data","Input",fsep=.Platform$file.sep)
objdir=file.path(home,"1_Data","Objects",fsep=.Platform$file.sep)
funcdir<-file.path(home,"2_Scripts","Functions",fsep=.Platform$file.sep)

#Read geoloc data
geo.tox<-readRDS(file.path(objdir,"geotox.rds"))
geo.aer<-readRDS(file.path(objdir,"geoaer.rds"))
geo.trap<-readRDS(file.path(objdir,"geotrap.rds"))

# Format data for period-level analysis ----------------------------------------

# * Get summaries for each pig/period ------------------------------------------

geo.tox.sums=
  geo.tox %>% 
  dplyr::group_by(animalid,removal.period.akdecalc) %>% 
  dplyr::summarise(num.locs=n(),
                   data_from=first(data_from),
                   strt.date=min(date_only),
                   end.date=max(date_only))
geo.tox.sums$difftime=as.Date(geo.tox.sums$end.date)-as.Date(geo.tox.sums$strt.date)

geo.trap.sums=
  geo.trap %>% 
  dplyr::group_by(animalid,removal.period.akdecalc) %>% 
  dplyr::summarise(num.locs=n(),
                   data_from=first(data_from),
                   strt.date=min(date_only),
                   end.date=max(date_only))
geo.trap.sums$difftime=as.Date(geo.trap.sums$end.date)-as.Date(geo.trap.sums$strt.date)

geo.aer.sums=
  geo.aer %>% 
  dplyr::group_by(animalid,removal.period.akdecalc) %>% 
  dplyr::summarise(num.locs=n(),
                   data_from=first(data_from),
                   strt.date=min(date_only),
                   end.date=max(date_only))
geo.aer.sums$difftime=as.Date(geo.aer.sums$end.date)-as.Date(geo.aer.sums$strt.date)

# * Trim pigs with incomplete periods ------------------------------------------

# Trim based on above summaries
#remove any pigs that were tracked less than half the number of days of other pigs in any period
datestooshort.trap=geo.trap.sums[geo.trap.sums$difftime<(trap.len/2),]$animalid
geo.trap=geo.trap[!geo.trap$animalid%in%datestooshort.trap,]

# Remove geolocations with very short minimum time intervals
thresh=as.difftime(5,units="mins")
Remove.Geo.Intervals<-function(geo.rem,thresh){
  geo.rem$dtlag=dplyr::lag(geo.rem$datetime)
  geo.rem$difft=geo.rem$datetime-geo.rem$dtlag
  if(!is.na(any(geo.rem$difft<=thresh))){
    geo.out=geo.rem[-which(geo.rem$difft<=thresh),]
    geo.out=geo.out[,-c(which(colnames(geo.out)=="dtlag"),which(colnames(geo.out)=="difft"))]
  } else{
    geo.out=geo.rem
  }
  return(geo.out)
}

# Make functions to model CTMMS ----------------------------------------------

#Function to convert geolocs to telemetry format
Convert.Telemetry<-function(geolocs){
  id.n=which(colnames(geolocs)=="animalid")
  id.t=which(colnames(geolocs)=="datetime")
  id.lon=which(colnames(geolocs)=="latitude")
  id.lat=which(colnames(geolocs)=="longitude")
  tk=geolocs[,c(id.n,id.t,id.lon,id.lat)]
  colnames(tk)<-c("ID","timestamp","longitude","latitude")
  crs_str="+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  te1 <- ctmm::as.telemetry(tk,projection=crs_str)
  return(te1)
}

#AKDE function for periods
GetAKDE_AC<-function(te1,akde.outdir,pigID,period.str,removal.str){
  out.df=data.frame(matrix(nrow=1,ncol=12))
  colnames(out.df)=c("animalid","period",
                     "area.CI.low","area.est","area.CI.high",
                     "ctr.x","ctr.y",
                     "rem.overlaps.CI.low","rem.overlaps.est","rem.overlaps.CI.high",
                     "Neff","akde.hrarea.units")

  out.df$animalid=pigID
  out.df$period=period.str

  if(nrow(te1)>1){
    #Fit ctmm model
    print("finding good start point for model fitting")
    GUESS1 <- ctmm.guess(te1, interactive = FALSE)
    
    print("fitting ctmm model")
    FIT1_pHREML <- ctmm.select(te1, GUESS1, method = 'pHREML')
    Neff=summary(FIT1_pHREML)$DOF["area"]
    if(!dir.exists(paste0(outdir))){dir.create(paste0(outdir))}
    if(!dir.exists(paste0(outdir,pigID))){dir.create(paste0(outdir,pigID))}
    saveRDS(FIT1_pHREML,paste0(outdir,pigID,"/",pigID,"_",period.str,".rds"))
    
    #Get area estimation
    print("getting akde area estimation")
    UD1_pHREML <- akde(te1, FIT1_pHREML)
    area.est=summary(UD1_pHREML)$CI
    area.units=rownames(summary(UD1_pHREML)$CI)
    
    #Get centerpoint
    print("getting centerpoint of home range akde")
    ctr.dim=which(UD1_pHREML$CDF==min(UD1_pHREML$CDF),arr.ind=TRUE)
    ctr.x=UD1_pHREML$r$x[ctr.dim[1]]
    ctr.y=UD1_pHREML$r$y[ctr.dim[2]]
    
    #is period during?
    if(period.str=="during"){
      #if so
      #Get sf polygon of UD area, 95%
      UDshp.95=as.sf(UD1_pHREML,level.UD=0.95)
      #ask which is not na
      if(removal.str=="aer"|removal.str=="ctrl"){ #if aer.index, get aerial overlap
        #overlap with fp.chulls
        rem.overlaps.area=c(NA,NA,NA)
      }
      
      if(removal.str=="tox"){ #if tox.index, get tox overlap
        int=st_intersection(UDshp.95,tox.chull)
        #get area of overlap
        #tox.overlaps.area=st_area(int)
        
        if(nrow(int)<3){
          rem.overlaps.area=c(0,0,0)
          
          if(length(grep("low",int$name))>0){
            rem.overlaps.area[1]<-st_area(int[grep("low",int$name),])
          }
          if(length(grep("est",int$name))>0){
            rem.overlaps.area[2]<-st_area(int[grep("est",int$name),])
          }
          if(length(grep("low",int$name))>0){
            rem.overlaps.area[3]<-st_area(int[grep("high",int$name),])
          }
        } else{
          rem.overlaps.area=st_area(int)
        }
        
      }
      
      if(removal.str=="trap"){ #if trap.index, get trap overlap
        int=st_intersection(UDshp.95,trap.chull)
        #get area of overlap
        
        if(nrow(int)<3){
          rem.overlaps.area=c(0,0,0)
          if(length(grep("low",int$name))>0){
            rem.overlaps.area[1]<-st_area(int[grep("low",int$name),])
          }
          if(length(grep("est",int$name))>0){
            rem.overlaps.area[2]<-st_area(int[grep("est",int$name),])
          }
          if(length(grep("low",int$name))>0){
            rem.overlaps.area[3]<-st_area(int[grep("high",int$name),])
          }
        } else{
          rem.overlaps.area=st_area(int)
        }
        
      } 
      
      
    } else{
      rem.overlaps.area=c(NA,NA,NA)
    }
    
    out.df[,3:5]=area.est
    out.df[,6:7]=c(ctr.x,ctr.y)
    out.df[8:10]=rem.overlaps.area
    out.df[11]=Neff
    out.df[12]=area.units
    
  } 
  
  return(out.df)
}

#Loop function for each removal type
#loop through each pig ID, getting akde for each period
#and saving CTMM model object
Loop.AKDE<-function(geo.rem,periods,outdir,test){
  geo.rem$datetime<-Neat.Dates.POSIXct(geo.rem$datetime,tz="UTC")
  IDs=unique(geo.rem$animalid)
  #if(test){IDs=IDs[1:3]}
  #object 'rem.overlaps.area' not found
  for(i in 1:length(IDs)){
    print(paste0("pig ",IDs[i],": ",i," of ",length(IDs)))
    #geo.pig=geo.rem[geo.rem$animalid=="85440_E2_E2",]
    geo.pig=geo.rem[geo.rem$animalid==IDs[i],]
    removal.str=geo.pig$Removal.Type[1]
    for(p in 1:length(periods)){
      print(paste0("period ",periods[p]))
      period.str=periods[p]
      geo.pig.per=geo.pig[geo.pig$removal.period.akdecalc==periods[p],]
      if(nrow(geo.pig.per)>1){
        #geo.pig.per$datetime<-Neat.Dates.POSIXct(geo.pig.per$datetime,tz="UTC")
      #geo.pig.per=Remove.Geo.Intervals(geo.pig.per,thresh)
      te1=Convert.Telemetry(geo.pig.per)
      out.df.per=GetAKDE_AC(te1,outdir,IDs[i],periods[p],removal.str)
      } else{ #else just make empty row with animalid and period
        out.df.per=data.frame(matrix(NA,nrow=1,ncol=12))
        colnames(out.df.per)=c("animalid","period",
                           "area.CI.low","area.est","area.CI.high",
                           "ctr.x","ctr.y",
                           "rem.overlaps.CI.low","rem.overlaps.est","rem.overlaps.CI.high",
                           "Neff","akde.hrarea.units")
        
        out.df.per$animalid=IDs[i]
        out.df.per$period=period.str
      }
      
      if(i==1&p==1){out.df.total=out.df.per
      } else{
        out.df.total=rbind(out.df.total,out.df.per)
      }
    }
    
  }
  beep()
  return(out.df.total)
  
}

# Run CTMMs for pigs in each removal type --------------------------------------

#Started trap period akde rerun at 1150AM
outdir=file.path(home,"1_Data","Objects","AKDE_trap")
out.df.trap=Loop.AKDE(geo.trap,c("before","during","after"),outdir,FALSE)
saveRDS(out.df.trap,"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/outdf_akde_trap.rds")

outdir=file.path(home,"1_Data","Objects","AKDE_tox")
out.df.tox=Loop.AKDE(geo.tox,c("before","during","after"),outdir,FALSE)
saveRDS(out.df.tox,"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/outdf_akde_tox.rds")

outdir=file.path(home,"1_Data","Objects","AKDE_aer")
out.df.aer=Loop.AKDE(geo.aer,c("before","after"),outdir,FALSE)
saveRDS(out.df.aer,"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/outdf_akde_aerial.rds")
