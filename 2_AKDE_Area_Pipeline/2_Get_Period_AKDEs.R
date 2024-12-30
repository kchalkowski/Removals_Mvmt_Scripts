#
library(amt)
library(sf)
library(stringr)
library(dplyr)
library(ctmm)
library(beepr)

####save objects needed to run on linux workstation
#ws.folder="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Workstation_Run"
#saveRDS(fp.chulls,paste0(ws.folder,"/fp.chulls"))

#read objects from input
input="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/1_Data/Input/"
sitelocs=readRDS(paste0(input,"trapntox_sites.raw.rds")) #trap/tox site locs
fp.chulls=readRDS(paste0(input,"fp.chulls.rds")) #flight path convex hulls
activities.raw=readRDS(paste0(input,"activities.rds"))

#read geolocation data
home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/1_Data/"
geo=read.csv(paste0(home,"geo_remtyp.csv"))
geo=geo[,-1]
geo=geo[!is.na(geo$Removal.Type),]
geo$date_only<-as.Date(geo$date_only)

#get ID refs to join later
geo.id=unique(geo[,c(1:4,14)])

#Get tox/trap chulls
tox.sites=sitelocs[sitelocs$activity=="toxic",]
trap.sites=sitelocs[sitelocs$activity=="trap",]

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

#Source needed functions
funcdir<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/2_Scripts/Functions"
func.list=list.files(funcdir,full.names=TRUE)
for(f in 1:length(func.list)){
  source(func.list[f])
}

{
  #trap dates
  trap.act=activities.raw[activities.raw$activity=="trap",]
  #trap.act$trap_datetime=Neat.Dates.POSIXct(trap.act$datetime,tz="UTC")
  trap.start.date=as.Date(min(trap.act$datetime))
  trap.end.date=as.Date(max(trap.act$datetime))
  
  #tox dates
  tox.act=activities.raw[activities.raw$activity=="toxic",]
  #tox.act$tox_datetime=Neat.Dates.POSIXct(tox.act$tox_datetime,tz="UTC")
  tox.start.date=as.Date(min(tox.act$datetime))
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

#verification
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

#View(geo.trap.sums)

######Trim based on above summaries
#remove any pigs that were tracked less than half the number of days of other pigs in any period
datestooshort.trap=geo.trap.sums[geo.trap.sums$difftime<(trap.len/2),]$animalid
geo.trap=geo.trap[!geo.trap$animalid%in%datestooshort.trap,]

#datestooshort.tox=geo.tox.sums[geo.tox.sums$difftime<(tox.len/2),]$animalid
datestooshort.tox=c("86070_H2_H2")
geo.tox=geo.tox[!geo.tox$animalid%in%datestooshort.tox,]

######

###Remove geolocations with very short minimum time intervals
#geo.rem<-geo.trap
thresh=as.difftime(5,units="mins")
#geo.rem<-geo.rem[geo.rem$animalid=="85428_B2_B2",]
#geo.rem$datetime<-Neat.Dates.POSIXct(geo.rem$datetime,tz="UTC")
#testing=Remove.Geo.Intervals(geo.rem,thresh)
#geo.rem=geo.pig.per
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

#Set AKDE functions
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
#pigID="48449_4A_4A"
#Revise AKDE function for periods instead of weeks
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
    #if(length(grep("square meters",rownames(summary(UD1_pHREML)$CI)))>0){
    #  area.est=area.est/1000000
    #}
    
    #if(length(grep("hectares",rownames(summary(UD1_pHREML)$CI)))>0){
    #  area.est=area.est/1000000
    #}
    
    #Get centerpoint
    print("getting centerpoint of home range akde")
    ctr.dim=which(UD1_pHREML$CDF==min(UD1_pHREML$CDF),arr.ind=TRUE)
    ctr.x=UD1_pHREML$r$x[ctr.dim[1]]
    ctr.y=UD1_pHREML$r$y[ctr.dim[2]]
    
    #is period during?
    #if(any(!(is.na(removal.key[removal.key$period==period.str,c(2,3,4)])))){
    if(period.str=="during"){
      #if so
      #Get sf polygon of UD area, 95%
      UDshp.95=as.sf(UD1_pHREML,level.UD=0.95)
      #ask which is not na
      #removals.apply=which(!(is.na(removal.key[removal.key$week==weeknum,])))
      
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
#enter geo.rem, other stuff
#loop through each pig ID, getting akde for each period
#geo.rem is subset made above
#periods is string like 'c("before","after")-- 
  #doing this way bc dont want during for aer

#geo.rem=geo.trap
#periods=c("before","during","after")
#outdir="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/AKDE_objects_test/"
#test=FALSE
#"48449_4A_4A"
#i=5
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

#Run looping function for each removal type

#Started trap period akde rerun at 1150AM
outdir="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/1_Data/AKDE_trap/"
out.df.trap=Loop.AKDE(geo.trap,c("before","during","after"),outdir,FALSE)
#saveRDS(out.df.trap,"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/outdf_akde_trap.rds")

#outdir="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/AKDE_tox/"
#out.df.tox=Loop.AKDE(geo.tox,c("before","during","after"),outdir,FALSE)
#saveRDS(out.df.tox,"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/outdf_akde_tox.rds")

#outdir="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/AKDE_aerial/"
#out.df.aer=Loop.AKDE(geo.aer,c("before","after"),outdir,FALSE)
#saveRDS(out.df.aer,"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/outdf_akde_aerial.rds")


