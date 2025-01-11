### Title: Summarizing weekly individual speed and distance moved from ctmm models
### Author: Abbey Feuka
### Date: 26NOV24
### uses ctmms fit by Kayleigh Chalkowski 

# In/out --------------------------------

#input: georem type refrence file, ctmm trajectories
#output: weekly summaries of speed and distance moved (i.e. ./1_Data/Objects/pig_weekly_speed_ctmm.rds and ./1_Data/Objects/pig_weekly_distance_ctmm.rds)

# Setup ---------------------------------

#Load libraries
library(tidyverse)
library(sf)
library(ctmm)
library(Rcpp)
library(RcppArmadillo)

#Set directories
# homedir <- "//aapcoftc3fp13/Projects/MUDD/ASF_NIFA/Pipelines/Removals_Mvmt"
# homedir <- "C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Contact Analysis/Removals_Mvmt"
 homedir <- "/Users/kayleigh.chalkowski/OneDrive/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline"
ctmm_dir <- file.path(homedir,"1_Data","Objects","ctmm_Predictions",fsep=.Platform$file.sep)
objdir=file.path(homedir,"1_Data","Objects",fsep=.Platform$file.sep)

#Load geolocation data
# georem <- read.csv("./1_Data/Objects/geo_remtyp_period.csv")
geo.aer<-readRDS(file.path(objdir,"geoaer.rds"))
geo.tox<-readRDS(file.path(objdir,"geotox.rds"))
geo.trap<-readRDS(file.path(objdir,"geotrap.rds"))
# * ! Note, combining geo_rem dfs 
# be very careful when combining any of these-- I see it's just used below to add sex back in
# any time geo_rem dfs are combined, needs to have a column for trt_ctrl-- otherwise no way to tell which controls go to which removal types
# all removal types had same control pigs, but the actual trajectories of those control pigs vary based on the time scales relevant for each removal period
geo.all <- rbind(geo.aer,geo.tox,geo.trap)

# georem$date_only <- as.Date(georem$date_only)
# georem$datetime <- as.POSIXct(georem$datetime,format="%Y-%m-%d %H:%M:%S", tz="UTC")

#Create summary of pig animalids, removal types, and periods
pigs_trt <- geo.all %>% 
  group_by(Removal.Type,removal.period.akdecalc) %>% 
  reframe(animalid=unique(animalid))%>% 
  filter(!(Removal.Type=="aer" & removal.period.akdecalc=="during"))

# get animalids for each treatment/ctrl type
ctrl_id <- pigs_trt %>% filter(Removal.Type=="ctrl")
aer_id <- pigs_trt %>% filter(Removal.Type=="aer")
tox_id <- pigs_trt %>% filter(Removal.Type=="tox")
trap_id <- pigs_trt %>% filter(Removal.Type=="trap")

#summarize dates by week for aerial removal trt/ctrl
trt_aer <- geo.aer %>%
  group_by(removal.period.akdecalc,Removal.Type,week) %>%
  summarise(week_start = min(date_only),
            week_end = max(date_only))

#summarize dates by week for tox removal trt/ctrl
trt_tox <- geo.tox %>%
  group_by(removal.period.akdecalc,Removal.Type,week) %>%
  summarise(week_start = min(date_only),
            week_end = max(date_only))

#summarize dates by week for trap removal trt/ctrl
trt_trap <- geo.trap %>%
  group_by(removal.period.akdecalc,Removal.Type,week) %>%
  summarise(week_start = min(date_only),
            week_end = max(date_only))

#get vectors for periods/removal type for looping
periods <- unique(geo.all$removal.period.akdecalc)
rem_typ <- unique(geo.all$Removal.Type)[-1] #no control group

#source rcpp function (below syntax should work on all platforms)
cpp_dir = file.path(homedir,"2_Scripts","Functions", fsep = .Platform$file.sep)
Rcpp::sourceCpp(file.path(cpp_dir,"TestPairwise.cpp", fsep = .Platform$file.sep), verbose=TRUE)

# Get weekly distances --------------------------------

#Initialize output lists
pig_dist_per_rem <- pig_speed_per_rem <- list()

#Set pig ids to remove from further analysis if they show up
#Remove aerial pigs that were accidentally culled:
aer_gone=c("48476_2_4Y_4Y","85401_2_U_U","85440_E2_E2")
#Remove tox pig with no data in 'before' period
tox_gone=c("86070_H2_H2")

#start loop through folders, get predictions
for(i in 1:length(rem_typ)){
  
  print(paste0("starting rem type ", rem_typ[i]))
  pig_id <- list.files(file.path(ctmm_dir,rem_typ[i], fsep = .Platform$file.sep))
  pig_id <- pig_id[!pig_id%in%rem_typ]
  
  #Remove gone pigs
  pig_id=pig_id[!(pig_id%in%aer_gone|pig_id%in%tox_gone)]
  
  #initialize empty lists
  pig_dist_per_rem[[i]] <- list()
  pig_speed_per_rem[[i]] <- list()
  
  #loop through periods
  for(k in 1:length(periods)){
    
    #Start sub-lists for period subsets
    pig_dist_per_rem[[i]][[k]] <- NA
    pig_speed_per_rem[[i]][[k]] <- NA
    
    if(!(rem_typ[i]=="aer" & periods[k]=="during")){
      
      #loop through pig ids
      for(j in 1:length(pig_id)){
        
        filename <- file.path(ctmm_dir,
                           rem_typ[i],
                           pig_id[j],
                           paste0(pig_id[j],
                           "_",periods[k],"_one_min_pred.rds"),
                           fsep = .Platform$file.sep)
        
        if(file.exists(filename)){
          traj <- readRDS(filename)
          
          lag_dist <- TestPairwise(traj$x,traj$y,
                                   lag(traj$x),lag(traj$y),NA)
          traj$dist <- lag_dist
          traj$date <- as.POSIXct(traj$t,tz="UTC")
          
          traj_xy <- cbind.data.frame(x=traj$x,y=traj$y,date=traj$date,dist=traj$dist)
          
          speed_mat <- traj_xy %>% 
            group_by(hour=floor_date(date,"hour")) %>% 
            summarise(hr_dist = sum(dist,na.rm=T),
                      mn_x=mean(x,na.rm=T),
                      mn_y=mean(y,na.rm=T)) %>% 
            mutate(animalid=pig_id[j],
                   rem_typ=rem_typ[i],
                   period=periods[k])
          
          dist_mat <- traj_xy %>% 
            group_by(hour=floor_date(date,"hour")) %>% 
            summarise(hr_dist = sum(dist,na.rm=T),
                      mn_x=mean(x,na.rm=T),
                      mn_y=mean(y,na.rm=T)) %>% 
            group_by(day=floor_date(hour,"day")) %>% 
            summarise(daily_dist_m=sum(hr_dist,na.rm=T),
                      mn_x=mean(mn_x,na.rm=T),
                      mn_y=mean(mn_y,na.rm=T)
                      # mn_m_per_hr=mean(hr_dist),
                      # min_m_per_hr=min(hr_dist),
                      # max_m_per_hr=max(hr_dist),
                      # md_m_per_hr=median(hr_dist)
                      ) %>% 
            mutate(animalid=pig_id[j],
                   rem_typ=rem_typ[i],
                   period=periods[k])
          
          if(rem_typ[i]=="aer"){
            trt_week <- trt_aer
          } else if(rem_typ[i]=="tox"){
            trt_week <- trt_tox
          } else {
            trt_week <- trt_trap
          }

          #Get weeks
          speed_mat$week <- sapply(1:nrow(speed_mat),function(x){
            wk <- unique(trt_week$week[trt_week$week_start<=as.Date(speed_mat$hour[x]) &
                                       trt_week$week_end>=as.Date(speed_mat$hour[x]) &
                                       trt_week$removal.period.akdecalc==periods[k] &
                                       trt_week$Removal.Type==rem_typ[i]])
            if(sum(wk)==0){
              # wk<-NA
              wk <- unique(trt_week$week[trt_week$week_start<=as.Date(speed_mat$hour[x]) &
                                         trt_week$week_end>=as.Date(speed_mat$hour[x]) &
                                         trt_week$Removal.Type=="ctrl"])
            }
            return(wk)
          })
          
          dist_mat$week <- sapply(1:nrow(dist_mat),function(x){
            wk <- unique(trt_week$week[trt_week$week_start<=as.Date(dist_mat$day[x]) &
                                       trt_week$week_end>=as.Date(dist_mat$day[x]) &
                                       trt_week$removal.period.akdecalc==periods[k] &
                                       trt_week$Removal.Type==rem_typ[i]])
            if(sum(wk)==0){
              # wk<-NA
              wk <- unique(trt_week$week[trt_week$week_start<=as.Date(dist_mat$day[x]) &
                                         trt_week$week_end>=as.Date(dist_mat$day[x]) &
                                         trt_week$Removal.Type=="ctrl"])
            }
            return(wk)
          })
          
          #filter out weeks that have less than 6 days of data
          dist_mat=
            dist_mat %>% 
            dplyr::group_by(week) %>%
            dplyr::mutate(nd=n()) %>%
            dplyr::filter(nd>=6) %>%
            dplyr::select(!nd)
          
          speed_mat=
            speed_mat %>% 
            dplyr::group_by(week) %>%
            dplyr::mutate(nd=n()) %>%
            dplyr::filter(nd>=6) %>%
            dplyr::select(!nd)
          
          #combine into list
          if(!is.data.frame(pig_dist_per_rem[[i]][[k]])){
            
            pig_dist_per_rem[[i]][[k]] <- dist_mat
            
          } else {
            
            pig_dist_per_rem[[i]][[k]] <- rbind.data.frame(pig_dist_per_rem[[i]][[k]],
                                                           dist_mat)
          }
          
          if(!is.data.frame(pig_speed_per_rem[[i]][[k]])){
            
            pig_speed_per_rem[[i]][[k]] <- speed_mat
            
          } else {
            
            pig_speed_per_rem[[i]][[k]] <- rbind.data.frame(pig_speed_per_rem[[i]][[k]],
                                                           speed_mat)
          }
        } #ctmm file exists
        
      } #pig id loop
      
    } #skip aerial during
    
  } # treatment period loop
  
  pig_dist_per_rem[[i]] <- do.call("rbind.data.frame",pig_dist_per_rem[[i]])
  pig_speed_per_rem[[i]] <- do.call("rbind.data.frame",pig_speed_per_rem[[i]])

  #remove any cols with all NA (result of missing period, like in aerial)
  pig_dist_per_rem[[i]]=pig_dist_per_rem[[i]][complete.cases(pig_dist_per_rem[[i]]),]
  pig_speed_per_rem[[i]]=pig_speed_per_rem[[i]][complete.cases(pig_speed_per_rem[[i]]),]
  
  #add column for trt/control
  pig_dist_per_rem[[i]]$trt_ctrl="trt"
  pig_speed_per_rem[[i]]$trt_ctrl="trt"
  
  #set control pigs in trt_ctrl-- if animalidid in ctrl_id
  pig_dist_per_rem[[i]]$trt_ctrl[pig_dist_per_rem[[i]]$animalid%in%as.character(ctrl_id$animalid)]<-"ctrl"
  pig_speed_per_rem[[i]]$trt_ctrl[pig_speed_per_rem[[i]]$animalid%in%as.character(ctrl_id$animalid)]<-"ctrl"

  } # removal type loop

#rbind everything
pig_dist_all <- do.call("rbind.data.frame",pig_dist_per_rem)           
pig_speed_all <- do.call("rbind.data.frame",pig_speed_per_rem)   

#Adjust week numbering
#Make another function to adjust week intervals
adjust_intvals<-function(geo){
  
  for(rem in 1:length(unique(geo$rem_typ))){
  geo.rem=geo[geo$rem_typ==unique(geo$rem_typ)[rem],]
  ids=unique(geo.rem$animalid)
  for(i in 1:length(ids)){
    geo.rem_i=geo.rem[geo.rem$animalid==ids[i],]
    key1=unique(geo.rem_i$week)
    key2=1:length(key1)
    
    for(k in 1:length(key1)){
      geo.rem[geo.rem$animalid==ids[i]&
                geo.rem$week==key1[k],]$week<-key2[k]
    }
    
  } #for animalid
  if(rem==1){
    geo.out=geo.rem
  } else{
    geo.out=rbind(geo.out,geo.rem)
  }
  } #for remtype
    return(geo.out)
}

pig_dist_all=adjust_intvals(pig_dist_all)
pig_speed_all=adjust_intvals(pig_speed_all)

# weekly summaries-----------------
## distance ----------------

#add sex
pig_dist_all <- pig_dist_all %>% dplyr::left_join(geo.all %>% dplyr::select(sex,animalid) %>% distinct())

#per individual per week
pig_dist_wk <- pig_dist_all %>% 
  group_by(animalid,rem_typ,trt_ctrl,period,sex,week) %>% 
  summarise(weekly_dist_m=sum(daily_dist_m),
            weekly_dist_km=weekly_dist_m/1000,
            mn_x=mean(mn_x),
            mn_y=mean(mn_y)) %>% 
  filter(!is.na(rem_typ)) %>% 
  dplyr::select(-weekly_dist_m) %>% 
  rename(Removal.Type=rem_typ,
         removal.period.akdecalc=period,
         mX=mn_x,
         mY=mn_y)

#formatting 
#offset to prevent 0's
pig_dist_wk$weekly_dist_km <- pig_dist_wk$weekly_dist_km +0.0001
pig_dist_wk$animalid <- factor(pig_dist_wk$animalid)
pig_dist_wk$sex <- factor(pig_dist_wk$sex)

#trt_ctrl levels
pig_dist_wk$trt_ctrl <- factor(pig_dist_wk$trt_ctrl,levels=c('ctrl','trt'))

saveRDS(pig_dist_wk,paste0(objdir,"/pig_weekly_distance_ctmm.rds"))

## speed -------------------------

#add sex
pig_speed_all <- pig_speed_all %>% left_join(geo.all %>% dplyr::select(sex,animalid) %>% distinct())

#per individual per week
pig_speed_wk <- pig_speed_all %>% 
  mutate(hr_dist_km=hr_dist/1000) %>% 
  group_by(animalid,rem_typ,trt_ctrl,period,sex,week) %>% 
  summarise(weekly_md_km_hr=median(hr_dist_km),
            mX=mean(mn_x),
            mY=mean(mn_y)) %>% 
 filter(!is.na(rem_typ)) %>% 
  rename(Removal.Type=rem_typ,
         removal.period.akdecalc=period)


#offset by very small amount
pig_speed_wk$weekly_md_km_hr <- pig_speed_wk$weekly_md_km_hr + 0.0001

#set factors
pig_speed_wk$animalid <- factor(pig_speed_wk$animalid)
pig_speed_wk$sex <- factor(pig_speed_wk$sex)
pig_speed_wk$trt_ctrl <- factor(pig_speed_wk$trt_ctrl,levels=c('ctrl','trt'))

#save output
saveRDS(pig_speed_wk,paste0(objdir,"/pig_weekly_speed_ctmm.rds"))

