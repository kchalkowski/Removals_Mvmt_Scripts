### Title: Summarizing weekly individual speed and distance moved from ctmm models
### Author: Abbey Feuka
### Date: 26NOV24
### uses ctmms fit by Kayleigh Chalkowski 

#input: georem type refrence file, ctmm trajectories
#output: weekly summaries of speed and distance moved (i.e. ./1_Data/Objects/pig_weekly_speed_ctmm.rds and ./1_Data/Objects/pig_weekly_distance_ctmm.rds)

library(tidyverse)
library(sf)
library(ctmm)
library(Rcpp)
library(RcppArmadillo)

homedir <- "//aapcoftc3fp13/Projects/MUDD/ASF_NIFA/Pipelines/Removals_Mvmt"
# homedir <- "C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Contact Analysis/Removals_Mvmt"

ctmm_dir <- paste0(homedir,"./1_Data/Objects/ctmm_Predictions/")
out_dir <- paste0(homedir,"./1_Data/Objects/")
objdir=file.path(homedir,"1_Data","Objects",fsep=.Platform$file.sep)

# georem <- read.csv("./1_Data/Objects/geo_remtyp_period.csv")
geo.aer<-readRDS(file.path(objdir,"geoaer.rds"))
geo.tox<-readRDS(file.path(objdir,"geotox.rds"))
geo.trap<-readRDS(file.path(objdir,"geotrap.rds"))
geo.all <- rbind(geo.aer,geo.tox,geo.trap)

# georem$date_only <- as.Date(georem$date_only)
# georem$datetime <- as.POSIXct(georem$datetime,format="%Y-%m-%d %H:%M:%S", tz="UTC")

pigs_trt <- geo.all %>% 
  group_by(Removal.Type,removal.period.akdecalc) %>% 
  reframe(animalid=unique(animalid))%>% 
  filter(!(Removal.Type=="aer" & removal.period.akdecalc=="during"))

ctrl_id <- pigs_trt %>% filter(Removal.Type=="ctrl")
aer_id <- pigs_trt %>% filter(Removal.Type=="aer")
tox_id <- pigs_trt %>% filter(Removal.Type=="tox")
trap_id <- pigs_trt %>% filter(Removal.Type=="trap")

trt_aer <- geo.aer %>%
  group_by(removal.period.akdecalc,Removal.Type,week) %>%
  summarise(week_start = min(date_only),
            week_end = max(date_only))

trt_tox <- geo.tox %>%
  group_by(removal.period.akdecalc,Removal.Type,week) %>%
  summarise(week_start = min(date_only),
            week_end = max(date_only))
trt_trap <- geo.trap %>%
  group_by(removal.period.akdecalc,Removal.Type,week) %>%
  summarise(week_start = min(date_only),
            week_end = max(date_only))

periods <- unique(geo.all$removal.period.akdecalc)
rem_typ <- unique(geo.all$Removal.Type)[-1] #no control group

#source rcpp function
cpp_dir <- "./Functions/"
Rcpp::sourceCpp(paste0(cpp_dir,"TestPairwise.cpp"), verbose=TRUE)

pig_dist_per_rem <- pig_speed_per_rem <- list()
# i<-1
# k<-1
# j<-2
geo.aer %>% filter(animalid==pig_id[j]) %>% select(date_only) %>% reframe(range(date_only))
pig_id[j]
for(i in 1:length(rem_typ)){
  
  pig_id <- list.files(paste0(ctmm_dir,rem_typ[i]))
  pig_id <- pig_id[!pig_id%in%rem_typ]
  
  pig_dist_per_rem[[i]] <- list()
  pig_speed_per_rem[[i]] <- list()
  
  for(k in 1:length(periods)){
    
    pig_dist_per_rem[[i]][[k]] <- NA
    pig_speed_per_rem[[i]][[k]] <- NA
    
    if(!(rem_typ[i]=="aer" & periods[k]=="during")){
      
      for(j in 1:length(pig_id)){
        
        filename <- paste0(ctmm_dir,
                           rem_typ[i],
                           "/",pig_id[j],
                           "/",pig_id[j],
                           "_",periods[k],"_one_min_pred.rds")
        
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
} # removal type loop

pig_dist_all <- do.call("rbind.data.frame",pig_dist_per_rem)           
pig_speed_all <- do.call("rbind.data.frame",pig_speed_per_rem)   

# saveRDS(pig_dist_all,file=paste0(out_dir,"daily_distance_nifa.rds"))
# saveRDS(pig_speed_all,file=paste0(out_dir,"hourly_speed_nifa.rds"))


# weekly summaries-----------------
## distance ----------------
# pig_dist_all <- readRDS(paste0(out_dir,"daily_distance_nifa.rds"))

#add controls in 
ctrl_ids <- georem %>% filter(Removal.Type=="ctrl") %>%
  reframe(animalid=unique(animalid))
pig_dist_all$rem_typ[pig_dist_all$animalid%in%ctrl_ids$animalid] <- "ctrl"

#add sex
pig_dist_all <- pig_dist_all %>% left_join(georem %>% select(sex,animalid) %>% distinct())

#per individual per week
pig_dist_wk <- pig_dist_all %>% 
  group_by(animalid,rem_typ,period,sex,week) %>% 
  summarise(weekly_dist_m=sum(daily_dist_m),
            weekly_dist_km=weekly_dist_m/1000,
            mn_x=mean(mn_x),
            mn_y=mean(mn_y)) %>% 
  filter(!is.na(rem_typ)) %>% 
  select(-weekly_dist_m) %>% 
  rename(Removal.Type=rem_typ,
         removal.period.akdecalc=period,
         mX=mn_x,
         mY=mn_y)

saveRDS(pig_dist_wk,paste0(out_dir,"pig_weekly_distance_ctmm.rds"))

#entire treatment group per week -----------
pig_dist_wk_trt <- pig_dist_wk %>% 
  group_by(week,Removal.Type,removal.period.akdecalc) %>% 
  summarise(md=median(weekly_dist_km),
            sd=sqrt(var(weekly_dist_km)),
            mn=mean(weekly_dist_km),
            lci=mn - 1.96*(sd/sqrt(n())),
            uci=mn + 1.96*(sd/sqrt(n())),
            min=min(weekly_dist_km),
            max=max(weekly_dist_km))

## speed -------------------------
# pig_speed_all <- readRDS(paste0(out_dir,"hourly_speed_nifa.rds"))

pig_speed_all$rem_typ[pig_speed_all$animalid%in%ctrl_ids$animalid] <- "ctrl"

#add sex
pig_speed_all <- pig_speed_all %>% left_join(georem %>% select(sex,animalid) %>% distinct())

#per individual per week
pig_speed_wk <- pig_speed_all %>% 
  mutate(hr_dist_km=hr_dist/1000) %>% 
  group_by(animalid,rem_typ,period,sex,week) %>% 
  summarise(weekly_md_km_hr=median(hr_dist_km),
            mX=mean(mn_x),
            mY=mean(mn_y)) %>% 
 filter(!is.na(rem_typ)) %>% 
  rename(Removal.Type=rem_typ,
         removal.period.akdecalc=period)

saveRDS(pig_speed_wk,paste0(out_dir,"pig_weekly_speed_ctmm.rds"))

