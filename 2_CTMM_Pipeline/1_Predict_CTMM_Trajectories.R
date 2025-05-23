### Title: Creating paths from NIFA ctmm outputs
### Author: Abbey Feuka
### Date: 30OCT24
### Notes: ctmms fit by Kayliegh Chalkowski

#input: georem_type reference file, AKDE files in "./1_Data/Objects/"
#output: ctmm trajectories in ctmm_Predictions e.g. "./1_Data/Objects/ctmm_Predictions/aer/33336_1T_1T/33336_1T_1T_before_one_min_preds.rds

library(ctmm)
library(tidyverse)
library(move)

# homedir <- "//aapcoftc3fp13/Projects/MUDD/ASF_NIFA/Pipelines/Removals_Mvmt"
# homedir <- "C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Contact Analysis/Removals_Mvmt"
# homedir <- "/Users/kayleigh.chalkowski/OneDrive/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline"

akde_dir <- paste0(homedir,"/1_Data/Objects/")
ctmm_dir <- paste0(homedir,"/1_Data/Objects/ctmm_Predictions/")
objdir=file.path(homedir,"1_Data","Objects",fsep=.Platform$file.sep)

#load data -------
#georem <- read.csv("./1_Data/Objects/georemtyp_period.csv")
geo.aer<-readRDS(file.path(objdir,"geoaer.rds"))
geo.tox<-readRDS(file.path(objdir,"geotox.rds"))
geo.trap<-readRDS(file.path(objdir,"geotrap.rds"))
#####^ needs be separated geo.rem, since new files contain already split periods and assigned weeks

periods <- unique(geo.tox$removal.period.akdecalc)
rem_typ <- c("aer","tox","trap")

#predict one-minute trajectories using original data
# per removal type and treatment period
# k<-1
# i<-2
# j<-1

for(k in 1:length(rem_typ)){ #removal type loop
  pig_mod <- pig_pred <-list()
  pig_files <- list.files(paste0(akde_dir,"AKDE_",rem_typ[k]))
  
  #create removal type directory
  if(!file.exists(paste0(ctmm_dir,"/",rem_typ[k]))){
    dir.create(file.path(ctmm_dir,"/",rem_typ[k]))
  }

  for(i in 1:length(pig_files)){ #individual pig loop
    
    #create individual directory within removal type dir
    if(!file.exists(paste0(ctmm_dir,rem_typ[k],"/",pig_files[i]))){
      dir.create(file.path(paste0(ctmm_dir,rem_typ[k],"/"),pig_files[i]))
    }

    for(j in 1:length(periods)){ # treatment period loop
      if(!(rem_typ[k]=="aerial" & j==2)){ #no during for aerial ops
        
        #skip missing AKDE files
        if(file.exists(paste0(akde_dir,"AKDE_",rem_typ[k],"/",
                              pig_files[i],"/",pig_files[i],"_",
                              periods[j],".rds"))
        # & 
        #    #skip already made predictions
        #    !file.exists(paste0(ctmm_dir,
        #           rem_typ[k],"/",
        #           pig_files[i],"/",
        #           pig_files[i],"_",
        #           periods[j],
        #           "_one_min_pred.rds"))
           ){
          #grab pig ctmm model
          pig_mod[[i]] <- readRDS(paste0(akde_dir,"AKDE_",rem_typ[k],"/",
                                         pig_files[i],"/",pig_files[i],"_",
                                         periods[j],".rds"))

          if(rem_typ[k]=="aer"){
            georem <- geo.aer
          } else if(rem_typ[k]=="tox"){
            georem <- geo.tox
          } else {
            georem <- geo.trap
          }
          
          pig_dat <- georem %>% filter(animalid==pig_files[i] & 
                                          removal.period.akdecalc==periods[j]) %>% 
            arrange(datetime)

          pig_dat <- pig_dat[!duplicated(pig_dat$datetime),]

          if(nrow(pig_dat)>0){ #pigs with akdes
            pig_dat$datetime <- as.POSIXct(pig_dat$datetime,tz="UTC")
            
            #use UTM coordinates only
            pig_dat <- pig_dat %>% dplyr::select(-c(longitude,latitude))
            
            #create move object
            pig_dat_move <- move(x=pig_dat$X,
                                 y=pig_dat$Y,
                                 time=pig_dat$datetime,
                                 proj="+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",#"epsg:32614",
                                 sensor="gps",
                                 animal=unique(pig_dat$animalid))
            
            #convert move object to telemetry object
            pig_dat_telem <- list(as.telemetry(pig_dat_move,
                                               timeformat="%Y-%m-%d %H:%M:%S",
                                               timezone = "UTC",
                                               projection="+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
            )
            
            #identify pig
            pig_name <- as.character(unique(pig_dat$animalid))
            names(pig_dat_telem) <- pig_name
            attr(pig_dat_telem[[1]],"info")[1]$identity <- pig_name
            
            #predict trajectory at one minute intervals
            pig_pred[[i]] <- ctmm::predict(pig_mod[[i]],
                                           data=pig_dat_telem[[1]],
                                           t=seq(from=min(pig_dat_telem[[1]]$timestamp),
                                                 to=max(pig_dat_telem[[1]]$timestamp),
                                                 by="min")
            )
            
            #save to directory
            saveRDS(pig_pred[[i]],
                    file=paste0(ctmm_dir,
                                rem_typ[k],"/",
                                pig_name,"/",
                                pig_name,"_",
                                periods[j],
                                "_one_min_pred.rds")
            )
          }
        } #skip missing AKDEs
      }
    }
  }
}

