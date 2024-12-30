### Title: Contact analysis for NIFA pigs - parallelized 
### Author: Abbey Feuka
### Date: 19NOV24
### Notes: ctmm predictions from predict_ctmms_abf.R
### uses Kayleigh Chalkowski's TestPairwise.cpp function (faster than ctmm::distances()

#home dir is //aapcoftc3fp13/Projects/MUDD/ASF_NIFA/Pipelines/Removals_Mvmt

#input: georem type reference file, ctmm trajectories, AKDE files
#output: list of dataframes for each animal's pairwise contacts and degree within three distance thresholds

# first level is removal/treatment type
## second level is period (before/during/after)
### each data frame has all individual combinations per week with # of contacts and # of individuals contacted (Degree) (i.e. ./1_Data/Objects/pairwise_contacts.RDS)

library(tidyverse)
library(snow)

georem <- read.csv("./1_Data/Objects/geo_remtyp_period.csv")

georem$date_only <- as.Date(georem$date_only)
georem$datetime <- as.POSIXct(georem$datetime,format="%Y-%m-%d %H:%M:%S", tz="UTC")

pigs_trt <- georem %>% 
  group_by(Removal.Type,removal.period.akdecalc) %>% 
  reframe(animalid=unique(animalid))%>% 
  filter(!(Removal.Type=="aer" & removal.period.akdecalc=="during"))

trt_week_orig <- georem %>% 
  group_by(removal.period.akdecalc,Removal.Type,week,) %>% 
  summarise(week_start = min(date_only),
            week_end = max(date_only))

periods <- unique(georem$removal.period.akdecalc)
rem_typ <- unique(georem$Removal.Type)[-1] #no control group
cdist <- c(10,30,50)

#parallel function ----------
contact_from_ctmm <- function(i){ # parallel over rem types (3)
  
  require(tidyverse)
  require(ctmm)
  require(sf)
  require(Rcpp)
  require(RcppArmadillo)
  
  #source rcpp function
  cpp_dir <- "./2_Scripts/Functions/"
  Rcpp::sourceCpp(paste0(cpp_dir,"TestPairwise.cpp"), verbose=TRUE)
  
  out_dir <-"./1_Data/Objects/"
  
    # if(!dir.exists(paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Contact Analysis/ctmm Contacts/",rem_typ[i]))){ #create
    #   dir.create(paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Contact Analysis/ctmm Contacts/",rem_typ[i]))
    # }
    
    dist_mat <- list()
    
    for(j in 1:length(periods)){ #treatment period sub directories
      if(!(rem_typ[i]=="aer" & periods[j]=="during")){ #no during aerial 
        
          pigs_trt_temp <- pigs_trt %>% filter(Removal.Type==rem_typ[i] &
                                                 removal.period.akdecalc==periods[j])
          trt_week <- trt_week_orig %>% 
            filter(removal.period.akdecalc==periods[j] & Removal.Type==rem_typ[i])
          
          #distance matrix for number of contacts and number of individuals contacted (degree)
          #for each removal type and treatment period combo
          dist_mat[[j]] <- expand_grid(rem_typ=unique(pigs_trt_temp$Removal.Type),
                                            period=unique(pigs_trt_temp$removal.period.akdecalc),
                                            animalid=unique(pigs_trt_temp$animalid),
                                            week=unique(trt_week$week),
                                            dist=cdist,
                                            num_contacts=NA,
                                            num_indivs=NA)
          for(k in 1:nrow(pigs_trt_temp)){ # individual pig in treatment
            
            main_pig_id <- pigs_trt_temp$animalid[k]
            
            if(file.exists(paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Contact Analysis/ctmm Predicitions/",
                                  rem_typ[i],"/",main_pig_id,"/",main_pig_id,"_",periods[j],"_one_min_pred.rds"))){
              
              main_pig_telem <- readRDS(paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Contact Analysis/ctmm Predicitions/",
                                               rem_typ[i],"/",main_pig_id,"/",main_pig_id,"_",periods[j],"_one_min_pred.rds"))
              main_pig_ctmm <- readRDS(paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Contact Analysis/Data/AKDE_",
                                              rem_typ[i],"/",main_pig_id,"/",main_pig_id,"_",
                                              periods[j],".rds"))
              
              for(w in trt_week$week){ #for each week in treatment period
                
                week_t <- as.POSIXct(main_pig_telem$t,tz="UTC")
                
                #main pigs dates within week
                week_t <- week_t[as.Date(week_t)>=trt_week$week_start[trt_week$week==w] & 
                                   as.Date(week_t)<=trt_week$week_end[trt_week$week==w]]
                week_t_n <- as.numeric(week_t)
                
                #clip main pig's trajectory to week
                main_pig_telem_wk <- main_pig_telem[main_pig_telem$t%in%week_t_n,]
                
                #matrix of all pig pairs, over individual week w
                temp_contacts <- expand_grid(dist=cdist,
                                             main_pig=unique(pigs_trt_temp$animalid),
                                             pair_pig=unique(pigs_trt_temp$animalid),
                                             week=w,
                                             n_contacts=NA)

                for(l in 1:nrow(pigs_trt_temp)){ # all pairs for indiv pig
                  
                  pair_pig_id <- pigs_trt_temp$animalid[l]
                  
                  if(pair_pig_id!=main_pig_id){ #not same pig if statement
                    
                    pair_pig_telem <- readRDS(paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Contact Analysis/ctmm Predicitions/",
                                                     rem_typ[i],"/",pair_pig_id,"/",pair_pig_id,"_",periods[j],"_one_min_pred.rds"))
                    
                    pair_pig_ctmm <- readRDS(paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Contact Analysis/Data/AKDE_",
                                                    rem_typ[i],"/",pair_pig_id,"/",pair_pig_id,"_",
                                                    periods[j],".rds"))
                    
                    pair_pig_pred <- ctmm::predict(pair_pig_ctmm,
                                                   data=pair_pig_telem,
                                                   t=main_pig_telem_wk$t)
                    
                    main_pair_dist <- TestPairwise(main_pig_telem_wk$x,main_pig_telem_wk$y,
                                                   pair_pig_pred$x,pair_pig_pred$y,NA)
                    
                    for(d in 1:length(cdist)){#distance thresholds
                      temp_contacts$n_contacts[temp_contacts$main_pig==main_pig_id &
                                                 temp_contacts$pair_pig==pair_pig_id &
                                                 temp_contacts$dist==cdist[d]] <- sum(main_pair_dist<cdist[d])
                    } #distance thresholds
                    # temp_contacts %>% filter(is.na(n_contacts)) %>% filter(main_pig==main_pig_id &
                    #                                                        pair_pig==pair_pig_id) #check
                  }  #not same pig if statement
                }  # all pairs for indiv pig
                
                # temp_contacts %>% filter(is.na(n_contacts)) %>% filter(main_pig==main_pig_id) #check
                # temp_contacts %>% filter(n_contacts>0) %>% filter(main_pig==main_pig_id) #check
                
                for(d in 1:length(cdist)){ # distance - summing contacts and indivs per pig
                  #number of contacts per week at distance
                  dist_mat[[j]]$num_contacts[dist_mat[[j]]$animalid==main_pig_id &
                                                    dist_mat[[j]]$week==w &
                                                    dist_mat[[j]]$dist==cdist[d]] <-
                    sum(temp_contacts$n_contacts[temp_contacts$main_pig==main_pig_id &
                                                   temp_contacts$dist==cdist[d]],
                        na.rm=T)
                  
                  #degree per week at distance
                  dist_mat[[j]][dist_mat[[j]]$animalid==main_pig_id &
                                       dist_mat[[j]]$week==w &
                                       dist_mat[[j]]$dist==cdist[d],]
                  
                  dist_mat[[j]]$num_indivs[dist_mat[[j]]$animalid==main_pig_id &
                                                  dist_mat[[j]]$week==w &
                                                  dist_mat[[j]]$dist==cdist[d]] <-
                    length(temp_contacts$pair_pig[temp_contacts$main_pig==main_pig_id &
                                                    temp_contacts$dist==cdist[d] &
                                                    temp_contacts$n_contacts>0 &
                                                    !is.na(temp_contacts$n_contacts)])
                } # distance - summing contacts and indivs per pig
                
                # dist_mat[[j]] %>% filter(week==w & animalid==main_pig_id) #check
              }# week in treatment period
              
              # dist_mat[[j]] %>% filter(animalid==main_pig_id) %>% filter(is.na(num_contacts))#check
            }#if ctmm exists statement
          }#individual pig 
          # dist_mat[[j]] %>% filter(num_contacts>0)
          # dist_mat_all %>% filter(num_indivs>0)
      }# aerial during if statement
    }#treatment period loop
    dist_mat_all <- do.call("rbind.data.frame",dist_mat)
    
    return(dist_mat_all)
}

cl <- makeCluster(length(rem_typ), "SOCK")
clusterExport(cl, list("periods","rem_typ","cdist","pigs_trt",
                       'trt_week_orig','contact_from_ctmm'))

system.time(
  parContacts <- clusterApply(cl, 1:length(rem_typ), contact_from_ctmm)
)
stopCluster(cl)

saveRDS(parContacts,file=paste0(out_dir,"pairwise_contacts.RDS"))
