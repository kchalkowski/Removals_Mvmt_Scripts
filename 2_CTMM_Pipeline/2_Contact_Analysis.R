### Title: Contact analysis for NIFA pigs - parallelized 
### Author: Abbey Feuka
### Date: 19NOV24
### Notes: ctmm predictions from predict_ctmms_abf.R
### uses Kayleigh Chalkowski's TestPairwise.cpp function (faster than ctmm::distances()

#input: georem type reference file, ctmm trajectories, AKDE files
#output: list of dataframes for each animal's pairwise contacts and degree within three distance thresholds

# first level is removal/treatment type
## second level is period (before/during/after)
### each data frame has all individual combinations per week with # of contacts and # of individuals contacted (Degree) (i.e. ./1_Data/Objects/pairwise_contacts.RDS)

library(tidyverse)
library(snow)

# homedir <- "//aapcoftc3fp13/Projects/MUDD/ASF_NIFA/Pipelines/Removals_Mvmt" #nifa server
homedir <- "C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Contact Analysis/Removals_Mvmt" #abbey local
# homedir <- "/cm/shared/NFS/Projects/AbigailFeuka/Removals_Mvmt" #hpc
# homedir <- "/home/abigail.feuka/files/Removals_Mvmt" #hpc (win scp)

cpp_dir <- "./Functions/"
ctmm_dir <- paste0(homedir,"/1_Data/Objects/ctmm_Predictions/")
objdir <- file.path(homedir,"1_Data","Objects",fsep=.Platform$file.sep)

#load data -------
georem <- read.csv(paste0(objdir,"/geo_remtyp.csv"))

geo.aer<-readRDS(file.path(objdir,"geoaer.rds"))
geo.tox<-readRDS(file.path(objdir,"geotox.rds"))
geo.trap<-readRDS(file.path(objdir,"geotrap.rds"))

periods <- unique(geo.trap$removal.period.akdecalc)
rem_typ <- c("ctrl","aer","tox","trap")
rem_typ_folders <- rem_typ[-which(rem_typ=="ctrl")]
cdist <- 10

#all pairs within 95% MCP HR
hr_pairs <- readRDS(paste0(objdir,"/pairs.rds"))
hr_pairs$area <- round(as.numeric(hr_pairs$area),0)

#filter out satellite animals
sat_pigs <- georem %>% 
  group_by(animalid) %>% 
  summarise(data_from=unique(data_from)) %>% 
  filter(data_from=='satellite')
  
# rem_idx<-2
# per_idx <- 1
# w<-2
# k<-1
# l<-2
# d<-1


#parallel function ----------
contact_from_ctmm <- function(rem_idx,per_idx){ # parallel over rem type FOLDERS (3)
  
  require(tidyverse)
  require(ctmm)
  require(Rcpp)
  require(RcppArmadillo)

  #source rcpp function
  Rcpp::sourceCpp(paste0(cpp_dir,"TestPairwise.cpp"), verbose=TRUE)

    dist_mat <- list()
    
      if(!(rem_typ_folders[rem_idx]=="aer" & periods[per_idx]=="during")){ #no during aerial 
          
        # ref_trt <- akde_files %>% filter(rem_typ_folder==rem_typ_folders[rem_idx])
        # folder_trmts <- unique(ref_trt$rem_typ)
        folder_trmts <- c("ctrl","trt")
        
        hr_pairs_rem <- hr_pairs %>% filter(rem==rem_typ_folders[rem_idx])

          for(w in 1:length(folder_trmts)){# ctrl vs trt

            hr_pairs_rem_trt <- hr_pairs_rem %>% filter(trt_ctrl_1==folder_trmts[w])
            main_pigs_rem_trt <- unique(c(hr_pairs_rem_trt$id_1,hr_pairs_rem_trt$id_2))
            
            # pigs_trt_temp <- akde_files %>%
            #   filter(rem_typ_folder==rem_typ_folders[rem_idx]) %>%
            #   filter(rem_typ==folder_trmts[w])
            
            # distance matrix for number of contacts and number of individuals contacted (degree)
            # for each removal type and treatment period combo
            dist_mat[[w]] <- expand_grid(rem_typ=rem_typ_folders[rem_idx],
                                         trt_typ=folder_trmts[w],
                                         period=periods[per_idx],
                                         animalid=main_pigs_rem_trt,
                                         dist=cdist,
                                         num_contacts=NA,
                                         num_indivs=NA)
            
            #pairs version
            # for(d in 1:length(cdist)){
            #   if(d==1){
            #     dist_mat[[w]] <- hr_pairs_rem_trt %>% 
            #       dplyr::select(rem,id_1,id_2,trt_ctrl_1) %>% 
            #       rename(rem_typ=rem,
            #              trt_typ=trt_ctrl_1) %>% 
            #       mutate(period=periods[per_idx],
            #              dist=cdist[d],
            #              num_contacts=NA,
            #              num_indivs=NA)
            #   } else {
            #     temp <- hr_pairs_rem_trt %>% 
            #       dplyr::select(rem,id_1,id_2,trt_ctrl_1) %>% 
            #       rename(rem_typ=rem,
            #              trt_typ=trt_ctrl_1) %>% 
            #       mutate(period=periods[per_idx],
            #              dist=cdist[d],
            #              num_contacts=NA,
            #              num_indivs=NA)
            #     
            #     dist_mat[[w]] <- rbind(dist_mat[[w]],temp)
            #   }
            # }
              
            # for(k in 1:2){
            for(k in 1:length(main_pigs_rem_trt)){ # each indiv in treatment/ctrl group
              
              main_pig_id <- main_pigs_rem_trt[k]
              
              main_pig_pairs <- hr_pairs_rem_trt %>% 
                filter(id_1==main_pig_id |
                         id_2==main_pig_id) %>% 
                filter(!duplicated(area))
              
              if(file.exists(paste0(ctmm_dir,
                                    rem_typ_folders[rem_idx],
                                    "/",main_pig_id,
                                    "/",main_pig_id,
                                    "_",periods[per_idx],"_one_min_pred.rds")) & #file exists
                 !(main_pig_id%in%sat_pigs) #data not from satellite
                 ){
                
                main_pig_telem <- readRDS(paste0(ctmm_dir,
                                                 rem_typ_folders[rem_idx],
                                                 "/",main_pig_id,
                                                 "/",main_pig_id,
                                                 "_",periods[per_idx],
                                                 "_one_min_pred.rds"))
                
                main_pig_ctmm <- readRDS(paste0(objdir,"/AKDE_",
                                                rem_typ_folders[rem_idx],
                                                "/",main_pig_id,
                                                "/",main_pig_id,"_",
                                                periods[per_idx],".rds"))
                
                num_contacts_main <- matrix(NA,nrow=length(cdist),ncol=nrow(main_pig_pairs))
                
                pairs_tmp <- unique(c(unique(main_pig_pairs$id_1),unique(main_pig_pairs$id_2)))
                pairs_tmp <- pairs_tmp[which(pairs_tmp!=main_pig_id)]
                
                # for(l in 1:2){
                for(l in 1:length(pairs_tmp)){ # all pairs for indiv pig
                  
                  pair_pig_id <- pairs_tmp[l]
                    
                  if(file.exists(paste0(ctmm_dir,
                                        rem_typ_folders[rem_idx],
                                        "/",pair_pig_id,
                                        "/",pair_pig_id,
                                        "_",periods[per_idx],"_one_min_pred.rds")) & #file exists
                     !(pair_pig_id%in%sat_pigs) #data not from satellite
                  ){
                    
                    pair_pig_telem <- readRDS(paste0(ctmm_dir,
                                                     rem_typ_folders[rem_idx],
                                                     "/",pair_pig_id,
                                                     "/",pair_pig_id,
                                                     "_",periods[per_idx],
                                                     "_one_min_pred.rds"))
                    
                    pair_pig_ctmm <- readRDS(paste0(objdir,"/AKDE_",
                                                    rem_typ_folders[rem_idx],
                                                    "/",pair_pig_id,
                                                    "/",pair_pig_id,"_",
                                                    periods[per_idx],".rds"))
                    
                    pair_pig_pred <- ctmm::predict(pair_pig_ctmm,
                                                   data=pair_pig_telem,
                                                   t=main_pig_telem$t)
                    
                    main_pair_dist <- TestPairwise(main_pig_telem$x,main_pig_telem$y,
                                                   pair_pig_pred$x,pair_pig_pred$y,NA)
                    
                    for(d in 1:length(cdist)){#distance thresholds
                      
                      num_contacts_main[d,l] <- sum(main_pair_dist<cdist[d])
                      
                    } #distance thresholds
                    
                  } #if pair pig exists  
                    
                }  # all pairs for indiv pig
                
                for(d in 1:length(cdist)){ # distance - summing contacts and indivs per pig
                  
                  #number of contacts per week at distance
                  dist_mat[[w]]$num_contacts[dist_mat[[w]]$animalid==main_pig_id &
                                               dist_mat[[w]]$dist==cdist[d]] <-
                    sum(num_contacts_main[d,],na.rm=T)
                  
                  #degree per week at distance
                  cd <- num_contacts_main[d,]
                  cd <- cd[!is.na(cd)]
                  dist_mat[[w]]$num_indivs[dist_mat[[w]]$animalid==main_pig_id &
                                             dist_mat[[w]]$dist==cdist[d]] <-
                    sum(cd>0)
                  
                } # distance - summing contacts and indivs per pig
                
              }#if ctmm exists statement
              
            }#individual pig 
            
          } #ctrl vs trt within folder
        
        dist_mat_all<- do.call("rbind.data.frame",dist_mat)
        dist_mat_all <- as.data.frame(dist_mat_all)
        
      }# aerial during if statement
    
    list(dist_mat_all=dist_mat_all)
}

#trmt/period combinations
perms <- expand_grid(rem_idx=1:length(rem_typ_folders),per_idx=1:length(periods))
perms$rem_typ_folder <- rem_typ_folders[perms$rem_idx]
perms$period <- periods[perms$per_idx]
perms <- perms %>% filter(!(rem_typ_folder=="aer" & period=="during"))
perms


cl <- makeCluster(nrow(perms))
clusterExport(cl, list("periods","cdist",
                       "rem_typ_folders","perms",#"akde_files",
                       "sat_pigs",
                       "cpp_dir","objdir","homedir","ctmm_dir",
                       "hr_pairs",'contact_from_ctmm'))

system.time(
  parContacts <- clusterMap(cl, fun=contact_from_ctmm, rem_idx=perms$rem_idx, per_idx=perms$per_idx)
)
stopCluster(cl)

parContacts1 <- lapply(1:length(parContacts),function(i)parContacts[[i]]$dist_mat_all)
parContacts1 <- do.call("rbind.data.frame",parContacts1)


# split by rem typ and link to sex
contacts_aer <- parContacts1 %>% filter(rem_typ=="aer")
contacts_tox <- parContacts1 %>% filter(rem_typ=="tox")
contacts_trap <- parContacts1 %>% filter(rem_typ=="trap")

contacts_aer <- contacts_aer %>% left_join(geo.aer %>% dplyr::select(animalid,sex) %>% distinct())
contacts_trap <- contacts_trap %>% left_join(geo.trap %>% dplyr::select(animalid,sex) %>% distinct())
contacts_tox <- contacts_tox %>% left_join(geo.tox %>% dplyr::select(animalid,sex) %>% distinct())

#standardize by number of days in period
aer_per_days <- geo.aer %>% 
  group_by(removal.period.akdecalc) %>% 
  summarise(start=min(date_only),
            end=max(date_only),
            period_days=as.numeric(end-start))

trap_per_days <- geo.trap %>% 
  group_by(removal.period.akdecalc) %>% 
  summarise(start=min(date_only),
            end=max(date_only),
            period_days=as.numeric(end-start))

tox_per_days <- geo.tox %>% 
  group_by(removal.period.akdecalc) %>% 
  summarise(start=min(date_only),
            end=max(date_only),
            period_days=as.numeric(end-start))

contacts_aer <- contacts_aer %>% rename(removal.period.akdecalc=period) %>% 
  left_join(aer_per_days %>% select(-c(start,end))) %>%
  mutate(contacts_per_day=num_contacts/period_days,
         indivs_per_day=num_indivs/period_days)
contacts_trap<- contacts_trap %>% rename(removal.period.akdecalc=period) %>% 
  left_join(trap_per_days %>% select(-c(start,end)))%>%
  mutate(contacts_per_day=num_contacts/period_days,
         indivs_per_day=num_indivs/period_days)
contacts_tox <- contacts_tox %>% rename(removal.period.akdecalc=period) %>% 
  left_join(tox_per_days %>% select(-c(start,end)))%>%
  mutate(contacts_per_day=num_contacts/period_days,
         indivs_per_day=num_indivs/period_days)

head(contacts_aer)
head(contacts_trap)
head(contacts_tox)

#save ouptut
saveRDS(contacts_aer,file=paste0(objdir,"/pairwise_contacts_aer.rds"))
saveRDS(contacts_tox,file=paste0(objdir,"/pairwise_contacts_tox.rds"))
saveRDS(contacts_trap,file=paste0(objdir,"/pairwise_contacts_trap.rds"))

#check
contacts_all <- rbind.data.frame(contacts_aer,contacts_trap,contacts_tox)
contacts_all$removal.period.akdecalc <- factor(contacts_all$removal.period.akdecalc,levels=c('before','during','after'))


ggplot(contacts_all)+
  geom_boxplot(aes(x=removal.period.akdecalc,y=contacts_per_day,fill=trt_typ))+
  facet_grid(dist~rem_typ)+
  guides(col="none")+
  theme(text=element_text(size=15))

ggplot(contacts_all)+
  geom_boxplot(aes(x=removal.period.akdecalc,y=indivs_per_day,fill=trt_typ))+
  facet_grid(dist~rem_typ)+
  guides(col="none")+
  theme(text=element_text(size=15))
