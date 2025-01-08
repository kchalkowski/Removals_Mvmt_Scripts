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
georem <- read.csv(paste0(objdir,"/geo_remtyp_period.csv"))

geo.aer<-readRDS(file.path(objdir,"geoaer.rds"))
geo.tox<-readRDS(file.path(objdir,"geotox.rds"))
geo.trap<-readRDS(file.path(objdir,"geotrap.rds"))
geo.all <-rbind(geo.aer,geo.tox,geo.trap)

# pigs_trt <- geo.all %>%
#   group_by(Removal.Type,removal.period.akdecalc) %>%
#   reframe(animalid=unique(animalid))%>%
#   filter(!(Removal.Type=="aer" & removal.period.akdecalc=="during"))
# saveRDS(pigs_trt,"C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Contact Analysis/Removals_Mvmt/1_Data/Input/pigs_trt.rds")

pigs_trt <- readRDS(paste0(homedir,"/1_Data/Input/pigs_trt.rds"))

periods <- unique(geo.all$removal.period.akdecalc)
rem_typ <- unique(geo.all$Removal.Type) #no control group
rem_typ_folders <- rem_typ[-which(rem_typ=="ctrl")]
cdist <- c(10,30,50)

#identify controls in each folder 
akde_files<-list()
for(rem in c("aer","tox","trap")){
  
  i<-which(c('aer','tox','trap')==rem)
  
  dir_files <- list.files(paste0(objdir,"/AKDE_",rem))
  if(rem%in%dir_files){
    dir_files <- dir_files[-which(dir_files==rem)]
  }

  akde_files[[i]]  <- data.frame(animalid=dir_files,
                           rem_typ=NA,
                           rem_typ_folder=rem)
  
  pigs_trt_temp <- pigs_trt %>% filter(Removal.Type==rem) %>% 
    dplyr::select(animalid) %>% distinct()
  
  akde_files[[i]] $rem_typ[!akde_files[[i]] $animalid%in%unique(pigs_trt_temp$animalid)] <- "ctrl"
  akde_files[[i]] $rem_typ[akde_files[[i]] $animalid%in%unique(pigs_trt_temp$animalid)] <- rem
}
akde_files <- do.call("rbind.data.frame",akde_files)

akde_files %>% ungroup() %>%  dplyr::group_by(factor(animalid)) %>% 
  dplyr::summarise(n=length(unique(rem_typ))) %>% 
  dplyr::summarise(length(unique(n)))

hr_pairs <- readRDS(paste0(objdir,"/pairs.rds"))
  
# rem_idx<-2
# per_idx<-2
# k<-1
# l<-2
# d<-1
# w<-1

#parallel function ----------
contact_from_ctmm <- function(rem_idx,per_idx){ # parallel over rem type FOLDERS (3)
  
  require(tidyverse)
  require(ctmm)
  # require(sf)
  require(Rcpp)
  require(RcppArmadillo)
  
  #source rcpp function
  Rcpp::sourceCpp(paste0(cpp_dir,"TestPairwise.cpp"), verbose=TRUE)
  
    # if(!dir.exists(paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Contact Analysis/ctmm Contacts/",rem_typ_folders[rem_idx]))){ #create
    #   dir.create(paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Contact Analysis/ctmm Contacts/",rem_typ_folders[rem_idx]))
    # }
    
    dist_mat <- list()
    
      if(!(rem_typ_folders[rem_idx]=="aer" & periods[per_idx]=="during")){ #no during aerial 
          
        ref_trt <- akde_files %>% filter(rem_typ_folder==rem_typ_folders[rem_idx])
        folder_trmts <- unique(ref_trt$rem_typ)

          for(w in 1:length(folder_trmts)){# ctrl vs trtmt
            
            pigs_trt_temp <- akde_files %>% 
              filter(rem_typ_folder==rem_typ_folders[rem_idx]) %>% 
              filter(rem_typ==folder_trmts[w])
            
            #distance matrix for number of contacts and number of individuals contacted (degree)
            #for each removal type and treatment period combo
            dist_mat[[w]] <- expand_grid(rem_typ=unique(pigs_trt_temp$rem_typ),
                                         period=periods[per_idx],
                                         animalid=unique(pigs_trt_temp$animalid),
                                         # week=unique(trt_week$week),
                                         dist=cdist,
                                         num_contacts=NA,
                                         num_indivs=NA)
            # for(k in 1:2){
            for(k in 1:nrow(pigs_trt_temp)){ # individual pig in treatment/ctrl group
              
              main_pig_id <- pigs_trt_temp$animalid[k]
              
              if(file.exists(paste0(ctmm_dir,
                                    rem_typ_folders[rem_idx],
                                    "/",main_pig_id,
                                    "/",main_pig_id,
                                    "_",periods[per_idx],"_one_min_pred.rds"))
                 ){
                
                main_pig_telem <- readRDS(paste0(ctmm_dir,
                                                 rem_typ_folders[rem_idx],"/",main_pig_id,"/",main_pig_id,"_",periods[per_idx],"_one_min_pred.rds"))
                main_pig_ctmm <- readRDS(paste0(objdir,"/AKDE_",
                                                rem_typ_folders[rem_idx],"/",main_pig_id,"/",main_pig_id,"_",
                                                periods[per_idx],".rds"))
                
                num_contacts_main <- matrix(NA,nrow=length(cdist),ncol=nrow(pigs_trt_temp))
                
                # for(l in 1:2){
                for(l in 1:nrow(pigs_trt_temp)){ # all pairs for indiv pig
                  
                  hr_trt_pairs <- hr_pairs %>% filter(id_1==main_pig_id & 
                                     rem==rem_typ_folders[rem_idx])
                  
                  pair_pig_id <- pigs_trt_temp$animalid[l]
                  
                  if(pair_pig_id!=main_pig_id &  #not same pig and
                     pair_pig_id%in%hr_trt_pairs$id_2){ #at least 1% HR overlap
                    
                    pair_pig_telem <- readRDS(paste0(ctmm_dir,
                                                     rem_typ_folders[rem_idx],
                                                     "/",pair_pig_id,
                                                     "/",pair_pig_id,
                                                     "_",periods[per_idx],"_one_min_pred.rds"))
                    
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
                  }  #not same pig / hr overlap if statement
                }  # all pairs for indiv pig
                
                for(d in 1:length(cdist)){ # distance - summing contacts and indivs per pig
                  
                  #number of contacts per week at distance
                  dist_mat[[w]]$num_contacts[dist_mat[[w]]$animalid==main_pig_id &
                                               # dist_mat[[per_idx]]$week==w &
                                               dist_mat[[w]]$dist==cdist[d]] <-
                    sum(num_contacts_main[d,],na.rm=T)
                  
                  #degree per week at distance
                  cd <- num_contacts_main[d,]
                  cd <- cd[!is.na(cd)]
                  dist_mat[[w]]$num_indivs[dist_mat[[w]]$animalid==main_pig_id &
                                             # dist_mat[[per_idx]]$week==w &
                                             dist_mat[[w]]$dist==cdist[d]] <-
                    sum(cd>0)
                  
                } # distance - summing contacts and indivs per pig

                # }# week in treatment period
                
              }#if ctmm exists statement
              
            }#individual pig 
            
          } #ctrl vs trtm within folder
        
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

cl <- makeCluster(nrow(perms), "SOCK")
clusterExport(cl, list("periods","rem_typ","cdist","pigs_trt",
                       "rem_typ_folders","perms","akde_files",
                       "cpp_dir","objdir","homedir","ctmm_dir",
                       "hr_pairs",'contact_from_ctmm'))

system.time(
  parContacts <- clusterMap(cl, fun=contact_from_ctmm, rem_idx=perms$rem_idx, per_idx=perms$per_idx)
)
stopCluster(cl)

parContacts[[4]]$dist_mat_all

contacts_all <- do.call("rbind.data.frame",lapply(1:nrow(perms),function(i)parContacts[[i]]$dist_mat_all))
contacts_all %>% filter(!is.na(num_contacts))
contacts_all

pig_sex_age <- georem %>% 
  group_by(animalid) %>% 
  dplyr::summarise(age=unique(age),
            sex=unique(sex))

contacts_all <- contacts_all %>% left_join(pig_sex_age) %>% 
  filter(is.na(sex))

saveRDS(contacts_all,file=paste0(objdir,"/pairwise_contacts.rds"))


#check
# parContacts <- readRDS(paste0(objdir,"/pairwise_contacts.rds"))
# con <- do.call("rbind.data.frame",parContacts)
# con$period <- factor(con$period,levels=c('before','during','after'))
# 
# pig_info <- geo.all %>% 
#   dplyr::select(animalid,Removal.Type,sex,age) %>% 
#   distinct()
# 
# ctrl <- pig_info %>% filter(Removal.Type=="ctrl")
# con$rem_typ[con$animalid%in%ctrl$animalid]
# 
# 
# ggplot(con)+
#   geom_jitter(aes(x=period,y=num_contacts,col=animalid),alpha=0.5,
#               width=0.1)+
#   facet_grid(dist~rem_typ)+
#   guides(col="none")+
#   theme(text=element_text(size=15))
# 
# ggplot(con)+
#   geom_jitter(aes(x=period,y=num_indivs,col=animalid),alpha=0.5,
#               width=0.1)+
#   facet_grid(dist~rem_typ)+
#   guides(col="none")+
#   theme(text=element_text(size=15))
