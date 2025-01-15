### Title: Fitting GLMMs to step length, speed, and contact NIFA metrics
### Author: Abbey Feuka
### Date: 16DEC24
### Notes: 

#input: georem_typ reference file, toxicant akde dataframe (outdf_akde_tox_corrected.rds),
# weekly distance and speed summaries (e.g. pig_weekly_distance_ctmm.rds),
# contact summaries (pairwise_contacts.RDS)

#output: formatted data used to fit models (e.g. distaer, speedtox),
# model outputs (e.g. res_distance_rp_aer.rds, res_speed_rps_tox.rds),
# summary plots (e.g. ./Disp_GLM_Results/res_distance_rp_aer_fig.png)

#setup ----------------

#load libraries
  library(glmmTMB)
  library(DHARMa)
  library(ctmm)
  library(fitdistrplus)
  library(ggplot2)
  library(plyr)
  library(dplyr)
  library(hrbrthemes)
  library(sdmTMB)

#sneaky homedir substitution
switch(Sys.info()[['sysname']],
       Windows= {homedir <- "C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Contact Analysis/Removals_Mvmt"
        },
       Linux  = {print("I'm a penguin.")},
       Darwin = {
         homedir <- "/Users/kayleigh.chalkowski/OneDrive/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline"
        }
       )

  #set directories
  #homedir <- "/Users/kayleigh.chalkowski/OneDrive/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline"
  # homedir <- "//aapcoftc3fp13/Projects/MUDD/ASF_NIFA/Pipelines/Removals_Mvmt"
  #homedir <- "C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Contact Analysis/Removals_Mvmt"
  objdir=file.path(homedir,"1_Data","Objects",fsep=.Platform$file.sep)
  results_dir <- file.path(homedir,"/3_Output/",fsep=.Platform$file.sep)
  
  #pull geo_rem objects
  # georem <- read.csv("./1_Data/Objects/geo_remtyp_period.csv")
  geo.aer<-readRDS(file.path(objdir,"geoaer.rds"))
  geo.tox<-readRDS(file.path(objdir,"geotox.rds"))
  geo.trap<-readRDS(file.path(objdir,"geotrap.rds"))
  geo.all <- rbind(geo.aer,geo.tox,geo.trap)
  
  #identify animals that died in toxicant treatment
  aktox=readRDS(paste0(objdir,"/outdf_akde_tox_corrected.rds"))
  #Remove NA after periods for tox-killed pigs
  aktox=aktox[!is.na(aktox$area.CI.low),]
  
  aktox.sumperiods=aktox %>% group_by(animalid,period) %>% dplyr::summarise(n()) %>% tidyr::pivot_wider(names_from="period",values_from=`n()`) %>% as.data.frame()
  died.tox=aktox.sumperiods[which(is.na(aktox.sumperiods$after)),1]
  results_dir <- paste0(homedir,"/1_Data/Objects/Models/")
  
  #set mesh cutoff for all spatial models
  mesh_cutoff<-1
  spatial_res <- 1000

# distance -------------------
## formatting distance dfs ----------
  #read distance
  dist<- readRDS(paste0(objdir,"/pig_weekly_distance_ctmm.rds"))

  #subset dfs
  distaer <- dist %>% filter(Removal.Type%in%c('aer')) %>% filter(removal.period.akdecalc!='during') 
  disttox <- dist %>% filter(Removal.Type%in%c('tox')) 
  disttrap <- dist %>% filter(Removal.Type%in%c('trap')) 
  
  #relevel when during isn't recorded
  distaer$removal.period.akdecalc <- factor(distaer$removal.period.akdecalc,
                                            levels=c('before','after'))
  disttox$removal.period.akdecalc <- factor(disttox$removal.period.akdecalc,
                                            levels=c('before','during','after'))
  disttrap$removal.period.akdecalc <- factor(disttrap$removal.period.akdecalc,
                                             levels=c('before','during','after'))
  

## temporal autocorrelation test --------------
  res_aer=glmmTMB(weekly_dist_km~(1|animalid)+trt_ctrl*removal.period.akdecalc,
                data=distaer,family=Gamma(link="log"))
simout_aer <- simulateResiduals(fittedModel = res_aer, plot = F)
res_aer2 = recalculateResiduals(simout_aer, group = distaer$week,rotation="estimated")
testTemporalAutocorrelation(res_aer2, time = unique(distaer$week))
#no temporal autocorrelation, p=0.054

#sex 
res_aer=glmmTMB(weekly_dist_km~(1|animalid)+trt_ctrl*removal.period.akdecalc*sex,
                data=distaer,family=Gamma(link="log"))
simout_aer <- simulateResiduals(fittedModel = res_aer, plot = F)
res_aer2 = recalculateResiduals(simout_aer, group = distaer$week,rotation="estimated")
testTemporalAutocorrelation(res_aer2, time = unique(distaer$week))
#no temporal autocorrelation, p=0.055

res_trap=glmmTMB(weekly_dist_km~(1|animalid)+trt_ctrl*removal.period.akdecalc,
                 data=disttrap,family=Gamma(link="log"))
simout_trap <- simulateResiduals(fittedModel = res_trap, plot = F)
res_trap2 = recalculateResiduals(simout_trap, group = disttrap$week,rotation="estimated")
testTemporalAutocorrelation(res_trap2, time = unique(disttrap$week))
#no temporal autocorrelation, p=0.29

#sex
res_trap=glmmTMB(weekly_dist_km~(1|animalid)+trt_ctrl*removal.period.akdecalc*sex,
                 data=disttrap,family=Gamma(link="log"))
simout_trap <- simulateResiduals(fittedModel = res_trap, plot = F)
res_trap2 = recalculateResiduals(simout_trap, group = disttrap$week,rotation="estimated")
testTemporalAutocorrelation(res_trap2, time = unique(disttrap$week))
#no temporal autocorrelation, p=0.29

res_tox=glmmTMB(weekly_dist_km~(1|animalid)+trt_ctrl*removal.period.akdecalc,
                data=disttox,family=Gamma(link="log"))
simout_tox <- simulateResiduals(fittedModel = res_tox, plot = F)
res_tox2 = recalculateResiduals(simout_tox, group = disttox$week,rotation="estimated")
testTemporalAutocorrelation(res_tox2, time = unique(disttox$week))
#no temporal autocorrelation, p=0.58

res_tox=glmmTMB(weekly_dist_km~(1|animalid)+trt_ctrl*removal.period.akdecalc*sex,
                data=disttox,family=Gamma(link="log"))
simout_tox <- simulateResiduals(fittedModel = res_tox, plot = F)
res_tox2 = recalculateResiduals(simout_tox, group = disttox$week,rotation="estimated")
testTemporalAutocorrelation(res_tox2, time = unique(disttox$week))
#no temporal autocorrelation, p=0.27

## spatial autocorrelation test ------------------
distaer$X=floor(distaer$mX/5000)
distaer$Y=floor(distaer$mY/5000)
aer_test=glmmTMB(weekly_dist_km~(1|animalid)+trt_ctrl*removal.period.akdecalc,
                 data=distaer,family=Gamma(link="log"))
aer_test_resid <- simulateResiduals(aer_test)
aer_groupLocations = aggregate(distaer[,c("mX","mY")], list(distaer$animalid), mean)
aer_test_resid2 = recalculateResiduals(aer_test_resid, group = distaer$animalid, rotation="estimated")
testSpatialAutocorrelation(aer_test_resid2,aer_groupLocations$mX, aer_groupLocations$mY)
# spatial autocorrelation, p=0.012

#sex
aer_test=glmmTMB(weekly_dist_km~(1|animalid)+trt_ctrl*removal.period.akdecalc*sex,
                 data=distaer,family=Gamma(link="log"))
aer_test_resid <- simulateResiduals(aer_test)
aer_test_resid2 = recalculateResiduals(aer_test_resid, group = distaer$animalid, rotation="estimated")
aer_groupLocations = aggregate(distaer[,c("mX","mY")], list(distaer$animalid, distaer$sex), mean)
DHARMa::testSpatialAutocorrelation(aer_test_resid2,aer_groupLocations$mX, aer_groupLocations$mY)
# no spatial autocorrelation, p=0.99

trap_test=glmmTMB(weekly_dist_km~(1|animalid)+trt_ctrl*removal.period.akdecalc,
                  data=disttrap,family=Gamma(link="log"))
trap_test_resid <- simulateResiduals(trap_test)
trap_groupLocations = aggregate(disttrap[,c("mX","mY")], list(disttrap$animalid), mean)
trap_test_resid2 = recalculateResiduals(trap_test_resid, group = disttrap$animalid, rotation="estimated")
testSpatialAutocorrelation(trap_test_resid2,trap_groupLocations$mX, trap_groupLocations$mY)
#no spatial autocorrelation, p=0.15

#sex 
trap_test=glmmTMB(weekly_dist_km~(1|animalid)+trt_ctrl*removal.period.akdecalc*sex,
                  data=disttrap,family=Gamma(link="log"))
trap_test_resid <- simulateResiduals(trap_test)
trap_test_resid2 = recalculateResiduals(trap_test_resid, group = disttrap$animalid, rotation="estimated")
testSpatialAutocorrelation(trap_test_resid2,trap_groupLocations$mX, trap_groupLocations$mY)
#no spatial autocorrelation, p=0.17

tox_test=glmmTMB(weekly_dist_km~(1|animalid)+trt_ctrl*removal.period.akdecalc,
                 data=disttox,family=Gamma(link="log"))
tox_test_resid <- simulateResiduals(tox_test)
tox_groupLocations = aggregate(disttox[,c("mX","mY")], list(disttox$animalid), mean)
tox_test_resid2 = recalculateResiduals(tox_test_resid, group = disttox$animalid, rotation="estimated")
testSpatialAutocorrelation(tox_test_resid2,tox_groupLocations$mX, tox_groupLocations$mY)
#spatial autocorrelation, p=0.0006

#sex
tox_test=glmmTMB(weekly_dist_km~(1|animalid)+trt_ctrl*removal.period.akdecalc*sex,
                 data=disttox,family=Gamma(link="log"))
tox_test_resid <- simulateResiduals(tox_test)
tox_test_resid2 = recalculateResiduals(tox_test_resid, group = disttox$animalid, rotation="estimated")
testSpatialAutocorrelation(tox_test_resid2,tox_groupLocations$mX, tox_groupLocations$mY)
#no spatial autocorrelation, p=0.49

## distance autocorrelation conclusions ----------------
  #Temporal autocorrelation: none
  #Spatial autocorrelation:
    #aerial removal*period model
    #tox removal*period model

## fit models ------------------

### aerial -----------------

#### removal type * period ------
spatial_res <- 100
mesh_cutoff=1
distaer$mX_sc <- floor(distaer$mX/spatial_res)
distaer$mY_sc <- floor(distaer$mY/spatial_res)
meshtox_sp <- make_mesh(distaer,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)
res_distance_rp_aer=sdmTMB(weekly_dist_km ~ (1|animalid) + trt_ctrl*removal.period.akdecalc,
                        data=distaer,
                        mesh=meshtox_sp,
                        spatial='on',
                        family=Gamma(link='log'))

sanity(res_distance_rp_aer)
tox_res_rp <- simulate(res_distance_rp_aer, nsim = 544, type = "mle-mvn") %>% 
  dharma_residuals(res_distance_rp_aer, return_DHARMa = TRUE)
tox_res_rp2 = recalculateResiduals(tox_res_rp, group = as.factor(distaer$animalid),rotation="estimated")
groupLocations = aggregate(distaer[,c("mX","mY")], list(distaer$animalid), mean)
testSpatialAutocorrelation(tox_res_rp2,groupLocations$mX,groupLocations$mY)
#removed spatial autocorrelation, p=0.931

#saveRDS(res_dist_rp_aer,paste0(results_dir,"res_distance_rp_aer.rds"))

#### removal type * period * sex ------

res_distance_rps_aer <-glmmTMB(weekly_dist_km ~ (1|animalid) +
                                 trt_ctrl*removal.period.akdecalc*sex,
                               data=distaer,family=Gamma(link="log"))

#saveRDS(res_distance_rps_aer,paste0(results_dir,"res_distance_rps_aer.rds"))

### trap -----------------

#### removal type * period ------

res_distance_rp_trap=glmmTMB(weekly_dist_km ~ trt_ctrl*removal.period.akdecalc+
                               (1|animalid),
                             data=disttrap,family=Gamma(link="log"))

#saveRDS(res_distance_rp_trap,paste0(results_dir,"res_distance_rp_trap.rds"))

#### removal type * period * sex ------
res_distance_rps_trap=glmmTMB(weekly_dist_km ~ trt_ctrl*removal.period.akdecalc*sex+
                                (1|animalid),
                              data=disttrap,family=Gamma(link="log"))

#saveRDS(res_distance_rps_trap,paste0(results_dir,"res_distance_rps_trap.rds"))

### toxicant -----------------

#### removal type * period ------
spatial_res <- 100
mesh_cutoff=2
disttox$mX_sc <- floor(disttox$mX/spatial_res)
disttox$mY_sc <- floor(disttox$mY/spatial_res)
meshtox_sp <- make_mesh(disttox,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)
res_distance_rp_tox=sdmTMB(weekly_dist_km ~ (1|animalid) + trt_ctrl*removal.period.akdecalc,
                       data=disttox,
                       mesh=meshtox_sp,
                       spatial='on',
                       family=Gamma(link='log'))

sanity(res_distance_rp_tox)
tox_res_rp <- simulate(res_distance_rp_tox, nsim = 544, type = "mle-mvn") %>% 
  dharma_residuals(res_distance_rp_tox, return_DHARMa = TRUE)
tox_res_rp2 = recalculateResiduals(tox_res_rp, group = as.factor(disttox$animalid),rotation="estimated")
groupLocations = aggregate(disttox[,c("mX","mY")], list(disttox$animalid), mean)
testSpatialAutocorrelation(tox_res_rp2,groupLocations$mX,groupLocations$mY)
#removed spatial autocorrelation, p=0.106

#saveRDS(res_distance_rp_tox,paste0(results_dir,"res_distance_rp_tox.rds"))

#### removal type * period * sex ------
res_distance_rps_tox=glmmTMB(weekly_dist_km ~ trt_ctrl*removal.period.akdecalc*sex+
                               (1|animalid),
                             data=disttox,family=Gamma(link="log"))

#saveRDS(res_distance_rps_tox,paste0(results_dir,"res_distance_rps_tox.rds"))

#speed -------------------

  speed<- readRDS(paste0(objdir,"/pig_weekly_speed_ctmm.rds"))
  
  speedaer <- speed %>% filter(Removal.Type%in%c('aer')) %>% 
    filter(removal.period.akdecalc!='during') 
  speedtox <- speed %>% filter(Removal.Type%in%c('tox')) 
  speedtrap <- speed %>% filter(Removal.Type%in%c('trap')) 
  
  #relevel when during isn't recorded
  speedaer$removal.period.akdecalc <- factor(speedaer$removal.period.akdecalc,
                                             levels=c('before','after'))
  speedtox$removal.period.akdecalc <- factor(speedtox$removal.period.akdecalc,
                                             levels=c('before','during','after'))
  speedtrap$removal.period.akdecalc <- factor(speedtrap$removal.period.akdecalc,
                                              levels=c('before','during','after'))

## spatial autocorrelation test ------------------
## aerial
aer_tests=glmmTMB(weekly_md_km_hr~(1|animalid)+trt_ctrl*removal.period.akdecalc,
                  data=speedaer,family=Gamma(link='log'))
aer_test_resids <- simulateResiduals(aer_tests)
aer_groupLocations_sp = aggregate(speedaer[,c("mX","mY")], list(speedaer$animalid), mean)
aer_test_resids2 = recalculateResiduals(aer_test_resids, group = speedaer$animalid, rotation="estimated")
testSpatialAutocorrelation(aer_test_resids2,aer_groupLocations_sp$mX, aer_groupLocations_sp$mY)
# no spatial autocorrelation

#sex
aer_tests=glmmTMB(weekly_md_km_hr~(1|animalid)+trt_ctrl*removal.period.akdecalc*sex,
                  data=speedaer,family=Gamma(link='log'))
aer_test_resids <- simulateResiduals(aer_tests)
aer_test_resids2 = recalculateResiduals(aer_test_resids, group = speedaer$animalid, rotation="estimated")
testSpatialAutocorrelation(aer_test_resids2,aer_groupLocations_sp$mX, aer_groupLocations_sp$mY)
# no spatial autocorrelation, p=0.93

trap_tests=glmmTMB(weekly_md_km_hr~(1|animalid)+trt_ctrl*removal.period.akdecalc,
                   data=speedtrap,family=Gamma(link='log'))
trap_test_resids <- simulateResiduals(trap_tests)
trap_groupLocations_sp = aggregate(speedtrap[,c("mX","mY")], list(speedtrap$animalid), mean)
trap_test_resids2 = recalculateResiduals(trap_test_resids, group = speedtrap$animalid, rotation="estimated")
testSpatialAutocorrelation(trap_test_resids2,trap_groupLocations_sp$mX, trap_groupLocations_sp$mY)
#no spatial autocorrelation, p=0.69

#sex
trap_tests=glmmTMB(weekly_md_km_hr~(1|animalid)+trt_ctrl*removal.period.akdecalc*sex,
                   data=speedtrap,family=Gamma(link='log'))
trap_test_resids <- simulateResiduals(trap_tests)
trap_test_resids2 = recalculateResiduals(trap_test_resids, group = speedtrap$animalid, rotation="estimated")
testSpatialAutocorrelation(trap_test_resids2,trap_groupLocations_sp$mX, trap_groupLocations_sp$mY)
#no spatial autocorrelation, p=0.86

tox_tests=glmmTMB(weekly_md_km_hr~(1|animalid)+trt_ctrl*removal.period.akdecalc,
                  data=speedtox,family=Gamma(link='log'))
tox_test_resids <- simulateResiduals(tox_tests)
tox_groupLocations_sp = aggregate(speedtox[,c("mX","mY")], list(speedtox$animalid), mean)
tox_test_resids2 = recalculateResiduals(tox_test_resids, group = speedtox$animalid, rotation="estimated")
testSpatialAutocorrelation(tox_test_resids2,tox_groupLocations_sp$mX, tox_groupLocations_sp$mY)
# spatial autocorrelation, p=1.079e-07

tox_tests=glmmTMB(weekly_md_km_hr~(1|animalid)+trt_ctrl*removal.period.akdecalc*sex,
                  data=speedtox,family=Gamma(link='log'))
tox_test_resids <- simulateResiduals(tox_tests)
tox_test_resids2 = recalculateResiduals(tox_test_resids, group = speedtox$animalid, rotation="estimated")
testSpatialAutocorrelation(tox_test_resids2,tox_groupLocations_sp$mX, tox_groupLocations_sp$mY)
# spatial autocorrelation, p=0.0001852

## speed autocorrelation conclusions --------
  #Spatial autocorrelation:
    #tox removal*period model
    #tox removal*period*sex model

## fit glmms ------------------

### aerial ----------------

#### removal type * period  ------
res_speed_rp_aer <- glmmTMB(weekly_md_km_hr ~ trt_ctrl*removal.period.akdecalc+
                                 (1|animalid),
                               data=speedaer,
                               family=Gamma(link="log"),
                               verbose=TRUE)

#saveRDS(res_speed_rp_aer,paste0(results_dir,"res_speed_rp_aer.rds"))

#### removal type * period * sex ------
res_speed_rps_aer <- glmmTMB(weekly_md_km_hr ~ trt_ctrl*removal.period.akdecalc*sex+
                              (1|animalid),
                            data=speedaer,
                            family=Gamma(link="log"),
                            verbose=TRUE)

#saveRDS(res_speed_rps_aer,paste0(results_dir,"res_speed_rps_aer.rds"))

### trap ----------------
#### removal type * period ------
res_speed_rp_trap=glmmTMB(weekly_md_km_hr~ 
                            (1|animalid)+
                            trt_ctrl*removal.period.akdecalc,
                          data=speedtrap,family=Gamma(link='log'))

#saveRDS(res_speed_rp_trap,paste0(results_dir,"res_speed_rp_trap.rds"))

#### removal type * period * sex ------
res_speed_rps_trap=glmmTMB(weekly_md_km_hr~ 
                             (1|animalid)+
                             trt_ctrl*removal.period.akdecalc*sex,
                           data=speedtrap,family=Gamma(link='log'))

#saveRDS(res_speed_rps_trap,paste0(results_dir,"res_speed_rps_trap.rds"))

### toxicant ----------------
#set mesh cutoff for spatial models
spatial_res <- 100
mesh_cutoff=1
speedtox$mX_sc <- floor(speedtox$mX/spatial_res)
speedtox$mY_sc <- floor(speedtox$mY/spatial_res)

#### removal type * period ------
meshtox_sp <- make_mesh(speedtox,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)
res_speed_rp_tox=sdmTMB(weekly_md_km_hr ~ (1|animalid) + trt_ctrl*removal.period.akdecalc,
                        data=speedtox,
                       mesh=meshtox_sp,
                        spatial='on',
                        family=Gamma(link='log'))

sanity(res_speed_rp_tox)
tox_res_rp <- simulate(res_speed_rp_tox, nsim = 544, type = "mle-mvn") %>% 
  dharma_residuals(res_speed_rp_tox, return_DHARMa = TRUE)
tox_res_rp2 = recalculateResiduals(tox_res_rp, group = as.factor(speedtox$animalid),rotation="estimated")
groupLocations = aggregate(speedtox[,c("mX","mY")], list(speedtox$animalid), mean)
testSpatialAutocorrelation(tox_res_rp2,groupLocations$mX,groupLocations$mY)
#removed spatial autocorrelation, p=0.5807

#saveRDS(res_speed_rp_tox,paste0(results_dir,"res_speed_rp_tox.rds"))

#### removal type * period * sex ------
meshtox_sp <- make_mesh(speedtox,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)
res_speed_rps_tox=sdmTMB(weekly_md_km_hr ~ (1|animalid) + trt_ctrl*removal.period.akdecalc*sex,
                        data=speedtox,
                        mesh=meshtox_sp,
                        spatial='on',
                        family=Gamma(link='log'))

sanity(res_speed_rp_tox)
tox_res_rp <- simulate(res_speed_rp_tox, nsim = 544, type = "mle-mvn") %>% 
  dharma_residuals(res_speed_rp_tox, return_DHARMa = TRUE)
tox_res_rp2 = recalculateResiduals(tox_res_rp, group = as.factor(speedtox$animalid),rotation="estimated")
groupLocations = aggregate(speedtox[,c("mX","mY")], list(speedtox$animalid), mean)
testSpatialAutocorrelation(tox_res_rp2,groupLocations$mX,groupLocations$mY)
#removed spatial autocorrelation, p=0.6789

#saveRDS(res_speed_rps_tox,paste0(results_dir,"res_speed_rps_tox.rds"))

# Format model info -----------------------------------------------------------
mods<-list(res_distance_rp_aer,
           res_distance_rps_aer,
           res_distance_rp_tox,
           res_distance_rps_tox,
           res_distance_rp_trap,
           res_distance_rps_trap,
           res_speed_rp_aer,
           res_speed_rps_aer,
           res_speed_rp_tox,
           res_speed_rps_tox,
           res_speed_rp_trap,
           res_speed_rps_trap)

models=c("res_distance_rp_aer",
            "res_distance_rps_aer",
            "res_distance_rp_tox",
            "res_distance_rps_tox",
            "res_distance_rp_trap",
            "res_distance_rps_trap",
            "res_speed_rp_aer",
            "res_speed_rps_aer",
            "res_speed_rp_tox",
            "res_speed_rps_tox",
            "res_speed_rp_trap",
            "res_speed_rps_trap")

# * Format model prediction df -----------
for(i in 1:length(mods)){
  #i=1
  #i=8
  #if(class(mods[[i]])=="glmmTMB"){
  #coefs=summary(mods[[i]])$coef$cond
  if(length(grep("rps",models[i]))==0){
  tmp=as.data.frame(ggeffects::predict_response(mods[[i]], terms=c("trt_ctrl","removal.period.akdecalc")))
    tmp$facet=NA
  } else{
    tmp=as.data.frame(ggeffects::predict_response(mods[[i]], terms=c("trt_ctrl","removal.period.akdecalc","sex")))
  }
  
  #}
  
  tmp$model=models[i]
  
  if(i==1){
    preds=tmp
  } else{
    preds=rbind(preds,tmp)
  }
  
}

#make colnames more informative
colnames(preds)[c(1,6,7)]<-c("trt","per","sex")

#Make rem column based on model
preds$rem=NA
preds$rem[grep("aer",preds$model)]<-"aer"
preds$rem[grep("trap",preds$model)]<-"trap"
preds$rem[grep("tox",preds$model)]<-"tox"

#make response column
preds$response=NA
preds$response[grep("speed",preds$model)]<-"speed"
preds$response[grep("distance",preds$model)]<-"distance"

# Pull interactions table ------------------------------------------------------

for(i in 1:length(mods)){
  
  if(class(mods[[i]])=="sdmTMB"){
  df=broom.mixed::tidy(mods[[i]]) %>% as.data.frame()
  z=df$estimate/df$std.error
  df$p=exp(-0.717*z-0.416*z^2)
  df$model=models[i]
  df=df[,c(2,3,4,5,1)]
  df=as.data.frame(df[grep(":",df$term),])
  colnames(df)=c("Estimate","Std. Error","Pr(>|z|)","model","effect")
  
  } else{
    #Pull coefficients
    coefs=summary(mods[[i]])$coef$cond
    
    #Make table with just interaction effect, std.error, p value
    coefs2=as.data.frame(coefs[grep(":",rownames(coefs)),c(1,2,4),drop=FALSE])
    if(nrow(coefs2)!=0){
      
      if(nrow(coefs2)>1){
        coefs2=as.data.frame(coefs2[grep("trt_ctrl",rownames(coefs2)),,drop=FALSE])
        coefs2=as.data.frame(coefs2[grep("period",rownames(coefs2)),,drop=FALSE])
      }
      
      
      coefs2$model=models[i]
      coefs2$effect=rownames(coefs2)
      rownames(coefs2)=NULL
      
      #rename to rbind with sdmTMB results
      df=coefs2
    }
  }
  
  if(i==1){
    allc=df
  } else{
    allc=rbind(allc,df)
  }
}

#check that all captured
all(models%in%allc$model)
#which(!(models%in%allc$model))

#period
allc$period=NA
allc$period[grep("after",allc$effect)]<-"after"
allc$period[grep("during",allc$effect)]<-"during"

#remove NA for period, means is intxn between treatment and sex alone
allc=allc[!is.na(allc$period),]

#sex
allc$sex=NA
allc$sex[grep("rps",allc$model)]<-"female"
allc$sex[grep("sexMale",allc$effect)]<-"male"

#get treatment type
allc$trt=NA
allc$trt[grep("aer",allc$model)]<-"aer"
allc$trt[grep("trap",allc$model)]<-"trap"
allc$trt[grep("tox",allc$model)]<-"tox"

#get response type
allc$response=NA
allc$response[grep("speed",allc$model)]<-"speed"
allc$response[grep("distance",allc$model)]<-"distance"

#Remove ugly effect column
allc=allc[,-which(colnames(allc)=="effect")]

#change sex NA
allc$sex[is.na(allc$sex)]<-"whole"

#subset to only significant interactions
allc$alpha=0.1
allc$alpha[allc$`Pr(>|z|)`<0.05]<-1

# Combine full parameter tables -----------------------------------------------------------
for(i in 1:length(mods)){
  if(class(mods[[i]])=="sdmTMB"){
  df=broom.mixed::tidy(mods[[i]])
  df$statistic=df$estimate/df$std.error
  df$p.value=exp(-0.717*df$statistic-0.416*df$statistic^2)
  df<-as.data.frame(df)
  }
  
  #add statistic, p value
  if(class(mods[[i]])=="glmmTMB"){
    df=broom.mixed::tidy(mods[[i]])
    #df$model=models[i]
    df=df[df$effect=="fixed",]
    #term, estimate, std error, statistic, p value
    df=df[,c("term","estimate","std.error","statistic","p.value")]
    df<-as.data.frame(df)
    
  }
  
  df$model=models[i]
  
  if(i==1){
    parms=df
  } else{
    parms=rbind(parms,df)
  }
  
}

# Save tidied model output ----------------
outdir=file.path(homedir,"3_Output",fsep=.Platform$file.sep)
if(!dir.exists(file.path(outdir,"Model_Output"))){dir.create(file.path(outdir,"Model_Output"))}
saveRDS(preds,file.path(outdir,"Model_Output","spdist_preds.rds"))
saveRDS(allc,file.path(outdir,"Model_Output","spdist_intxns.rds"))
saveRDS(parms,file.path(outdir,"Model_Output","spdist_full_param.rds"))

# contacts - number --------
## aerial ------------
conaer <- readRDS(paste0(objdir,"/pairwise_contacts_aer.rds"))

conaer$Removal.Type <- 
  factor(conaer$trt_typ,levels=c('ctrl','trt'))
conaer$removal.period.akdecalc <- 
  factor(conaer$removal.period.akdecalc,levels=c('before','after'))

#distance thresholds
ggplot(conaer)+
  geom_histogram(aes(x=contacts_per_day))+
  facet_wrap(.~dist)

conaer <- conaer %>% filter(dist==10)

conaer <- conaer %>% left_join(
  geo.aer %>% 
    group_by(animalid,removal.period.akdecalc) %>% 
    summarise(mX=mean(X),
              mY=mean(Y)))

conaer$animalid <- 
  factor(conaer$animalid)

conaer$contacts_per_day_offset <- conaer$contacts_per_day+0.0001

conaer %>% group_by(contacts_per_day==0) %>% summarise(n=n())

## removal type * period -----
res_ncon_rp_aer=glmmTMB(contacts_per_day_offset~(1|animalid)+
                          Removal.Type*removal.period.akdecalc,
                        data=conaer,
                        family=Gamma(link="log")
)

aer_res_ncon_rp <- simulateResiduals(res_ncon_rp_aer)
aer_res_ncon_rp2 = recalculateResiduals(aer_res_ncon_rp, 
                                        group = as.factor(conaer$animalid), 
                                        rotation="estimated")
plot(aer_res_ncon_rp2)
plot(aer_res_ncon_rp2$simulatedResponse)
descdist(aer_res_ncon_rp2$fittedResiduals)
DHARMa::testDispersion(aer_res_ncon_rp2)

#test spatial autocorrelation
aer_groupLocations_con = aggregate(conaer[,c("mX","mY")],
                                    list(conaer$animalid),
                                    mean)
testSpatialAutocorrelation(aer_res_ncon_rp2,
                           aer_groupLocations_con$mX, 
                           aer_groupLocations_con$mY)
#spatially autocorrelated p-value = 0.002235

#correcting spatial autocorrelation
spatial_res <- 1000
mesh_cutoff=0.5
conaer$mX_sc <- floor(conaer$mX/spatial_res)
conaer$mY_sc <- floor(conaer$mY/spatial_res)
meshaer_con <- make_mesh(conaer,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)

res_ncon_rp_aer_sf = sdmTMB(contacts_per_day_offset~
                              (1|animalid)+
                               Removal.Type*removal.period.akdecalc,
                             data=conaer,
                             family=Gamma(link="log"),
                             mesh=meshaer_con,
                             spatial="on")

sanity(res_ncon_rp_aer_sf)
aer_res_rp_sf <- simulate(res_ncon_rp_aer_sf, 
                          nsim = 544,type = 'mle-mvn') %>%
  dharma_residuals(res_ncon_rp_aer_sf, return_DHARMa = TRUE)
aer_res_rp_sf2 = recalculateResiduals(aer_res_rp_sf,
                                      group = as.factor(conaer$animalid),
                                      rotation="estimated")
testSpatialAutocorrelation(aer_res_rp_sf2,
                           aer_groupLocations_con$mX,
                           aer_groupLocations_con$mY)

saveRDS(res_ncon_rp_aer_sf,paste0(results_dir,"res_ncon_rp_aer.rds"))

## removal type * period *sex -----
res_ncon_rps_aer=glmmTMB(contacts_per_day_offset~
                           (1|animalid)+
                           Removal.Type*removal.period.akdecalc*sex,
                         data=conaer,
                         family=Gamma(link="log")
)

aer_res_ncon_rps <- simulateResiduals(res_ncon_rps_aer)
aer_res_ncon_rps2 = recalculateResiduals(aer_res_ncon_rps, 
                                        group = as.factor(conaer$animalid), 
                                        rotation="estimated")
plot(aer_res_ncon_rps2)
plot(aer_res_ncon_rps2$simulatedResponse)
descdist(aer_res_ncon_rps2$fittedResiduals)
DHARMa::testDispersion(aer_res_ncon_rps2)

#test spatial autocorrelation
testSpatialAutocorrelation(aer_res_ncon_rps2,
                           aer_groupLocations_con$mX, 
                           aer_groupLocations_con$mY)
# spatial autocorrelation p=0.02

#correcting for spatial autocorrelation
### NOTE: model overparameterized, very poor fit
spatial_res <- 1000
mesh_cutoff=4
conaer$mX_sc <- floor(conaer$mX/spatial_res)
conaer$mY_sc <- floor(conaer$mY/spatial_res)
meshaer_con <- make_mesh(conaer,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)

res_ncon_rps_aer_sf = sdmTMB(contacts_per_day_offset~
                               (1|animalid)+
                              Removal.Type*removal.period.akdecalc*sex,
                            data=conaer,
                            family=Gamma(link="log"),
                            mesh=meshaer_con,
                            spatial="on")

sanity(res_ncon_rps_aer_sf)
aer_res_rps_sf <- simulate(res_ncon_rps_aer_sf,
                          nsim = 544,type = 'mle-mvn') %>%
  dharma_residuals(res_ncon_rp_aer_sf, return_DHARMa = TRUE)
aer_res_rps_sf2 = recalculateResiduals(aer_res_rps_sf,
                                      group = as.factor(conaer$animalid),
                                      rotation="estimated")
testSpatialAutocorrelation(aer_res_rps_sf2,
                           aer_groupLocations_con$mX,
                           aer_groupLocations_con$mY)

saveRDS(res_ncon_rps_aer,paste0(results_dir,"res_ncon_rps_aer.rds"))

#trap ------
contrap <- readRDS(paste0(objdir,"/pairwise_contacts_trap.rds"))

contrap$Removal.Type <- 
  factor(contrap$trt_typ,levels=c('ctrl','trt'))
contrap$removal.period.akdecalc <- 
  factor(contrap$removal.period.akdecalc,levels=c('before','during','after'))

#distance thresholds
ggplot(contrap)+
  geom_histogram(aes(x=contacts_per_day))+
  facet_wrap(.~dist)

contrap <- contrap %>% filter(dist==10)

contrap <- contrap %>% left_join(
  geo.trap %>% 
    group_by(animalid) %>% 
    summarise(mX=mean(X),
              mY=mean(Y)))

contrap$animalid <- 
  factor(contrap$animalid)

contrap$contacts_per_day_offset <- contrap$contacts_per_day+0.0001

contrap %>% group_by(contacts_per_day==0) %>% summarise(n=n())

## removal type * period -----
res_ncon_rp_trap=glmmTMB(contacts_per_day_offset~(1|animalid)+
                           Removal.Type*removal.period.akdecalc,
                         data=contrap,
                         family=Gamma(link="log")
                         )
trap_res_ncon_rp <- simulateResiduals(res_ncon_rp_trap)
trap_res_ncon_rp2 = recalculateResiduals(trap_res_ncon_rp, 
                                        group = as.factor(contrap$animalid), 
                                        rotation="estimated")
plot(trap_res_ncon_rp2)
plot(trap_res_ncon_rp2$simulatedResponse)
descdist(trap_res_ncon_rp2$fittedResiduals)
DHARMa::testDispersion(trap_res_ncon_rp2)

#test spatial autocorrelation
trap_groupLocations_con = aggregate(contrap[,c("mX","mY")], 
                                    list(contrap$animalid), 
                                    mean)
testSpatialAutocorrelation(trap_res_ncon_rp2,
                           trap_groupLocations_con$mX, 
                           trap_groupLocations_con$mY)
#spatial autocorrelation p=0.02

#correcting for spatial autocorrelation
spatial_res <- 1000
mesh_cutoff=0.5
contrap$mX_sc <- floor(contrap$mX/spatial_res)
contrap$mY_sc <- floor(contrap$mY/spatial_res)
meshtrap_con <- make_mesh(contrap,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)

res_ncon_rp_trap_sf = sdmTMB(contacts_per_day_offset~(1|animalid)+
                              Removal.Type*removal.period.akdecalc*sex,
                            data=contrap,
                            family=Gamma(link="log"),
                            mesh=meshtrap_con,
                            spatial="on")

sanity(res_ncon_rp_trap_sf)
trap_res_rp_sf <- simulate(res_ncon_rp_trap_sf,
                          nsim = 544,type = 'mle-mvn') %>%
  dharma_residuals(res_ncon_rp_trap_sf, return_DHARMa = TRUE)
trap_res_rp_sf2 = recalculateResiduals(trap_res_rp_sf,
                                      group = as.factor(contrap$animalid),
                                      rotation="estimated")
testSpatialAutocorrelation(trap_res_rp_sf2,
                           trap_groupLocations_con$mX,
                           trap_groupLocations_con$mY)
#spatial autocorrelation p=0.19

saveRDS(res_ncon_rp_trap_sf,paste0(results_dir,"res_ncon_rp_trap.rds"))

## removal type * period *sex -----
res_ncon_rps_trap=glmmTMB(contacts_per_day_offset~
                            (1|animalid)+
                            Removal.Type*removal.period.akdecalc*sex,
                          data=contrap,
                          family=Gamma(link="log")
                          )
trap_res_ncon_rps <- simulateResiduals(res_ncon_rps_trap)
trap_res_ncon_rps2 = recalculateResiduals(trap_res_ncon_rps, 
                                         group = as.factor(contrap$animalid), 
                                         rotation="estimated")
plot(trap_res_ncon_rps2)
plot(trap_res_ncon_rps2$simulatedResponse)
descdist(trap_res_ncon_rps2$fittedResiduals)
DHARMa::testDispersion(trap_res_ncon_rps2)

#test spatial autocorrelation
testSpatialAutocorrelation(trap_res_ncon_rps2,
                           trap_groupLocations_con$mX, 
                           trap_groupLocations_con$mY)
#no spatial autocorrelation p=0.08

saveRDS(res_ncon_rps_trap,paste0(results_dir,"res_ncon_rps_trap.rds"))

#tox ------
contox <- readRDS(paste0(objdir,"/pairwise_contacts_tox.rds"))

contox$Removal.Type <- 
  factor(contox$trt_typ,levels=c('ctrl','trt'))
contox$removal.period.akdecalc <- 
  factor(contox$removal.period.akdecalc,levels=c('before','during','after'))

#distance thresholds
ggplot(contox)+
  geom_histogram(aes(x=contacts_per_day))+
  facet_wrap(.~dist)

contox <- contox %>% filter(dist==10)

contox <- contox %>% left_join(
  geo.tox %>% 
    group_by(animalid) %>% 
    summarise(mX=mean(X),
              mY=mean(Y)))

contox$animalid <- 
  factor(contox$animalid)

contox$contacts_per_day_offset <- contox$contacts_per_day+0.0001

#remove animals that died
contox <- contox %>% filter(!is.na(num_contacts))

contox %>% group_by(contacts_per_day==0) %>% summarise(n=n())

## removal type * period -----
res_ncon_rp_tox=glmmTMB(contacts_per_day_offset~
                          (1|animalid)+
                          Removal.Type*removal.period.akdecalc,
                        data=contox,
                        family=Gamma(link="log")
                        )
tox_res_ncon_rp <- simulateResiduals(res_ncon_rp_tox)
tox_res_ncon_rp2 = recalculateResiduals(tox_res_ncon_rp, 
                                        group = as.factor(contox$animalid), 
                                        rotation="estimated")
plot(tox_res_ncon_rp2)
plot(tox_res_ncon_rp2$simulatedResponse)
descdist(tox_res_ncon_rp2$fittedResiduals)
DHARMa::testDispersion(tox_res_ncon_rp2)

#test spatial autocorrelation
tox_groupLocations_con = aggregate(contox[,c("mX","mY")], 
                                   list(contox$animalid), 
                                   mean)
testSpatialAutocorrelation(tox_res_ncon_rp2,
                           tox_groupLocations_con$mX, 
                           tox_groupLocations_con$mY)
#no spatial autocorrelation p=0.77

saveRDS(res_ncon_rp_tox,paste0(results_dir,"res_ncon_rp_tox.rds"))

## removal type * period *sex -----
res_ncon_rps_tox=glmmTMB(contacts_per_day_offset~
                           (1|animalid)+
                           Removal.Type*removal.period.akdecalc,
                         data=contox,
                         family=Gamma(link="log")
                         )

tox_res_ncon_rps <- simulateResiduals(res_ncon_rps_tox)
tox_res_ncon_rps2 = recalculateResiduals(tox_res_ncon_rps, 
                                         group = as.factor(contox$animalid), 
                                         rotation="estimated")
plot(tox_res_ncon_rps2)
plot(tox_res_ncon_rps2$simulatedResponse)
descdist(tox_res_ncon_rps2$fittedResiduals)
DHARMa::testDispersion(tox_res_ncon_rps2)

#test spatial autocorrelation
testSpatialAutocorrelation(tox_res_ncon_rps2,
                           tox_groupLocations_con$mX, 
                           tox_groupLocations_con$mY)
#no spatial autocorrelation p=0.77

saveRDS(res_ncon_rps_tox,paste0(results_dir,"res_ncon_rps_tox.rds"))


# contacts - degree --------
ggplot(conaer)+
  geom_histogram(aes(x=indivs_per_day))+
  facet_wrap(.~dist)

conaer$indivs_per_day_offset <- conaer$indivs_per_day+0.0001

#aerial ------
## removal type * period -----
res_nind_rp_aer=glmmTMB(indivs_per_day_offset~
                          (1|animalid)+
                          Removal.Type*removal.period.akdecalc,
                        data=conaer,
                        family=Gamma(link="log")
                        )

aer_res_nind_rp <- simulateResiduals(res_nind_rp_aer)
aer_res_nind_rp2 = recalculateResiduals(aer_res_nind_rp, 
                                        group = as.factor(conaer$animalid), 
                                        rotation="estimated")
plot(aer_res_nind_rp2)
plot(aer_res_nind_rp2$simulatedResponse)
descdist(aer_res_nind_rp2$fittedResiduals)
DHARMa::testDispersion(aer_res_nind_rp2)

#test spatial autocorrelation
testSpatialAutocorrelation(aer_res_nind_rp2,
                           aer_groupLocations_con$mX, 
                           aer_groupLocations_con$mY)
#no spatial autocorrelation p=0.09

saveRDS(res_nind_rp_aer,paste0(results_dir,"res_nind_rp_aer.rds"))

## removal type * period *sex -----
res_nind_rps_aer=glmmTMB(indivs_per_day_offset~
                           (1|animalid)+
                           Removal.Type*removal.period.akdecalc*sex,
                         data=conaer,
                         family=Gamma(link="log")
                         )

aer_res_nind_rps <- simulateResiduals(res_nind_rps_aer)
aer_res_nind_rps2 = recalculateResiduals(aer_res_nind_rps, 
                                         group = as.factor(conaer$animalid), 
                                         rotation="estimated")
plot(aer_res_nind_rps2)
plot(aer_res_nind_rps2$simulatedResponse)
descdist(aer_res_nind_rps2$fittedResiduals)
DHARMa::testDispersion(aer_res_nind_rps2)

#test spatial autocorrelation
testSpatialAutocorrelation(aer_res_nind_rps2,
                           aer_groupLocations_con$mX, 
                           aer_groupLocations_con$mY)
#no spatial autocorrelation p=0.06

saveRDS(res_nind_rps_aer,paste0(results_dir,"res_nind_rps_aer.rds"))

#trap ------
ggplot(contrap)+
  geom_histogram(aes(x=indivs_per_day))+
  facet_wrap(.~dist)

contrap$indivs_per_day_offset <- contrap$indivs_per_day+0.0001

## removal type * period -----
res_nind_rp_trap=glmmTMB(indivs_per_day_offset~
                           (1|animalid)+
                           Removal.Type*removal.period.akdecalc,
                         data=contrap,
                         family=Gamma(link="log")
                         )

trap_res_nind_rp <- simulateResiduals(res_nind_rp_trap)
trap_res_nind_rp2 = recalculateResiduals(trap_res_nind_rp, 
                                         group = as.factor(contrap$animalid), 
                                         rotation="estimated")
plot(trap_res_nind_rp2)
plot(trap_res_nind_rp2$simulatedResponse)
descdist(trap_res_nind_rp2$fittedResiduals)
DHARMa::testDispersion(trap_res_nind_rp2)

#test spatial autocorrelation
testSpatialAutocorrelation(trap_res_nind_rp2,
                           trap_groupLocations_con$mX, 
                           trap_groupLocations_con$mY)
#spatial autocorrelation p=0.02

#correcting for spatial autocorrelation
spatial_res <- 1000
mesh_cutoff=0.5
contrap$mX_sc <- floor(contrap$mX/spatial_res)
contrap$mY_sc <- floor(contrap$mY/spatial_res)
meshtrap_con <- make_mesh(contrap,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)

res_nind_rp_trap_sf = sdmTMB(indivs_per_day_offset~
                               (1|animalid)+
                               Removal.Type*removal.period.akdecalc,
                             data=contrap,
                             family=Gamma(link="log"),
                             mesh=meshtrap_con,
                             spatial="on")

sanity(res_nind_rp_trap_sf)
trap_res_nind_rp_sf <- simulate(res_nind_rp_trap_sf,
                           nsim = 544,type = 'mle-mvn') %>%
  dharma_residuals(res_nind_rp_trap_sf, return_DHARMa = TRUE)
trap_res_nind_rp_sf2 = recalculateResiduals(trap_res_nind_rp_sf,
                                       group = as.factor(contrap$animalid),
                                       rotation="estimated")
testSpatialAutocorrelation(trap_res_nind_rp_sf2,
                           trap_groupLocations_con$mX,
                           trap_groupLocations_con$mY)
#no spatial autocorrelation p=0.08

saveRDS(res_nind_rp_trap_sf,paste0(results_dir,"res_nind_rp_trap.rds"))

## removal type * period *sex -----
res_nind_rps_trap=glmmTMB(indivs_per_day_offset~
                            (1|animalid)+
                            Removal.Type*removal.period.akdecalc*sex,
                          data=contrap,
                          family=Gamma(link="log")
                          )

trap_res_nind_rps <- simulateResiduals(res_nind_rps_trap)
trap_res_nind_rps2 = recalculateResiduals(trap_res_nind_rps, 
                                          group = as.factor(contrap$animalid), 
                                          rotation="estimated")
plot(trap_res_nind_rps2)
plot(trap_res_nind_rps2$simulatedResponse)
descdist(trap_res_nind_rps2$fittedResiduals)
DHARMa::testDispersion(trap_res_nind_rps2)

#test spatial autocorrelation
testSpatialAutocorrelation(trap_res_nind_rps2,
                           trap_groupLocations_con$mX, 
                           trap_groupLocations_con$mY)
#no patial autocorrelation p=0.06

saveRDS(res_nind_rps_trap,paste0(results_dir,"res_nind_rps_trap.rds"))

#tox ------
ggplot(contox)+
  geom_histogram(aes(x=indivs_per_day))+
  facet_wrap(.~dist)

contox$indivs_per_day_offset <- contox$indivs_per_day+0.0001

## removal type * period -----
res_nind_rp_tox=glmmTMB(indivs_per_day_offset~
                          (1|animalid)+
                          Removal.Type*removal.period.akdecalc,
                        data=contox,
                        family=Gamma(link="log")
                        )
tox_res_nind_rp <- simulateResiduals(res_nind_rp_tox)
tox_res_nind_rp2 = recalculateResiduals(tox_res_nind_rp, 
                                        group = as.factor(contox$animalid), 
                                        rotation="estimated")
plot(tox_res_nind_rp2)
plot(tox_res_nind_rp2$simulatedResponse)
descdist(tox_res_nind_rp2$fittedResiduals)
DHARMa::testDispersion(tox_res_nind_rp2)

#test spatial autocorrelation
testSpatialAutocorrelation(tox_res_nind_rp2,
                           tox_groupLocations_con$mX, 
                           tox_groupLocations_con$mY)
#spatial autocorrelation p=0.000005

#correcting for spatial autocorrelation
spatial_res <- 1000
mesh_cutoff=1
contox$mX_sc <- floor(contox$mX/spatial_res)
contox$mY_sc <- floor(contox$mY/spatial_res)
meshtox_con <- make_mesh(contox,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)

res_nind_rp_tox_sf = sdmTMB(indivs_per_day_offset~
                               (1|animalid)+
                               Removal.Type*removal.period.akdecalc,
                             data=contox,
                             family=Gamma(link="log"),
                             mesh=meshtox_con,
                             spatial="on")

sanity(res_nind_rp_tox_sf)
tox_res_nind_rp_sf <- simulate(res_nind_rp_tox_sf,
                                nsim = 544,type = 'mle-mvn') %>%
  dharma_residuals(res_nind_rp_tox_sf, return_DHARMa = TRUE)
tox_res_nind_rp_sf2 = recalculateResiduals(tox_res_nind_rp_sf,
                                            group = as.factor(contox$animalid),
                                            rotation="estimated")
testSpatialAutocorrelation(tox_res_nind_rp_sf2,
                           tox_groupLocations_con$mX,
                           tox_groupLocations_con$mY)
#no spatial autocorrelation p=0.13

saveRDS(res_nind_rp_tox_sf,paste0(results_dir,"res_nind_rp_tox.rds"))

## removal type * period *sex -----
res_nind_rps_tox=glmmTMB(indivs_per_day_offset~
                           (1|animalid)+
                           Removal.Type*removal.period.akdecalc*sex,
                         data=contox,
                         family=Gamma(link="log")
                         )

tox_res_nind_rps <- simulateResiduals(res_nind_rps_tox)
tox_res_nind_rps2 = recalculateResiduals(tox_res_nind_rps, 
                                         group = as.factor(contox$animalid), 
                                         rotation="estimated")
plot(tox_res_nind_rps2)
plot(tox_res_nind_rps2$simulatedResponse)
descdist(tox_res_nind_rps2$fittedResiduals)
DHARMa::testDispersion(tox_res_nind_rps2)

#test spatial autocorrelation
testSpatialAutocorrelation(tox_res_nind_rps2,
                           tox_groupLocations_con$mX, 
                           tox_groupLocations_con$mY)
#spatial autocorrelation p=0.00007

#correcting for spatial autocorrelation
# spatial_res <- 1000
# mesh_cutoff=1
# contox$mX_sc <- floor(contox$mX/spatial_res)
# contox$mY_sc <- floor(contox$mY/spatial_res)
# meshtox_con <- make_mesh(contox,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)

res_nind_rps_tox_sf = sdmTMB(indivs_per_day_offset~
                              (1|animalid)+
                              Removal.Type*removal.period.akdecalc*sex,
                            data=contox,
                            family=Gamma(link="log"),
                            mesh=meshtox_con,
                            spatial="on")

sanity(res_nind_rps_tox_sf)
tox_res_nind_rps_sf <- simulate(res_nind_rps_tox_sf,
                               nsim = 544,type = 'mle-mvn') %>%
  dharma_residuals(res_nind_rps_tox_sf, return_DHARMa = TRUE)
tox_res_nind_rps_sf2 = recalculateResiduals(tox_res_nind_rps_sf,
                                           group = as.factor(contox$animalid),
                                           rotation="estimated")
testSpatialAutocorrelation(tox_res_nind_rps_sf2,
                           tox_groupLocations_con$mX,
                           tox_groupLocations_con$mY)
#no spatial autocorrelation p=0.08

saveRDS(res_nind_rps_tox_sf,paste0(results_dir,"res_nind_rps_tox.rds"))

