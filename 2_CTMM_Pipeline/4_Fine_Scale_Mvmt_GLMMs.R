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

  #set directories
   homedir <- "/Users/kayleigh.chalkowski/OneDrive/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline"
  # homedir <- "//aapcoftc3fp13/Projects/MUDD/ASF_NIFA/Pipelines/Removals_Mvmt"
  # homedir <- "C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Contact Analysis/Removals_Mvmt"
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
res_dist_rp_aer=sdmTMB(weekly_dist_km ~ (1|animalid) + trt_ctrl*removal.period.akdecalc,
                        data=distaer,
                        mesh=meshtox_sp,
                        spatial='on',
                        family=Gamma(link='log'))

sanity(res_dist_rp_aer)
tox_res_rp <- simulate(res_dist_rp_aer, nsim = 544, type = "mle-mvn") %>% 
  dharma_residuals(res_dist_rp_aer, return_DHARMa = TRUE)
tox_res_rp2 = recalculateResiduals(tox_res_rp, group = as.factor(distaer$animalid),rotation="estimated")
groupLocations = aggregate(distaer[,c("mX","mY")], list(distaer$animalid), mean)
testSpatialAutocorrelation(tox_res_rp2,groupLocations$mX,groupLocations$mY)
#removed spatial autocorrelation, p=0.931

saveRDS(res_dist_rp_aer,paste0(results_dir,"res_distance_rp_aer.rds"))

#### removal type * period * sex ------

res_distance_rps_aer <-glmmTMB(weekly_dist_km ~ (1|animalid) +
                                 trt_ctrl*removal.period.akdecalc*sex,
                               data=distaer,family=Gamma(link="log"))

saveRDS(res_distance_rps_aer,paste0(results_dir,"res_distance_rps_aer.rds"))

### trap -----------------

#### removal type * period ------

res_distance_rp_trap=glmmTMB(weekly_dist_km ~ trt_ctrl*removal.period.akdecalc+
                               (1|animalid),
                             data=disttrap,family=Gamma(link="log"))

saveRDS(res_distance_rp_trap,paste0(results_dir,"res_distance_rp_trap.rds"))

#### removal type * period * sex ------
res_distance_rps_trap=glmmTMB(weekly_dist_km ~ trt_ctrl*removal.period.akdecalc*sex+
                                (1|animalid),
                              data=disttrap,family=Gamma(link="log"))

saveRDS(res_distance_rps_trap,paste0(results_dir,"res_distance_rps_trap.rds"))

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

sanity(res_dist_rp_aer)
tox_res_rp <- simulate(res_distance_rp_tox, nsim = 544, type = "mle-mvn") %>% 
  dharma_residuals(res_distance_rp_tox, return_DHARMa = TRUE)
tox_res_rp2 = recalculateResiduals(tox_res_rp, group = as.factor(disttox$animalid),rotation="estimated")
groupLocations = aggregate(disttox[,c("mX","mY")], list(disttox$animalid), mean)
testSpatialAutocorrelation(tox_res_rp2,groupLocations$mX,groupLocations$mY)
#removed spatial autocorrelation, p=0.106

saveRDS(res_distance_rp_tox,paste0(results_dir,"res_distance_rp_tox.rds"))

#### removal type * period * sex ------
res_distance_rps_tox=glmmTMB(weekly_dist_km ~ trt_ctrl*removal.period.akdecalc+
                               (1|animalid),
                             data=disttox,family=Gamma(link="log"))

saveRDS(res_distance_rps_tox,paste0(results_dir,"res_distance_rps_tox.rds"))

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

saveRDS(res_speed_rp_aer,paste0(results_dir,"res_speed_rp_aer.rds"))

#### removal type * period * sex ------
res_speed_rps_aer <- glmmTMB(weekly_md_km_hr ~ trt_ctrl*removal.period.akdecalc*sex+
                              (1|animalid),
                            data=speedaer,
                            family=Gamma(link="log"),
                            verbose=TRUE)

saveRDS(res_speed_rps_aer,paste0(results_dir,"res_speed_rps_aer.rds"))

### trap ----------------
#### removal type * period ------
res_speed_rp_trap=glmmTMB(weekly_md_km_hr~ 
                            (1|animalid)+
                            trt_ctrl*removal.period.akdecalc,
                          data=speedtrap,family=Gamma(link='log'))

saveRDS(res_speed_rp_trap,paste0(results_dir,"res_speed_rp_trap.rds"))

#### removal type * period * sex ------
res_speed_rps_trap=glmmTMB(weekly_md_km_hr~ 
                             (1|animalid)+
                             trt_ctrl*removal.period.akdecalc*sex,
                           data=speedtrap,family=Gamma(link='log'))

saveRDS(res_speed_rps_trap,paste0(results_dir,"res_speed_rps_trap.rds"))

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

saveRDS(res_speed_rp_tox,paste0(results_dir,"res_speed_rp_tox.rds"))

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

saveRDS(res_speed_rps_tox,paste0(results_dir,"res_speed_rps_tox.rds"))

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

#offset for 0's 
conaer$contacts_per_day<- conaer$contacts_per_day+0.0001

## removal type * period -----
res_ncon_rp_aer=glmmTMB(contacts_per_day~(1|animalid)+
                        Removal.Type*removal.period.akdecalc,
                         data=conaer,
                         family=Gamma(link="log")
                        # family=nbinom12(link="log")
                        )

aer_res_ncon_rp <- simulateResiduals(res_ncon_rp_aer)
aer_res_ncon_rp2 = recalculateResiduals(aer_res_ncon_rp, 
                                        group = as.factor(conaer$animalid), 
                                        rotation="estimated")
plot(aer_res_ncon_rp2)
plot(aer_res_ncon_rp2$simulatedResponse)
descdist(aer_res_ncon_rp2$fittedResiduals)
DHARMa::testDispersion(aer_res_ncon_rp2)

saveRDS(res_ncon_rp_aer,paste0(results_dir,"res_ncon_rp_aer.rds"))

## removal type * period *sex -----
res_ncon_rps_aer=glmmTMB(num_contacts~(1|animalid)+
                          Removal.Type*removal.period.akdecalc*sex,
                        data=conaer,
                        family=poisson(link='log'))

aer_res_ncon_rps <- simulateResiduals(res_ncon_rps_aer)
aer_res_ncon_rps2 = recalculateResiduals(aer_res_ncon_rps, 
                                        group = as.factor(conaer$animalid), 
                                        rotation="estimated")
plot(aer_res_ncon_rps2)
plot(aer_res_ncon_rps2$simulatedResponse)
descdist(aer_res_ncon_rps2$fittedResiduals)
DHARMa::testDispersion(aer_res_ncon_rps2)

saveRDS(res_ncon_rps_aer,paste0(results_dir,"res_ncon_rps_aer.rds"))

#trap ------
## removal type * period -----
res_ncon_rp_trap=glmmTMB(num_contacts~(1|animalid)+
                          Removal.Type*removal.period.akdecalc,
                        data=contrap,
                        family=poisson(link='log'))

trap_res_ncon_rp <- simulateResiduals(res_ncon_rp_trap)
trap_res_ncon_rp2 = recalculateResiduals(trap_res_ncon_rp, 
                                        group = as.factor(contrap$animalid), 
                                        rotation="estimated")
plot(trap_res_ncon_rp2)
plot(trap_res_ncon_rp2$simulatedResponse)
descdist(trap_res_ncon_rp2$fittedResiduals)
DHARMa::testDispersion(trap_res_ncon_rp2)

saveRDS(res_ncon_rp_trap,paste0(results_dir,"res_ncon_rp_trap.rds"))

## removal type * period *sex -----
res_ncon_rps_trap=glmmTMB(num_contacts~(1|animalid)+
                           Removal.Type*removal.period.akdecalc*sex,
                         data=contrap,
                         family=poisson(link='log'))

trap_res_ncon_rps <- simulateResiduals(res_ncon_rps_trap)
trap_res_ncon_rps2 = recalculateResiduals(trap_res_ncon_rps, 
                                         group = as.factor(contrap$animalid), 
                                         rotation="estimated")
plot(trap_res_ncon_rps2)
plot(trap_res_ncon_rps2$simulatedResponse)
descdist(trap_res_ncon_rps2$fittedResiduals)
DHARMa::testDispersion(trap_res_ncon_rps2)

saveRDS(res_ncon_rps_trap,paste0(results_dir,"res_ncon_rps_trap.rds"))

#tox ------
## removal type * period -----
res_ncon_rp_tox=glmmTMB(num_contacts~(1|animalid)+
                          Removal.Type*removal.period.akdecalc,
                        data=contox,
                        family=poisson(link='log'))

tox_res_ncon_rp <- simulateResiduals(res_ncon_rp_tox)
tox_res_ncon_rp2 = recalculateResiduals(tox_res_ncon_rp, 
                                        group = as.factor(contox$animalid), 
                                        rotation="estimated")
plot(tox_res_ncon_rp2)
plot(tox_res_ncon_rp2$simulatedResponse)
descdist(tox_res_ncon_rp2$fittedResiduals)
DHARMa::testDispersion(tox_res_ncon_rp2)

saveRDS(res_ncon_rp_tox,paste0(results_dir,"res_ncon_rp_tox.rds"))

## removal type * period *sex -----
res_ncon_rps_tox=glmmTMB(num_contacts~(1|animalid)+
                           Removal.Type*removal.period.akdecalc*sex,
                         data=contox,
                         family=poisson(link='log'))

tox_res_ncon_rps <- simulateResiduals(res_ncon_rps_tox)
tox_res_ncon_rps2 = recalculateResiduals(tox_res_ncon_rps, 
                                         group = as.factor(contox$animalid), 
                                         rotation="estimated")
plot(tox_res_ncon_rps2)
plot(tox_res_ncon_rps2$simulatedResponse)
descdist(tox_res_ncon_rps2$fittedResiduals)
DHARMa::testDispersion(tox_res_ncon_rps2)

saveRDS(res_ncon_rps_tox,paste0(results_dir,"res_ncon_rps_tox.rds"))


# contacts - degree --------

#aerial ------
## removal type * period -----
res_nind_rp_aer=glmmTMB(num_indivs~(1|animalid)+
                          Removal.Type*removal.period.akdecalc,
                        data=conaer,
                        family=poisson(link='log'))

aer_res_nind_rp <- simulateResiduals(res_nind_rp_aer)
aer_res_nind_rp2 = recalculateResiduals(aer_res_nind_rp, 
                                        group = as.factor(conaer$animalid), 
                                        rotation="estimated")
plot(aer_res_nind_rp2)
plot(aer_res_nind_rp2$simulatedResponse)
descdist(aer_res_nind_rp2$fittedResiduals)
DHARMa::testDispersion(aer_res_nind_rp2)

saveRDS(res_nind_rp_aer,paste0(results_dir,"res_nind_rp_aer.rds"))

## removal type * period *sex -----
res_nind_rps_aer=glmmTMB(num_indivs~(1|animalid)+
                           Removal.Type*removal.period.akdecalc*sex,
                         data=conaer,
                         family=poisson(link='log'))

aer_res_nind_rps <- simulateResiduals(res_nind_rps_aer)
aer_res_nind_rps2 = recalculateResiduals(aer_res_nind_rps, 
                                         group = as.factor(conaer$animalid), 
                                         rotation="estimated")
plot(aer_res_nind_rps2)
plot(aer_res_nind_rps2$simulatedResponse)
descdist(aer_res_nind_rps2$fittedResiduals)
DHARMa::testDispersion(aer_res_nind_rps2)

saveRDS(res_nind_rps_aer,paste0(results_dir,"res_nind_rps_aer.rds"))

#trap ------
## removal type * period -----
res_nind_rp_trap=glmmTMB(num_indivs~(1|animalid)+
                           Removal.Type*removal.period.akdecalc,
                         data=contrap,
                         family=poisson(link='log'))

trap_res_nind_rp <- simulateResiduals(res_nind_rp_trap)
trap_res_nind_rp2 = recalculateResiduals(trap_res_nind_rp, 
                                         group = as.factor(contrap$animalid), 
                                         rotation="estimated")
plot(trap_res_nind_rp2)
plot(trap_res_nind_rp2$simulatedResponse)
descdist(trap_res_nind_rp2$fittedResiduals)
DHARMa::testDispersion(trap_res_nind_rp2)

saveRDS(res_nind_rp_trap,paste0(results_dir,"res_nind_rp_trap.rds"))

## removal type * period *sex -----
res_nind_rps_trap=glmmTMB(num_indivs~(1|animalid)+
                            Removal.Type*removal.period.akdecalc*sex,
                          data=contrap,
                          family=poisson(link='log'))

trap_res_nind_rps <- simulateResiduals(res_nind_rps_trap)
trap_res_nind_rps2 = recalculateResiduals(trap_res_nind_rps, 
                                          group = as.factor(contrap$animalid), 
                                          rotation="estimated")
plot(trap_res_nind_rps2)
plot(trap_res_nind_rps2$simulatedResponse)
descdist(trap_res_nind_rps2$fittedResiduals)
DHARMa::testDispersion(trap_res_nind_rps2)

saveRDS(res_nind_rps_trap,paste0(results_dir,"res_nind_rps_trap.rds"))

#tox ------
## removal type * period -----
res_nind_rp_tox=glmmTMB(num_indivs~(1|animalid)+
                          Removal.Type*removal.period.akdecalc,
                        data=contox,
                        family=poisson(link='log'))

tox_res_nind_rp <- simulateResiduals(res_nind_rp_tox)
tox_res_nind_rp2 = recalculateResiduals(tox_res_nind_rp, 
                                        group = as.factor(contox$animalid), 
                                        rotation="estimated")
plot(tox_res_nind_rp2)
plot(tox_res_nind_rp2$simulatedResponse)
descdist(tox_res_nind_rp2$fittedResiduals)
DHARMa::testDispersion(tox_res_nind_rp2)

saveRDS(res_nind_rp_tox,paste0(results_dir,"res_nind_rp_tox.rds"))

## removal type * period *sex -----
res_nind_rps_tox=glmmTMB(num_indivs~(1|animalid)+
                           Removal.Type*removal.period.akdecalc*sex,
                         data=contox,
                         family=poisson(link='log'))

tox_res_nind_rps <- simulateResiduals(res_nind_rps_tox)
tox_res_nind_rps2 = recalculateResiduals(tox_res_nind_rps, 
                                         group = as.factor(contox$animalid), 
                                         rotation="estimated")
plot(tox_res_nind_rps2)
plot(tox_res_nind_rps2$simulatedResponse)
descdist(tox_res_nind_rps2$fittedResiduals)
DHARMa::testDispersion(tox_res_nind_rps2)

saveRDS(res_nind_rps_tox,paste0(results_dir,"res_nind_rps_tox.rds"))

