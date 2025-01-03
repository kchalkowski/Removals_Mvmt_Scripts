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

{
  library(glmmTMB)
  library(DHARMa)
  library(ctmm)
  library(fitdistrplus)
  library(ggplot2)
  library(plyr)
  library(dplyr)
  library(hrbrthemes)
  library(sdmTMB)
}

#setup ----------------
{
  
  homedir <- "//aapcoftc3fp13/Projects/MUDD/ASF_NIFA/Pipelines/Removals_Mvmt"
  homedir <- "C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Contact Analysis/Removals_Mvmt"
  
  objdir=file.path(homedir,"1_Data","Objects",fsep=.Platform$file.sep)
  
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
  results_dir <- paste0(homedir,"/3_Output/")
  
  #set mesh cutoff for all spatial models
  mesh_cutoff<-3
  spatial_res <- 1000
}

#distance -------------------
{
  dist<- readRDS(paste0(objdir,"/pig_weekly_distance_ctmm.rds"))
  
  #offset to prevent 0's
  dist$weekly_dist_km <- dist$weekly_dist_km +0.0001
  dist$animalid <- factor(dist$animalid)
  
  dist$sex <- factor(dist$sex)
  
  distaer <- dist %>% filter(Removal.Type%in%c('aer','ctrl')) %>% 
    filter(removal.period.akdecalc!='during') 
  disttox <- dist %>% filter(Removal.Type%in%c('tox','ctrl')) 
  disttrap <- dist %>% filter(Removal.Type%in%c('trap','ctrl')) 
  
  distaer$Removal.Type <- factor(distaer$Removal.Type,levels=c('ctrl','aer'))
  disttox$Removal.Type <- factor(disttox$Removal.Type,levels=c('ctrl','tox'))
  disttrap$Removal.Type <- factor(disttrap$Removal.Type,levels=c('ctrl','trap'))
  
  #relevel when during isn't recorded
  distaer$removal.period.akdecalc <- factor(distaer$removal.period.akdecalc,
                                            levels=c('before','after'))
  disttox$removal.period.akdecalc <- factor(disttox$removal.period.akdecalc,
                                            levels=c('before','during','after'))
  disttrap$removal.period.akdecalc <- factor(disttrap$removal.period.akdecalc,
                                             levels=c('before','during','after'))
  
  #id toxicant pigs that died 
  # disttox$died_tox <- NA
  # disttox$died_tox <- ifelse(disttox$animalid%in%died.tox,"0","1")
  # disttox$died_tox <- factor(disttox$died_tox)
  
  saveRDS(distaer,paste0(results_dir,"distaer.rds"))
  saveRDS(disttox,paste0(results_dir,"disttox.rds"))
  saveRDS(disttrap,paste0(results_dir,"distrap.rds"))
}

## temporal autocorrelation test --------------
res_aer=glmmTMB(weekly_dist_km~(1|animalid)+Removal.Type*removal.period.akdecalc,
                data=distaer,family=Gamma(link="log"))
simout_aer <- simulateResiduals(fittedModel = res_aer, plot = F)
res_aer2 = recalculateResiduals(simout_aer, group = distaer$week,rotation="estimated")
testTemporalAutocorrelation(res_aer2, time = unique(distaer$week))
#no temporal autocorrelation

#sex 
res_aer=glmmTMB(weekly_dist_km~(1|animalid)+Removal.Type*removal.period.akdecalc*sex,
                data=distaer,family=Gamma(link="log"))
simout_aer <- simulateResiduals(fittedModel = res_aer, plot = F)
res_aer2 = recalculateResiduals(simout_aer, group = distaer$week,rotation="estimated")
testTemporalAutocorrelation(res_aer2, time = unique(distaer$week))
#no temporal autocorrelation

res_trap=glmmTMB(weekly_dist_km~(1|animalid)+Removal.Type*removal.period.akdecalc,
                 data=disttrap,family=Gamma(link="log"))
simout_trap <- simulateResiduals(fittedModel = res_trap, plot = F)
res_trap2 = recalculateResiduals(simout_trap, group = disttrap$week,rotation="estimated")
testTemporalAutocorrelation(res_trap2, time = unique(disttrap$week))
#temporal autocorrelation

#sex
res_trap=glmmTMB(weekly_dist_km~(1|animalid)+Removal.Type*removal.period.akdecalc*sex,
                 data=disttrap,family=Gamma(link="log"))
simout_trap <- simulateResiduals(fittedModel = res_trap, plot = F)
res_trap2 = recalculateResiduals(simout_trap, group = disttrap$week,rotation="estimated")
testTemporalAutocorrelation(res_trap2, time = unique(disttrap$week))
#temporal autocorrelation

res_tox=glmmTMB(weekly_dist_km~(1|animalid)+Removal.Type*removal.period.akdecalc,
                data=disttox,family=Gamma(link="log"))
simout_tox <- simulateResiduals(fittedModel = res_tox, plot = F)
res_tox2 = recalculateResiduals(simout_tox, group = disttox$week,rotation="estimated")
testTemporalAutocorrelation(res_tox2, time = unique(disttox$week))
#temporal autocorrelation

res_tox=glmmTMB(weekly_dist_km~(1|animalid)+Removal.Type*removal.period.akdecalc*sex,
                data=disttox,family=Gamma(link="log"))
simout_tox <- simulateResiduals(fittedModel = res_tox, plot = F)
res_tox2 = recalculateResiduals(simout_tox, group = disttox$week,rotation="estimated")
testTemporalAutocorrelation(res_tox2, time = unique(disttox$week))
#temporal autocorrelation

## spatial autocorrelation test ------------------
aer_test=glmmTMB(weekly_dist_km~(1|animalid)+Removal.Type*removal.period.akdecalc,
                 data=distaer,family=Gamma(link="log"))
aer_test_resid <- simulateResiduals(aer_test)
aer_groupLocations = aggregate(distaer[,c("mX","mY")], list(distaer$animalid), mean)
aer_test_resid2 = recalculateResiduals(aer_test_resid, group = distaer$animalid, rotation="estimated")
testSpatialAutocorrelation(aer_test_resid2,aer_groupLocations$mX, aer_groupLocations$mY)
# spatial autocorrelation

#sex
aer_test=glmmTMB(weekly_dist_km~(1|animalid)+Removal.Type*removal.period.akdecalc*sex,
                 data=distaer,family=Gamma(link="log"))
aer_test_resid <- simulateResiduals(aer_test)
aer_test_resid2 = recalculateResiduals(aer_test_resid, group = distaer$animalid, rotation="estimated")
testSpatialAutocorrelation(aer_test_resid2,aer_groupLocations$mX, aer_groupLocations$mY)
# no spatial autocorrelation

trap_test=glmmTMB(weekly_dist_km~(1|animalid)+Removal.Type*removal.period.akdecalc,
                  data=disttrap,family=Gamma(link="log"))
trap_test_resid <- simulateResiduals(trap_test)
trap_groupLocations = aggregate(disttrap[,c("mX","mY")], list(disttrap$animalid), mean)
trap_test_resid2 = recalculateResiduals(trap_test_resid, group = disttrap$animalid, rotation="estimated")
testSpatialAutocorrelation(trap_test_resid2,trap_groupLocations$mX, trap_groupLocations$mY)
#no spatial autocorrelation

#sex 
trap_test=glmmTMB(weekly_dist_km~(1|animalid)+Removal.Type*removal.period.akdecalc*sex,
                  data=disttrap,family=Gamma(link="log"))
trap_test_resid <- simulateResiduals(trap_test)
trap_test_resid2 = recalculateResiduals(trap_test_resid, group = disttrap$animalid, rotation="estimated")
testSpatialAutocorrelation(trap_test_resid2,trap_groupLocations$mX, trap_groupLocations$mY)
#no spatial autocorrelation

tox_test=glmmTMB(weekly_dist_km~(1|animalid)+Removal.Type*removal.period.akdecalc,
                 data=disttox,family=Gamma(link="log"))
tox_test_resid <- simulateResiduals(tox_test)
tox_groupLocations = aggregate(disttox[,c("mX","mY")], list(disttox$animalid), mean)
tox_test_resid2 = recalculateResiduals(tox_test_resid, group = disttox$animalid, rotation="estimated")
testSpatialAutocorrelation(tox_test_resid2,tox_groupLocations$mX, tox_groupLocations$mY)
#no spatial autocorrelation (p=0.056)

#sex
tox_test=glmmTMB(weekly_dist_km~(1|animalid)+Removal.Type*removal.period.akdecalc*sex,
                 data=disttox,family=Gamma(link="log"))
tox_test_resid <- simulateResiduals(tox_test)
tox_test_resid2 = recalculateResiduals(tox_test_resid, group = disttox$animalid, rotation="estimated")
testSpatialAutocorrelation(tox_test_resid2,tox_groupLocations$mX, tox_groupLocations$mY)
#no spatial autocorrelation

#removal type * period 
### aer spatial
### trap temporal 
### tox temporal

#removal type * period * sex
### aer 
### trap temporal 
### tox temporal

## fit models ------------------
### aerial -----------------
distaer$mX_sc <- distaer$mX/spatial_res
distaer$mY_sc <- distaer$mY/spatial_res
meshaer <- make_mesh(distaer,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)
distaer$pos <- numFactor(distaer$mX_sc, distaer$mY_sc) 

#### removal type * period ------
# res_distance_rp_aer <- glmmTMB(weekly_dist_km ~ #ar1(as.factor(week)+0|animalid) + 
#                                  exp(pos + 0 | animalid) +
#                                  Removal.Type*removal.period.akdecalc,
#                                data=distaer,
#                                family=Gamma(link="log"),
#                                verbose=TRUE)

res_distance_rp_aer <- sdmTMB(weekly_dist_km ~ (1|animalid) +
                              Removal.Type*removal.period.akdecalc,
                              data=distaer,
                              mesh=meshaer,
                              spatial='on',
                              family=Gamma(link='log'))

sanity(res_distance_rp_aer)

aer_res <- simulate(res_distance_rp_aer, nsim = 250, type = "mle-mvn") %>%
  dharma_residuals(res_distance_rp_aer, return_DHARMa = TRUE,rotation="estimated")
aer_res2 = recalculateResiduals(aer_res, group = as.factor(distaer$animalid),rotation="estimated")
plot(aer_res2)
plot(aer_res2$simulatedResponse)
DHARMa::testDispersion(aer_res2)

testSpatialAutocorrelation(aer_res2,aer_groupLocations$mX,aer_groupLocations$mY)
#took care of spatial autocorrelation
saveRDS(res_distance_rp_aer,paste0(results_dir,"res_distance_rp_aer.rds"))

#### removal type * period * sex ------
# res_distance_rps_aer=sdmTMB(weekly_dist_km ~ (1|animalid) + 
#                             Removal.Type*removal.period.akdecalc*sex,
#                             data=distaer,
#                             mesh=meshaer,
#                             time='week',
#                             spatial='off',
#                             spatiotemporal = "ar1",
#                             family=Gamma(link='log'))
# sanity(res_distance_rps_aer)

# aer_res_rps <- simulate(res_distance_rps_aer, nsim = 250, type = "mle-mvn") %>% 
#   dharma_residuals(res_distance_rps_aer, return_DHARMa = TRUE)

# dharma_residuals(aer_res_rps,res_distance_rps_aer)
# aer_res_rps2 = recalculateResiduals(aer_res_rps, group = as.factor(distaer$animalid),rotation="estimated")
# plot(aer_res_rps2)
# plot(aer_res_rps2$simulatedResponse)
# DHARMa::testDispersion(aer_res_rps2)

res_distance_rps_aer <-glmmTMB(weekly_dist_km ~ #ar1(as.factor(week)+0|animalid) + 
                                 (1|animalid)+
                                 Removal.Type*removal.period.akdecalc*sex,
                               data=distaer,family=Gamma(link="log"))
aer_res_rps <- simulateResiduals(res_distance_rps_aer)
aer_res_rps2 = recalculateResiduals(aer_res_rps, group = as.factor(distaer$animalid),rotation="estimated")
plot(aer_res_rps2)
plot(aer_res_rps2$simulatedResponse)
descdist(aer_res_rps2$fittedResiduals)
DHARMa::testDispersion(aer_res_rps2) 
aer_res_rsp_week <- recalculateResiduals(aer_res_rps, group = as.factor(distaer$week),rotation="estimated")
testTemporalAutocorrelation(aer_res_rsp_week,time=unique(distaer$week))
#took care of temporal autocorrelation
saveRDS(res_distance_rps_aer,paste0(results_dir,"res_distance_rps_aer.rds"))

### trap -----------------
#### removal type * period ------
res_distance_rp_trap=glmmTMB(weekly_dist_km ~ ar1(as.factor(week)+0|animalid) + 
                               Removal.Type*removal.period.akdecalc,
                             data=disttrap,family=Gamma(link="log"))
trap_res <- simulateResiduals(res_distance_rp_trap)
trap_res2 = recalculateResiduals(trap_res, group = as.factor(disttrap$animalid),rotation="estimated")
plot(trap_res2)
plot(trap_res2$simulatedResponse)
descdist(trap_res2$fittedResiduals)
DHARMa::testDispersion(trap_res2) 
trap_res_week <- recalculateResiduals(trap_res2, group = as.factor(disttrap$week),rotation="estimated")
testTemporalAutocorrelation(trap_res_week,time=unique(disttrap$week))
saveRDS(res_distance_rp_trap,paste0(results_dir,"res_distance_rp_trap.rds"))

#### removal type * period * sex ------
res_distance_rps_trap=glmmTMB(weekly_dist_km ~ ar1(as.factor(week)+0|animalid) + 
                                Removal.Type*removal.period.akdecalc*sex,
                              data=disttrap,family=Gamma(link="log"))
trap_res_rps <- simulateResiduals(res_distance_rps_trap)
trap_res_rps2 = recalculateResiduals(trap_res_rps, group = as.factor(disttrap$animalid),rotation="estimated")
plot(trap_res_rps2)
plot(trap_res_rps2$simulatedResponse)
descdist(trap_res_rps2$fittedResiduals)
DHARMa::testDispersion(trap_res_rps2)
trap_res_rps_week <- recalculateResiduals(trap_res_rps2, group = as.factor(disttrap$week),rotation="estimated")
testTemporalAutocorrelation(trap_res_rps_week,time=unique(disttrap$week))
#took care of spatial autocorrelation
saveRDS(res_distance_rps_trap,paste0(results_dir,"res_distance_rps_trap.rds"))

### toxicant -----------------
#### removal type * period ------
res_distance_rp_tox=glmmTMB(weekly_dist_km ~ ar1(as.factor(week)+0|animalid) + 
                              Removal.Type*removal.period.akdecalc,
                            data=disttox,family=Gamma(link="log"))
tox_res <- simulateResiduals(res_distance_rp_tox)
tox_res2 = recalculateResiduals(tox_res, group = as.factor(disttox$animalid),rotation="estimated")
plot(tox_res2)
plot(tox_res2$simulatedResponse)
descdist(tox_res2$fittedResiduals)
DHARMa::testDispersion(tox_res2)

tox_res_week <- recalculateResiduals(tox_res2, group = as.factor(disttox$week),rotation="estimated")
testTemporalAutocorrelation(tox_res_week,time=unique(disttox$week))
#took care of temporal autocorrelation
saveRDS(res_distance_rp_tox,paste0(results_dir,"res_distance_rp_tox.rds"))

#### removal type * period * sex ------
res_distance_rps_tox=glmmTMB(weekly_dist_km ~ ar1(as.factor(week)+0|animalid) + 
                               Removal.Type*removal.period.akdecalc*sex,
                             data=disttox,family=Gamma(link="log"))
tox_res_rps <- simulateResiduals(res_distance_rps_tox)
tox_res_rps2 = recalculateResiduals(tox_res_rps, group = as.factor(disttox$animalid),rotation="estimated")
plot(tox_res_rps2)
plot(tox_res_rps2$simulatedResponse)
descdist(tox_res_rps2$fittedResiduals)
DHARMa::testDispersion(tox_res_rps2)

tox_res_rps_week <- recalculateResiduals(tox_res_rps, group = as.factor(disttox$week),rotation="estimated")
testTemporalAutocorrelation(tox_res_rps_week,time=unique(disttox$week))
#still/more temporally autocorrelated
saveRDS(res_distance_rps_tox,paste0(results_dir,"res_distance_rps_tox.rds"))


#speed -------------------
{
  speed<- readRDS(paste0(objdir,"/pig_weekly_speed_ctmm.rds"))
  
  #offset by very small amount
  speed$weekly_md_km_hr <- speed$weekly_md_km_hr + 0.0001
  speed$animalid <- factor(speed$animalid)
  
  speed$sex <- factor(speed$sex)
  
  speedaer <- speed %>% filter(Removal.Type%in%c('aer','ctrl')) %>% 
    filter(removal.period.akdecalc!='during') 
  speedtox <- speed %>% filter(Removal.Type%in%c('tox','ctrl')) 
  speedtrap <- speed %>% filter(Removal.Type%in%c('trap','ctrl')) 
  
  speedaer$Removal.Type <- factor(speedaer$Removal.Type,levels=c('ctrl','aer'))
  speedtox$Removal.Type <- factor(speedtox$Removal.Type,levels=c('ctrl','tox'))
  speedtrap$Removal.Type <- factor(speedtrap$Removal.Type,levels=c('ctrl','trap'))
  
  #relevel when during isn't recorded
  #relevel when during isn't recorded
  speedaer$removal.period.akdecalc <- factor(speedaer$removal.period.akdecalc,
                                             levels=c('before','after'))
  speedtox$removal.period.akdecalc <- factor(speedtox$removal.period.akdecalc,
                                             levels=c('before','during','after'))
  speedtrap$removal.period.akdecalc <- factor(speedtrap$removal.period.akdecalc,
                                              levels=c('before','during','after'))
  
  #id toxicant pigs that died 
  # speedtox$died_tox <- NA
  # speedtox$died_tox <- ifelse(speedtox$animalid%in%died.tox,"0","1")
  # speedtox$died_tox <- factor(speedtox$died_tox)
  
  saveRDS(speedaer,paste0(results_dir,"speedaer.rds"))
  saveRDS(speedtox,paste0(results_dir,"speedtox.rds"))
  saveRDS(speedtrap,paste0(results_dir,"speedtrap.rds"))
}

## temporal autocorrelation test --------------
res_aers=glmmTMB(weekly_md_km_hr~(1|animalid)+Removal.Type*removal.period.akdecalc,
                 data=speedaer,family=Gamma(link='log'))
simout_aers <- simulateResiduals(fittedModel = res_aers, plot = F)
res_aers2 = recalculateResiduals(simout_aers, group = speedaer$week,rotation="estimated")
testTemporalAutocorrelation(res_aers2, time = unique(speedaer$week))
# no temporal autocorrelation

#sex
res_aers=glmmTMB(weekly_md_km_hr~(1|animalid)+Removal.Type*removal.period.akdecalc*sex,
                 data=speedaer,family=Gamma(link='log'))
simout_aers <- simulateResiduals(fittedModel = res_aers, plot = F)
res_aers2 = recalculateResiduals(simout_aers, group = speedaer$week,rotation="estimated")
testTemporalAutocorrelation(res_aers2, time = unique(speedaer$week))
# no temporal autocorrelation

res_traps=glmmTMB(weekly_md_km_hr~(1|animalid)+Removal.Type*removal.period.akdecalc,
                  data=speedtrap,family=Gamma(link='log'))
simout_traps <- simulateResiduals(fittedModel = res_traps, plot = F)
res_traps2 = recalculateResiduals(simout_traps, group = speedtrap$week,rotation="estimated")
testTemporalAutocorrelation(res_traps2, time = unique(speedtrap$week))
#no temporal autocorrelation

#sex
res_traps=glmmTMB(weekly_md_km_hr~(1|animalid)+Removal.Type*removal.period.akdecalc*sex,
                  data=speedtrap,family=Gamma(link='log'))
simout_traps <- simulateResiduals(fittedModel = res_traps, plot = F)
res_traps2 = recalculateResiduals(simout_traps, group = speedtrap$week,rotation="estimated")
testTemporalAutocorrelation(res_traps2, time = unique(speedtrap$week))
#no temporal autocorrelation

res_toxs=glmmTMB(weekly_md_km_hr~(1|animalid)+Removal.Type*removal.period.akdecalc,
                 data=speedtox,family=Gamma(link='log'))
simout_toxs <- simulateResiduals(fittedModel = res_toxs, plot = F)
res_toxs2 = recalculateResiduals(simout_toxs, group = speedtox$week,rotation="estimated")
testTemporalAutocorrelation(res_toxs2, time = unique(speedtox$week))
# no temporal autocorrelation (p=0.06)

#sex
res_toxs=glmmTMB(weekly_md_km_hr~(1|animalid)+Removal.Type*removal.period.akdecalc*sex,
                 data=speedtox,family=Gamma(link='log'))
simout_toxs <- simulateResiduals(fittedModel = res_toxs, plot = F)
res_toxs2 = recalculateResiduals(simout_toxs, group = speedtox$week,rotation="estimated")
testTemporalAutocorrelation(res_toxs2, time = unique(speedtox$week))
#no temporal autocorrelation (p=0.05)

## spatial autocorrelation test ------------------
aer_tests=glmmTMB(weekly_md_km_hr~(1|animalid)+Removal.Type*removal.period.akdecalc,
                  data=speedaer,family=Gamma(link='log'))
aer_test_resids <- simulateResiduals(aer_tests)
aer_groupLocations_sp = aggregate(speedaer[,c("mX","mY")], list(speedaer$animalid), mean)
aer_test_resids2 = recalculateResiduals(aer_test_resids, group = speedaer$animalid, rotation="estimated")
testSpatialAutocorrelation(aer_test_resids2,aer_groupLocations_sp$mX, aer_groupLocations_sp$mY)
# spatial autocorrelation

#sex
aer_tests=glmmTMB(weekly_md_km_hr~(1|animalid)+Removal.Type*removal.period.akdecalc*sex,
                  data=speedaer,family=Gamma(link='log'))
aer_test_resids <- simulateResiduals(aer_tests)
aer_test_resids2 = recalculateResiduals(aer_test_resids, group = speedaer$animalid, rotation="estimated")
testSpatialAutocorrelation(aer_test_resids2,aer_groupLocations_sp$mX, aer_groupLocations_sp$mY)
# no spatial autocorrelation


trap_tests=glmmTMB(weekly_md_km_hr~(1|animalid)+Removal.Type*removal.period.akdecalc,
                   data=speedtrap,family=Gamma(link='log'))
trap_test_resids <- simulateResiduals(trap_tests)
trap_groupLocations_sp = aggregate(speedtrap[,c("mX","mY")], list(speedtrap$animalid), mean)
trap_test_resids2 = recalculateResiduals(trap_test_resids, group = speedtrap$animalid, rotation="estimated")
testSpatialAutocorrelation(trap_test_resids2,trap_groupLocations_sp$mX, trap_groupLocations_sp$mY)
#no spatial autocorrelation

#sex
trap_tests=glmmTMB(weekly_md_km_hr~(1|animalid)+Removal.Type*removal.period.akdecalc*sex,
                   data=speedtrap,family=Gamma(link='log'))
trap_test_resids <- simulateResiduals(trap_tests)
trap_test_resids2 = recalculateResiduals(trap_test_resids, group = speedtrap$animalid, rotation="estimated")
testSpatialAutocorrelation(trap_test_resids2,trap_groupLocations_sp$mX, trap_groupLocations_sp$mY)
#no spatial autocorrelation

tox_tests=glmmTMB(weekly_md_km_hr~(1|animalid)+Removal.Type*removal.period.akdecalc,
                  data=speedtox,family=Gamma(link='log'))
tox_test_resids <- simulateResiduals(tox_tests)
tox_groupLocations_sp = aggregate(speedtox[,c("mX","mY")], list(speedtox$animalid), mean)
tox_test_resids2 = recalculateResiduals(tox_test_resids, group = speedtox$animalid, rotation="estimated")
testSpatialAutocorrelation(tox_test_resids2,tox_groupLocations_sp$mX, tox_groupLocations_sp$mY)
# spatial autocorrelation

tox_tests=glmmTMB(weekly_md_km_hr~(1|animalid)+Removal.Type*removal.period.akdecalc*sex,
                  data=speedtox,family=Gamma(link='log'))
tox_test_resids <- simulateResiduals(tox_tests)
tox_test_resids2 = recalculateResiduals(tox_test_resids, group = speedtox$animalid, rotation="estimated")
testSpatialAutocorrelation(tox_test_resids2,tox_groupLocations_sp$mX, tox_groupLocations_sp$mY)
# spatial autocorrelation

#removal type * period
# aer: spatial
# trap: 
# tox:  spatial

#removal type * period *sex
# aer: 
# trap: 
# tox: spatial

## fit glmms ------------------
### aerial ----------------
speedaer$mX_sc <- speedaer$mX/spatial_res
speedaer$mY_sc <- speedaer$mY/spatial_res
meshaer_sp <- make_mesh(speedaer,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)

#### removal type * period  ------
res_speed_rp_aer=sdmTMB(weekly_md_km_hr ~ (1|animalid)+
                        Removal.Type*removal.period.akdecalc,
                        data=speedaer,
                        mesh=meshaer_sp,
                        # time='week',
                        spatial='on',
                        family=Gamma(link='log'))
sanity(res_speed_rp_aer)

aer_res_sp_rp <- simulate(res_speed_rp_aer, nsim = 250, type = "mle-mvn") %>% 
  dharma_residuals(res_speed_rp_aer, return_DHARMa = TRUE)
aer_res_sp_rp2 = recalculateResiduals(aer_res_sp_rp, group = as.factor(speedaer$animalid),rotation="estimated")
plot(aer_res_sp_rp2)
plot(aer_res_sp_rp2$simulatedResponse)
DHARMa::testDispersion(aer_res_sp_rp2)

testSpatialAutocorrelation(aer_res_sp_rp2,aer_groupLocations_sp$mX,aer_groupLocations_sp$mY)
#took care of spatial autocorrelation 

saveRDS(res_speed_rp_aer,paste0(results_dir,"res_speed_rp_aer.rds"))

#### removal type * period * sex ------
res_speed_rps_aer=glmmTMB(weekly_md_km_hr~ #ar1(as.factor(week)+0 | animalid) + 
                            (1|animalid)+
                            Removal.Type*removal.period.akdecalc*sex,
                          data=speedaer,family=Gamma(link='log'))
aer_res_rps_sp <- simulateResiduals(res_speed_rps_aer)
aer_res_rps_sp2 = recalculateResiduals(aer_res_rps_sp, group = as.factor(speedaer$animalid), rotation="estimated")
plot(aer_res_rps_sp2)
plot(aer_res_rps_sp2$simulatedResponse)
descdist(aer_res_rps_sp2$fittedResiduals)
DHARMa::testDispersion(aer_res_rps_sp2)

# aer_res_sp_rps_week = recalculateResiduals(aer_res_rps_sp2, group = as.factor(speedaer$week),rotation="estimated")
# testTemporalAutocorrelation(aer_res_sp_rps_week,time=unique(speedaer$week))
# #took care of temporal autocorrelation

saveRDS(res_speed_rps_aer,paste0(results_dir,"res_speed_rps_aer.rds"))


### trap ----------------
#### removal type * period ------
res_speed_rp_trap=glmmTMB(weekly_md_km_hr~ #ar1(as.factor(week)+0 | animalid) + 
                            (1|animalid)+
                            Removal.Type*removal.period.akdecalc,
                          data=speedtrap,family=Gamma(link='log'))
trap_res_sp <- simulateResiduals(res_speed_rp_trap)
trap_res_sp2 = recalculateResiduals(trap_res_sp, group = as.factor(speedtrap$animalid), rotation="estimated")
plot(trap_res_sp2)
plot(trap_res_sp2$simulatedResponse)
descdist(trap_res_sp2$fittedResiduals)
DHARMa::testDispersion(trap_res_sp2)

# trap_res_sp_week = recalculateResiduals(trap_res_sp2, group = as.factor(speedtrap$week),rotation="estimated")
# testTemporalAutocorrelation(trap_res_sp_week,time=unique(speedtrap$week))
# #took care of temporal autocorrelation

saveRDS(res_speed_rp_trap,paste0(results_dir,"res_speed_rp_trap.rds"))

#### removal type * period * sex ------
res_speed_rps_trap=glmmTMB(weekly_md_km_hr~#ar1(as.factor(week)+0 | animalid) + 
                             (1|animalid)+
                             Removal.Type*removal.period.akdecalc*sex,
                           data=speedtrap,family=Gamma(link='log'))

trap_res_sp_rps <- simulateResiduals(res_speed_rps_trap)
trap_res_sp_rps2 = recalculateResiduals(trap_res_sp_rps, group = as.factor(speedtrap$animalid), rotation="estimated")
plot(trap_res_sp_rps2)
plot(trap_res_sp_rps2$simulatedResponse)
descdist(trap_res_sp_rps2$fittedResiduals)
DHARMa::testDispersion(trap_res_sp_rps2)

# trap_res_sp_rps_week = recalculateResiduals(trap_res_sp_rps2, group = as.factor(speedtrap$week),rotation="estimated")
# testTemporalAutocorrelation(trap_res_sp_rps_week,time=unique(speedtrap$week))
# #took care of temporal autocorrelation

saveRDS(res_speed_rps_trap,paste0(results_dir,"res_speed_rps_trap.rds"))


### toxcicant ----------------
speedtox$mX_sc <- speedtox$mX/spatial_res
speedtox$mY_sc <- speedtox$mY/spatial_res
meshtox_sp <- make_mesh(speedtox,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)

#### removal type * period ------
res_speed_rp_tox=sdmTMB(weekly_md_km_hr ~ (1|animalid) + 
                        Removal.Type*removal.period.akdecalc,
                        data=speedtox,
                        mesh=meshtox_sp,
                        # time='week',
                        spatial='on',
                        # spatiotemporal = "ar1",
                        family=Gamma(link='log'))
sanity(res_speed_rp_tox)

tox_res_rp <- simulate(res_speed_rp_tox, nsim = 250, type = "mle-mvn") %>% 
  dharma_residuals(res_speed_rp_tox, return_DHARMa = TRUE)
tox_res_rp2 = recalculateResiduals(tox_res_rp, group = as.factor(distaer$animalid),rotation="estimated")
plot(tox_res_rp2)
plot(tox_res_rp2$simulatedResponse)
DHARMa::testDispersion(tox_res_rp2)

testSpatialAutocorrelation(tox_res_rp2,tox_groupLocations_sp$mX,tox_groupLocations_sp$mY)
#took care of spatial autocorrelation

tox_res_sp_week <- recalculateResiduals(tox_res_rp,group=factor(speedtox$week))
testTemporalAutocorrelation(tox_res_sp_week,time=unique(speedtox$week))
#took care of temporal autocorrealtion

saveRDS(res_speed_rp_tox,paste0(results_dir,"res_speed_rp_tox.rds"))

#### removal type * period * sex ------
res_speed_rps_tox=sdmTMB(weekly_md_km_hr ~ (1|animalid) + 
                         Removal.Type*removal.period.akdecalc*sex,
                         data=speedtox,
                         mesh=meshtox_sp,
                         # time='week',
                         spatial='on',
                         # spatiotemporal = "ar1",
                         family=Gamma(link='log'))
sanity(res_speed_rps_tox)

tox_res_rps <- simulate(res_speed_rps_tox, nsim = 250, type = "mle-mvn") %>% 
  dharma_residuals(res_speed_rps_tox, return_DHARMa = TRUE)
tox_res_rps2 = recalculateResiduals(tox_res_rps, group = as.factor(distaer$animalid),rotation="estimated")
plot(tox_res_rps2)
plot(tox_res_rps2$simulatedResponse)
DHARMa::testDispersion(tox_res_rps2)

testSpatialAutocorrelation(tox_res_rps2,tox_groupLocations_sp$mX,tox_groupLocations_sp$mY)
#took care of spatial autocorrelation

tox_res_sp_rps_week <- recalculateResiduals(tox_res_rps,group=factor(speedtox$week))
testTemporalAutocorrelation(tox_res_sp_rps_week,time=unique(speedtox$week))
#took care of temporal autocorrealtion

saveRDS(res_speed_rps_tox,paste0(results_dir,"res_speed_rps_tox.rds"))
