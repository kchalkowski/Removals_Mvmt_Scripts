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

# Setup ----------------

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
  library(gt)
  library(gtsummary)

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
  results_dir <- file.path(homedir,"3_Output",fsep=.Platform$file.sep)

  #pull geo_rem objects
  # georem <- read.csv("./1_Data/Objects/geo_remtyp_period.csv")
  geo.aer<-readRDS(file.path(objdir,"geoaer.rds"))
  geo.tox<-readRDS(file.path(objdir,"geotox.rds"))
  geo.trap<-readRDS(file.path(objdir,"geotrap.rds"))
  geo.all <- rbind(geo.aer,geo.tox,geo.trap)
  
  #change colnames to match all responses
  colnames(geo.aer)[12]<-"period"
  colnames(geo.tox)[12]<-"period"
  colnames(geo.trap)[12]<-"period"
  colnames(geo.all)[12]<-"period"
  
  #identify animals that died in toxicant treatment
  aktox=readRDS(paste0(objdir,"/outdf_akde_tox_corrected.rds"))
  #Remove NA after periods for tox-killed pigs
  aktox=aktox[!is.na(aktox$area.CI.low),]
  
  aktox.sumperiods=aktox %>% group_by(animalid,period) %>% dplyr::summarise(n()) %>% tidyr::pivot_wider(names_from="period",values_from=`n()`) %>% as.data.frame()
  died.tox=aktox.sumperiods[which(is.na(aktox.sumperiods$after)),1]
  #results_dir <- paste0(homedir,"/1_Data/Objects/Models/")
  
  #set mesh cutoff for all spatial models
  mesh_cutoff<-1
  spatial_res <- 1000

# Distance -------------------
## formatting distance dfs ----------
  #read distance
  dist<- readRDS(paste0(objdir,"/pig_weekly_distance_ctmm.rds"))
  
  #rename period
  colnames(dist)[4]<-"period"
  
  #subset dfs
  distaer <- dist %>% filter(Removal.Type%in%c('aer')) %>% filter(period!='during') 
  disttox <- dist %>% filter(Removal.Type%in%c('tox')) 
  disttrap <- dist %>% filter(Removal.Type%in%c('trap')) 
  
  #relevel when during isn't recorded
  distaer$period <- factor(distaer$period,
                                            levels=c('before','after'))
  disttox$period <- factor(disttox$period,
                                            levels=c('before','during','after'))
  disttrap$period <- factor(disttrap$period,
                                             levels=c('before','during','after'))
  
  #relevel treatment/ctrl-- ctrl needs be ref
  distaer$trt_ctrl<-forcats::fct_relevel(distaer$trt_ctrl,c("ctrl","trt"))
  disttrap$trt_ctrl<-forcats::fct_relevel(disttrap$trt_ctrl,c("ctrl","trt"))
  disttox$trt_ctrl<-forcats::fct_relevel(disttox$trt_ctrl,c("ctrl","trt"))
  
  #temp releveling to get significance of intxn when male is ref
  #distaer$sex<-forcats::fct_relevel(distaer$sex,c("Male","Female"))
  #disttrap$sex<-forcats::fct_relevel(disttrap$sex,c("Male","Female"))
  #disttox$sex<-forcats::fct_relevel(disttox$sex,c("Male","Female"))
  

## temporal autocorrelation test --------------
  res_aer=glmmTMB(weekly_dist_km~(1|animalid)+trt_ctrl*period,
                data=distaer,family=Gamma(link="log"))
simout_aer <- simulateResiduals(fittedModel = res_aer, plot = F)
res_aer2 = recalculateResiduals(simout_aer, group = distaer$week,rotation="estimated")
testTemporalAutocorrelation(res_aer2, time = unique(distaer$week))
#no temporal autocorrelation, p=0.054

#sex 
res_aer=glmmTMB(weekly_dist_km~(1|animalid)+trt_ctrl*period*sex,
                data=distaer,family=Gamma(link="log"))
simout_aer <- simulateResiduals(fittedModel = res_aer, plot = F)
res_aer2 = recalculateResiduals(simout_aer, group = distaer$week,rotation="estimated")
testTemporalAutocorrelation(res_aer2, time = unique(distaer$week))
#no temporal autocorrelation, p=0.055

res_trap=glmmTMB(weekly_dist_km~(1|animalid)+trt_ctrl*period,
                 data=disttrap,family=Gamma(link="log"))

simout_trap <- simulateResiduals(fittedModel = res_trap, plot = F)
res_trap2 = recalculateResiduals(simout_trap, group = disttrap$week,rotation="estimated")
testTemporalAutocorrelation(res_trap2, time = unique(disttrap$week))
#no temporal autocorrelation, p=0.29

#sex
res_trap=glmmTMB(weekly_dist_km~(1|animalid)+trt_ctrl*period*sex,
                 data=disttrap,family=Gamma(link="log"))
simout_trap <- simulateResiduals(fittedModel = res_trap, plot = F)
res_trap2 = recalculateResiduals(simout_trap, group = disttrap$week,rotation="estimated")
testTemporalAutocorrelation(res_trap2, time = unique(disttrap$week))
#no temporal autocorrelation, p=0.29

res_tox=glmmTMB(weekly_dist_km~(1|animalid)+trt_ctrl*period,
                data=disttox,family=Gamma(link="log"))
simout_tox <- simulateResiduals(fittedModel = res_tox, plot = F)
res_tox2 = recalculateResiduals(simout_tox, group = disttox$week,rotation="estimated")
testTemporalAutocorrelation(res_tox2, time = unique(disttox$week))
#no temporal autocorrelation, p=0.58

res_tox=glmmTMB(weekly_dist_km~(1|animalid)+trt_ctrl*period*sex,
                data=disttox,family=Gamma(link="log"))
simout_tox <- simulateResiduals(fittedModel = res_tox, plot = F)
res_tox2 = recalculateResiduals(simout_tox, group = disttox$week,rotation="estimated")
testTemporalAutocorrelation(res_tox2, time = unique(disttox$week))
#no temporal autocorrelation, p=0.27

## spatial autocorrelation test ------------------
distaer$X=floor(distaer$mX/5000)
distaer$Y=floor(distaer$mY/5000)
aer_test=glmmTMB(weekly_dist_km~(1|animalid)+trt_ctrl*period,
                 data=distaer,family=Gamma(link="log"))
aer_test_resid <- simulateResiduals(aer_test)
aer_groupLocations = aggregate(distaer[,c("mX","mY")], list(distaer$animalid), mean)
aer_test_resid2 = recalculateResiduals(aer_test_resid, group = distaer$animalid, rotation="estimated")
testSpatialAutocorrelation(aer_test_resid2,aer_groupLocations$mX, aer_groupLocations$mY)
# spatial autocorrelation, p=0.012

#sex
aer_test=glmmTMB(weekly_dist_km~(1|animalid)+trt_ctrl*period*sex,
                 data=distaer,family=Gamma(link="log"))
aer_test_resid <- simulateResiduals(aer_test)
aer_test_resid2 = recalculateResiduals(aer_test_resid, group = distaer$animalid, rotation="estimated")
aer_groupLocations = aggregate(distaer[,c("mX","mY")], list(distaer$animalid, distaer$sex), mean)
DHARMa::testSpatialAutocorrelation(aer_test_resid2,aer_groupLocations$mX, aer_groupLocations$mY)
# no spatial autocorrelation, p=0.99

trap_test=glmmTMB(weekly_dist_km~(1|animalid)+trt_ctrl*period,
                  data=disttrap,family=Gamma(link="log"))
trap_test_resid <- simulateResiduals(trap_test)
trap_groupLocations = aggregate(disttrap[,c("mX","mY")], list(disttrap$animalid), mean)
trap_test_resid2 = recalculateResiduals(trap_test_resid, group = disttrap$animalid, rotation="estimated")
testSpatialAutocorrelation(trap_test_resid2,trap_groupLocations$mX, trap_groupLocations$mY)
#no spatial autocorrelation, p=0.15

#sex 
trap_test=glmmTMB(weekly_dist_km~(1|animalid)+trt_ctrl*period*sex,
                  data=disttrap,family=Gamma(link="log"))
trap_test_resid <- simulateResiduals(trap_test)
trap_test_resid2 = recalculateResiduals(trap_test_resid, group = disttrap$animalid, rotation="estimated")
testSpatialAutocorrelation(trap_test_resid2,trap_groupLocations$mX, trap_groupLocations$mY)
#no spatial autocorrelation, p=0.17

tox_test=glmmTMB(weekly_dist_km~(1|animalid)+trt_ctrl*period,
                 data=disttox,family=Gamma(link="log"))
tox_test_resid <- simulateResiduals(tox_test)
tox_groupLocations = aggregate(disttox[,c("mX","mY")], list(disttox$animalid), mean)
tox_test_resid2 = recalculateResiduals(tox_test_resid, group = disttox$animalid, rotation="estimated")
testSpatialAutocorrelation(tox_test_resid2,tox_groupLocations$mX, tox_groupLocations$mY)
#spatial autocorrelation, p=0.0006

#sex
tox_test=glmmTMB(weekly_dist_km~(1|animalid)+trt_ctrl*period*sex,
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

## Fitting spatial autocorrelation ------------------

### aerial -----------------

#### removal type * period ------
spatial_res <- 100
mesh_cutoff=1
distaer$mX_sc <- floor(distaer$mX/spatial_res)
distaer$mY_sc <- floor(distaer$mY/spatial_res)
meshaer_sp <- make_mesh(distaer,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)
res_distance_rp_aer=sdmTMB(weekly_dist_km ~ (1|animalid) + trt_ctrl*period,
                        data=distaer,
                        mesh=meshaer_sp,
                        spatial='on',
                        family=Gamma(link='log'))
#finalmod

sanity(res_distance_rp_aer)
tox_res_rp <- simulate(res_distance_rp_aer, nsim = 544, type = "mle-mvn") %>% 
  dharma_residuals(res_distance_rp_aer, return_DHARMa = TRUE)
tox_res_rp2 = recalculateResiduals(tox_res_rp, group = as.factor(distaer$animalid),rotation="estimated")
groupLocations = aggregate(distaer[,c("mX","mY")], list(distaer$animalid), mean)
testSpatialAutocorrelation(tox_res_rp2,groupLocations$mX,groupLocations$mY)
#removed spatial autocorrelation, p=0.931

### toxicant -----------------

#### removal type * period ------
spatial_res <- 100
mesh_cutoff=2
disttox$mX_sc <- floor(disttox$mX/spatial_res)
disttox$mY_sc <- floor(disttox$mY/spatial_res)
meshtox_sp <- make_mesh(disttox,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)
res_distance_rp_tox=sdmTMB(weekly_dist_km ~ (1|animalid) + trt_ctrl*period,
                       data=disttox,
                       mesh=meshtox_sp,
                       spatial='on',
                       family=Gamma(link='log'))
#finalmod

sanity(res_distance_rp_tox)
tox_res_rp <- simulate(res_distance_rp_tox, nsim = 544, type = "mle-mvn") %>% 
  dharma_residuals(res_distance_rp_tox, return_DHARMa = TRUE)
tox_res_rp2 = recalculateResiduals(tox_res_rp, group = as.factor(disttox$animalid),rotation="estimated")
groupLocations = aggregate(disttox[,c("mX","mY")], list(disttox$animalid), mean)
testSpatialAutocorrelation(tox_res_rp2,groupLocations$mX,groupLocations$mY)
#removed spatial autocorrelation, p=0.106

# Speed -------------------

  speed<- readRDS(paste0(objdir,"/pig_weekly_speed_ctmm.rds"))
  
  #rename cols
  colnames(speed)[4]<-"period"

  speedaer <- speed %>% filter(Removal.Type%in%c('aer')) %>% 
    filter(period!='during') 
  speedtox <- speed %>% filter(Removal.Type%in%c('tox')) 
  speedtrap <- speed %>% filter(Removal.Type%in%c('trap')) 
  
  #relevel when during isn't recorded
  speedaer$period <- factor(speedaer$period,
                                             levels=c('before','after'))
  speedtox$period <- factor(speedtox$period,
                                             levels=c('before','during','after'))
  speedtrap$period <- factor(speedtrap$period,
                                              levels=c('before','during','after'))

  #relevel treatment/ctrl-- ctrl needs be ref
  speedaer$trt_ctrl<-forcats::fct_relevel(speedaer$trt_ctrl,c("ctrl","trt"))
  speedtrap$trt_ctrl<-forcats::fct_relevel(speedtrap$trt_ctrl,c("ctrl","trt"))
  speedtox$trt_ctrl<-forcats::fct_relevel(speedtox$trt_ctrl,c("ctrl","trt"))
  
  #temp releveling to get significance of intxn when male is ref
  #speedaer$sex<-forcats::fct_relevel(speedaer$sex,c("Male","Female"))
  #speedtrap$sex<-forcats::fct_relevel(speedtrap$sex,c("Male","Female"))
  #speedtox$sex<-forcats::fct_relevel(speedtox$sex,c("Male","Female"))
  
  
## spatial autocorrelation testing ------------------
## aerial
aer_tests=glmmTMB(weekly_md_km_hr~(1|animalid)+trt_ctrl*period,
                  data=speedaer,family=Gamma(link='log'))
aer_test_resids <- simulateResiduals(aer_tests)
aer_groupLocations_sp = aggregate(speedaer[,c("mX","mY")], list(speedaer$animalid), mean)
aer_test_resids2 = recalculateResiduals(aer_test_resids, group = speedaer$animalid, rotation="estimated")
testSpatialAutocorrelation(aer_test_resids2,aer_groupLocations_sp$mX, aer_groupLocations_sp$mY)
# no spatial autocorrelation

#sex
aer_tests=glmmTMB(weekly_md_km_hr~(1|animalid)+trt_ctrl*period*sex,
                  data=speedaer,family=Gamma(link='log'))
aer_test_resids <- simulateResiduals(aer_tests)
aer_test_resids2 = recalculateResiduals(aer_test_resids, group = speedaer$animalid, rotation="estimated")
testSpatialAutocorrelation(aer_test_resids2,aer_groupLocations_sp$mX, aer_groupLocations_sp$mY)
# no spatial autocorrelation, p=0.93

trap_tests=glmmTMB(weekly_md_km_hr~(1|animalid)+trt_ctrl*period,
                   data=speedtrap,family=Gamma(link='log'))
trap_test_resids <- simulateResiduals(trap_tests)
trap_groupLocations_sp = aggregate(speedtrap[,c("mX","mY")], list(speedtrap$animalid), mean)
trap_test_resids2 = recalculateResiduals(trap_test_resids, group = speedtrap$animalid, rotation="estimated")
testSpatialAutocorrelation(trap_test_resids2,trap_groupLocations_sp$mX, trap_groupLocations_sp$mY)
#no spatial autocorrelation, p=0.69

#sex
trap_tests=glmmTMB(weekly_md_km_hr~(1|animalid)+trt_ctrl*period*sex,
                   data=speedtrap,family=Gamma(link='log'))
trap_test_resids <- simulateResiduals(trap_tests)
trap_test_resids2 = recalculateResiduals(trap_test_resids, group = speedtrap$animalid, rotation="estimated")
testSpatialAutocorrelation(trap_test_resids2,trap_groupLocations_sp$mX, trap_groupLocations_sp$mY)
#no spatial autocorrelation, p=0.86

tox_tests=glmmTMB(weekly_md_km_hr~(1|animalid)+trt_ctrl*period,
                  data=speedtox,family=Gamma(link='log'))
tox_test_resids <- simulateResiduals(tox_tests)
tox_groupLocations_sp = aggregate(speedtox[,c("mX","mY")], list(speedtox$animalid), mean)
tox_test_resids2 = recalculateResiduals(tox_test_resids, group = speedtox$animalid, rotation="estimated")
testSpatialAutocorrelation(tox_test_resids2,tox_groupLocations_sp$mX, tox_groupLocations_sp$mY)
# spatial autocorrelation, p=1.079e-07

tox_tests=glmmTMB(weekly_md_km_hr~(1|animalid)+trt_ctrl*period*sex,
                  data=speedtox,family=Gamma(link='log'))
tox_test_resids <- simulateResiduals(tox_tests)
tox_test_resids2 = recalculateResiduals(tox_test_resids, group = speedtox$animalid, rotation="estimated")
testSpatialAutocorrelation(tox_test_resids2,tox_groupLocations_sp$mX, tox_groupLocations_sp$mY)
# spatial autocorrelation, p=0.0001852

## speed autocorrelation conclusions --------
  #Spatial autocorrelation:
    #tox removal*period model
    #tox removal*period*sex model

## Fitting spatial autocorrelation ------------------

### toxicant ----------------
#set mesh cutoff for spatial models
spatial_res <- 100
mesh_cutoff=1
speedtox$mX_sc <- floor(speedtox$mX/spatial_res)
speedtox$mY_sc <- floor(speedtox$mY/spatial_res)

#### removal type * period ------
meshtox_sp <- make_mesh(speedtox,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)
res_speed_rp_tox=sdmTMB(weekly_md_km_hr ~ (1|animalid) + trt_ctrl*period,
                        data=speedtox,
                       mesh=meshtox_sp,
                        spatial='on',
                        family=Gamma(link='log'))
#finalmod

sanity(res_speed_rp_tox)
tox_res_rp <- simulate(res_speed_rp_tox, nsim = 544, type = "mle-mvn") %>% 
  dharma_residuals(res_speed_rp_tox, return_DHARMa = TRUE)
tox_res_rp2 = recalculateResiduals(tox_res_rp, group = as.factor(speedtox$animalid),rotation="estimated")
groupLocations = aggregate(speedtox[,c("mX","mY")], list(speedtox$animalid), mean)
testSpatialAutocorrelation(tox_res_rp2,groupLocations$mX,groupLocations$mY)
#removed spatial autocorrelation, p=0.5807

#### removal type * period * sex ------
spatial_res <- 100
mesh_cutoff=1
speedtox$mX_sc <- floor(speedtox$mX/spatial_res)
speedtox$mY_sc <- floor(speedtox$mY/spatial_res)
meshtox_sp <- make_mesh(speedtox,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)
res_speed_rps_tox=sdmTMB(weekly_md_km_hr ~ (1|animalid) + 
                           trt_ctrl*period*sex,
                        data=speedtox,
                        mesh=meshtox_sp,
                        spatial='on',
                        family=Gamma(link='log'))
#finalmod

sanity(res_speed_rp_tox)
tox_res_rp <- simulate(res_speed_rp_tox, nsim = 544, type = "mle-mvn") %>% 
  dharma_residuals(res_speed_rp_tox, return_DHARMa = TRUE)
tox_res_rp2 = recalculateResiduals(tox_res_rp, group = as.factor(speedtox$animalid),rotation="estimated")
groupLocations = aggregate(speedtox[,c("mX","mY")], list(speedtox$animalid), mean)
testSpatialAutocorrelation(tox_res_rp2,groupLocations$mX,groupLocations$mY)
#removed spatial autocorrelation, p=0.6789

# Contacts - number --------
## Formatting ncon dfs -------------------------------------
### aerial ------------
conaer <- readRDS(paste0(objdir,"/pairwise_contacts_aer.rds"))

#Rename period to match all responses
colnames(conaer)[c(3,2)]<-c("period","trt_ctrl")

conaer$trt_ctrl <- 
  factor(conaer$trt_ctrl,levels=c('ctrl','trt'))
conaer$period <- 
  factor(conaer$period,levels=c('before','after'))

#set to 10m
conaer <- conaer %>% filter(dist==10)
# 
# geo.aer$period <- 
#   factor(geo.aer$period,levels=c('before','after'))

#join mean X/Y locs to conaer
conaer <- conaer %>% left_join(
  geo.aer %>% 
    dplyr::group_by(animalid,period) %>% 
    dplyr::summarise(mX=mean(X),
                     mY=mean(Y)))

conaer$animalid <- 
  factor(conaer$animalid)
conaer$period <- 
  factor(conaer$period,levels=c('before','after'))


#set negligible offset to get gamma to work
conaer$contacts_per_day_offset <- conaer$contacts_per_day+0.0001

#get summary of 0 contacts
conaer %>% dplyr::group_by(contacts_per_day==0) %>% dplyr::summarise(n=n())

###trap ------
contrap <- readRDS(paste0(objdir,"/pairwise_contacts_trap.rds"))
colnames(contrap)[c(3,2)]<-c("period","trt_ctrl")

contrap$trt_ctrl <- 
  factor(contrap$trt_ctrl,levels=c('ctrl','trt'))
contrap$period <- 
  factor(contrap$period,levels=c('before','during','after'))

contrap <- contrap %>% filter(dist==10)

contrap <- contrap %>% left_join(
  geo.trap %>% 
    dplyr::group_by(animalid) %>% 
    dplyr::summarise(mX=mean(X),
                     mY=mean(Y)))

contrap$trt_ctrl <- 
  factor(contrap$trt_ctrl,levels=c('ctrl','trt'))
contrap$animalid <- 
  factor(contrap$animalid)
contrap$period <- 
  factor(contrap$period,levels=c('before','during','after'))


contrap$contacts_per_day_offset <- contrap$contacts_per_day+0.0001

contrap %>%  dplyr::group_by(contacts_per_day==0) %>%  dplyr::summarise(n=n())


###tox ------
contox <- readRDS(paste0(objdir,"/pairwise_contacts_tox.rds"))
colnames(contox)[c(3,2)]<-c("period","trt_ctrl")

contox$trt_ctrl <- 
  factor(contox$trt_ctrl,levels=c('ctrl','trt'))
contox$period <- 
  factor(contox$period,levels=c('before','during','after'))

contox <- contox %>% filter(dist==10)

contox <- contox %>% left_join(
  geo.tox %>% 
    dplyr::group_by(animalid) %>% 
    dplyr::summarise(mX=mean(X),
                     mY=mean(Y)))

contox$animalid <- 
  factor(contox$animalid)
contox$trt_ctrl <- 
  factor(contox$trt_ctrl,levels=c('ctrl','trt'))
contox$period <- 
  factor(contox$period,levels=c('before','during','after'))
contox$contacts_per_day_offset <- contox$contacts_per_day+0.0001

#remove animals that died
contox <- contox %>% filter(!is.na(num_contacts))

contox %>% dplyr::group_by(contacts_per_day==0) %>% dplyr::summarise(n=n())

## Testing spatial autocorrelation ------------------

### aerial ----
#### removal type * period -----
res_ncon_rp_aer=glmmTMB(contacts_per_day_offset~(1|animalid)+
                          trt_ctrl*period,
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

res_ncon_rp_aer = sdmTMB(contacts_per_day_offset~
                              (1|animalid)+
                               trt_ctrl*period,
                             data=conaer,
                             family=Gamma(link="log"),
                             mesh=meshaer_con,
                             spatial="on")
#finalmod

sanity(res_ncon_rp_aer)
aer_res_rp <- simulate(res_ncon_rp_aer, 
                          nsim = 544,type = 'mle-mvn') %>%
  dharma_residuals(res_ncon_rp_aer, return_DHARMa = TRUE)
aer_res_rp2 = recalculateResiduals(aer_res_rp,
                                      group = as.factor(conaer$animalid),
                                      rotation="estimated")
testSpatialAutocorrelation(aer_res_rp2,
                           aer_groupLocations_con$mX,
                           aer_groupLocations_con$mY)


#### removal type * period *sex -----
#### model not run - spatially autocorrelated by spatial model won't converge ####
# 
# res_ncon_rps_aer=glmmTMB(contacts_per_day_offset~
#                            (1|animalid)+
#                            trt_ctrl*period*sex,
#                          data=conaer,
#                          family=Gamma(link="log")
# )
# 
# aer_res_ncon_rps <- simulateResiduals(res_ncon_rps_aer)
# aer_res_ncon_rps2 = recalculateResiduals(aer_res_ncon_rps,
#                                         group = as.factor(conaer$animalid),
#                                         rotation="estimated")
# plot(aer_res_ncon_rps2)
# plot(aer_res_ncon_rps2$simulatedResponse)
# descdist(aer_res_ncon_rps2$fittedResiduals)
# DHARMa::testDispersion(aer_res_ncon_rps2)
# 
# #test spatial autocorrelation
# testSpatialAutocorrelation(aer_res_ncon_rps2,
#                            aer_groupLocations_con$mX,
#                            aer_groupLocations_con$mY)
# spatial autocorrelation p=0.02

#correcting for spatial autocorrelation
### NOTE: model overparameterized, very poor fit
# spatial_res <- 1000
# mesh_cutoff=0.5
# conaer$mX_sc <- floor(conaer$mX/spatial_res)
# conaer$mY_sc <- floor(conaer$mY/spatial_res)
# meshaer_con <- make_mesh(conaer,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)
# 
# res_ncon_rps_aer = sdmTMB(contacts_per_day_offset~
#                                (1|animalid)+
#                               trt_ctrl*period*sex,
#                             data=conaer,
#                             family=Gamma(link="log"),
#                             mesh=meshaer_con,
#                             spatial="on")
# 
# sanity(res_ncon_rps_aer)
# aer_res_rps <- simulate(res_ncon_rps_aer,
#                           nsim = 544,type = 'mle-mvn') %>%
#   dharma_residuals(res_ncon_rp_aer, return_DHARMa = TRUE)
# aer_res_rps2 = recalculateResiduals(aer_res_rps,
#                                       group = as.factor(conaer$animalid),
#                                       rotation="estimated")
# testSpatialAutocorrelation(aer_res_rps2,
#                            aer_groupLocations_con$mX,
#                            aer_groupLocations_con$mY)

### trap ----

#### removal type * period -----
res_ncon_rp_trap=glmmTMB(contacts_per_day_offset~(1|animalid)+
                           trt_ctrl*period,
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

res_ncon_rp_trap = sdmTMB(contacts_per_day_offset~(1|animalid)+
                              trt_ctrl*period,
                            data=contrap,
                            family=Gamma(link="log"),
                            mesh=meshtrap_con,
                            spatial="on")
#finalmod

sanity(res_ncon_rp_trap)
trap_res_rp <- simulate(res_ncon_rp_trap,
                          nsim = 544,type = 'mle-mvn') %>%
  dharma_residuals(res_ncon_rp_trap, return_DHARMa = TRUE)
trap_res_rp2 = recalculateResiduals(trap_res_rp,
                                      group = as.factor(contrap$animalid),
                                      rotation="estimated")
testSpatialAutocorrelation(trap_res_rp2,
                           trap_groupLocations_con$mX,
                           trap_groupLocations_con$mY)
#spatial autocorrelation p=0.85

#### removal type * period *sex -----
res_ncon_rps_trap=glmmTMB(contacts_per_day_offset~
                            (1|animalid)+
                            trt_ctrl*period*sex,
                          data=contrap,
                          family=Gamma(link="log")
                          )
#finalmod

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

### tox ----

#### removal type * period -----
res_ncon_rp_tox=glmmTMB(contacts_per_day_offset~
                          (1|animalid)+
                          trt_ctrl*period,
                        data=contox,
                        family=Gamma(link="log")
                        )
#finalmod

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
#no spatial autocorrelation p=0.09

#### removal type * period *sex -----
res_ncon_rps_tox=glmmTMB(contacts_per_day_offset~
                           (1|animalid)+
                           trt_ctrl*period*sex,
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
#spatial autocorrelation p=0.001

#correcting for spatial autocorrelation
spatial_res <- 1000
mesh_cutoff=0.5
contox$mX_sc <- floor(contox$mX/spatial_res)
contox$mY_sc <- floor(contox$mY/spatial_res)
meshtox_con <- make_mesh(contox,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)

res_ncon_rps_tox = sdmTMB(contacts_per_day_offset~(1|animalid)+
                             trt_ctrl*period*sex,
                           data=contox,
                           family=Gamma(link="log"),
                           mesh=meshtox_con,
                           spatial="on")
#finalmod

sanity(res_ncon_rps_tox)
tox_res_rps <- simulate(res_ncon_rps_tox,
                        nsim = 544,type = 'mle-mvn') %>%
  dharma_residuals(res_ncon_rps_tox, return_DHARMa = TRUE)
tox_res_rps2 = recalculateResiduals(tox_res_rps,
                                    group = as.factor(contox$animalid),
                                    rotation="estimated")
testSpatialAutocorrelation(tox_res_rps2,
                           tox_groupLocations_con$mX,
                           tox_groupLocations_con$mY)
#spatial autocorrelation p=0.53

#Set model list to pull into model info dfs at end of script
cn_model_list=list(res_ncon_rp_tox,
               res_ncon_rps_tox,
               res_ncon_rps_trap,
               res_ncon_rp_trap,
               # res_ncon_rps_aer, #model had spatal autocorr but couldn't fit correction
               res_ncon_rp_aer
               )

# Contacts - degree --------
## Formating degree contacts df --------
conaer$indivs_per_day_offset <- conaer$indivs_per_day+0.0001
contrap$indivs_per_day_offset <- contrap$indivs_per_day+0.0001
contox$indivs_per_day_offset <- contox$indivs_per_day+0.0001

## Testing spatial autocorrelation ------------
### aerial ------
#### removal type * period -----
res_nind_rp_aer=glmmTMB(indivs_per_day_offset~
                          (1|animalid)+
                          trt_ctrl*period,
                        data=conaer,
                        family=Gamma(link="log")
                        )
#finalmod

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


#### removal type * period *sex -----
res_nind_rps_aer=glmmTMB(indivs_per_day_offset~
                           (1|animalid)+
                           trt_ctrl*period*sex,
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
# spatial autocorrelation p=0.02

#correcting for spatial autocorrelation
spatial_res <- 1000
mesh_cutoff=2
conaer$mY_sc <- floor(conaer$mY/spatial_res)
meshtaer_con <- make_mesh(conaer,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)

res_nind_rps_aer = sdmTMB(indivs_per_day_offset~
                            (1|animalid)+
                            trt_ctrl*period*sex,
                          data=conaer,
                          family=Gamma(link="log"),
                          mesh=meshtaer_con,
                          spatial="on")
#finalmod

sanity(res_nind_rps_aer)
aer_res_nind_rps <- simulate(res_nind_rps_aer,
                             nsim = 544,type = 'mle-mvn') %>%
  dharma_residuals(res_nind_rps_aer, return_DHARMa = TRUE)
aer_res_nind_rps2 = recalculateResiduals(aer_res_nind_rps,
                                         group = as.factor(conaer$animalid),
                                         rotation="estimated")
testSpatialAutocorrelation(aer_res_nind_rps2,
                           aer_groupLocations_con$mX,
                           aer_groupLocations_con$mY)
#no spatial autocorrelation p=0.08


###trap ------

#### removal type * period -----
res_nind_rp_trap=glmmTMB(indivs_per_day_offset~
                           (1|animalid)+
                           trt_ctrl*period,
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

res_nind_rp_trap = sdmTMB(indivs_per_day_offset~
                               (1|animalid)+
                               trt_ctrl*period,
                             data=contrap,
                             family=Gamma(link="log"),
                             mesh=meshtrap_con,
                             spatial="on")
#finalmod

sanity(res_nind_rp_trap)
trap_res_nind_rp <- simulate(res_nind_rp_trap,
                           nsim = 544,type = 'mle-mvn') %>%
  dharma_residuals(res_nind_rp_trap, return_DHARMa = TRUE)
trap_res_nind_rp2 = recalculateResiduals(trap_res_nind_rp,
                                       group = as.factor(contrap$animalid),
                                       rotation="estimated")
testSpatialAutocorrelation(trap_res_nind_rp2,
                           trap_groupLocations_con$mX,
                           trap_groupLocations_con$mY)
#no spatial autocorrelation p=0.33

#### removal type * period *sex -----
res_nind_rps_trap=glmmTMB(indivs_per_day_offset~
                            (1|animalid)+
                            trt_ctrl*period*sex,
                          data=contrap,
                          family=Gamma(link="log")
                          )
#finalmod

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


###tox ------

#### removal type * period -----
res_nind_rp_tox=glmmTMB(indivs_per_day_offset~
                          (1|animalid)+
                          trt_ctrl*period,
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
#spatial autocorrelation p=0.0003

#correcting for spatial autocorrelation
spatial_res <- 1000
mesh_cutoff=0.1
contox$mX_sc <- floor(contox$mX/spatial_res)
contox$mY_sc <- floor(contox$mY/spatial_res)
meshtox_con <- make_mesh(contox,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)

res_nind_rp_tox = sdmTMB(indivs_per_day_offset~
                               (1|animalid)+
                               trt_ctrl*period,
                             data=contox,
                             family=Gamma(link="log"),
                             mesh=meshtox_con,
                             spatial="on")
#finalmod

sanity(res_nind_rp_tox)
tox_res_nind_rp <- simulate(res_nind_rp_tox,
                                nsim = 544,type = 'mle-mvn') %>%
  dharma_residuals(res_nind_rp_tox, return_DHARMa = TRUE)
tox_res_nind_rp2 = recalculateResiduals(tox_res_nind_rp,
                                            group = as.factor(contox$animalid),
                                            rotation="estimated")
testSpatialAutocorrelation(tox_res_nind_rp2,
                           tox_groupLocations_con$mX,
                           tox_groupLocations_con$mY)
#no spatial autocorrelation p=0.11

#### removal type * period *sex -----
res_nind_rps_tox=glmmTMB(indivs_per_day_offset~
                           (1|animalid)+
                           trt_ctrl*period*sex,
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
#spatial autocorrelation p=0.0004

#correcting for spatial autocorrelation
spatial_res <- 1000
mesh_cutoff=0.1
contox$mX_sc <- floor(contox$mX/spatial_res)
contox$mY_sc <- floor(contox$mY/spatial_res)
meshtox_con <- make_mesh(contox,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)

res_nind_rps_tox = sdmTMB(indivs_per_day_offset~
                              (1|animalid)+
                              trt_ctrl*period*sex,
                            data=contox,
                            family=Gamma(link="log"),
                            mesh=meshtox_con,
                            spatial="on")
#finalmod

sanity(res_nind_rps_tox)
tox_res_nind_rps <- simulate(res_nind_rps_tox,
                               nsim = 544,type = 'mle-mvn') %>%
  dharma_residuals(res_nind_rps_tox, return_DHARMa = TRUE)
tox_res_nind_rps2 = recalculateResiduals(tox_res_nind_rps,
                                           group = as.factor(contox$animalid),
                                           rotation="estimated")
testSpatialAutocorrelation(tox_res_nind_rps2,
                           tox_groupLocations_con$mX,
                           tox_groupLocations_con$mY)
#no spatial autocorrelation p=0.18

# Run all finalmods together ------------------------------------------

## Distance GLMMs ----------------------------------------

### Aerial distance rp ----------------------------------------
spatial_res <- 100
mesh_cutoff=1
distaer$mX_sc <- floor(distaer$mX/spatial_res)
distaer$mY_sc <- floor(distaer$mY/spatial_res)
meshaer_sp <- make_mesh(distaer,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)
res_distance_rp_aer=sdmTMB(weekly_dist_km ~ (1|animalid) + trt_ctrl*period,
                           data=distaer,
                           mesh=meshaer_sp,
                           spatial='on',
                           family=Gamma(link='log'))

### Aerial distance rps ----------------------------------------
res_distance_rps_aer <-glmmTMB(weekly_dist_km ~ (1|animalid) +
                                 trt_ctrl*period*sex,
                               data=distaer,family=Gamma(link="log"))
### Trap distance rp ----------------------------------------
res_distance_rp_trap=glmmTMB(weekly_dist_km ~ trt_ctrl*period+
                               (1|animalid),
                             data=disttrap,family=Gamma(link="log"))
### Trap distance rps ----------------------------------------
res_distance_rps_trap=glmmTMB(weekly_dist_km ~ trt_ctrl*period*sex+
                                (1|animalid),
                              data=disttrap,family=Gamma(link="log"))
### Tox distance rp ----------------------------------------
spatial_res <- 100
mesh_cutoff=2
disttox$mX_sc <- floor(disttox$mX/spatial_res)
disttox$mY_sc <- floor(disttox$mY/spatial_res)
meshtox_sp <- make_mesh(disttox,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)
res_distance_rp_tox=sdmTMB(weekly_dist_km ~ (1|animalid) + trt_ctrl*period,
                           data=disttox,
                           mesh=meshtox_sp,
                           spatial='on',
                           family=Gamma(link='log'))
### Tox distance rps ----------------------------------------
res_distance_rps_tox=glmmTMB(weekly_dist_km ~ trt_ctrl*period*sex+
                               (1|animalid),
                             data=disttox,family=Gamma(link="log"))

## Speed GLMMs ----------------------------------------

### Aer speed rp ----------------------------------------

res_speed_rp_aer <- glmmTMB(weekly_md_km_hr ~ trt_ctrl*period+
                              (1|animalid),
                            data=speedaer,
                            family=Gamma(link="log"),
                            verbose=TRUE)
### Aer speed rps ----------------------------------------

res_speed_rps_aer <- glmmTMB(weekly_md_km_hr ~ trt_ctrl*period*sex+
                               (1|animalid),
                             data=speedaer,
                             family=Gamma(link="log"),
                             verbose=TRUE)
### Trap speed rp ----------------------------------------
res_speed_rp_trap=glmmTMB(weekly_md_km_hr~ 
                            (1|animalid)+
                            trt_ctrl*period,
                          data=speedtrap,family=Gamma(link='log'))

### Trap speed rps ----------------------------------------
res_speed_rps_trap=glmmTMB(weekly_md_km_hr~ 
                             (1|animalid)+
                             trt_ctrl*period*sex,
                           data=speedtrap,family=Gamma(link='log'))
### Tox speed rp ----------------------------------------
spatial_res <- 100
mesh_cutoff=1
speedtox$mX_sc <- floor(speedtox$mX/spatial_res)
speedtox$mY_sc <- floor(speedtox$mY/spatial_res)
meshtox_sp <- make_mesh(speedtox,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)
res_speed_rp_tox=sdmTMB(weekly_md_km_hr ~ (1|animalid) + trt_ctrl*period,
                        data=speedtox,
                        mesh=meshtox_sp,
                        spatial='on',
                        family=Gamma(link='log'))
### Tox speed rps ----------------------------------------
spatial_res <- 100
mesh_cutoff=1
speedtox$mX_sc <- floor(speedtox$mX/spatial_res)
speedtox$mY_sc <- floor(speedtox$mY/spatial_res)
meshtox_sp <- make_mesh(speedtox,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)
res_speed_rps_tox=sdmTMB(weekly_md_km_hr ~ (1|animalid) + 
                           trt_ctrl*period*sex,
                         data=speedtox,
                         mesh=meshtox_sp,
                         spatial='on',
                         family=Gamma(link='log'))

## Ncon GLMMs ----------------------------------------

### Aer ncon rp ----------------------------------------
spatial_res <- 1000
mesh_cutoff=0.5
conaer$mX_sc <- floor(conaer$mX/spatial_res)
conaer$mY_sc <- floor(conaer$mY/spatial_res)
meshaer_con <- make_mesh(conaer,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)
res_ncon_rp_aer = sdmTMB(contacts_per_day_offset~
                           (1|animalid)+
                           trt_ctrl*period,
                         data=conaer,
                         family=Gamma(link="log"),
                         mesh=meshaer_con,
                         spatial="on")
### Aer ncon rps ----------------------------------------
#Removed, spatially autocorrelated but wouldn't converge

### Trap ncon rp ----------------------------------------
spatial_res <- 1000
mesh_cutoff=0.5
contrap$mX_sc <- floor(contrap$mX/spatial_res)
contrap$mY_sc <- floor(contrap$mY/spatial_res)
meshtrap_con <- make_mesh(contrap,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)
res_ncon_rp_trap = sdmTMB(contacts_per_day_offset~(1|animalid)+
                            trt_ctrl*period,
                          data=contrap,
                          family=Gamma(link="log"),
                          mesh=meshtrap_con,
                          spatial="on")
### Trap ncon rps ----------------------------------------
res_ncon_rps_trap=glmmTMB(contacts_per_day_offset~
                            (1|animalid)+
                            trt_ctrl*period*sex,
                          data=contrap,
                          family=Gamma(link="log")
)
### Tox ncon rp ----------------------------------------
res_ncon_rp_tox=glmmTMB(contacts_per_day_offset~
                          (1|animalid)+
                          trt_ctrl*period,
                        data=contox,
                        family=Gamma(link="log")
)
### Tox ncon rps ----------------------------------------
spatial_res <- 1000
mesh_cutoff=0.5
contox$mX_sc <- floor(contox$mX/spatial_res)
contox$mY_sc <- floor(contox$mY/spatial_res)
meshtox_con <- make_mesh(contox,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)
res_ncon_rps_tox = sdmTMB(contacts_per_day_offset~(1|animalid)+
                            trt_ctrl*period*sex,
                          data=contox,
                          family=Gamma(link="log"),
                          mesh=meshtox_con,
                          spatial="on")

## Nind GLMMs ----------------------------------------

### Aer nind rp ----------------------------------------
res_nind_rp_aer=glmmTMB(indivs_per_day_offset~
                          (1|animalid)+
                          trt_ctrl*period,
                        data=conaer,
                        family=Gamma(link="log")
)
### Aer nind rps ----------------------------------------
spatial_res <- 1000
mesh_cutoff=2
conaer$mY_sc <- floor(conaer$mY/spatial_res)
meshtaer_con <- make_mesh(conaer,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)
res_nind_rps_aer = sdmTMB(indivs_per_day_offset~
                            (1|animalid)+
                            trt_ctrl*period*sex,
                          data=conaer,
                          family=Gamma(link="log"),
                          mesh=meshtaer_con,
                          spatial="on")
### Trap nind rp ----------------------------------------
spatial_res <- 1000
mesh_cutoff=0.5
contrap$mX_sc <- floor(contrap$mX/spatial_res)
contrap$mY_sc <- floor(contrap$mY/spatial_res)
meshtrap_con <- make_mesh(contrap,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)
res_nind_rp_trap = sdmTMB(indivs_per_day_offset~
                            (1|animalid)+
                            trt_ctrl*period,
                          data=contrap,
                          family=Gamma(link="log"),
                          mesh=meshtrap_con,
                          spatial="on")
### Trap nind rps ----------------------------------------
res_nind_rps_trap=glmmTMB(indivs_per_day_offset~
                            (1|animalid)+
                            trt_ctrl*period*sex,
                          data=contrap,
                          family=Gamma(link="log")
)
### Tox nind rp ----------------------------------------
spatial_res <- 1000
mesh_cutoff=0.1
contox$mX_sc <- floor(contox$mX/spatial_res)
contox$mY_sc <- floor(contox$mY/spatial_res)
meshtox_con <- make_mesh(contox,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)
res_nind_rp_tox = sdmTMB(indivs_per_day_offset~
                           (1|animalid)+
                           trt_ctrl*period,
                         data=contox,
                         family=Gamma(link="log"),
                         mesh=meshtox_con,
                         spatial="on")
### Tox nind rps ----------------------------------------
spatial_res <- 1000
mesh_cutoff=0.1
contox$mX_sc <- floor(contox$mX/spatial_res)
contox$mY_sc <- floor(contox$mY/spatial_res)
meshtox_con <- make_mesh(contox,c("mX_sc","mY_sc"),cutoff=mesh_cutoff)
res_nind_rps_tox = sdmTMB(indivs_per_day_offset~
                            (1|animalid)+
                            trt_ctrl*period*sex,
                          data=contox,
                          family=Gamma(link="log"),
                          mesh=meshtox_con,
                          spatial="on")

# Format model info -----------------------------------------------------------

ci_model_list=list(res_nind_rp_aer,
                   res_nind_rps_aer,
                   res_nind_rp_trap,
                   res_nind_rps_trap,
                   res_nind_rp_tox,
                   res_nind_rps_tox)

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

cn_mods=cn_model_list
cn_models=c("res_ncon_rp_tox",
                "res_ncon_rps_tox",
                "res_ncon_rps_trap",
                "res_ncon_rp_trap",
                # "res_ncon_rps_aer", #taking out of running, wouldnt converge w spat autocorr
                "res_ncon_rp_aer")

ci_mods=ci_model_list
ci_models=c("res_nind_rp_aer",
               "res_nind_rps_aer",
               "res_nind_rp_trap",
               "res_nind_rps_trap",
               "res_nind_rp_tox",
               "res_nind_rps_tox")

#pull in contact models
mods=c(mods,cn_mods,ci_mods)
models=c(models,cn_models,ci_models)

##  Make gt summary tables----------------------

#source funcs to make gt tables for sdmTMB models
source(file.path(homedir,"2_Scripts","Functions","Table_Convert_Funcs.R",fsep=.Platform$file.sep))

### distance gt summary tables----------------------
aer_tbl<-make_sdmTMB_gt(res_distance_rp_aer,distaer)
trap_tbl <- tbl_regression(res_distance_rp_trap, exponentiate = TRUE)
tox_tbl<-make_sdmTMB_gt(res_distance_rp_tox,disttox)

dist_tbl=tbl_merge(
  tbls = list(aer_tbl,
              tox_tbl, 
              trap_tbl),
  tab_spanner = c("aerial","tox","trap")
) 

saveRDS(dist_tbl,file.path(results_dir,"Model_Output","dist_parm_gt.rds",fsep=.Platform$file.sep))

aer_tbl_s <- tbl_regression(res_distance_rps_aer, exponentiate = TRUE)
trap_tbl_s <- tbl_regression(res_distance_rps_trap, exponentiate = TRUE)
tox_tbl_s <- tbl_regression(res_distance_rps_tox, exponentiate = TRUE)

dist_tbl_s=tbl_merge(
  tbls = list(aer_tbl_s,
              tox_tbl_s,
              trap_tbl_s),
  tab_spanner = c("aerial","tox","trap")
) 

saveRDS(dist_tbl_s,file.path(results_dir,"Model_Output","dist_parm_gt_s.rds",fsep=.Platform$file.sep))

### speed gt summary tables----------------------
aer_tbl <- tbl_regression(res_speed_rp_aer, exponentiate = TRUE)
trap_tbl <- tbl_regression(res_speed_rp_trap, exponentiate = TRUE)
tox_tbl<-make_sdmTMB_gt(res_speed_rp_tox,speedtox)

speed_tbl=tbl_merge(
  tbls = list(aer_tbl, 
              tox_tbl,
              trap_tbl),
  tab_spanner = c("aerial","tox","trap")
) 

saveRDS(speed_tbl,file.path(results_dir,"Model_Output","speed_parm_gt.rds",fsep=.Platform$file.sep))


aer_tbl_s <- tbl_regression(res_speed_rps_aer, exponentiate = TRUE)
trap_tbl_s <- tbl_regression(res_speed_rps_trap, exponentiate = TRUE)
tox_tbl<-make_sdmTMB_gt(res_speed_rps_tox,speedtox)

speed_tbl_s=tbl_merge(
  tbls = list(aer_tbl_s, 
              tox_tbl_s,
              trap_tbl_s),
  tab_spanner = c("aerial","tox","trap")
) 

saveRDS(speed_tbl_s,file.path(results_dir,"Model_Output","speed_parm_gt_s.rds",fsep=.Platform$file.sep))


### ncon gt summary tables----------------------
aer_tbl <- make_sdmTMB_gt(res_ncon_rp_aer,conaer)
trap_tbl <- make_sdmTMB_gt(res_ncon_rp_trap, contrap)
tox_tbl <- tbl_regression(res_ncon_rp_tox, exponentiate = TRUE)

ncon_tbl=tbl_merge(
  tbls = list(tox_tbl,
              trap_tbl, 
              aer_tbl),
  tab_spanner = c("tox","trap","aer")
) 

saveRDS(ncon_tbl,file.path(results_dir,"Model_Output","ncon_parm_gt.rds",fsep=.Platform$file.sep))

#aer_tbl_s <- make_sdmTMB_gt(res_ncon_rps_aer, conaer)
trap_tbl_s <- tbl_regression(res_ncon_rps_trap, exponentiate = TRUE)
tox_tbl_s <- make_sdmTMB_gt(res_ncon_rps_tox, contox)

ncon_tbl_s=tbl_merge(
  tbls = list(#aer_tbl,
              tox_tbl_s, 
              trap_tbl_s),
  tab_spanner = c(#"aerial",
                  "tox","trap")
) 

saveRDS(ncon_tbl_s,file.path(results_dir,"Model_Output","ncon_parm_gt_s.rds",fsep=.Platform$file.sep))

### nind gt summary tables----------------------
aer_tbl <- tbl_regression(res_nind_rp_aer, exponentiate = TRUE)
trap_tbl <- make_sdmTMB_gt(res_nind_rp_trap, contrap)
tox_tbl <- make_sdmTMB_gt(res_nind_rp_tox, contox)

nind_tbl=tbl_merge(
  tbls = list(tox_tbl, 
              trap_tbl,
              aer_tbl),
  tab_spanner = c("tox","trap","aer")
) 

saveRDS(nind_tbl,file.path(results_dir,"Model_Output","nind_parm_gt.rds",fsep=.Platform$file.sep))


aer_tbl_s <- make_sdmTMB_gt(res_nind_rps_aer, conaer)
trap_tbl_s <- tbl_regression(res_nind_rps_trap, exponentiate = TRUE)
tox_tbl_s <- make_sdmTMB_gt(res_nind_rps_tox, contox)

nind_tbl_s=tbl_merge(
  tbls = list(tox_tbl_s,
              trap_tbl_s,
              aer_tbl_s),
  tab_spanner = c("tox","trap","aer")
) 

saveRDS(nind_tbl_s,file.path(results_dir,"Model_Output","nind_parm_gt_s.rds",fsep=.Platform$file.sep))

# * Format model prediction df -----------------------------------------------
for(i in 1:length(mods)){

  if(length(grep("rps",models[i]))==0){
    #add if statement here and inside else-- if ncon or nind in model name, trt_ctrl is trt_ctrl
    if(length(grep("ncon",models[i]))!=0|
       length(grep("nind",models[i]))!=0){
      tmp=as.data.frame(ggeffects::predict_response(mods[[i]], terms=c("trt_ctrl","period")))
    } else{
      tmp=as.data.frame(ggeffects::predict_response(mods[[i]], terms=c("trt_ctrl","period")))
    }
    tmp$facet=NA
  } else{
    if(length(grep("ncon",models[i]))!=0|
       length(grep("nind",models[i]))!=0){
      tmp=as.data.frame(ggeffects::predict_response(mods[[i]], terms=c("trt_ctrl","period","sex")))
      } else{
        tmp=as.data.frame(ggeffects::predict_response(mods[[i]], terms=c("trt_ctrl","period","sex")))
      }
  }
  
  
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
preds$response[grep("nind",preds$model)]<-"nind"
preds$response[grep("ncon",preds$model)]<-"ncon"


# * Pull interactions table ------------------------------------------------------
for(i in 1:length(mods)){
  
  if(class(mods[[i]])=="sdmTMB"){
    df=summary(mods[[i]]$sd_report,p.value=TRUE) %>% as.data.frame()
    df=df[grep("b_",rownames(df)),]
    names = broom.mixed::tidy(mods[[i]]) %>% as.data.frame() %>% select(term) 
    df$effect = names$term
    df$model=models[i]
    df=df[,c(1,2,4,6,5)]
    df=df[grep(":",df$effect),]
    colnames(df)=c("Estimate","Std. Error","Pr(>|z|)","model","effect")

  } else{
    #Pull coefficients
    coefs=summary(mods[[i]])$coef$cond
    
    #Make table with just interaction effect, std.error, p value
    coefs2=as.data.frame(coefs[grep(":",rownames(coefs)),c(1,2,4),drop=FALSE])
    if(nrow(coefs2)!=0){
      
      if(nrow(coefs2)>1){
        if(length(grep("ncon",models[i]))!=0|
           length(grep("nind",models[i]))!=0){
          coefs2=as.data.frame(coefs2[grep("trt_ctrl",rownames(coefs2)),,drop=FALSE])
          coefs2=as.data.frame(coefs2[grep("period",rownames(coefs2)),,drop=FALSE])
        } else{
          coefs2=as.data.frame(coefs2[grep("trt_ctrl",rownames(coefs2)),,drop=FALSE])
          coefs2=as.data.frame(coefs2[grep("period",rownames(coefs2)),,drop=FALSE])
          
        }
        
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
allc$response[grep("ncon",allc$model)]<-"ncon"
allc$response[grep("nind",allc$model)]<-"nind"

#Remove ugly effect column
allc=allc[,-which(colnames(allc)=="effect")]

#change sex NA
allc$sex[is.na(allc$sex)]<-"whole"

#subset to only significant interactions
allc$alpha=0.1
allc$alpha[allc$`Pr(>|z|)`<0.05]<-1

# * Combine full parameter tables -----------------------------------------------------------

for(i in 1:length(mods)){
  if(class(mods[[i]])=="sdmTMB"){
    df=summary(mods[[i]]$sd_report,p.value=TRUE) %>% as.data.frame()
    df=df[grep("b_",rownames(df)),]
    names = broom.mixed::tidy(mods[[i]]) %>% as.data.frame() %>% select(term) 
    df$term = names$term
    df$model=models[i]
    df=df[,c("term","Estimate","Std. Error","z value","Pr(>|z^2|)")]
    colnames(df) <- c("term","estimate","std.error","statistic","p.value")
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

#check to make sure p vals calc'd correctly
range(parms$p.value)

# * Save tidied model output ----------------
#outdir=file.path(homedir,"3_Output",fsep=.Platform$file.sep)
if(!dir.exists(file.path(results_dir,"Model_Output"))){dir.create(file.path(outdir,"Model_Output"))}
saveRDS(preds,file.path(results_dir,"Model_Output","spdist_preds.rds"))
saveRDS(allc,file.path(results_dir,"Model_Output","spdist_intxns.rds"))
saveRDS(parms,file.path(results_dir,"Model_Output","spdist_full_param.rds"))

