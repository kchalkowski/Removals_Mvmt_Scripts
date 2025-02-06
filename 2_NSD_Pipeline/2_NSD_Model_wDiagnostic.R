#set home dir of pipeline
home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline"

# Purpose ----------------------------------------------------------------------
#The purpose of this script is to formulate and run NSD glmms

# Process ----------------------------------------------------------------------

#Check for temporal autocorr
#Check for spatial autocorr
#GLMMs

# In/Out -----------------------------------------------------------------------

#inputs: geo_aerd_wk.rds, geo_toxd_wk.rds, geo_trapd_wk.rds

# Setup ------------------------------------------------------------------------

#Set dirs
input=file.path(home,"1_Data","Input",fsep=.Platform$file.sep)
objdir=file.path(home,"1_Data","Objects",fsep=.Platform$file.sep)
outdir=file.path(home,"3_Output",fsep=.Platform$file.sep)

#load libraries
library(ctmm)
library(ggplot2)
library(segclust2d)
library(NIFApackagev1.1)
library(stringr)
library(amt)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggeffects)
library(hrbrthemes)
library(glmmTMB)
library(DHARMa)
library(sdmTMB)
library(gt)
library(gtsummary)

#Read in data to use for NSD analysis:
geo.aerd.wk=readRDS(file.path(objdir,"NSDgeoaer.rds",fsep=.Platform$file.sep))
geo.trapd.wk=readRDS(file.path(objdir,"NSDgeotrap.rds",fsep=.Platform$file.sep))
geo.toxd.wk=readRDS(file.path(objdir,"NSDgeotox.rds",fsep=.Platform$file.sep))

#change colnames/trt naming
colnames(geo.aerd.wk)[c(2,3)]<-c("trt_ctrl","period")
colnames(geo.trapd.wk)[c(2,3)]<-c("trt_ctrl","period")
colnames(geo.toxd.wk)[c(2,3)]<-c("trt_ctrl","period")

geo.aerd.wk$trt_ctrl<-as.character(geo.aerd.wk$trt_ctrl)
geo.trapd.wk$trt_ctrl<-as.character(geo.trapd.wk$trt_ctrl)
geo.toxd.wk$trt_ctrl<-as.character(geo.toxd.wk$trt_ctrl)

geo.aerd.wk$trt_ctrl[geo.aerd.wk$trt_ctrl!="ctrl"]<-"trt"
geo.trapd.wk$trt_ctrl[geo.trapd.wk$trt_ctrl!="ctrl"]<-"trt"
geo.toxd.wk$trt_ctrl[geo.toxd.wk$trt_ctrl!="ctrl"]<-"trt"

geo.aerd.wk$trt_ctrl<-as.factor(geo.aerd.wk$trt_ctrl)
geo.trapd.wk$trt_ctrl<-as.factor(geo.trapd.wk$trt_ctrl)
geo.toxd.wk$trt_ctrl<-as.factor(geo.toxd.wk$trt_ctrl)

geo.aerd.wk$trt_ctrl<-forcats::fct_relevel(geo.aerd.wk$trt_ctrl,c("ctrl","trt"))
geo.trapd.wk$trt_ctrl<-forcats::fct_relevel(geo.trapd.wk$trt_ctrl,c("ctrl","trt"))
geo.toxd.wk$trt_ctrl<-forcats::fct_relevel(geo.toxd.wk$trt_ctrl,"ctrl","trt")


# Check for spatial autocorrelation --------------------------------------------

#Within each week, across pigs, do we see spatial autocorrelation of displacement values?

# * Aerial spatial autocorr ----------------------------------------------------

res1.ac.sac=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + trt_ctrl*period, data=geo.aerd.wk,family=Gamma(link=log))

#Test for spatial autocorrelation with dharma package
res <- simulateResiduals(res1.ac.sac)
groupLocations = aggregate(geo.aerd.wk[, 7:8], list(as.factor(geo.aerd.wk$animalid)), mean)
res2 = recalculateResiduals(res, group = as.factor(geo.aerd.wk$animalid), rotation="estimated")
testSpatialAutocorrelation(res2,groupLocations$mX, groupLocations$mY)
#No spatial autocorrelation, p=0.3733

res1.ac.sac=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + trt_ctrl*period*sex, data=geo.aerd.wk,family=Gamma(link=log))

#Test for spatial autocorrelation with dharma package
res <- simulateResiduals(res1.ac.sac)
groupLocations = aggregate(geo.aerd.wk[, 7:8], list(as.factor(geo.aerd.wk$animalid)), mean)
res2 = recalculateResiduals(res, group = as.factor(geo.aerd.wk$animalid), rotation="estimated")
testSpatialAutocorrelation(res2,groupLocations$mX, groupLocations$mY)
#No spatial autocorrelation, p=0.47

# * Trap spatial autocorr ----------------------------------------------------

#Run model
res1.ac.sac=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + trt_ctrl*period, data=geo.trapd.wk,family=Gamma(link=log))

#Test for spatial autocorrelation with dharma package
res <- simulateResiduals(res1.ac.sac) #mle-mvn?
groupLocations = aggregate(geo.trapd.wk[, 7:8], list(as.factor(geo.trapd.wk$animalid)), mean)
res2 = recalculateResiduals(res, group = as.factor(geo.trapd.wk$animalid), rotation="estimated")
testSpatialAutocorrelation(res2,groupLocations$mX, groupLocations$mY)
#is spatial autocorrelation for trap, p=9.262e-05

#Adjusting for spat autocor fit:
#adjust spat autocor type: mat, exp, gau (matern, exponential, gaussian)
#adjust grain (e.g., /100=100m resolution)

#Run model with spatial autocor structure
geo.trapd.wk$X=round(geo.trapd.wk$mX/1000)
geo.trapd.wk$Y=round(geo.trapd.wk$mY/1000)
geo.trapd.wk$animalid<-factor(geo.trapd.wk$animalid)
geo.trapd.wk$pos <- numFactor(geo.trapd.wk$X, geo.trapd.wk$Y)
res.spat=glmmTMB(mNSD~trt_ctrl*period+
                   exp(pos + 0 | animalid)+ar1(as.factor(week) + 0 | animalid),
                 data=geo.trapd.wk,
                 family=Gamma(link="log"),
                 verbose=TRUE)

#Test for spatial autocorrelation with dharma package
res <- simulateResiduals(res.spat) #mle-mvn?
groupLocations = aggregate(geo.trapd.wk[, 7:8], list(as.factor(geo.trapd.wk$animalid)), mean)
res2 = recalculateResiduals(res, group = as.factor(geo.trapd.wk$animalid), rotation="estimated")
testSpatialAutocorrelation(res2,groupLocations$mX, groupLocations$mY)
#0.05652 with /1000 and exp, took care of spatial autocorrelation

#Test trap model with sex interaction for spat autocorr
#Run model
res1.ac.sac=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + trt_ctrl*period*sex, data=geo.trapd.wk,family=Gamma(link=log))

#Test for spatial autocorrelation with dharma package
res <- simulateResiduals(res1.ac.sac)
groupLocations = aggregate(geo.trapd.wk[, 7:8], list(as.factor(geo.trapd.wk$animalid)), mean)
res2 = recalculateResiduals(res, group = as.factor(geo.trapd.wk$animalid), rotation="estimated")
testSpatialAutocorrelation(res2,groupLocations$mX, groupLocations$mY)
#is spatial autocorrelation for trap with sex interaction, p=.02253

geo.trapd.wk$X=round(geo.trapd.wk$mX/1000)
geo.trapd.wk$Y=round(geo.trapd.wk$mY/1000)
geo.trapd.wk$animalid<-factor(geo.trapd.wk$animalid)
geo.trapd.wk$pos <- numFactor(geo.trapd.wk$Y, geo.trapd.wk$Y)
res.spat=glmmTMB(mNSD~trt_ctrl*period*sex+
                   exp(pos + 0 | animalid)+ar1(as.factor(week) + 0 | animalid),
                 data=geo.trapd.wk,
                 family=Gamma(link="log"),
                 verbose=TRUE)

geo.trapd.wk$key=paste0(geo.trapd.wk$animalid,geo.trapd.wk$week)
res <- simulateResiduals(res.spat) 
groupLocations = aggregate(geo.trapd.wk[, 7:8], list(as.factor(geo.trapd.wk$animalid)), mean)
res2 = recalculateResiduals(res, group = as.factor(geo.trapd.wk$animalid), rotation="estimated")
testSpatialAutocorrelation(res2,groupLocations$mX, groupLocations$mY)
#/1000, Y-only, exp, p=0.16

# * Tox spatial autocorr ----------------------------------------------------

#Run model
res1.ac.sac=glmmTMB(mNSD ~ trt_ctrl*period, data=geo.toxd.wk,family=Gamma(link=log))

#Test for spatial autocorrelation with dharma package
res <- simulateResiduals(res1.ac.sac)
groupLocations = aggregate(geo.toxd.wk[, 7:8], list(as.factor(geo.toxd.wk$animalid)), mean)
res2 = recalculateResiduals(res, group = as.factor(geo.toxd.wk$animalid), rotation="estimated")
testSpatialAutocorrelation(res2,groupLocations$mX, groupLocations$mY)
#no spatial autocorrelation, p=0.3964

#Run model with sex interaction
res1.ac.sac=glmmTMB(mNSD ~ trt_ctrl*period*sex, data=geo.toxd.wk,family=Gamma(link=log))

#Test for spatial autocorrelation with dharma package
res <- simulateResiduals(res1.ac.sac)
groupLocations = aggregate(geo.toxd.wk[, 7:8], list(as.factor(geo.toxd.wk$animalid)), mean)
res2 = recalculateResiduals(res, group = as.factor(geo.toxd.wk$animalid), rotation="estimated")
testSpatialAutocorrelation(res2,groupLocations$mX, groupLocations$mY)
#is spatial autocorrelation, p=0.002

geo.toxd.wk$X=round(geo.toxd.wk$mX/1000)
geo.toxd.wk$Y=round(geo.toxd.wk$mY/1000)
geo.toxd.wk$animalid<-factor(geo.toxd.wk$animalid)
geo.toxd.wk$pos <- numFactor(geo.toxd.wk$X, geo.toxd.wk$Y)
res.spat=glmmTMB(mNSD~trt_ctrl*period*sex+
                   exp(pos + 0 | animalid),
                 data=geo.toxd.wk,
                 family=Gamma(link="log"),
                 verbose=TRUE)

res <- simulateResiduals(res.spat)
groupLocations = aggregate(geo.toxd.wk[, 7:8], list(as.factor(geo.toxd.wk$animalid)), mean)
res2 = recalculateResiduals(res, group = as.factor(geo.toxd.wk$animalid), rotation="estimated")
testSpatialAutocorrelation(res2,groupLocations$mX, groupLocations$mY)
#/1000 and exp removed spatial autocorr, p=0.18

# Autocorrelation summary -------------------

#Temporal autocorrelation:
  #all aerial models, all trap models
#Spatial autocorrelation:
  #trap spatial autocorr:
    #rp: 0.05652 with /1000 and exp
    #rps: /1000, Y-only, exp
  #tox autocorr:
    #rps: #/1000 and exp 

#Aerial glmms: all temporal only
#Trap glmms: all spatial and temporal
#Tox glmms: rps model needs spatial autocorr

# Run GLMMs --------------------------------------------------------------------

#set relevel for sex, check models for significance when male is ref
#geo.aerd.wk$sex<-forcats::fct_relevel(geo.aerd.wk$sex,c("Male","Female"))
#geo.trapd.wk$sex<-forcats::fct_relevel(geo.trapd.wk$sex,c("Male","Female"))
#geo.toxd.wk$sex<-forcats::fct_relevel(geo.toxd.wk$sex,c("Male","Female"))


# * Aerial GLMMs --------------------------------------------------------------------
geo.wkdf=geo.aerd.wk
geo.wkdf$period<-forcats::fct_relevel(geo.wkdf$period,c("before","after"))
#geo.wkdf$sex<-forcats::fct_relevel(geo.wkdf$sex,c("Female","Male"))
res.rp_aer=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + trt_ctrl*period, data=geo.wkdf,family=Gamma(link=log))
res.rps_aer=glmmTMB(mNSD ~ (1|animalid/week) + trt_ctrl*period*sex, data=geo.wkdf,family=Gamma(link=log))

#aerial rps model with temporal autocorrelation wouldn't converge

# * Trap GLMMs --------------------------------------------------------------------
#rp: 0.05652 with /1000 and exp
#rps: /1000, Y-only, exp
geo.wkdf=geo.trapd.wk
geo.wkdf$X=round(geo.wkdf$mX/1000)
geo.wkdf$Y=round(geo.wkdf$mY/1000)
geo.wkdf$animalid<-factor(geo.wkdf$animalid)
geo.wkdf$pos1 <- numFactor(geo.wkdf$X, geo.wkdf$Y)
geo.wkdf$pos2 <- numFactor(geo.wkdf$Y, geo.wkdf$Y)
geo.wkdf$trt_ctrl<-forcats::fct_relevel(geo.wkdf$trt_ctrl,c("ctrl","trap"))
geo.wkdf$period<-forcats::fct_relevel(geo.wkdf$period,c("before","during","after"))
#geo.wkdf$sex<-forcats::fct_relevel(geo.wkdf$sex,c("Female","Male"))
res.rp_trap=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + exp(pos1 + 0 | animalid) + trt_ctrl*period, 
                    data=geo.wkdf,
                    family=Gamma(link=log),
                    verbose=TRUE)
res.rps_trap=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + exp(pos2 + 0 | animalid) + trt_ctrl*period*sex, 
                     data=geo.wkdf,
                     family=Gamma(link=log),
                     verbose=TRUE)

# * Toxicant GLMMs --------------------------------------------------------------------
geo.wkdf=geo.toxd.wk
geo.wkdf$X=round(geo.wkdf$mX/1000)
geo.wkdf$Y=round(geo.wkdf$mY/1000)
geo.wkdf$animalid<-factor(geo.wkdf$animalid)
geo.wkdf$pos <- numFactor(geo.wkdf$X, geo.wkdf$Y)
#geo.wkdf$sex<-forcats::fct_relevel(geo.wkdf$sex,c("Female","Male"))
res.rp_tox=glmmTMB(mNSD ~ trt_ctrl*period, 
                   data=geo.wkdf,
                   family=Gamma(link=log),
                   verbose=TRUE)
res.rps_tox=glmmTMB(mNSD ~ exp(pos + 0 | animalid) + trt_ctrl*period*sex, 
                    data=geo.wkdf,
                    family=Gamma(link=log),
                    verbose=TRUE)


# Format model info -----------------------------------------------------------

mods<-list(res.rp_aer,
           res.rps_aer,
           res.rp_trap,
           res.rps_trap,
           res.rp_tox,
           res.rps_tox)
models=c("res.rp_aer",
            "res.rps_aer",
            "res.rp_trap",
            "res.rps_trap",
            "res.rp_tox",
            "res.rps_tox")

# * Make gt summary table----------------------

aer_tbl <- tbl_regression(res.rp_aer, exponentiate = TRUE)
trap_tbl <- tbl_regression(res.rp_trap, exponentiate = TRUE)
tox_tbl <- tbl_regression(res.rp_tox, exponentiate = TRUE)

nsd_tbl=tbl_merge(
  tbls = list(aer_tbl, 
              tox_tbl,
              trap_tbl),
  tab_spanner = c("aerial","tox","trap")
) 

saveRDS(nsd_tbl,file.path(outdir,"Model_Output","nsd_parm_gt.rds",fsep=.Platform$file.sep))

aer_tbl_s <- tbl_regression(res.rps_aer, exponentiate = TRUE)
trap_tbl_s <- tbl_regression(res.rps_trap, exponentiate = TRUE)
tox_tbl_s <- tbl_regression(res.rps_tox, exponentiate = TRUE)

nsd_tbl_s=tbl_merge(
  tbls = list(aer_tbl_s, 
              tox_tbl_s,
              trap_tbl_s),
  tab_spanner = c("aerial","tox","trap")
) 

saveRDS(nsd_tbl_s,file.path(outdir,"Model_Output","nsd_parm_gt_s.rds",fsep=.Platform$file.sep))


# * Format model prediction df -----------
for(i in 1:length(mods)){
  #i=1
  if(length(grep("rps",models[i]))==0){
    tmp=as.data.frame(ggeffects::predict_response(mods[[i]], terms=c("trt_ctrl","period")))
    tmp$facet=NA
  } else{
    tmp=as.data.frame(ggeffects::predict_response(mods[[i]], terms=c("trt_ctrl","period","sex")))
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
preds$response="nsd"

# Pull interactions table ------------------------------------------------------

#Loop through model objects
for(i in 1:length(mods)){
  
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
      
      if(i==1){
        allc=coefs2
      } else{
        allc=rbind(allc,coefs2)
      }
      

  
} #if nrow coef not 0
} #for loop

#period
allc$period=NA
allc$period[grep("after",allc$effect)]<-"after"
allc$period[grep("during",allc$effect)]<-"during"

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
allc$response="nsd"

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

# Save model output -----------------------------------------------------------

if(!dir.exists(file.path(outdir,"Model_Output"))){dir.create(file.path(outdir,"Model_Output"))}
saveRDS(preds,file.path(outdir,"Model_Output","nsd_preds.rds"))
saveRDS(allc,file.path(outdir,"Model_Output","nsd_intxns.rds"))
saveRDS(parms,file.path(outdir,"Model_Output","nsd_full_param.rds"))

