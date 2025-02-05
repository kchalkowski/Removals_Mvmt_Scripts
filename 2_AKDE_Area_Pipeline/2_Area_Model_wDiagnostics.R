#set home dir of pipeline
home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline"

# Purpose ----------------------------------------------------------------------
#The purpose of this script is run autocorrelation diagnostics and glmms for akde home range size estimates

# In/Out -----------------------------------------------------------------------

#inputs:
#outdf_akde_aer_corrected_f.rds, outdf_akde_tox_corrected_f.rds, outdf_akde_trap_corrected_f.rds

#outputs:

# Setup ------------------------------------------------------------------------

#load libraries
library(glmmTMB)
library(DHARMa)
library(ctmm)
library(fitdistrplus)
library(ggplot2)
library(plyr)
library(dplyr)
library(hrbrthemes)
library(sf)
library(mapview)
library(gt)
library(gtsummary)

#Set dirs
input=file.path(home,"1_Data","Input",fsep=.Platform$file.sep)
objdir=file.path(home,"1_Data","Objects",fsep=.Platform$file.sep)
if(!dir.exists(file.path(home,"3_Output","Area_GLM_Results",fsep=.Platform$file.sep))){
  dir.create(file.path(home,"3_Output","Area_GLM_Results",fsep=.Platform$file.sep))
}
outdir=file.path(home,"3_Output",fsep=.Platform$file.sep)

#Read input objects
akaer=readRDS(file.path(objdir,"outdf_akde_aer_corrected_f.rds",fsep=.Platform$file.sep))
aktrap=readRDS(file.path(objdir,"outdf_akde_trap_corrected_f.rds",fsep=.Platform$file.sep))
aktox=readRDS(file.path(objdir,"outdf_akde_aktox_corrected_f.rds",fsep=.Platform$file.sep))

colnames(akaer)[19]<-"trt_ctrl"
colnames(aktrap)[19]<-"trt_ctrl"
colnames(aktox)[19]<-"trt_ctrl"


#source needed functions
funcdir<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/2_Scripts/Functions"
func.list=list.files(funcdir,full.names=TRUE)
for(f in 1:(length(func.list)-1)){
  source(func.list[f])
}

# Area df object data formatting --------------------------------------------------

#fix units-- sq meters/hectares-- want all to be sq km
Fix.Area.Units<-function(akrem){
  akrem[grep("square meters",akrem$akde.hrarea.units),][3:5]=akrem[grep("square meters",akrem$akde.hrarea.units),][3:5]/1e6
  akrem[grep("hectares",akrem$akde.hrarea.units),][3:5]=akrem[grep("hectares",akrem$akde.hrarea.units),][3:5]/100
  return(akrem)
}

aktox=Fix.Area.Units(aktox)
aktrap=Fix.Area.Units(aktrap)
akaer=Fix.Area.Units(akaer)

#Remove pigs that were killed during aerial
killed.aer=akaer[is.na(akaer$area.CI.low),]$animalid
akaer=akaer[!(akaer$animalid%in%killed.aer),]

#Remove NA 'after' periods for tox-killed pigs
aktox=aktox[!is.na(aktox$area.CI.low),]

#Relevel periods
akaer$period<-forcats::fct_relevel(akaer$period,c("before","after"))
aktrap$period<-forcats::fct_relevel(aktrap$period,c("before","during","after"))
aktox$period<-forcats::fct_relevel(aktox$period,c("before","after"))

#change levels in trt_ctrl
akaer$trt_ctrl<-as.character(akaer$trt_ctrl)
aktrap$trt_ctrl<-as.character(aktrap$trt_ctrl)
aktox$trt_ctrl<-as.character(aktox$trt_ctrl)

akaer$trt_ctrl[akaer$trt_ctrl!="ctrl"]<-"trt"
aktrap$trt_ctrl[aktrap$trt_ctrl!="ctrl"]<-"trt"
aktox$trt_ctrl[aktox$trt_ctrl!="ctrl"]<-"trt"

#Relevel removal types
akaer$trt_ctrl<-forcats::fct_relevel(akaer$trt_ctrl,c("ctrl","trt"))
aktrap$trt_ctrl<-forcats::fct_relevel(aktrap$trt_ctrl,c("ctrl","trt"))
aktox$trt_ctrl<-forcats::fct_relevel(aktox$trt_ctrl,c("ctrl","trt"))

#Relevel sex, check male significance level for results figure
#akaer$sex<-forcats::fct_relevel(akaer$sex,c("Male","Female"))
#aktrap$sex<-forcats::fct_relevel(aktrap$sex,c("Male","Female"))
#aktox$sex<-forcats::fct_relevel(aktox$sex,c("Male","Female"))

# Check aerial data spatial autocorrelation -----------------------------------
{
res=glmmTMB(area.est~(1|animalid)+trt_ctrl*period,data=akaer,family=Gamma(link="log"))
ress <- simulateResiduals(res)
groupLocations = aggregate(akaer[, 6:7], list(akaer$animalid), mean)
ress2 = recalculateResiduals(ress, group = akaer$animalid)
#testSpatialAutocorrelation(ress2,groupLocations$ctr.x, groupLocations$ctr.y)
#no spatial autocorrelation

#Test model fit
simulationOutput <- simulateResiduals(fittedModel = res, plot = F)
so2=recalculateResiduals(simulationOutput, group = akaer$animalid)
#plot(so2)
#descdist(so2$fittedResiduals)
}

# Check trap data spatial autocorrelation -----------------------------------
{
#Look at distribution
#hist(aktrap$area.est)
#descdist(aktrap$area.est)
#Gamma seems like best option
res=glmmTMB(area.est~(1|animalid)+trt_ctrl*period,data=aktrap,family=Gamma(link="log"))
ress <- simulateResiduals(res)
groupLocations = aggregate(aktrap[, 6:7], list(aktrap$animalid), mean)
ress2 = recalculateResiduals(ress, group = aktrap$animalid, rotation="estimated")
#testSpatialAutocorrelation(ress2,groupLocations$ctr.x, groupLocations$ctr.y)
#is spatial autocorrelation

#Correct for spatial autocorrelation, try gaussian autocor structure
aktrap$pos <- numFactor(aktrap$ctr.x, aktrap$ctr.y)
res.sa=glmmTMB(area.est~trt_ctrl*period+exp(pos + 0 | animalid),data=aktrap,family=Gamma(link="log"))
#*Note, gaussian autocor structure wouldnt converge, going with exp
#*Note, inclusion of animalid as random intercept AND as group in spatial autocorr structure also will not converge. Removing animalid random intercept.

#Need correct for spatial autocorrelation
simulationOutput <- simulateResiduals(fittedModel = res.sa, plot = F)
so2=recalculateResiduals(simulationOutput, group = aktrap$animalid, rotation="estimated")
#testSpatialAutocorrelation(so2,groupLocations$ctr.x, groupLocations$ctr.y)
#spatial autocorrelation no longer significant, exp structure takes care of it

#descdist(so2$fittedResiduals)
#gamma still seems like best option
}

#if just running models, run this to get GLMMs:
{
#Correct for spatial autocorrelation, try gaussian autocor structure
aktrap$pos <- numFactor(aktrap$ctr.x, aktrap$ctr.y)
res.sa=glmmTMB(area.est~trt_ctrl*period+exp(pos + 0 | animalid),data=aktrap,family=Gamma(link="log"))
#*Note, gaussian autocor structure wouldnt converge, going with exp
#*Note, inclusion of animalid as random intercept AND as group in spatial autocorr structure also will not converge. Removing animalid random intercept.
}

# Check tox data spatial autocorrelation -----------------------------------
{
#####Do some summaries of home range areas used in analysis
aktox.sumperiods=aktox %>% group_by(animalid,period) %>% dplyr::summarise(n()) %>% tidyr::pivot_wider(names_from="period",values_from=`n()`) %>% as.data.frame()
died.tox=aktox.sumperiods[which(is.na(aktox.sumperiods$after)),1]
aktox$died.tox=0
aktox[aktox$animalid%in%died.tox,]$died.tox=1

#Look at distribution
#hist(aktox$area.est)
#descdist(aktox$area.est)
#Gamma seems like best option

#make sure died tox is factor
aktox$died.tox<-as.factor(aktox$died.tox)
aktox$died.tox<-forcats::fct_relevel(aktox$died.tox, c("0","1"))
aktox2=aktox[aktox$died.tox==0,]

#Check for spatial autocorrelation
res=glmmTMB(area.est~(1|animalid)+trt_ctrl*period*died.tox,data=aktox,family=Gamma(link="log"))
ress <- simulateResiduals(res)
groupLocations = aggregate(aktox[, 6:7], list(aktox$animalid), mean)
ress2 = recalculateResiduals(ress, group = aktox$animalid, rotation="estimated")
#testSpatialAutocorrelation(ress2,groupLocations$ctr.x, groupLocations$ctr.y)
#No spatial autocorrelation when died.tox is included-- opt to not include spatial autocorr in model, is assoc. with died. tox
}

# Run GLMMs ------------------------------------------------------------

# Aerial
res.rp_aer=glmmTMB(area.est ~ trt_ctrl*period + (1|animalid), data=akaer,family=Gamma(link=log))
res.rps_aer=glmmTMB(area.est ~ trt_ctrl*period*sex + (1|animalid), data=akaer,family=Gamma(link=log))

# Trap
res.rp_trap=glmmTMB(area.est ~ trt_ctrl*period + exp(pos + 0 | animalid), data=aktrap,family=Gamma(link=log))
#Note, model with trt_ctrl*period+sex would not converge, so excluding from model list

# Tox
res.rp_tox=glmmTMB(area.est ~ trt_ctrl*period + (1|animalid), data=aktox,family=Gamma(link=log))
res.rps_tox=glmmTMB(area.est ~ trt_ctrl*period*sex + (1|animalid), data=aktox,family=Gamma(link=log))


# Format model info -----------------------------------------------------------

mods<-list(res.rp_aer,
           res.rps_aer,
           res.rp_trap,
           res.rp_tox,
           res.rps_tox)
models=c("res.rp_aer",
         "res.rps_aer",
         "res.rp_trap",
         "res.rp_tox",
         "res.rps_tox")

## Make gt summary table -----------

aer_tbl <- tbl_regression(res.rp_aer, exponentiate = TRUE)
trap_tbl <- tbl_regression(res.rp_trap, exponentiate = TRUE)
tox_tbl <- tbl_regression(res.rp_tox, exponentiate = TRUE)

area_tbl=tbl_merge(
  tbls = list(aer_tbl, 
              trap_tbl,
              tox_tbl),
  tab_spanner = c("aerial","trap","tox")
) 

saveRDS(area_tbl,file.path(outdir,"Model_Output","area_parm_gt.rds",fsep=.Platform$file.sep))

aer_tbl_s <- tbl_regression(res.rps_aer, exponentiate = TRUE)
#trap_tbl_s <- tbl_regression(res.rps_trap, exponentiate = TRUE)
tox_tbl_s <- tbl_regression(res.rps_tox, exponentiate = TRUE)

area_tbl_s=tbl_merge(
  tbls = list(aer_tbl_s,
              tox_tbl_s),
  tab_spanner = c("aerial","tox")
) 

saveRDS(area_tbl_s,file.path(outdir,"Model_Output","area_parm_gt_s.rds",fsep=.Platform$file.sep))

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
preds$response="area"

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
allc$trt[grep("aer",allc$effect)]<-"aer"
allc$trt[grep("trap",allc$effect)]<-"trap"
allc$trt[grep("tox",allc$effect)]<-"tox"

#get response type
allc$response="area"

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
  #df$response="area"
  
  if(i==1){
    parms=df
  } else{
    parms=rbind(parms,df)
  }
  
}

if(!dir.exists(file.path(outdir,"Model_Output"))){dir.create(file.path(outdir,"Model_Output"))}
saveRDS(preds,file.path(outdir,"Model_Output","area_preds.rds"))
saveRDS(allc,file.path(outdir,"Model_Output","area_intxns.rds"))
saveRDS(parms,file.path(outdir,"Model_Output","area_full_param.rds"))


