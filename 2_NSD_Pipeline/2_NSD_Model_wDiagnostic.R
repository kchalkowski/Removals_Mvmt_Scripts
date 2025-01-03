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

#Read in data to use for NSD analysis:
geo.aerd.wk=readRDS(file.path(objdir,"NSDgeoaer.rds",fsep=.Platform$file.sep))
geo.trapd.wk=readRDS(file.path(objdir,"NSDgeotrap.rds",fsep=.Platform$file.sep))
geo.toxd.wk=readRDS(file.path(objdir,"NSDgeotox.rds",fsep=.Platform$file.sep))

# Temporal (weekly) autocorrelation checking -----------------------------------

#Functions
Do.PACF<-function(ys,len){
  
  #### PACF
  PACF = 0          # Starting up an empty storage vector.
  for(j in 2:len){   # Picked up 25 lag points to parallel R `pacf()` output.
    cols = j        
    rows = length(ys) - j + 1 # To end up with equal length vectors we clip.
    
    lag = matrix(0, rows, j)    # The storage matrix for different groups of lagged vectors.
    
    for(i in 1:cols){
      lag[ ,i] = ys[i : (i + rows - 1)]  #Clipping progressively to get lagged ts's.
    }
    lag = as.data.frame(lag)
    fit = lm(lag$V1 ~ . - 1, data = lag) # Running an OLS for every group.
    PACF[j] = coef(fit)[j - 1]           # Getting the slope for the last lagged ts.
  }
  
  PACF[1]=1 #first should be 1.. since it will be zero lag, corr. with itself
  
  #Make dataframe showing lags explicitly
  PACF.df=data.frame(lag=seq(0,(length(PACF)-1)),PACF=PACF)
  return(PACF.df)
}

# * Check partial autocorrelation plots ----------------------------------------

#Requires PACF.df, a dataframe with two cols
#lag=intervals to plot along, starting with zero through length of the data
#PACF is the PACF results
Plot.PACF<-function(PACF.df){
  q <- ggplot(data = PACF.df, mapping = aes(x = lag, y = PACF)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0))+
    geom_hline(aes(yintercept=1.96/sqrt(nrow(PACF.df))),linetype="dashed")+
    geom_hline(aes(yintercept=(-1.96/sqrt(nrow(PACF.df)))),linetype="dashed")+
    theme_ipsum()
}

Do.PACF.Exploratory<-function(geo.remd.wk,outdir,removal.str){
  #Want a PACF for each animalid.... 
  #subset to only needed cols
  #animalid, week, responses
  id.col=which(colnames(geo.remd.wk)=="animalid")
  wk.col=which(colnames(geo.remd.wk)=="week")
  response=which(colnames(geo.remd.wk)=="mNSD")
  ys.ID=geo.remd.wk[,c(id.col,wk.col,response)] 
  ys.ID=ys.ID %>% nest_by(animalid) #nest within animalid
  
  #Start with areas_km
  #do offset by 10 for each
  PACF.df.list<-vector(mode="list",length=nrow(ys.ID))
  for(i in 1:nrow(ys.ID)){
    print(i)
    #if(nrow(ys.ID$data[[i]])>10){
    ID.len=min(10,(floor(length(ys.ID$data[[i]]$mNSD)/2)))
    PACF.df.list[[i]]=Do.PACF(ys.ID$data[[i]]$mNSD,len=ID.len)
    PACF.df.list[[i]]$animalid=ys.ID$animalid[[i]]
    #} else{PACF.df.list[[i]]<-NA}
  }
  
  #PACF.df.list
  #which(is.na(PACF.df.list)) #check
  
  PACF.df=dplyr::bind_rows(PACF.df.list)
  IDs=unique(PACF.df$animalid)
  for(i in 1:length(IDs)){
    PACF.df.ID=PACF.df[PACF.df$animalid==IDs[i],]
    ggplot(data = PACF.df.ID, mapping = aes(x = lag, y = PACF)) +
      geom_hline(aes(yintercept = 0)) +
      geom_segment(mapping = aes(xend = lag, yend = 0))+
      geom_hline(aes(yintercept=1.96/sqrt(length(ys.ID$data[[i]]$mNSD))),linetype="dashed")+
      geom_hline(aes(yintercept=(-1.96/sqrt(length(ys.ID$data[[i]]$mNSD)))),linetype="dashed")+
      theme_ipsum()+theme(legend.position="none")+
      labs(title=ys.ID$animalid[i])
    if(!exists(paste0(outdir,removal.str))){dir.create(paste0(outdir,removal.str))}
    filename=paste0(paste0(outdir,removal.str),"/",IDs[i],".png")
    ggsave(filename,bg="white")
  }
  
}

outdir.pacf=file.path(outdir,"Exploratory_Outputs","PACF_animalid_plots_akde_displacement")
Do.PACF.Exploratory(geo.aerd.wk,outdir.pacf,"aer")
Do.PACF.Exploratory(geo.toxd.wk,outdir.pacf,"tox")
Do.PACF.Exploratory(geo.trapd.wk,outdir.pacf,"trap")

#Looked at plots, does appear to be some first-order autocorrelation in displacement between weeks
#Try adding autocorrelation structure to model, see if residuals become less autocorrelated across weeks

# * Check residual plots -------------------------------------------------------

#Make function to run different autocorrelation structures, check residual plots
Check.Week.Autocorr<-function(geo.aerd.wk,outdir,removal.str){
  
  geo.aerd.wk$weekfac=as.factor(geo.aerd.wk$week)
  
  res0=glmmTMB(mNSD ~ 1, data=geo.aerd.wk,family=Gamma(link=log))
  res1=glmmTMB(mNSD ~ ar1(weekfac + 0 | animalid), data=geo.aerd.wk,family=Gamma(link=log))
  
  df.plot=data.frame("animalid"=geo.aerd.wk$animalid,"week"=geo.aerd.wk$week,"res"=residuals(res0))
  IDs=unique(df.plot$animalid)
  for(i in 1:length(unique(df.plot$animalid))){
    ggplot(data=df.plot[df.plot$animalid==unique(df.plot$animalid)[i],],aes(x=week,y=res,group=1))+
      geom_line()
    if(!exists(paste0(outdir,"week_residuals_null/",removal.str,"/"))){dir.create(paste0(outdir,"week_residuals_null/",removal.str,"/"))}
    filename=paste0(paste0(outdir,"week_residuals_null/",removal.str,"/"),IDs[i],".png")
    ggsave(filename,bg="white")
  }
  
  df.plot=data.frame("animalid"=geo.aerd.wk$animalid,"week"=geo.aerd.wk$week,"res"=residuals(res1))
  IDs=unique(df.plot$animalid)
  for(i in 1:length(unique(df.plot$animalid))){
    ggplot(data=df.plot[df.plot$animalid==unique(df.plot$animalid)[i],],aes(x=week,y=res,group=1))+
      geom_line()
    if(!exists(paste0(outdir,"week_residuals_ac1/",removal.str,"/"))){dir.create(paste0(outdir,"week_residuals_ac1/",removal.str,"/"))}
    filename=paste0(paste0(outdir,"week_residuals_ac1/",removal.str,"/"),IDs[i],".png")
    ggsave(filename,bg="white")
  }
  
}

outdir.autocor=file.path(outdir,"Exploratory_Outputs")
Check.Week.Autocorr(geo.aerd.wk,outdir.autocor,"aer")
Check.Week.Autocorr(geo.trapd.wk,outdir.autocor,"trap")
Check.Week.Autocorr(geo.toxd.wk,outdir.autocor,"tox")

#Try models with ar1 structure
res.ac.aer=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid), data=geo.aerd.wk,family=Gamma(link=log))
res.ac.trap=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid), data=geo.trapd.wk,family=Gamma(link=log))
res.ac.tox=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid), data=geo.toxd.wk,family=Gamma(link=log))

#Test models with ar1 structure
simulationOutput <- simulateResiduals(fittedModel = res.ac.aer, plot = F)
so2=recalculateResiduals(simulationOutput, group = geo.aerd.wk$animalid)
plot(so2$simulatedResponse)
descdist(so2$fittedResiduals)
DHARMa::testDispersion(so2)
plot(so2)

simulationOutput <- simulateResiduals(fittedModel = res.ac.trap, plot = F)
so2=recalculateResiduals(simulationOutput, group = geo.trapd.wk$animalid)
plot(so2$simulatedResponse)
descdist(so2$fittedResiduals)
DHARMa::testDispersion(so2)
plot(so2)

simulationOutput <- simulateResiduals(fittedModel = res.ac.tox, plot = F)
so2=recalculateResiduals(simulationOutput, group = geo.toxd.wk$animalid)
plot(so2$simulatedResponse)
descdist(so2$fittedResiduals)
DHARMa::testDispersion(so2)
plot(so2)

#Including autocor structure does seem to reduce autocor of residuals considerably.

# Check distribution fits ------------------------------------------------------

# * Check histograms ----------------------------------------

#Check histos to determine appropriate model structures
Do.Displacement.Histos<-function(wak2,outdir.hist,removal.str){
  IDs=unique(wak2$animalid)
  for(i in 1:length(unique(wak2$animalid))){
    ggplot(data=wak2[wak2$animalid==unique(wak2$animalid)[i],],aes(x=mNSD))+
      geom_histogram()
    if(!exists(paste0(outdir,removal.str,"/"))){dir.create(paste0(outdir,removal.str,"/"))}
    filename=paste0(paste0(outdir,removal.str,"/"),IDs[i],".png")
    ggsave(filename,bg="white")
  }
}

outdir.hist=file.path(outdir,"Exploratory_Outputs","weekly_disp_hist")
Do.Displacement.Histos(geo.aerd.wk,outdir.hist,"aer")
Do.Displacement.Histos(geo.toxd.wk,outdir.hist,"tox")
Do.Displacement.Histos(geo.trapd.wk,outdir.hist,"trap")

# Check for spatial autocorrelation --------------------------------------------

#Within each week, across pigs, do we see spatial autocorrelation of displacement values?

# * Aerial spatial autocorr ----------------------------------------------------

res1.ac.sac=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + Removal.Type*removal.period.akdecalc, data=geo.aerd.wk,family=Gamma(link=log))

#Test for spatial autocorrelation with dharma package
res <- simulateResiduals(res1.ac.sac)
groupLocations = aggregate(geo.aerd.wk[, 7:8], list(as.factor(geo.aerd.wk$animalid)), mean)
res2 = recalculateResiduals(res, group = as.factor(geo.aerd.wk$animalid), rotation="estimated")
testSpatialAutocorrelation(res2,groupLocations$mX, groupLocations$mY)
#Is spatial autocorr

# * Trap spatial autocorr ----------------------------------------------------

#Run model
res1.ac.sac=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + Removal.Type*removal.period.akdecalc, data=geo.trapd.wk,family=Gamma(link=log))

#Test for spatial autocorrelation with dharma package
res <- simulateResiduals(res1.ac.sac) #mle-mvn?
groupLocations = aggregate(geo.trapd.wk[, 7:8], list(as.factor(geo.trapd.wk$animalid)), mean)
res2 = recalculateResiduals(res, group = as.factor(geo.trapd.wk$animalid), rotation="estimated")
testSpatialAutocorrelation(res2,groupLocations$mX, groupLocations$mY)
#is spatial autocorrelation for trap, p=0.027

#Adjusting for spat autocor fit:
#adjust spat autocor type: mat, exp, gau (matern, exponential, gaussian)
#adjust grain (e.g., /100=100m resolution)

#Run model with spatial autocor structure
geo.trapd.wk$X=round(geo.trapd.wk$mX/1000)
geo.trapd.wk$Y=round(geo.trapd.wk$mY/1000)
geo.trapd.wk$animalid<-factor(geo.trapd.wk$animalid)
geo.trapd.wk$pos <- numFactor(geo.trapd.wk$X, geo.trapd.wk$Y)
res.spat=glmmTMB(mNSD~Removal.Type*removal.period.akdecalc+
                   exp(pos + 0 | animalid)+ar1(as.factor(week) + 0 | animalid),
                 data=geo.trapd.wk,
                 family=Gamma(link="log"),
                 verbose=TRUE)

#Test for spatial autocorrelation with dharma package
res <- simulateResiduals(res.spat) #mle-mvn?
groupLocations = aggregate(geo.trapd.wk[, 7:8], list(as.factor(geo.trapd.wk$animalid)), mean)
res2 = recalculateResiduals(res, group = as.factor(geo.trapd.wk$animalid), rotation="estimated")
testSpatialAutocorrelation(res2,groupLocations$mX, groupLocations$mY)
#0.2369 with /1000 and exp

#Test trap model with sex interaction for spat autocorr
#Run model
res1.ac.sac=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + Removal.Type*removal.period.akdecalc*sex, data=geo.trapd.wk,family=Gamma(link=log))

#Test for spatial autocorrelation with dharma package
res <- simulateResiduals(res1.ac.sac)
groupLocations = aggregate(geo.trapd.wk[, 7:8], list(as.factor(geo.trapd.wk$animalid)), mean)
res2 = recalculateResiduals(res, group = as.factor(geo.trapd.wk$animalid), rotation="estimated")
testSpatialAutocorrelation(res2,groupLocations$mX, groupLocations$mY)
#is spatial autocorrelation for trap with sex interaction, p=0.0006

geo.trapd.wk$X=round(geo.trapd.wk$mX/1000)
geo.trapd.wk$Y=round(geo.trapd.wk$mY/1000)
geo.trapd.wk$animalid<-factor(geo.trapd.wk$animalid)
geo.trapd.wk$pos <- numFactor(geo.trapd.wk$Y, geo.trapd.wk$Y)
res.spat=glmmTMB(mNSD~Removal.Type*removal.period.akdecalc*sex+
                   exp(pos + 0 | animalid)+ar1(as.factor(week) + 0 | animalid),
                 data=geo.trapd.wk,
                 family=Gamma(link="log"),
                 verbose=TRUE)

geo.trapd.wk$key=paste0(geo.trapd.wk$animalid,geo.trapd.wk$week)
res <- simulateResiduals(res.spat) 
groupLocations = aggregate(geo.trapd.wk[, 7:8], list(as.factor(geo.trapd.wk$animalid)), mean)
res2 = recalculateResiduals(res, group = as.factor(geo.trapd.wk$animalid), rotation="estimated")
testSpatialAutocorrelation(res2,groupLocations$mX, groupLocations$mY)

summary(res.spat)
#/1000, Y-only, exp, p=0.16
max(geo.trapd.wk$mY)-min(geo.trapd.wk$mY)
max(geo.trapd.wk$mX)-min(geo.trapd.wk$mX)
#slightly more distance along Y

# * Tox spatial autocorr ----------------------------------------------------

#Run model
res1.ac.sac=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + Removal.Type*removal.period.akdecalc, data=geo.toxd.wk,family=Gamma(link=log))

#Test for spatial autocorrelation with dharma package
res <- simulateResiduals(res1.ac.sac)
groupLocations = aggregate(geo.toxd.wk[, 7:8], list(as.factor(geo.toxd.wk$animalid)), mean)
res2 = recalculateResiduals(res, group = as.factor(geo.toxd.wk$animalid), rotation="estimated")
testSpatialAutocorrelation(res2,groupLocations$mX, groupLocations$mY)
#is spatial autocorrelation, p=0.0097

geo.toxd.wk$X=round(geo.toxd.wk$mX/1000)
geo.toxd.wk$Y=round(geo.toxd.wk$mY/1000)
geo.toxd.wk$animalid<-factor(geo.toxd.wk$animalid)
geo.toxd.wk$pos <- numFactor(geo.toxd.wk$X, geo.toxd.wk$Y)
res.spat=glmmTMB(mNSD~Removal.Type*removal.period.akdecalc+
                   exp(pos + 0 | animalid)+ar1(as.factor(week) + 0 | animalid),
                 data=geo.toxd.wk,
                 family=Gamma(link="log"),
                 verbose=TRUE)

res <- simulateResiduals(res.spat)
groupLocations = aggregate(geo.toxd.wk[, 7:8], list(as.factor(geo.toxd.wk$animalid)), mean)
res2 = recalculateResiduals(res, group = as.factor(geo.toxd.wk$animalid), rotation="estimated")
testSpatialAutocorrelation(res2,groupLocations$mX, groupLocations$mY)
#removed spatial autocorr, p=0.07655

#Run model with sex interaction
res1.ac.sac=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + Removal.Type*removal.period.akdecalc*sex, data=geo.toxd.wk,family=Gamma(link=log))

#Test for spatial autocorrelation with dharma package
res <- simulateResiduals(res1.ac.sac)
groupLocations = aggregate(geo.toxd.wk[, 7:8], list(as.factor(geo.toxd.wk$animalid)), mean)
res2 = recalculateResiduals(res, group = as.factor(geo.toxd.wk$animalid), rotation="estimated")
testSpatialAutocorrelation(res2,groupLocations$mX, groupLocations$mY)
#is spatial autocorrelation, p=0.005


geo.toxd.wk$X=round(geo.toxd.wk$mX/1000)
geo.toxd.wk$Y=round(geo.toxd.wk$mY/1000)
geo.toxd.wk$animalid<-factor(geo.toxd.wk$animalid)
geo.toxd.wk$pos <- numFactor(geo.toxd.wk$X, geo.toxd.wk$Y)
res.spat=glmmTMB(mNSD~Removal.Type*removal.period.akdecalc*sex+
                   exp(pos + 0 | animalid)+ar1(as.factor(week) + 0 | animalid),
                 data=geo.toxd.wk,
                 family=Gamma(link="log"),
                 verbose=TRUE)

res <- simulateResiduals(res.spat)
groupLocations = aggregate(geo.toxd.wk[, 7:8], list(as.factor(geo.toxd.wk$animalid)), mean)
res2 = recalculateResiduals(res, group = as.factor(geo.toxd.wk$animalid), rotation="estimated")
testSpatialAutocorrelation(res2,groupLocations$mX, groupLocations$mY)
#removed spatial autocorr, p=0.9

#Spatial autocorrelation summary:
#aer: no spat autocorr
#trap: spat autocor for both models
#tox: spat autocorr for both models

# Run GLMMs --------------------------------------------------------------------

# * Aerial GLMMs --------------------------------------------------------------------
geo.wkdf=geo.aerd.wk
#geo.wkdf$Removal.Type<-forcats::fct_relevel(geo.wkdf$Removal.Type,c("ctrl","aer"))
#geo.wkdf$removal.period.akdecalc<-forcats::fct_relevel(geo.wkdf$removal.period.akdecalc,c("before","after"))
geo.wkdf$sex<-forcats::fct_relevel(geo.wkdf$sex,c("Female","Male"))
res.rp_aer=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + Removal.Type*removal.period.akdecalc, data=geo.wkdf,family=Gamma(link=log))
res.rps_aer=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + Removal.Type*removal.period.akdecalc*sex, data=geo.wkdf,family=Gamma(link=log))

# * Trap GLMMs --------------------------------------------------------------------

geo.wkdf=geo.trapd.wk
geo.wkdf$X=round(geo.wkdf$mX/1000)
geo.wkdf$Y=round(geo.wkdf$mY/1000)
geo.wkdf$animalid<-factor(geo.wkdf$animalid)
geo.wkdf$pos1 <- numFactor(geo.wkdf$X, geo.wkdf$Y)
geo.wkdf$pos2 <- numFactor(geo.wkdf$Y, geo.wkdf$Y)
#geo.wkdf$Removal.Type<-forcats::fct_relevel(geo.wkdf$Removal.Type,c("ctrl","trap"))
#geo.wkdf$removal.period.akdecalc<-forcats::fct_relevel(geo.wkdf$removal.period.akdecalc,c("before","during","after"))
geo.wkdf$sex<-forcats::fct_relevel(geo.wkdf$sex,c("Female","Male"))
res.rp_trap=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + exp(pos1 + 0 | animalid) + Removal.Type*removal.period.akdecalc, 
                    data=geo.wkdf,
                    family=Gamma(link=log),
                    verbose=TRUE)
res.rps_trap=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + exp(pos2 + 0 | animalid) + Removal.Type*removal.period.akdecalc*sex, 
                     data=geo.wkdf,
                     family=Gamma(link=log),
                     verbose=TRUE)

# * Toxicant GLMMs --------------------------------------------------------------------
geo.wkdf=geo.toxd.wk
geo.wkdf$X=round(geo.wkdf$mX/1000)
geo.wkdf$Y=round(geo.wkdf$mY/1000)
geo.wkdf$animalid<-factor(geo.wkdf$animalid)
geo.wkdf$pos <- numFactor(geo.wkdf$X, geo.wkdf$Y)
geo.wkdf$sex<-forcats::fct_relevel(geo.wkdf$sex,c("Female","Male"))
res.rp_tox=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + exp(pos + 0 | animalid) + Removal.Type*removal.period.akdecalc, 
                   data=geo.wkdf,
                   family=Gamma(link=log),
                   verbose=TRUE)
res.rps_tox=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + exp(pos + 0 | animalid) + Removal.Type*removal.period.akdecalc*sex, 
                    data=geo.wkdf,
                    family=Gamma(link=log),
                    verbose=TRUE)

# Save model objects -----------------------------------------------------------

saveRDS(res.rp_aer,file.path(objdir,"Models","res_NSD_rp_aer.rds",fsep=.Platform$file.sep))
saveRDS(res.rps_aer,file.path(objdir,"Models","res_NSD_rps_aer.rds",fsep=.Platform$file.sep))
saveRDS(res.rp_trap,file.path(objdir,"Models","res_NSD_rp_trap.rds",fsep=.Platform$file.sep))
saveRDS(res.rps_trap,file.path(objdir,"Models","res_NSD_rps_trap.rds",fsep=.Platform$file.sep))
saveRDS(res.rp_tox,file.path(objdir,"Models","res_NSD_rp_tox.rds",fsep=.Platform$file.sep))
saveRDS(res.rps_tox,file.path(objdir,"Models","res_NSD_rps_tox.rds",fsep=.Platform$file.sep))

