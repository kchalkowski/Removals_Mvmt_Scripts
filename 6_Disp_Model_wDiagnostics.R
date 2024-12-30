#set home dir of pipeline
home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline"

# Purpose ----------------------------------------------------------------------
#The purpose of this script is to format outputs to use for analysis in removals
#responses pipeline

#1. AKDE output formatting
#relevel, remove missing periods, etc.
#2. Geolocation data formatting
#General tidying
#set periods
#3. Calculate NSD
#4. Summarize weekly median NSD

# In/Out -----------------------------------------------------------------------

#for NSD calcs:
#geo_remtype.csv --> geo_aerd_wk.rds, geo_toxd_wk.rds, geo_trapd_wk.rds

#for AKDE calcs
#outdf_akde_aer_corrected.rds --> outdf_akde_aer_corrected_f.rds
#outdf_akde_tox_corrected.rds --> outdf_akde_tox_corrected_f.rds
#outdf_akde_trap_corrected.rds --> outdf_akde_trap_corrected_f.rds

# Setup ------------------------------------------------------------------------

#Set dirs
input=file.path(home,"1_Data","Input",fsep=.Platform$file.sep)
objdir=file.path(home,"1_Data","Objects",fsep=.Platform$file.sep)



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

#Read in data to use for NSD analysis:
geo.aerd.wk=readRDS(file.path(objdir,"geo_aerd_wk.rds",fsep=.Platform$file.sep))
geo.trapd.wk=readRDS(file.path(objdir,"geo_trapd_wk.rds",fsep=.Platform$file.sep))
geo.toxd.wk=readRDS(file.path(objdir,"geo_toxd_wk.rds",fsep=.Platform$file.sep))

#Relevel periods
geo.aerd.wk$removal.period.akdecalc<-fct_relevel(geo.aerd.wk$removal.period.akdecalc,c("before","during"))
geo.trapd.wk$removal.period.akdecalc<-fct_relevel(geo.trapd.wk$removal.period.akdecalc,c("before","during"))
geo.toxd.wk$removal.period.akdecalc<-fct_relevel(geo.toxd.wk$removal.period.akdecalc,c("before","during"))

#Relevel removal types
geo.aerd.wk$Removal.Type<-fct_relevel(geo.aerd.wk$Removal.Type,"ctrl")
geo.trapd.wk$Removal.Type<-fct_relevel(geo.trapd.wk$Removal.Type,"ctrl")
geo.toxd.wk$Removal.Type<-fct_relevel(geo.toxd.wk$Removal.Type,"ctrl")

# Temporal (weekly) autocorrelation checking -----------------------------------
{
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

outdir="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Exploratory_Outputs/PACF_animalid_plots_akde_displacement/"
geo.remd.wk=geo.aerd.wk
removal.str="aer"
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

Do.PACF.Exploratory(geo.aerd.wk,outdir,"aer")
Do.PACF.Exploratory(geo.toxd.wk,outdir,"tox")
Do.PACF.Exploratory(geo.trapd.wk,outdir,"trap")

#Looked at plots, does appear to be some first-order autocorrelation in displacement between weeks
#Try adding autocorrelation structure to model, see if residuals become less autocorrelated across weeks
#library(fitdistrplus)
#Check histos to determine appropriate model structures
outdir="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Exploratory_Outputs/weekly_disp_hist/"
Do.Displacement.Histos<-function(wak2,outdir,removal.str){
  IDs=unique(wak2$animalid)
for(i in 1:length(unique(wak2$animalid))){
ggplot(data=wak2[wak2$animalid==unique(wak2$animalid)[i],],aes(x=mNSD))+
    geom_histogram()
  if(!exists(paste0(outdir,removal.str,"/"))){dir.create(paste0(outdir,removal.str,"/"))}
  filename=paste0(paste0(outdir,removal.str,"/"),IDs[i],".png")
  ggsave(filename,bg="white")
}
}

Do.Displacement.Histos(geo.aerd.wk,outdir,"aer")
Do.Displacement.Histos(geo.toxd.wk,outdir,"tox")
Do.Displacement.Histos(geo.trapd.wk,outdir,"trap")

#Make function to run different autocorrelation structures, check residual plots

removal.str="aer"
outdir="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Exploratory_Outputs/"

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

Check.Week.Autocorr(geo.aerd.wk,outdir,"aer")
Check.Week.Autocorr(geo.trapd.wk,outdir,"trap")
Check.Week.Autocorr(geo.toxd.wk,outdir,"tox")

#Including autocor structure does seem to reduce autocor of residuals considerably.
res.ac.aer=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid), data=geo.aerd.wk,family=Gamma(link=log))
res.ac.trap=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid), data=geo.trapd.wk,family=Gamma(link=log))
res.ac.tox=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid), data=geo.toxd.wk,family=Gamma(link=log))

}

# Check distribution fits ------------------------------------------------------
{
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
}

# Check for spatial autocorrelation --------------------------------------------
{
  #Within each week, across pigs, do we see spatial autocorrelation of displacement values?
  res1.ac.sac=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + Removal.Type*removal.period.akdecalc, data=geo.aerd.wk,family=Gamma(link=log))
  #Test for spatial autocorrelation with dharma package
  res <- simulateResiduals(res1.ac.sac)
  groupLocations = aggregate(geo.aerd.wk[, 7:8], list(as.factor(geo.aerd.wk$animalid)), mean)
  res2 = recalculateResiduals(res, group = as.factor(geo.aerd.wk$animalid), rotation="estimated")
  testSpatialAutocorrelation(res2,groupLocations$mX, groupLocations$mY)
  
  #No spatial autocorr
  res1.ac.sac=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + Removal.Type*removal.period.akdecalc, data=geo.trapd.wk,family=Gamma(link=log))
  #Test for spatial autocorrelation with dharma package
  res <- simulateResiduals(res1.ac.sac)
  groupLocations = aggregate(geo.trapd.wk[, 7:8], list(as.factor(geo.trapd.wk$animalid)), mean)
  res2 = recalculateResiduals(res, group = as.factor(geo.trapd.wk$animalid), rotation="estimated")
  testSpatialAutocorrelation(res2,groupLocations$mX, groupLocations$mY)
  #no spatial autocorrelation
  
  res1.ac.sac=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + Removal.Type*removal.period.akdecalc, data=geo.toxd.wk,family=Gamma(link=log))
  #Test for spatial autocorrelation with dharma package
  res <- simulateResiduals(res1.ac.sac)
  groupLocations = aggregate(geo.toxd.wk[, 7:8], list(as.factor(geo.toxd.wk$animalid)), mean)
  res2 = recalculateResiduals(res, group = as.factor(geo.toxd.wk$animalid), rotation="estimated")
  testSpatialAutocorrelation(res2,groupLocations$mX, groupLocations$mY)
  #spatial autocorrelation

  #####Get which IDs died in aer/tox
  grps=geo.toxd.wk %>% group_by(animalid) %>% dplyr::summarise(n=unique(removal.period.akdecalc))
  grps$vals=1
  grps=tidyr::pivot_wider(grps,names_from=n,values_from=vals)
  died.tox=grps[which(is.na(grps$after)),]$animalid
  #Separate data by who died during tox treatment
  geo.toxd.wk$died_tox=NA
  geo.toxd.wk[geo.toxd.wk$animalid%in%died.tox,]$died_tox=1
  geo.toxd.wk[is.na(geo.toxd.wk$died_tox),]$died_tox<-0
  geo.toxd.wk$died_tox<-as.factor(geo.toxd.wk$died_tox)
  
  res1.ac.sac=glmmTMB(mNSD ~ sex + ar1(as.factor(week) + 0 | animalid) + Removal.Type*removal.period.akdecalc, data=geo.toxd.wk,family=Gamma(link=log))
  #Test for spatial autocorrelation with dharma package
  res <- simulateResiduals(res1.ac.sac)
  groupLocations = aggregate(geo.toxd.wk[, c(which(colnames(geo.toxd.wk)=="mX"),which(colnames(geo.toxd.wk)=="mY"))], list(as.factor(geo.toxd.wk$animalid)), mean)
  res2 = recalculateResiduals(res, group = as.factor(geo.toxd.wk$animalid), rotation="estimated")
  testSpatialAutocorrelation(res2,groupLocations$mX, groupLocations$mY)
  #spatial autocorrelation
  
  #spatial autocor goes away when either include died_tox or sex as a var
} 

# Run GLMMs --------------------------------------------------------------------

# Aerial
  geo.wkdf=geo.aerd.wk
  #geo.wkdf$Removal.Type<-forcats::fct_relevel(geo.wkdf$Removal.Type,c("ctrl","aer"))
  #geo.wkdf$removal.period.akdecalc<-forcats::fct_relevel(geo.wkdf$removal.period.akdecalc,c("before","after"))
  geo.wkdf$sex<-forcats::fct_relevel(geo.wkdf$sex,c("Female","Male"))
  res.rp_aer=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + Removal.Type*removal.period.akdecalc, data=geo.wkdf,family=Gamma(link=log))
  res.rps_aer=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + Removal.Type*removal.period.akdecalc*sex, data=geo.wkdf,family=Gamma(link=log))
  
# Trap
  geo.wkdf=geo.trapd.wk
  #geo.wkdf$Removal.Type<-forcats::fct_relevel(geo.wkdf$Removal.Type,c("ctrl","trap"))
  #geo.wkdf$removal.period.akdecalc<-forcats::fct_relevel(geo.wkdf$removal.period.akdecalc,c("before","during","after"))
  geo.wkdf$sex<-forcats::fct_relevel(geo.wkdf$sex,c("Female","Male"))
  res.rp_trap=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + Removal.Type*removal.period.akdecalc, data=geo.wkdf,family=Gamma(link=log))
  res.rps_trap=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + Removal.Type*removal.period.akdecalc*sex, data=geo.wkdf,family=Gamma(link=log))

# Toxicant 
  geo.wkdf=geo.toxd.wk
  #geo.wkdf$Removal.Type<-forcats::fct_relevel(geo.wkdf$Removal.Type,c("tox","ctrl"))
  #geo.wkdf$removal.period.akdecalc<-forcats::fct_relevel(geo.wkdf$removal.period.akdecalc,c("before","during","after"))
  geo.wkdf$sex<-forcats::fct_relevel(geo.wkdf$sex,c("Female","Male"))
  res.rp_tox=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + Removal.Type*removal.period.akdecalc, data=geo.wkdf,family=Gamma(link=log))
  res.rps_tox=glmmTMB(mNSD ~ ar1(as.factor(week) + 0 | animalid) + Removal.Type*removal.period.akdecalc*sex, data=geo.wkdf,family=Gamma(link=log))
  
  #Save model objects
  saveRDS(res.rp_aer,file.path(objdir,"Models","res_NSD_rp_aer.rds",fsep=.Platform$file.sep))
  saveRDS(res.rps_aer,file.path(objdir,"Models","res_NSD_rps_aer.rds",fsep=.Platform$file.sep))
  saveRDS(res.rp_trap,file.path(objdir,"Models","res_NSD_rp_trap.rds",fsep=.Platform$file.sep))
  saveRDS(res.rps_trap,file.path(objdir,"Models","res_NSD_rps_trap.rds",fsep=.Platform$file.sep))
  saveRDS(res.rp_tox,file.path(objdir,"Models","res_NSD_rp_tox.rds",fsep=.Platform$file.sep))
  saveRDS(res.rps_tox,file.path(objdir,"Models","res_NSD_rps_tox.rds",fsep=.Platform$file.sep))
  
  
  
  
  
  
# Make results figures (need tidy/fix below) --------------------------------------------------------------------

  Tidy.Parms=function(topmodel,removal.str,model.num){
    topmodel=broom.mixed::tidy(topmodel) %>% as.data.frame()
    topmodel=topmodel[topmodel$effect=="fixed",]
    topmodel$removal=removal.str
    topmodel$model=model.num
    return(topmodel)
  }
  
  aerm1=Tidy.Parms(aer.NSD.topmodel1,"aer",1)
  aerm2=Tidy.Parms(aer.NSD.topmodel2,"aer",2)
  aerm3=Tidy.Parms(aer.NSD.topmodel3,"aer",3)
  trapm1=Tidy.Parms(trap.NSD.topmodel1,"trap",1)
  toxm1=Tidy.Parms(tox.NSD.topmodel1,"tox",1)
  toxm2=Tidy.Parms(tox.NSD.topmodel2,"tox",2)
  
  parm.tbls=rbind(aerm1,aerm2,aerm3,trapm1,toxm1,toxm2)
  
  results.dir="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Output/Disp_GLM_Results/"
  write.csv(parm.tbls,paste0(results.dir,"topmodels_parmtbl.csv"))
  #write.csv(toxmdt0,paste0(results.dir,"toxdt0parms.csv"))
  #write.csv(aerm2,paste0(results.dir,"aermodel2.csv"))
  
  
  AIC.out=rbind(AIC.aer.NSD,
        AIC.trap.NSD,
        AIC.tox.NSD)
  
  write.csv(AIC.out,paste0(results.dir,"AIC_NSD_all.csv"))
  
  
  ########################## Figures for model effects
  
  #aer... aer.NSD.topmodel2
  aer.tmp <- ggeffects::predict_response(aer.NSD.topmodel2, terms=c("Removal.Type","removal.period.akdecalc","sex"))
  plot(aer.tmp)
  aer.tmp.df=as.data.frame(aer.tmp)
  
  aer.tmp.df$x<-forcats::fct_relevel(aer.tmp.df$x,c("ctrl","aer"))
  aer.tmp.df$group<-forcats::fct_relevel(aer.tmp.df$group,c("before","after"))
  aer.tmp.df$facet<-forcats::fct_relevel(aer.tmp.df$facet,c("Female","Male"))
  
  aer.ND.plot=ggplot(aer.tmp.df,aes(x=group,yend=sqrt(conf.high),group=x,color=x))+
    geom_point(aes(y=sqrt(predicted)),size=6,position = position_dodge(width = 0.2))+
    geom_line(aes(y=sqrt(predicted),group=x),linewidth=2,position = position_dodge(width = 0.2))+
    geom_segment(aes(y=sqrt(conf.low)),linewidth=3,position = position_dodge(width = 0.2))+
    theme_ipsum(axis_title_size=15,axis_text_size=18,strip_text_size =18,base_size=15)+
    theme(legend.text=element_text(size=15))+
    scale_color_manual(name="Removal\ntreatment",
                       values=c("#FFD065","#822AFF"))+
    ylab("Med. wk. disp. (m)")+
    xlab("Removal period")+
    facet_wrap(~facet)
  #ggsave(paste0(results.dir,"aer_NSD_fig.png"),bg="white",height=4.65,width=8.43,units="in")
  
  ########################## TRAP
  #trap, rps2
  #trap.NSD.topmodel1
  trap.tmp <- ggeffects::predict_response(trap.NSD.topmodel1, terms=c("Removal.Type","removal.period.akdecalc","sex"))
  plot(trap.tmp)
  trap.tmp.df=as.data.frame(trap.tmp)
  
  trap.tmp.df$x<-forcats::fct_relevel(trap.tmp.df$x,c("ctrl","trap"))
  trap.tmp.df$group<-forcats::fct_relevel(trap.tmp.df$group,c("before","during","after"))
  trap.tmp.df$facet<-forcats::fct_relevel(trap.tmp.df$facet,c("Female","Male"))
  
  trap.ND.plot=ggplot(trap.tmp.df,aes(x=group,yend=sqrt(conf.high),group=x,color=x))+
    geom_point(aes(y=sqrt(predicted)),size=6,position = position_dodge(width = 0.2))+
    geom_line(aes(y=sqrt(predicted),group=x),linewidth=2,position = position_dodge(width = 0.2))+
    geom_segment(aes(y=sqrt(conf.low)),linewidth=3,position = position_dodge(width = 0.2))+
    theme_ipsum(axis_title_size=15,axis_text_size=18,strip_text_size =18,base_size=15)+
    theme(legend.text=element_text(size=15))+
    scale_color_manual(name="Removal\ntreatment",
                       values=c("#FFD065","#F74DD4"))+
    ylab("Med. wk. disp. (m)")+
    xlab("Removal period")+
    facet_wrap(~facet)
  #ggsave(paste0(results.dir,"trap_NSD_fig.png"))
  #ggsave(paste0(results.dir,"trap_NSD_fig.png"),bg="white",height=4.65,width=8.43,units="in")
  
  
  ########################## TOX
  #trap, rps2
  #trap.NSD.topmodel1
  tox.tmp <- ggeffects::predict_response(res.rp.dto2, terms=c("Removal.Type","removal.period.akdecalc","died_tox"))
  tox.tmp.df=as.data.frame(tox.tmp)
  
  tox.tmp.df$x<-forcats::fct_relevel(tox.tmp$x,c("ctrl","tox"))
  tox.tmp.df$group<-forcats::fct_relevel(tox.tmp$group,c("before","during","after"))
  tox.tmp.df=tox.tmp.df[!(tox.tmp.df$group=="after"&tox.tmp.df$facet==1&tox.tmp.df$x=="tox"),]
  tox.tmp.df=tox.tmp.df[!(tox.tmp.df$facet==1&tox.tmp.df$x=="ctrl"),]
  tox.tmp.df$x=as.character(tox.tmp.df$x)
  tox.tmp.df[tox.tmp.df$facet==1&tox.tmp.df$x=="tox",]$x="tox_fate.d"
  tox.tmp.df[tox.tmp.df$facet==0&tox.tmp.df$x=="tox",]$x="tox_fate.s"
  
  tox.tmp.df$x=as.factor(tox.tmp.df$x)
  
  tox.ND.plot=ggplot(tox.tmp.df,aes(x=group,yend=sqrt(conf.high),group=x,color=x))+
    geom_point(aes(y=sqrt(predicted)),size=6,position = position_dodge(width = 0.2))+
    geom_line(aes(y=sqrt(predicted),group=x),linewidth=2,position = position_dodge(width = 0.2))+
    geom_segment(aes(y=sqrt(conf.low)),linewidth=3,position = position_dodge(width = 0.2))+
    theme_ipsum(axis_title_size=15,axis_text_size=18,strip_text_size =18,base_size=15)+
    theme(legend.text=element_text(size=15))+
    scale_color_manual(name="Removal\ntreatment",
                       values=c("#FFD065","#ff5f00","#FF985A"))+
    ylab("Med. wk. disp. (m)")+
    xlab("Removal period")
  tox.ND.plot
  #ggsave(paste0(results.dir,"tox_NSD_fig.png"))
  
  

  #tox.tmp.dt0 <- ggeffects::predict_response(tox.dt0.NSD.topmodel1, terms=c("Removal.Type","removal.period.akdecalc"))
  #tox.tmp.df.dt0=as.data.frame(tox.tmp.dt0)
  
  #tox.tmp.df=tox.tmp.df.dt0
  #tox.tmp.df$x<-forcats::fct_relevel(tox.tmp$x,c("ctrl","tox"))
  #tox.tmp.df$group<-forcats::fct_relevel(tox.tmp$group,c("before","during","after"))
  
  #toxdt0.ND.plot=ggplot(tox.tmp.df,aes(x=group,yend=sqrt(conf.high),group=x,color=x))+
  #  geom_point(aes(y=sqrt(predicted)),size=6,position = position_dodge(width = 0.2))+
  #  geom_line(aes(y=sqrt(predicted),group=x),linewidth=2,position = position_dodge(width = 0.2))+
  #  geom_segment(aes(y=sqrt(conf.low)),linewidth=3,position = position_dodge(width = 0.2))+
  #  theme_ipsum(axis_title_size=15,axis_text_size=18,strip_text_size =18,base_size=15)+
  #  theme(legend.text=element_text(size=15))+
  #  scale_color_manual(name="Removal\ntreatment",
  #                     values=c("#FFD065","#FF7D2E"))+
  #  ylab("Med. wk. disp. (m)")+
  #  xlab("Removal period")
  ##ggsave(paste0(results.dir,"toxdt0_NSD_fig.png"))
  
  
  #aer.ND.plot
  #trap.ND.plot
  #tox.ND.plot
  #toxdt0.ND.plot
  
  cowplot::plot_grid(aer.ND.plot,trap.ND.plot,tox.ND.plot,ncol=1)
  ggsave(paste0(results.dir,"allplot.png"),bg="white",height=10,width=10,units="in")
  
  #ctrl: #FFD065
  #trap: #FF57A9
  #aer: #822AFF
  #tox: #FF985A

  
  #Get values for reporting results in manuscript
  Get.Pred.Diffs<-function(trap.tmp.df,removal.str){
  trap.tmp.df.t=trap.tmp.df
  trap.tmp.df.t[,2:5]=sqrt(trap.tmp.df.t[,2:5])
  
  removal.str="trap"
  #Control females
  trap.tmp.df.t.FC=trap.tmp.df.t[trap.tmp.df.t$facet=="Female"&
                                   trap.tmp.df.t$x=="ctrl",]
  trap.tmp.df.t.FC$pred.diff=trap.tmp.df.t.FC$predicted-trap.tmp.df.t.FC[trap.tmp.df.t.FC$group=="before",]$predicted
  trap.tmp.df.t.FC$conf.low.diff=trap.tmp.df.t.FC$conf.low-trap.tmp.df.t.FC[trap.tmp.df.t.FC$group=="before",]$conf.low
  trap.tmp.df.t.FC$conf.high.diff=trap.tmp.df.t.FC$conf.high-trap.tmp.df.t.FC[trap.tmp.df.t.FC$group=="before",]$conf.high
  
  #Control males
  trap.tmp.df.t.MC=trap.tmp.df.t[trap.tmp.df.t$facet=="Male"&
                                   trap.tmp.df.t$x=="ctrl",]
  trap.tmp.df.t.MC$pred.diff=trap.tmp.df.t.MC$predicted-trap.tmp.df.t.MC[trap.tmp.df.t.MC$group=="before",]$predicted
  trap.tmp.df.t.MC$conf.low.diff=trap.tmp.df.t.MC$conf.low-trap.tmp.df.t.MC[trap.tmp.df.t.MC$group=="before",]$conf.low
  trap.tmp.df.t.MC$conf.high.diff=trap.tmp.df.t.MC$conf.high-trap.tmp.df.t.MC[trap.tmp.df.t.MC$group=="before",]$conf.high
  
  #Treatment females
  trap.tmp.df.t.FT=trap.tmp.df.t[trap.tmp.df.t$facet=="Female"&
                                   trap.tmp.df.t$x==removal.str,]
  trap.tmp.df.t.FT$pred.diff=trap.tmp.df.t.FT$predicted-trap.tmp.df.t.FT[trap.tmp.df.t.FT$group=="before",]$predicted
  trap.tmp.df.t.FT$conf.low.diff=trap.tmp.df.t.FT$conf.low-trap.tmp.df.t.FT[trap.tmp.df.t.FT$group=="before",]$conf.low
  trap.tmp.df.t.FT$conf.high.diff=trap.tmp.df.t.FT$conf.high-trap.tmp.df.t.FT[trap.tmp.df.t.FT$group=="before",]$conf.high
  
  
  #Treatment males
  trap.tmp.df.t.MT=trap.tmp.df.t[trap.tmp.df.t$facet=="Male"&
                                   trap.tmp.df.t$x==removal.str,]
  trap.tmp.df.t.MT$pred.diff=trap.tmp.df.t.MT$predicted-trap.tmp.df.t.MT[trap.tmp.df.t.MT$group=="before",]$predicted
  trap.tmp.df.t.MT$conf.low.diff=trap.tmp.df.t.MT$conf.low-trap.tmp.df.t.MT[trap.tmp.df.t.MT$group=="before",]$conf.low
  trap.tmp.df.t.MT$conf.high.diff=trap.tmp.df.t.MT$conf.high-trap.tmp.df.t.MT[trap.tmp.df.t.MT$group=="before",]$conf.high
  
  
out.diffs=rbind(trap.tmp.df.t.FC,trap.tmp.df.t.MC,trap.tmp.df.t.FT,trap.tmp.df.t.MT)
  return(out.diffs)
  }

  trap.diffs=Get.Pred.Diffs(trap.tmp.df,"trap")
  trap.diffs[trap.diffs$x=="trap"&trap.diffs$facet=="Male",]
    #Control females displacement:
      #increased 393 m from before - during (95% C.L. 289-533m)
      #increased 558 m from before - after (95% C.L. 397-784m)
    #Control males displacement
      #increased 182 m from before - during (95% C.L. 130-254m)
      #increased 682 m from before -after (95% C.L. 455-1021m)
  
    #Trap females displacement:
      #increased 31m from before - during (95% C.L. 30-32m)
      #increased 19m from before - after (95% C.L. 30-34m)
    #Trap males displacement:
      #increased 1401m from before - during (95% C.L. 995-1973m)
      #increased 1459m from before - after (95% C.L. 1017-2095m)
  
  #summarize by week
  geo.aerd.wk=geo.aerd2 %>% group_by(animalid, Removal.Type, removal.period.akdecalc,sex, week) %>% dplyr::summarise(mNSD=median(NSD),NSD.25=quantile(NSD,0.25),NSD.75=quantile(NSD,0.75),mX=mean(X),mY=mean(Y)) %>% as.data.frame()
  geo.trapd.wk=geo.trapd2 %>% group_by(animalid, Removal.Type, removal.period.akdecalc,sex, week) %>% dplyr::summarise(mNSD=median(NSD),NSD.25=quantile(NSD,0.25),NSD.75=quantile(NSD,0.75),mX=mean(X),mY=mean(Y)) %>% as.data.frame()
  geo.toxd.wk=geo.toxd2 %>% group_by(animalid, Removal.Type, removal.period.akdecalc,sex, week) %>% dplyr::summarise(mNSD=median(NSD),NSD.25=quantile(NSD,0.25),NSD.75=quantile(NSD,0.75),mX=mean(X),mY=mean(Y)) %>% as.data.frame()

  aermeds=geo.aerd2 %>% group_by(Removal.Type,removal.period.akdecalc,sex,week) %>% dplyr::summarise(mNSD=median(NSD),NSD.25=quantile(NSD,0.25),NSD.75=quantile(NSD,0.75)) %>% as.data.frame()
  trapmeds=geo.trapd2 %>% group_by(Removal.Type,removal.period.akdecalc,sex,week) %>% dplyr::summarise(mNSD=median(NSD),NSD.25=quantile(NSD,0.25),NSD.75=quantile(NSD,0.75)) %>% as.data.frame()
  toxmeds=geo.toxd2 %>% group_by(Removal.Type,removal.period.akdecalc,sex,week) %>% dplyr::summarise(mNSD=median(NSD),NSD.25=quantile(NSD,0.25),NSD.75=quantile(NSD,0.75)) %>% as.data.frame()
  
  aer.aft.w=min(geo.aerd.wk[geo.aerd.wk$removal.period.akdecalc=="after",]$week)
  ggplot(aermeds,aes(group=sex, color=sex))+
    geom_vline(xintercept=6,
              color='#b7ffb5',linewidth=2)+
    geom_line(aes(y=mNSD,x=week))+
    facet_wrap(~Removal.Type)+
    theme_ipsum()+
    geom_ribbon(aes(x=week, ymin=NSD.25, ymax=NSD.75,fill=sex), linetype="blank",alpha=0.1)
    
  trap.dur.w=min(geo.trapd.wk[geo.trapd.wk$removal.period.akdecalc=="during",]$week)
  trap.aft.w=min(geo.trapd.wk[geo.trapd.wk$removal.period.akdecalc=="after",]$week)
  ggplot(trapmeds,aes(group=sex))+
    geom_vline(xintercept=trap.dur.w,
               color='#8fff8c',linewidth=2)+
    geom_vline(xintercept=trap.aft.w,
               color='#8fff8c',linewidth=2)+
    geom_line(aes(y=mNSD,x=week))+
    facet_wrap(~Removal.Type)+
    theme_ipsum()+
    geom_ribbon(aes(x=week, ymin=NSD.25, ymax=NSD.75,fill=sex), alpha=0.3)
  
  
  tox.dur.w=min(geo.toxd.wk[geo.toxd.wk$removal.period.akdecalc=="during",]$week)
  tox.aft.w=min(geo.toxd.wk[geo.toxd.wk$removal.period.akdecalc=="after",]$week)
  ggplot(toxmeds,aes(group=sex))+
    geom_line(aes(y=mNSD,x=week,color=sex))+
    geom_vline(xintercept=tox.dur.w,
               color='#8fff8c',linewidth=2)+
    geom_vline(xintercept=tox.aft.w,
               color='#8fff8c',linewidth=2)+
    facet_wrap(~Removal.Type)+
    theme_ipsum()+
    geom_ribbon(aes(x=week, ymin=NSD.25, ymax=NSD.75,fill=sex), linetype=2, alpha=0.1)
  
  
  #86101_J3_J3
  toxmeds=geo.toxd2[geo.toxd2$animalid!="86101_J3_J3",] %>% group_by(Removal.Type,removal.period.akdecalc,sex,week) %>% dplyr::summarise(mNSD=mean(NSD),NSD.25=quantile(NSD,0.25),NSD.75=quantile(NSD,0.75)) %>% as.data.frame()
  toxmeds=geo.toxd2 %>% group_by(Removal.Type,removal.period.akdecalc,sex,week) %>% dplyr::summarise(mNSD=mean(NSD),NSD.25=quantile(NSD,0.25),NSD.75=quantile(NSD,0.75)) %>% as.data.frame()
  
  tox.dur.w=min(geo.toxd.wk[geo.toxd.wk$removal.period.akdecalc=="during",]$week)
  tox.aft.w=min(geo.toxd.wk[geo.toxd.wk$removal.period.akdecalc=="after",]$week)
  ggplot(toxmeds,aes(group=sex))+
    geom_line(aes(y=mNSD,x=week,color=sex))+
    geom_vline(xintercept=tox.dur.w,
               color='#8fff8c',linewidth=2)+
    geom_vline(xintercept=tox.aft.w,
               color='#8fff8c',linewidth=2)+
    facet_wrap(~Removal.Type)+
    theme_ipsum()+
    geom_ribbon(aes(x=week, ymin=NSD.25, ymax=NSD.75,fill=sex), linetype=2, alpha=0.1)
  
  #86101_J3_J3
  
  library(sf)
  library(mapview)
  pig=geo.toxd[geo.toxd$animalid=="86101_J3_J3",]
  pigsf=st_as_sf(pig,coords=c(7,8),crs=st_crs(4326))  
  mapview(pigsf,zcol="NSD")  
  
  
  toxs=geo.toxd %>% group_by(animalid,Removal.Type,removal.period.akdecalc) %>% dplyr::summarise(maxNSD=max(NSD),mX=median(X),mY=median(Y))  
  traps=geo.trapd %>% group_by(animalid,Removal.Type,removal.period.akdecalc) %>% dplyr::summarise(maxNSD=max(NSD),mX=median(X),mY=median(Y))  
  aers=geo.aerd %>% group_by(animalid,Removal.Type,removal.period.akdecalc) %>% dplyr::summarise(maxNSD=max(NSD),mX=median(X),mY=median(Y))  
  
  toxsf=st_as_sf(toxs,coords=c(5,6),crs=st_crs(32614))
  trapsf=st_as_sf(traps,coords=c(5,6),crs=st_crs(32614))
  aersf=st_as_sf(aers,coords=c(5,6),crs=st_crs(32614))
  
  mapview(toxsf[toxsf$removal.period.akdecalc=="before",],zcol="maxNSD")
    mapview(trapsf,zcol="maxNSD")
    mapview(aersf,zcol="maxNSD")
    
    ggplot(toxs,aes(x=maxNSD,fill=Removal.Type))+
      geom_density()
    
    x=toxs[toxs$Removal.Type=="ctrl"&toxs$removal.period.akdecalc=="before",]$maxNSD
    y=toxs[toxs$Removal.Type=="tox"&toxs$removal.period.akdecalc=="before",]$maxNSD
    
    tox.ks.res=ks.test(x,y)
    
    
    
    
    ggplot2()
    
    
    ?ks.test
    
    hist(toxsf$maxNSD)
  #i 7 86101 j3 j3
  #i 9 86103 u2 u2  
    i=27
  IDs=unique(geo.toxd2[geo.toxd2$Removal.Type=="tox",]$animalid)
  pigex=geo.toxd2[geo.toxd2$animalid==IDs[i],]
ggplot(pigex,
       aes(x=as.numeric(datetime),
           y=NSD,color=removal.period.akdecalc))+
      geom_line()+
      labs(title=paste(IDs[i],pigex$Removal.Type[1],pigex$sex[1]))
    
  
  