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

#Set dirs
input=file.path(home,"1_Data","Input",fsep=.Platform$file.sep)
objdir=file.path(home,"1_Data","Objects",fsep=.Platform$file.sep)
if(!dir.exists(file.path(home,"3_Output","Area_GLM_Results",fsep=.Platform$file.sep))){
  dir.create(file.path(home,"3_Output","Area_GLM_Results",fsep=.Platform$file.sep))
}
outdir=file.path(home,"3_Output","Area_GLM_Results",fsep=.Platform$file.sep)

#Read input objects
akaer=readRDS(file.path(objdir,"outdf_akde_aer_corrected_f.rds",fsep=.Platform$file.sep))
aktrap=readRDS(file.path(objdir,"outdf_akde_trap_corrected_f.rds",fsep=.Platform$file.sep))
aktox=readRDS(file.path(objdir,"outdf_akde_aktox_corrected_f.rds",fsep=.Platform$file.sep))

#source needed functions
funcdir<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/2_Scripts/Functions"
func.list=list.files(funcdir,full.names=TRUE)
for(f in 1:length(func.list)){
  source(func.list[f])
}

# Check aerial data spatial autocorrelation -----------------------------------
{
res=glmmTMB(area.est~(1|animalid)+Removal.Type*period,data=akaer,family=Gamma(link="log"))
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
res=glmmTMB(area.est~(1|animalid)+Removal.Type*period,data=aktrap,family=Gamma(link="log"))
ress <- simulateResiduals(res)
groupLocations = aggregate(aktrap[, 6:7], list(aktrap$animalid), mean)
ress2 = recalculateResiduals(ress, group = aktrap$animalid, rotation="estimated")
#testSpatialAutocorrelation(ress2,groupLocations$ctr.x, groupLocations$ctr.y)
#is spatial autocorrelation

#Correct for spatial autocorrelation, try gaussian autocor structure
aktrap$pos <- numFactor(aktrap$ctr.x, aktrap$ctr.y)
res.sa=glmmTMB(area.est~Removal.Type*period+exp(pos + 0 | animalid),data=aktrap,family=Gamma(link="log"))
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
res.sa=glmmTMB(area.est~Removal.Type*period+exp(pos + 0 | animalid),data=aktrap,family=Gamma(link="log"))
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
res=glmmTMB(area.est~(1|animalid)+Removal.Type*period*died.tox,data=aktox,family=Gamma(link="log"))
ress <- simulateResiduals(res)
groupLocations = aggregate(aktox[, 6:7], list(aktox$animalid), mean)
ress2 = recalculateResiduals(ress, group = aktox$animalid, rotation="estimated")
#testSpatialAutocorrelation(ress2,groupLocations$ctr.x, groupLocations$ctr.y)
#No spatial autocorrelation when died.tox is included-- opt to not include spatial autocorr in model, is assoc. with died. tox
}

# Run GLMMs ------------------------------------------------------------

# Aerial
res.rp_aer=glmmTMB(area.est ~ Removal.Type*period + (1|animalid), data=akaer,family=Gamma(link=log))
res.rps_aer=glmmTMB(area.est ~ Removal.Type*period*sex + (1|animalid), data=akaer,family=Gamma(link=log))

# Trap
res.rp_trap=glmmTMB(area.est ~ Removal.Type*period + exp(pos + 0 | animalid), data=aktrap,family=Gamma(link=log))
#Note, model with removal.type*period+sex would not converge, so excluding from model list

# Tox
res.rp_tox=glmmTMB(area.est ~ Removal.Type*period + (1|animalid), data=aktox,family=Gamma(link=log))
res.rps_tox=glmmTMB(area.est ~ Removal.Type*period*sex + (1|animalid), data=aktox,family=Gamma(link=log))

#Save model objects
saveRDS(res.rp_aer,file.path(objdir,"Models","res_area_rp_aer.rds",fsep=.Platform$file.sep))
saveRDS(res.rps_aer,file.path(objdir,"Models","res_area_rps_aer.rds",fsep=.Platform$file.sep))
saveRDS(res.rp_trap,file.path(objdir,"Models","res_area_rp_trap.rds",fsep=.Platform$file.sep))
saveRDS(res.rp_tox,file.path(objdir,"Models","res_area_rp_tox.rds",fsep=.Platform$file.sep))
saveRDS(res.rps_tox,file.path(objdir,"Models","res_area_rps_tox.rds",fsep=.Platform$file.sep))





















# Make results figures ------------------------------------------------------------
# Below needs tidying with new structure

tox.tmp <- ggeffects::predict_response(tox.topmodel, terms=c("Removal.Type","period","sex"))
tox.tmp.df=as.data.frame(tox.tmp)

ggplot(tox.tmp.df,aes(x=x,yend=conf.high,group=group,color=group))+
  geom_point(aes(y=predicted),size=6,position = position_dodge(width = 0.2))+
  geom_segment(aes(y=sqrt(conf.low)),linewidth=3,position = position_dodge(width = 0.2))+
  #geom_line(aes(y=sqrt(predicted),group=x),linewidth=2,position = position_dodge(width = 0.2))+
  theme_ipsum(axis_title_size=15,axis_text_size=18,strip_text_size =18,base_size=15)+
  theme(legend.text=element_text(size=15))+
  scale_color_manual(name="Removal\ntreatment",
                     values=c("#FFD065","#FF7D2E"))+
  ylab("Home range area (km2)")+
  xlab("Removal period")+
  facet_wrap(~facet)
ggsave(paste0(outdir,"tox_prediction_figure.png"),bg="white",height=4,width=6,units="in")

#season
tox.tmpdt0 <- ggeffects::predict_response(tox.topmodeldt0, terms=c("Removal.Type","period","sex"))
tox.tmp.df.dt0=as.data.frame(tox.tmpdt0)

ggplot(tox.tmp.df.dt0,aes(x=x,yend=conf.high,group=group,color=group))+
  geom_point(aes(y=predicted),size=6,position = position_dodge(width = 0.2))+
  geom_segment(aes(y=sqrt(conf.low)),linewidth=3,position = position_dodge(width = 0.2))+
  #geom_line(aes(y=sqrt(predicted),group=x),linewidth=2,position = position_dodge(width = 0.2))+
  theme_ipsum(axis_title_size=15,axis_text_size=18,strip_text_size =18,base_size=15)+
  theme(legend.text=element_text(size=15))+
  scale_color_manual(name="Removal\ntreatment",
                     values=c("#FFD065","#FF7D2E"))+
  ylab("Home range area (km2)")+
  xlab("Removal period")+
  facet_wrap(~facet)

ggsave(paste0(outdir,"toxdt0_prediction_figure.png"),bg="white",height=4,width=6,units="in")

tidy.toxmodel=broom.mixed::tidy(tox.topmodel) %>% as.data.frame()
tidy.trapmodel=broom.mixed::tidy(trap.topmodel) %>% as.data.frame()
tidy.aermodel=broom.mixed::tidy(aer.topmodel) %>% as.data.frame()
tidy.toxdt0model=broom.mixed::tidy(tox.topmodeldt0) %>% as.data.frame()

tidy.toxmodel$removal="tox"
tidy.trapmodel$removal="trap"
tidy.aermodel$removal="aer"
tidy.toxdt0model$removal="toxdt0"

allparms=rbind(tidy.toxmodel, 
               tidy.trapmodel, 
               tidy.aermodel,
               tidy.toxdt0model)




# Format model data ------------------------------------------------------------

#Loop through model objects
for(i in 1:length(mods)){
  
  #Pull coefficients
  coefs=summary(mods[[i]])$coef$cond
  
  #Make table with just interaction effect, std.error, p value
  coefs2=as.data.frame(coefs[grep(":",rownames(coefs)),c(1,2,4),drop=FALSE])
  if(nrow(coefs2)>1){
    coefs2=as.data.frame(coefs2[grep("Removal.Type",rownames(coefs2)),,drop=FALSE])
    coefs2=as.data.frame(coefs2[grep("period",rownames(coefs2)),,drop=FALSE])
  }
  
  coefs2$model=names(mods)[i]
  coefs2$effect=rownames(coefs2)
  rownames(coefs2)=NULL
  
  if(i==1){
    allc=coefs2
  } else{
    allc=rbind(allc,coefs2)
  }
  
  
}

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
allc$response=NA
allc$response[grep("area",allc$model)]<-"area"
allc$response[grep("NSD",allc$model)]<-"NSD"

#Remove ugly effect column
allc=allc[,-which(colnames(allc)=="effect")]

# Format data for plotting -----------------------------------------------------
#Make col to adjust alpha transparency to p val
allc$alpha=0.1
allc$alpha[allc$`Pr(>|z|)`<0.05]<-1

#change sex NA
allc$sex[is.na(allc$sex)]<-"whole"

#Make keys
allc$rem_per_sex=paste0(allc$trt,allc$period,allc$sex)

#pivot wider for m/f estimates
#a1=allc %>% pivot_wider(id_cols=c(period,trt,response),names_from=sex,values_from=c(Estimate,alpha))

# Plot data ------------------------------------------------------------

#ta=allc[is.na(allc$sex),]
#ta=allc[allc$period=="after",]

df <- mtcars
df$col <- ifelse( df$cyl < 6, "cyl < 6", as.character(df$carb))
non_green <- unique(df$col[df$col != "cyl < 6"])
non_green <- non_green[order(non_green)]

aer_vals <- unique(allc$Estimate[allc$trt=="aer"])

ggplot(df, aes(x=qsec, y=hp,  colour= col )) +
  scale_color_manual(values = c(setNames(brewer.pal(length(non_green), name = "Reds"), non_green), "cyl < 6" = "green")) +
  geom_point() +
  theme_bw()

allc %>%
  ggplot(.,aes(x=response,y=rem_per_sex,fill=Estimate,alpha=alpha))+
  geom_tile()+
  theme_ipsum()+
  facet_wrap(~trt)

install.packages("ggnewscale")
library(ggnewscale)
aer_colors <- rev(c('#cdabff', "#822AFF", "#311854"))
trap_colors <- rev(c('#fcaed4', '#FF57A9', '#632242'))
tox_colors <- rev(c('#ffbf99', '#FF985A', '#8a4820'))


#dat1=allc[allc$sex=="whole",]
dat1=allc
ggplot() + 
  geom_tile(
    data = dat1 %>% filter(trt=="aer") %>% droplevels, 
    aes(response, rem_per_sex, fill=Estimate,alpha=alpha)
  ) + 
  scale_fill_gradientn(colors = aer_colors)+
  scale_alpha(range=c(0.25,1))+
  new_scale_fill() + 
  geom_tile(
    data = dat1 %>% filter(trt=="trap") %>% droplevels, 
    aes(response, rem_per_sex, fill=Estimate,alpha=alpha)
  )+
  scale_fill_gradientn(colors = trap_colors)+
  new_scale_fill() + 
  geom_tile(
    data = dat1 %>% filter(trt=="tox") %>% droplevels, 
    aes(response, rem_per_sex, fill=Estimate,alpha=alpha)
  )+
  scale_fill_gradientn(colors = tox_colors)+
  theme_ipsum()


dat2=allc[allc$sex!="whole",]
ggplot() + 
  geom_tile(
    data = dat2 %>% filter(trt=="aer") %>% droplevels, 
    aes(response, rem_per_sex, fill=Estimate,alpha=alpha)
  ) + 
  scale_fill_gradientn(colors = aer_colors)+
  scale_alpha(range=c(0.25,1))+
  new_scale_fill() + 
  geom_tile(
    data = dat2 %>% filter(trt=="trap") %>% droplevels, 
    aes(response, rem_per_sex, fill=Estimate,alpha=alpha)
  )+
  scale_fill_gradientn(colors = trap_colors)+
  new_scale_fill() + 
  geom_tile(
    data = dat2 %>% filter(trt=="tox") %>% droplevels, 
    aes(response, rem_per_sex, fill=Estimate,alpha=alpha)
  )+
  scale_fill_gradientn(colors = tox_colors)+
  theme_ipsum()




#current challenges for this heat map
#values scale needs to be different between each one
#heat map is not great when you have different units....
#because you'd have to have the color legends over on the right for each one
#as is currently, is just on one scale, which is also not appropriate for obvious reasons

#geom_col is way better for what we want to show here
#also note, bit misleading to show differences that weren't sig in the model
#moving forward, ONLY show removal types that had significant effects
#within those, use lower opacity/different color to indicate effects that were not significant

aer_colors <- rev(c('#cdabff', "#822AFF", "#311854"))
trap_colors <- rev(c('#fcaed4', '#FF57A9', '#632242'))
tox_colors <- rev(c('#ffbf99', '#FF985A', '#8a4820'))

#Next steps:
#1. determine which interactions in model are significant
# just do this manually, easier
#2. plot only significant effects-- including male/female intxn terms
#3. add labels to end of columns

p1=pdiffs5 %>% filter(rem=="tox", response=="area", sex=="whole") %>%
  ggplot()+
  geom_col(
    data = . %>% droplevels, 
    aes(x=key, y=pred, group=response),
    fill="#FF985A"
  )+theme_ipsum()

p2=pdiffs5 %>% filter(rem=="tox", response=="nsd", sex=="whole") %>%
  ggplot()+
  geom_col(
    data = . %>% droplevels, 
    aes(x=key, y=pred, group=response),
    fill="#FF985A"
  )+theme_ipsum()



