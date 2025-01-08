#set home dir of pipeline
home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline"

# Purpose ----------------------------------------------------------------------

#The purpose of this script is run autocorrelation diagnostics and glmms for akde home range size estimates

# In/Out -----------------------------------------------------------------------

#inputs:
#outdf_akde_aer_corrected_f.rds, outdf_akde_tox_corrected_f.rds, outdf_akde_trap_corrected_f.rds

#outputs: stitched results figure

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

#Read distance input objects
speedaer=readRDS(file.path(objdir,"speedaer.rds",fsep=.Platform$file.sep))
speedtox=readRDS(file.path(objdir,"speedtox.rds",fsep=.Platform$file.sep))
speedtrap=readRDS(file.path(objdir,"speedtrap.rds",fsep=.Platform$file.sep))

#Read speed input objects
distaer=readRDS(file.path(objdir,"distaer.rds",fsep=.Platform$file.sep))
disttox=readRDS(file.path(objdir,"disttox.rds",fsep=.Platform$file.sep))
disttrap=readRDS(file.path(objdir,"distrap.rds",fsep=.Platform$file.sep))

#Read area input objects
akaer=readRDS(file.path(objdir,"outdf_akde_aer_corrected_f.rds",fsep=.Platform$file.sep))
aktrap=readRDS(file.path(objdir,"outdf_akde_trap_corrected_f.rds",fsep=.Platform$file.sep))
aktox=readRDS(file.path(objdir,"outdf_akde_aktox_corrected_f.rds",fsep=.Platform$file.sep))

#Read in data to use for NSD analysis:
geo.aerd.wk=readRDS(file.path(objdir,"geo_aerd_wk.rds",fsep=.Platform$file.sep))
geo.trapd.wk=readRDS(file.path(objdir,"geo_trapd_wk.rds",fsep=.Platform$file.sep))
geo.toxd.wk=readRDS(file.path(objdir,"geo_toxd_wk.rds",fsep=.Platform$file.sep))

#Read in model objects
mods=lapply(list.files(file.path(objdir,"Models",fsep=.Platform$file.sep),full.names=TRUE),readRDS)
mod_names=list.files(file.path(objdir,"Models",fsep=.Platform$file.sep))
names(mods)=stringr::str_sub(mod_names,5L,-5L)

#model read-ins from Abbey are broken,
#probably lame version issue. 
#just re-eval calls, since can still pull those and have the data
rerun_seq=c(6:11,18:21)
for(i in rerun_seq){
  print(i)
  call=mods[[i]]$call
  print(call)
  #mods[[i]]=eval(call)
}


#source needed functions
funcdir<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/2_Scripts/Functions"
func.list=list.files(funcdir,full.names=TRUE)
for(f in 1:length(func.list)){
  source(func.list[f])
}

# Format model predictions ------------------------------------------------------------
library(ggeffects)
for(i in 1:length(mods)){
  
coefs=summary(mods[[i]])$coef$cond

if(length(grep("sex",rownames(coefs)))==0){
if(length(grep("area",names(mods)[i]))==0){
tmp=ggeffects::predict_response(mods[[i]], terms=c("Removal.Type","removal.period.akdecalc"))
} 
if(length(grep("area",names(mods)[i]))>0){
tmp=ggeffects::predict_response(mods[[i]], terms=c("Removal.Type","period"))
} 
} else{
  if(length(grep("area",names(mods)[i]))==0){
    tmp=ggeffects::predict_response(mods[[i]], terms=c("Removal.Type","removal.period.akdecalc","sex"))
  } 
  if(length(grep("area",names(mods)[i]))>0){
    tmp=ggeffects::predict_response(mods[[i]], terms=c("Removal.Type","period","sex"))
  }   
}

tmp.df=as.data.frame(tmp)
if(length(grep("facet",colnames(tmp.df)))==0){
  tmp.df$facet=NA
}
tmp.df$model=names(mods)[i]

if(i==1){
  preds=tmp.df
} else{
  preds=rbind(preds,tmp.df)
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
preds$response[grep("NSD",preds$model)]<-"nsd"
preds$response[grep("area",preds$model)]<-"area"
preds$response[grep("speed",preds$model)]<-"speed"
preds$response[grep("distance",preds$model)]<-"distance"

write.csv(preds,"~/Desktop/rr_preds.csv")
#This output has been formatted into tables in excel
#and saved as Predictions.xlsx in the Outputs folder
#Formatting done for NSD/area so far.. need add abbey's responses

# Pull interactions table ------------------------------------------------------

#Loop through model objects
for(i in 1:length(mods)){
  
  #Pull coefficients
  coefs=summary(mods[[i]])$coef$cond
  
  #Make table with just interaction effect, std.error, p value
  coefs2=as.data.frame(coefs[grep(":",rownames(coefs)),c(1,2,4),drop=FALSE])
  if(nrow(coefs2)!=0){

  if(nrow(coefs2)>1){
    
    coefs2=as.data.frame(coefs2[grep("Removal.Type",rownames(coefs2)),,drop=FALSE])
    coefs2=as.data.frame(coefs2[grep("period",rownames(coefs2)),,drop=FALSE])
  
    
  }
  
  if(nrow(coefs2)!=0){
  coefs2$model=names(mods)[i]
  coefs2$effect=rownames(coefs2)
  rownames(coefs2)=NULL
  
  if(i==1){
    allc=coefs2
  } else{
    allc=rbind(allc,coefs2)
  }
  
  }
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
allc$response[grep("speed",allc$model)]<-"speed"
allc$response[grep("distance",allc$model)]<-"distance"

#Remove ugly effect column
allc=allc[,-which(colnames(allc)=="effect")]

#change sex NA
allc$sex[is.na(allc$sex)]<-"whole"

#subset to only significant interactions
allc$alpha=0.1
allc$alpha[allc$`Pr(>|z|)`<0.05]<-1
sigintxn=allc[allc$alpha==1,]

# Filter out only significant effects for plotting -----------------------------

#reference sigintxn
unique(sigintxn$model)

#preds
View(preds)
predsg=preds[preds$model%in%unique(sigintxn$model),]

#pdiffs5
#remove aerial
#remove whole 
#temporarily remove anything with "distance_rp_" to get pipe to run
predsg=predsg[-grep("distance_rp_",predsg$model),]
predsg=predsg[-grep("distance_rps_tox",predsg$model),]

# Format pred diffs ------------------------------------------------------------

#Overall (sex is NA)
  #get differences for ctrl: btwn before-during and before-after
  #get differences for trt: btwn before-during and before-after
    #get differnces between ctrl/trt

#clean up
pdiffs=predsg[,c(1,2,6,7,9,10)]

pdiffs2=pivot_wider(pdiffs,id_cols=c(trt,sex,rem,response),values_from=predicted,names_from=per)


#period differences
pdiffs2$before_during=pdiffs2$during-pdiffs2$before
pdiffs2$before_after=pdiffs2$after-pdiffs2$before

head(pdiffs2)
pdiffs2$trt<-as.character(pdiffs2$trt)
pdiffs2$trt[pdiffs2$trt!="ctrl"]<-"trt"

#pivot wider to calc treatment diffs
pdiffs3<-pivot_wider(pdiffs2,id_cols=c(sex,rem,response),values_from=c(before_during,before_after),names_from=trt)

#calc treatment diffs
pdiffs3$during_trt_ctrl<-pdiffs3$before_during_trt-pdiffs3$before_during_ctrl
pdiffs3$after_trt_ctrl<-pdiffs3$before_after_trt-pdiffs3$before_after_ctrl

#trim unneeded cols
pdiffs4=pdiffs3[,c(1:3,8,9)]

#pivot_logner
pdiffs5=pivot_longer(pdiffs4,cols=4:5,names_to="per",values_to="pred")

pdiffs5$period=NA
pdiffs5$period[grep("during",pdiffs5$per)]<-"during"
pdiffs5$period[grep("after",pdiffs5$per)]<-"after"
unique(pdiffs$per)

#remove NA rows
pdiffs5=pdiffs5[!is.na(pdiffs5$pred),]

#replace NA in sex with whole model
pdiffs5$sex[is.na(pdiffs5$sex)]<-"whole"


#make key for preds to match with pdiffs easily
predsg$sex[is.na(predsg$sex)]<-"whole"
predsg$sm<-"sex"
predsg$sm[predsg$sex=="whole"]<-"whole"
pdiffs5$sex[is.na(pdiffs5$sex)]<-"whole"
pdiffs5$sm<-"sex"
pdiffs5$sm[pdiffs5$sex=="whole"]<-"whole"

#make key
pdiffs5$key=paste(pdiffs5$response,pdiffs5$rem,pdiffs5$sm,sep="_")
predsg$key=paste(predsg$response,predsg$rem,predsg$sm,sep="_")

# Plot results -----------------------------------------------------------------

#predsg for dot/whisker
#pdiffs5 for charts showing differences

keys=unique(predsg$key)

#choose to just show sex differences if significant
trt.keys=keys[c(1:4)]
trt.col="#FF7D2E"
pg.list=vector(mode="list",length=length(trt.keys))
preds_ylab_list=c(
  "Home range area (km2)",
  "NSD (m2)"
)
diffs_ylab_list=c(
  "Home range area (km2) difference",
  "NSD (m2) difference"
)



MakePairCharts=function(trt.keys,trt.col,preds_ylab_list,diffs_ylab_list){
  pg.list=vector(mode="list",length=length(trt.keys))
  
  for(k in 1:length(trt.keys)){
#subset by key
  
p1=predsg %>% filter(key==trt.keys[k]) %>%
  ggplot(aes(x=per,yend=conf.high,group=trt,color=trt))+
  geom_point(aes(y=predicted),size=6,position = position_dodge(width = 0.2))+
  geom_segment(aes(y=conf.low),linewidth=3,position = position_dodge(width = 0.2))+
  theme_ipsum(axis_title_size=15,axis_text_size=18,strip_text_size =18,base_size=15)+
  theme(legend.text=element_text(size=15))+
  scale_color_manual(name="Removal\ntreatment",
                     values=c("#FFD065",trt.col))+
  ylab(preds_ylab_list[k])+
  xlab("Removal period")+
  facet_wrap(~sex)

#barchart pair
maxp=ceiling(abs(max(pdiffs5[pdiffs5$key==trt.keys[k],]$pred)))
minp=ceiling(abs(min(pdiffs5[pdiffs5$key==trt.keys[k],]$pred)))
scale_max=max(maxp,minp)

if(scale_max==1){
  maxp=abs(max(pdiffs5[pdiffs5$key==trt.keys[k],]$pred))+0.1*abs(max(pdiffs5[pdiffs5$key==trt.keys[k],]$pred))
  minp=abs(min(pdiffs5[pdiffs5$key==trt.keys[k],]$pred))+0.1*abs(min(pdiffs5[pdiffs5$key==trt.keys[k],]$pred))
  scale_max=max(maxp,minp)
  }

p2=pdiffs5 %>% filter(key==trt.keys[k]) %>%
  ggplot()+
  geom_col(
    data = . %>% droplevels, 
    aes(x=period, y=pred, group=response),
    fill=trt.col
  )+ylim(-scale_max,scale_max)+
  theme_ipsum(axis_title_size=15,axis_text_size=18,strip_text_size =18,base_size=15)+
  ylab(diffs_ylab_list[k])+
  facet_wrap(~sex)
}
return(list(p1,p2))
}

#Run function to make paired plots

preds_ylab_list=c(
  "Home range area (km2)"
)
diffs_ylab_list=c(
  "Home range area (km2) difference"
  )
trt.keys=keys[c(2)]
trt.col="#FF7D2E"
plist=MakePairCharts(trt.keys,trt.col,preds_ylab_list,diffs_ylab_list)
k1=cowplot::plot_grid(plist[[1]],plist[[2]],nrow=1,rel_widths=c(1,0.8),labels=c("Predictions","Differences"))
k1

preds_ylab_list=c(
  "NSD (m2)"
)
diffs_ylab_list=c(
  "NSD (m2) difference"
)
trt.keys=keys[c(4)]
trt.col="#FF7D2E"
plist=MakePairCharts(trt.keys,trt.col,preds_ylab_list,diffs_ylab_list)
k2=cowplot::plot_grid(plist[[1]],plist[[2]],nrow=1,rel_widths=c(1,0.8))
k2

preds_ylab_list=c(
  "NSD (m2)"
)
diffs_ylab_list=c(
  "NSD (m2) difference"
)
trt.keys=keys[c(5)]
trt.col="#F74DD4"
plist=MakePairCharts(trt.keys,trt.col,preds_ylab_list,diffs_ylab_list)
k3=cowplot::plot_grid(plist[[1]],plist[[2]],nrow=1,rel_widths=c(1,0.8))
k3


preds_ylab_list=c(
  "speed (km/hr)"
)
diffs_ylab_list=c(
  "speed (km/hr) difference"
)
trt.keys=keys[c(6)]
trt.col="#F74DD4"
plist=MakePairCharts(trt.keys,trt.col,preds_ylab_list,diffs_ylab_list)
k4=cowplot::plot_grid(plist[[1]],plist[[2]],nrow=1,rel_widths=c(1,0.8))
k4

# Stitch plots together -----------------------------------------------------------------

out=cowplot::plot_grid(k1,
                       k2,
                       k3,
                       ncol=1)

ggsave("~/Desktop/plot_out.png",plot=out,height=10,width=14,units="in")

#2.5,000,000


# Scrap code -----------------------------------------------------------------
preds
ggplot(preds,aes(x=x,yend=conf.high,group=group,color=group))+
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






