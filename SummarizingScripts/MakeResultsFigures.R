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
library(sdmTMB)

#Set dirs
input=file.path(home,"1_Data","Input",fsep=.Platform$file.sep)
objdir=file.path(home,"1_Data","Objects",fsep=.Platform$file.sep)
if(!dir.exists(file.path(home,"3_Output","Area_GLM_Results",fsep=.Platform$file.sep))){
  dir.create(file.path(home,"3_Output","Area_GLM_Results",fsep=.Platform$file.sep))
}
outdir=file.path(home,"3_Output",fsep=.Platform$file.sep)

#Read tidied model output
fs.mods=list.files(file.path(outdir,"Model_Output"),full.names=TRUE)
mod.names=list.files(file.path(outdir,"Model_Output"),full.names=FALSE)
mods=lapply(fs.mods,readRDS)
mod.names
params=rbind(mods[[grep("full_param",mod.names)[1]]],
      mods[[grep("full_param",mod.names)[2]]],
      mods[[grep("full_param",mod.names)[3]]])

intxn=rbind(mods[[grep("intxn",mod.names)[1]]],
             mods[[grep("intxn",mod.names)[2]]],
            mods[[grep("intxn",mod.names)[3]]])

preds=rbind(mods[[grep("preds",mod.names)[1]]],
            mods[[grep("preds",mod.names)[2]]],
            mods[[grep("preds",mod.names)[3]]])

#relevel periods
preds$per<-forcats::fct_relevel(preds$per,c("before","during","after"))

# Format pred diffs ------------------------------------------------------------

#Overall (sex is NA)
  #get differences for ctrl: btwn before-during and before-after
  #get differences for trt: btwn before-during and before-after
    #get differnces between ctrl/trt

#clean up
#pdiffs=preds[,c(1,2,6,7,9,10)]
pdiffs=preds[,c("trt",
                "predicted",
                "per",
                "sex",
                "model",
                "rem",
                "response")]

pdiffs2=tidyr::pivot_wider(pdiffs,id_cols=c(trt,sex,rem,response),values_from=predicted,names_from=per)

#period differences
pdiffs2$before_during=pdiffs2$during-pdiffs2$before
pdiffs2$before_after=pdiffs2$after-pdiffs2$before

head(pdiffs2)
pdiffs2$trt<-as.character(pdiffs2$trt)
pdiffs2$trt[pdiffs2$trt!="ctrl"]<-"trt"

#pivot wider to calc treatment diffs
pdiffs3<-tidyr::pivot_wider(pdiffs2,id_cols=c(sex,rem,response),values_from=c(before_during,before_after),names_from=trt)

#calc treatment diffs
pdiffs3$during_trt_ctrl<-pdiffs3$before_during_trt-pdiffs3$before_during_ctrl
pdiffs3$after_trt_ctrl<-pdiffs3$before_after_trt-pdiffs3$before_after_ctrl

#trim unneeded cols, don't need diffs within control/treatment
pdiffs4=pdiffs3[,c("sex",
                   "rem",
                   "response",
                   "during_trt_ctrl",
                   "after_trt_ctrl")]

#pivot_logner
pdiffs5=tidyr::pivot_longer(pdiffs4,cols=4:5,names_to="per",values_to="pred")

pdiffs5$period=NA
pdiffs5$period[grep("during",pdiffs5$per)]<-"during"
pdiffs5$period[grep("after",pdiffs5$per)]<-"after"
pdiffs5$period<-forcats::fct_relevel(pdiffs5$period,c("during","after"))

#remove NA rows (aerial, no during period)
pdiffs5=pdiffs5[!is.na(pdiffs5$pred),]

#replace NA in sex with whole model
pdiffs5$sex[is.na(pdiffs5$sex)]<-"whole"

#make key for preds to match with pdiffs easily
preds$sex[is.na(preds$sex)]<-"whole"
preds$sm<-"sex"
preds$sm[preds$sex=="whole"]<-"whole"
pdiffs5$sex[is.na(pdiffs5$sex)]<-"whole"
pdiffs5$sm<-"sex"
pdiffs5$sm[pdiffs5$sex=="whole"]<-"whole"

#make key
pdiffs5$key=paste(pdiffs5$response,pdiffs5$rem,pdiffs5$sm,sep="_")
preds$key=paste(preds$response,preds$rem,preds$sm,sep="_")

# Plot results -----------------------------------------------------------------

#preds for dot/whisker
#pdiffs5 for charts showing differences

keys=unique(preds$key)
trt.keys=keys[grep("area_aer_sex",keys)]
trt.col=aer.hex
preds_ylab_list=""
xlabtitle=""
# Plot results -----------------------------------------------------------------

MakeDotWhisker=function(trt.keys,trt.col,preds_ylab_list,xlabtitle,xaxislabels,legend,facet_label,yaxislabels,pt_size,l_size,axtitlesize){
  
  if(xaxislabels==TRUE){
p1=preds %>% filter(key==trt.keys) %>%
  ggplot(aes(x=per,yend=conf.high,group=trt,color=trt))+
  geom_segment(aes(y=conf.low),linewidth=l_size,position = position_dodge(width = 0.2))+
  geom_point(aes(y=predicted),size=pt_size,position = position_dodge(width = 0.2))+
  theme_ipsum(axis_title_size=axtitlesize,
              axis_text_size=8,
              strip_text_size=axtitlesize,
              base_size=15,
              caption_margin=2,
              plot_margin = margin(3, 3, 3, 3))+
  theme(legend.text=element_text(size=15))+
  scale_color_manual(name="Removal\ntreatment",
                     values=c("#FFD065",trt.col))+
  ylab(preds_ylab_list)+
  xlab(xlabtitle)+
  facet_wrap(~sex)
    }
  
  if(xaxislabels==FALSE){
    p1=preds %>% filter(key==trt.keys) %>%
      ggplot(aes(x=per,yend=conf.high,group=trt,color=trt))+
      geom_segment(aes(y=conf.low),linewidth=l_size,position = position_dodge(width = 0.2))+
      geom_point(aes(y=predicted),size=pt_size,position = position_dodge(width = 0.2))+
      theme_ipsum(axis_title_size=axtitlesize,
                  axis_text_size=8,
                  strip_text_size =axtitlesize,
                  base_size=15,
                  caption_margin=2,
                  plot_margin = margin(3, 3, 3, 3))+
      theme(legend.text=element_text(size=15),
            axis.text.x = element_blank(),
            axis.title.x = element_blank())+
      scale_color_manual(name="Removal\ntreatment",
                         values=c("#FFD065",trt.col))+
      ylab(preds_ylab_list)+
      xlab(xlabtitle)+
      facet_wrap(~sex)
  }
  
  if(legend==FALSE){
  p1=p1 + theme(legend.position="none")
  }
  
  if(facet_label==FALSE){
    p1=p1 + theme(strip.text=element_blank())
  }
  
  if(yaxislabels==FALSE){
    p1=p1+theme(axis.title.y=element_blank())
  }

return(p1)
}

#Run function to make paired plots
tox.hex="#FF985A"
trap.hex="#FF57A9"
aer.hex="#822AFF"
ctrl.hex="#ffc400"

# Make results figures (whole) ------------------------------------------------------------

# NSD
a_nsd_w=MakeDotWhisker(keys[grep("nsd_aer_whole",keys)],aer.hex,"NSD (m^2)","Removal Period",FALSE,FALSE,FALSE,TRUE)
t_nsd_w=MakeDotWhisker(keys[grep("nsd_trap_whole",keys)],trap.hex,"NSD (m^2)","Removal Period",FALSE,FALSE,FALSE,FALSE)
x_nsd_w=MakeDotWhisker(keys[grep("nsd_tox_whole",keys)],tox.hex,"NSD (m^2)","Removal Period",FALSE,FALSE,FALSE,FALSE)

# distance
a_dis_w=MakeDotWhisker(keys[grep("distance_aer_whole",keys)],aer.hex,"distance (km)","Removal Period",FALSE,FALSE,FALSE,TRUE)
t_dis_w=MakeDotWhisker(keys[grep("distance_trap_whole",keys)],trap.hex,"distance (km)","Removal Period",FALSE,FALSE,FALSE,FALSE)
x_dis_w=MakeDotWhisker(keys[grep("distance_tox_whole",keys)],tox.hex,"distance (km)","Removal Period",FALSE,FALSE,FALSE,FALSE)

# speed
a_sp_w=MakeDotWhisker(keys[grep("speed_aer_whole",keys)],aer.hex,"speed (km/hr)","Removal Period",FALSE,FALSE,FALSE,TRUE)
t_sp_w=MakeDotWhisker(keys[grep("speed_trap_whole",keys)],trap.hex,"speed (km/hr)","Removal Period",FALSE,FALSE,FALSE,FALSE)
x_sp_w=MakeDotWhisker(keys[grep("speed_tox_whole",keys)],tox.hex,"speed (km/hr)","Removal Period",FALSE,FALSE,FALSE,FALSE)

# area
a_are_w=MakeDotWhisker(keys[grep("area_aer_whole",keys)],aer.hex,"area (km^2)","Removal Period",TRUE,FALSE,FALSE,TRUE)
t_are_w=MakeDotWhisker(keys[grep("area_trap_whole",keys)],trap.hex,"area (km^2)","Removal Period",TRUE,FALSE,FALSE,FALSE)
x_are_w=MakeDotWhisker(keys[grep("area_tox_whole",keys)],tox.hex,"area (km^2)","Removal Period",TRUE,FALSE,FALSE,FALSE)

#combine removal types
nsd_pg=cowplot::plot_grid(a_nsd_w,t_nsd_w,x_nsd_w,nrow=1,rel_widths=c(0.65,0.65,0.65))
dis_pg=cowplot::plot_grid(a_dis_w,t_dis_w,x_dis_w,nrow=1,rel_widths=c(0.65,0.65,0.65))
sp_pg=cowplot::plot_grid(a_sp_w,t_sp_w,x_sp_w,nrow=1,rel_widths=c(0.65,0.65,0.65))
are_pg=cowplot::plot_grid(a_are_w,t_are_w,x_are_w,nrow=1,rel_widths=c(0.65,0.65,0.65))

#combine responses
allp=cowplot::plot_grid(nsd_pg,dis_pg,sp_pg,are_pg,ncol=1)
ggsave("~/Downloads/allp.png",allp,width=9,height=6.5,units="in")

#put asterisk on significant interactions
#doing manually in preview
  intxn=intxn[,c(8,4,5,6,7,1,2,3,9)]
  intxn[intxn$sex=="whole"&intxn$alpha==1,]
    #area:
      #tox reduction than controls
    #nsd:
      #tox increase than controls
    #distance:
      #tox: 
      #trap
    #speed
      #aer

# Make results figures (sex) ------------------------------------------------------------

  # NSD
  #MakeDotWhisker=function(trt.keys,trt.col,preds_ylab_list,xlabtitle,xaxislabels,legend,facet_label,yaxislabels,pt_size,l_size){
  #trt.keys,trt.col,preds_ylab_list,xlabtitle,xaxislabels,legend,facet_label,yaxislabels,sex_model,pt_size,l_size
  a_nsd_w_s=MakeDotWhisker(keys[grep("nsd_aer_sex",keys)],aer.hex,"NSD (m^2)","Removal Period", FALSE,FALSE,TRUE,TRUE,4,2)
  t_nsd_w_s=MakeDotWhisker(keys[grep("nsd_trap_sex",keys)],trap.hex,"NSD (m^2)","Removal Period",FALSE,FALSE,TRUE,FALSE,4,2)
  x_nsd_w_s=MakeDotWhisker(keys[grep("nsd_tox_sex",keys)],tox.hex,"NSD (m^2)","Removal Period",FALSE,FALSE,TRUE,FALSE,4,2)
  
  # distance
  a_dis_w_s=MakeDotWhisker(keys[grep("distance_aer_sex",keys)],aer.hex,"distance (km)","Removal Period",FALSE,FALSE,FALSE,TRUE,4,2)
  t_dis_w_s=MakeDotWhisker(keys[grep("distance_trap_sex",keys)],trap.hex,"distance (km)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,2)
  x_dis_w_s=MakeDotWhisker(keys[grep("distance_tox_sex",keys)],tox.hex,"distance (km)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,2)
  
  # speed
  a_sp_w_s=MakeDotWhisker(keys[grep("speed_aer_sex",keys)],aer.hex,"speed (km/hr)","Removal Period",FALSE,FALSE,FALSE,TRUE,4,2)
  t_sp_w_s=MakeDotWhisker(keys[grep("speed_trap_sex",keys)],trap.hex,"speed (km/hr)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,2)
  x_sp_w_s=MakeDotWhisker(keys[grep("speed_tox_sex",keys)],tox.hex,"speed (km/hr)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,2)
    
  # area
  a_are_w_s=MakeDotWhisker(keys[grep("area_aer_sex",keys)],aer.hex,"area (km^2)","Removal Period",TRUE,FALSE,FALSE,TRUE,4,2)
  #t_are_w_s=MakeDotWhisker(keys[grep("area_trap_sex",keys)],trap.hex,"area (km^2)","Removal Period",TRUE,FALSE,TRUE,FALSE)
  x_are_w_s=MakeDotWhisker(keys[grep("area_tox_sex",keys)],tox.hex,"area (km^2)","Removal Period",TRUE,FALSE,FALSE,FALSE,4,2)
  
  #combine removal types
  nsd_pg_s=cowplot::plot_grid(a_nsd_w_s,t_nsd_w_s,x_nsd_w_s,nrow=1,rel_widths=c(0.65,0.65,0.65))
  dis_pg_s=cowplot::plot_grid(a_dis_w_s,t_dis_w_s,x_dis_w_s,nrow=1,rel_widths=c(0.65,0.65,0.65))
  sp_pg_s=cowplot::plot_grid(a_sp_w_s,t_sp_w_s,x_sp_w_s,nrow=1,rel_widths=c(0.65,0.65,0.65))
  are_pg_s=cowplot::plot_grid(a_are_w_s,NULL,x_are_w_s,nrow=1,rel_widths=c(0.65,0.65,0.65))
  
  #combine responses
  allp_s=cowplot::plot_grid(nsd_pg_s,dis_pg_s,sp_pg_s,are_pg_s,ncol=1,rel_heights=c(1,0.9,0.9,0.9))
  ggsave("~/Downloads/allps.png",allp_s,width=9,height=6.5,units="in")
  
  #put asterisk on significant interactions
  #doing manually in preview
  #intxn=intxn[,c(8,4,5,6,7,1,2,3,9)]
  intxn[intxn$sex!="whole"&intxn$alpha==1,]
  #area:
  #tox reduction than controls
  #nsd:
  #tox increase than controls
  #distance:
  #tox: 
  #trap
  #speed
  #aer
  
# Make powerpoint figures -----------------
  
  #MakeDotWhisker=function(trt.keys,trt.col,preds_ylab_list,xlabtitle,xaxislabels,legend,facet_label,yaxislabels,pt_size,l_size){

  # NSD
  a_nsd_w=MakeDotWhisker(keys[grep("nsd_aer_whole",keys)],aer.hex,"NSD (m^2)","Removal Period",TRUE,FALSE,FALSE,TRUE,4,1,14)
  t_nsd_w=MakeDotWhisker(keys[grep("nsd_trap_whole",keys)],trap.hex,"NSD (m^2)","Removal Period",TRUE,FALSE,FALSE,TRUE,4,1,14)
  x_nsd_w=MakeDotWhisker(keys[grep("nsd_tox_whole",keys)],tox.hex,"NSD (m^2)","Removal Period",TRUE,FALSE,FALSE,TRUE,4,1,14)
  
  # distance
  a_dis_w=MakeDotWhisker(keys[grep("distance_aer_whole",keys)],aer.hex,"distance (km)","Removal Period",TRUE,FALSE,FALSE,TRUE,4,1,14)
  t_dis_w=MakeDotWhisker(keys[grep("distance_trap_whole",keys)],trap.hex,"distance (km)","Removal Period",TRUE,FALSE,FALSE,TRUE,4,1,14)
  x_dis_w=MakeDotWhisker(keys[grep("distance_tox_whole",keys)],tox.hex,"distance (km)","Removal Period",TRUE,FALSE,FALSE,TRUE,4,1,14)
  
  # speed
  a_sp_w=MakeDotWhisker(keys[grep("speed_aer_whole",keys)],aer.hex,"speed (km/hr)","Removal Period",TRUE,FALSE,FALSE,TRUE,4,1,14)
  t_sp_w=MakeDotWhisker(keys[grep("speed_trap_whole",keys)],trap.hex,"speed (km/hr)","Removal Period",TRUE,FALSE,FALSE,TRUE,4,1,14)
  x_sp_w=MakeDotWhisker(keys[grep("speed_tox_whole",keys)],tox.hex,"speed (km/hr)","Removal Period",TRUE,FALSE,FALSE,TRUE,4,1,14)
  
  # area
  a_are_w=MakeDotWhisker(keys[grep("area_aer_whole",keys)],aer.hex,"area (km^2)","Removal Period",TRUE,FALSE,FALSE,TRUE,4,1,14)
  t_are_w=MakeDotWhisker(keys[grep("area_trap_whole",keys)],trap.hex,"area (km^2)","Removal Period",TRUE,FALSE,FALSE,TRUE,4,1,14)
  x_are_w=MakeDotWhisker(keys[grep("area_tox_whole",keys)],tox.hex,"area (km^2)","Removal Period",TRUE,FALSE,FALSE,TRUE,4,1,14)
  
  #combine removal types
  aer_ppt=cowplot::plot_grid(a_nsd_w,a_dis_w,a_sp_w,a_are_w,nrow=1,rel_widths=c(0.65,0.65,0.65,0.65))
  ggsave("~/Downloads/aerppt.png",aer_ppt,width=8,height=2,units="in")
  
  trap_ppt=cowplot::plot_grid(t_nsd_w,t_dis_w,t_sp_w,t_are_w,nrow=1,rel_widths=c(0.65,0.65,0.65,0.65))
  ggsave("~/Downloads/trapppt.png",trap_ppt,width=8,height=2,units="in")
  
  aer_ppt=cowplot::plot_grid(x_nsd_w,x_dis_w,x_sp_w,x_are_w,nrow=1,rel_widths=c(0.65,0.65,0.65,0.65))
  ggsave("~/Downloads/toxppt.png",aer_ppt,width=8,height=2,units="in")
  
# Sex interactions --------
  # NSD
  a_nsd_w=MakeDotWhisker(keys[grep("nsd_aer_sex",keys)],aer.hex,"NSD (m^2)","Removal Period",TRUE,FALSE,TRUE,TRUE,4,1,12)
  t_nsd_w=MakeDotWhisker(keys[grep("nsd_trap_sex",keys)],trap.hex,"NSD (m^2)","Removal Period",TRUE,FALSE,TRUE,TRUE,4,1,12)
  x_nsd_w=MakeDotWhisker(keys[grep("nsd_tox_sex",keys)],tox.hex,"NSD (m^2)","Removal Period",TRUE,FALSE,TRUE,TRUE,4,1,12)
  
  # distance
  a_dis_w=MakeDotWhisker(keys[grep("distance_aer_sex",keys)],aer.hex,"distance (km)","Removal Period",TRUE,FALSE,TRUE,TRUE,4,1,12)
  t_dis_w=MakeDotWhisker(keys[grep("distance_trap_sex",keys)],trap.hex,"distance (km)","Removal Period",TRUE,FALSE,TRUE,TRUE,4,1,12)
  x_dis_w=MakeDotWhisker(keys[grep("distance_tox_sex",keys)],tox.hex,"distance (km)","Removal Period",TRUE,FALSE,TRUE,TRUE,4,1,12)
  
  # speed
  a_sp_w=MakeDotWhisker(keys[grep("speed_aer_sex",keys)],aer.hex,"speed (km/hr)","Removal Period",TRUE,FALSE,TRUE,TRUE,4,1,12)
  t_sp_w=MakeDotWhisker(keys[grep("speed_trap_sex",keys)],trap.hex,"speed (km/hr)","Removal Period",TRUE,FALSE,TRUE,TRUE,4,1,12)
  x_sp_w=MakeDotWhisker(keys[grep("speed_tox_sex",keys)],tox.hex,"speed (km/hr)","Removal Period",TRUE,FALSE,TRUE,TRUE,4,1,12)
  
  # area
  a_are_w=MakeDotWhisker(keys[grep("area_aer_sex",keys)],aer.hex,"area (km^2)","Removal Period",TRUE,FALSE,TRUE,TRUE,4,1,12)
  #t_are_w=MakeDotWhisker(keys[grep("area_trap_sex",keys)],trap.hex,"area (km^2)","Removal Period",TRUE,FALSE,TRUE,TRUE,4,1,14)
  x_are_w=MakeDotWhisker(keys[grep("area_tox_sex",keys)],tox.hex,"area (km^2)","Removal Period",TRUE,FALSE,TRUE,TRUE,4,1,12)
  
  #combine removal types
  aer_ppt=cowplot::plot_grid(a_nsd_w,a_dis_w,a_sp_w,a_are_w,nrow=1,rel_widths=c(0.65,0.65,0.65,0.65))
  ggsave("~/Downloads/s_aerppt.png",aer_ppt,width=9,height=2,units="in")
  
  trap_ppt=cowplot::plot_grid(t_nsd_w,t_dis_w,t_sp_w,nrow=1,rel_widths=c(0.65,0.65,0.65))
  ggsave("~/Downloads/s_trapppt.png",trap_ppt,width=9,height=2,units="in")
  
  aer_ppt=cowplot::plot_grid(x_nsd_w,x_dis_w,x_sp_w,x_are_w,nrow=1,rel_widths=c(0.65,0.65,0.65,0.65))
  ggsave("~/Downloads/s_toxppt.png",aer_ppt,width=9,height=2,units="in")
  
