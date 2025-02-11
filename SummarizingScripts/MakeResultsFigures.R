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
library(gt)
library(gtsummary)
library(colorspace)

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

#add response to params
mods[[grep("full_param",mod.names)[[1]]]]$response="area"
mods[[grep("full_param",mod.names)[[2]]]]$response="nsd"
mods[[grep("full_param",mod.names)[[3]]]]$response="spdist"

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

# Pull in gt objects -----------------------------------------------------------

## Parameter table for movement removal*period models --------------------------

area_tbl=readRDS(file.path(outdir,"Model_Output","area_parm_gt.rds",fsep=.Platform$file.sep))
nsd_tbl=readRDS(file.path(outdir,"Model_Output","nsd_parm_gt.rds",fsep=.Platform$file.sep))
dist_tbl=readRDS(file.path(outdir,"Model_Output","dist_parm_gt.rds",fsep=.Platform$file.sep))
speed_tbl=readRDS(file.path(outdir,"Model_Output","speed_parm_gt.rds",fsep=.Platform$file.sep))

stk=tbl_stack(tbls=list(area_tbl,nsd_tbl,dist_tbl,speed_tbl),
              group_header=c("Area (km^2)","NSD (m^2)","distance (km)","speed (km/hr)"))

stk2=as_gt(stk) %>%             # convert gtsummary table to gt
  gt::tab_style(           # use gt::tab_style() to shade column
    style = list(gt::cell_fill(color = "white")),
    locations = gt::cells_body(columns = c(estimate_1,ci_1,p.value_1))
  ) %>%
  gt::tab_style(           # use gt::tab_style() to shade column
    style = list(gt::cell_fill(color = "white")),
    locations = gt::cells_body(columns = c(estimate_2,ci_2,p.value_2))
  ) %>%
  gt::tab_style(           # use gt::tab_style() to shade column
    style = list(gt::cell_fill(color = "white")),
    locations = gt::cells_body(columns = c(estimate_3,ci_3,p.value_3))
  ) %>%
  gt::tab_style(
    style=cell_text(size=11),
    locations=gt::cells_body(columns=c(estimate_1,ci_1,p.value_1,
                                       estimate_2,ci_2,p.value_2,
                                       estimate_3,ci_3,p.value_3),
                             rows=everything())
  ) %>%
  gt::tab_style(
    style=list(cell_text(size=11,weight="bold",color="black"),
               cell_fill(color="white"),
               cell_borders(sides="bottom",weight=px(3))),
    locations=gt::cells_row_groups()
  ) %>%
  gt::tab_style(
    style=list(cell_text(size=11,weight="bold",color="black"),
               cell_fill(color="white"),
               cell_borders(sides="bottom",weight=px(3),color="white")),
    locations=list(gt::cells_column_labels(),gt::cells_column_spanners())
  ) %>%
  gt::tab_style(
    style=list(cell_text(size=11,weight="bold"),
               cell_borders(sides="bottom",weight=px(1))),
    locations=list(gt::cells_body(columns=label,rows=c(1,4,8,11,14,18,21,24,28,31,34,38)))
  ) %>% gt::tab_style(
    style=list(cell_borders(sides="bottom",weight=px(1))),
    locations=list(gt::cells_body(columns=everything(),rows=c(which(stk$table_body$label=="period"),which(stk$table_body$label=="trt_ctrl"),which(stk$table_body$label=="trt_ctrl * period")))
    )
    ) %>%
  gt::tab_style(
    style=list(cell_text(size=11)),
    locations=list(gt::cells_body(columns=label,rows=c(which(stk$table_body$label=="period"),which(stk$table_body$label=="trt_ctrl"),which(stk$table_body$label=="trt_ctrl * period")))
   )) %>%
  tab_options(data_row.padding = px(1))

#For now, taking screenshot to put on sharepoint
stk2 |> gtsave(file.path(outdir,"Parameter_Table","mvmt_rp_parms.html"))

## Parameter table for movement removal*period*sex models ----------------------

#need change order of all tables to aer, tox, trap
area_tbl_s=readRDS(file.path(outdir,"Model_Output","area_parm_gt_s.rds",fsep=.Platform$file.sep))
nsd_tbl_s=readRDS(file.path(outdir,"Model_Output","nsd_parm_gt_s.rds",fsep=.Platform$file.sep))
dist_tbl_s=readRDS(file.path(outdir,"Model_Output","dist_parm_gt_s.rds",fsep=.Platform$file.sep))
speed_tbl_s=readRDS(file.path(outdir,"Model_Output","speed_parm_gt_s.rds",fsep=.Platform$file.sep))

stk=tbl_stack(tbls=list(nsd_tbl_s,area_tbl_s,dist_tbl_s,speed_tbl_s),
              group_header=c("NSD (m2)","Area (km2)","distance (km)","speed (km/hr)"))

stk2<-as_gt(stk) %>%             # convert gtsummary table to gt
  gt::tab_style(           # use gt::tab_style() to shade column
    style = list(gt::cell_fill(color = "white")),
    locations = gt::cells_body(columns = c(estimate_1,ci_1,p.value_1))
  ) %>%
  gt::tab_style(           # use gt::tab_style() to shade column
    style = list(gt::cell_fill(color = "white")),
    locations = gt::cells_body(columns = c(estimate_2,ci_2,p.value_2))
  ) %>%
  gt::tab_style(           # use gt::tab_style() to shade column
    style = list(gt::cell_fill(color = "white")),
    locations = gt::cells_body(columns = c(estimate_3,ci_3,p.value_3))
  ) %>%
  gt::tab_style(
    style=cell_text(size=11),
    locations=gt::cells_body(columns=c(estimate_1,ci_1,p.value_1,
                                       estimate_2,ci_2,p.value_2,
                                       estimate_3,ci_3,p.value_3),
                             rows=everything())
  ) %>%
  gt::tab_style(
    style=list(cell_text(size=11,weight="bold",color="black"),
               cell_fill(color="white"),
               cell_borders(sides="bottom",weight=px(3))),
    locations=gt::cells_row_groups()
  ) %>%
  gt::tab_style(
    style=list(cell_text(size=11,weight="bold",color="black"),
               cell_fill(color="white"),
               cell_borders(sides="bottom",weight=px(3),color="white")),
    locations=list(gt::cells_column_labels(),gt::cells_column_spanners())
  ) %>%
  gt::tab_style(
    style=list(cell_text(size=11,weight="bold"),
               cell_borders(sides="bottom",weight=px(2))),
    locations=list(gt::cells_body(columns=everything(),rows=c(which(stk$table_body$row_type=="label"))))
  ) %>% gt::tab_style(
    style=list(cell_borders(sides="bottom",weight=px(1))),
    locations=list(gt::cells_body(columns=everything(),rows=which(stk$table_body$row_type!="label")))) 

stk2 |> gtsave(file.path(outdir,"Parameter_Table","mvmt_rps_parms.html"))


## Parameter table for contact removal*period models ---------------------------

ncon_tbl=readRDS(file.path(outdir,"Model_Output","ncon_parm_gt.rds",fsep=.Platform$file.sep))
nind_tbl=readRDS(file.path(outdir,"Model_Output","nind_parm_gt.rds",fsep=.Platform$file.sep))
#nind_tbl has sex, shouldn't
stk=tbl_stack(tbls=list(ncon_tbl,nind_tbl),
              group_header=c("N total contacts",
                             "N unique contacts"))

stk2<-as_gt(stk) %>%             # convert gtsummary table to gt
  gt::tab_style(           # use gt::tab_style() to shade column
    style = list(gt::cell_fill(color = "white")),
    locations = gt::cells_body(columns = c(estimate_1,ci_1,p.value_1))
  ) %>%
  gt::tab_style(           # use gt::tab_style() to shade column
    style = list(gt::cell_fill(color = "white")),
    locations = gt::cells_body(columns = c(estimate_2,ci_2,p.value_2))
  ) %>%
  gt::tab_style(           # use gt::tab_style() to shade column
    style = list(gt::cell_fill(color = "white")),
    locations = gt::cells_body(columns = c(estimate_3,ci_3,p.value_3))
  ) %>%
  gt::tab_style(
    style=cell_text(size=11),
    locations=gt::cells_body(columns=c(estimate_1,ci_1,p.value_1,
                                       estimate_2,ci_2,p.value_2,
                                       estimate_3,ci_3,p.value_3),
                             rows=everything())
  ) %>%
  gt::tab_style(
    style=list(cell_text(size=11,weight="bold",color="black"),
               cell_fill(color="white"),
               cell_borders(sides="bottom",weight=px(3))),
    locations=gt::cells_row_groups()
  ) %>%
  gt::tab_style(
    style=list(cell_text(size=11,weight="bold",color="black"),
               cell_fill(color="white"),
               cell_borders(sides="bottom",weight=px(3),color="white")),
    locations=list(gt::cells_column_labels(),gt::cells_column_spanners())
  ) %>%
  gt::tab_style(
    style=list(cell_text(size=11,weight="bold"),
               cell_borders(sides="bottom",weight=px(2))),
    locations=list(gt::cells_body(columns=everything(),rows=c(which(stk$table_body$row_type=="label"))))
  ) %>% gt::tab_style(
    style=list(cell_borders(sides="bottom",weight=px(1))),
    locations=list(gt::cells_body(columns=everything(),rows=which(stk$table_body$row_type!="label")))) 

stk2 |> gtsave(file.path(outdir,"Parameter_Table","contact_rp_parms.html"))

## Parameter table for contact removal*period*sex models -----------------------

ncon_tbl_s=readRDS(file.path(outdir,"Model_Output","ncon_parm_gt_s.rds",fsep=.Platform$file.sep))
nind_tbl_s=readRDS(file.path(outdir,"Model_Output","nind_parm_gt_s.rds",fsep=.Platform$file.sep))

stk=tbl_stack(tbls=list(ncon_tbl_s,nind_tbl_s),
              group_header=c("N total contacts",
                             "N unique contacts"))

stk2<-as_gt(stk) %>%             # convert gtsummary table to gt
  gt::tab_style(           # use gt::tab_style() to shade column
    style = list(gt::cell_fill(color = "white")),
    locations = gt::cells_body(columns = c(estimate_1,ci_1,p.value_1))
  ) %>%
  gt::tab_style(           # use gt::tab_style() to shade column
    style = list(gt::cell_fill(color = "white")),
    locations = gt::cells_body(columns = c(estimate_2,ci_2,p.value_2))
  ) %>%
  gt::tab_style(           # use gt::tab_style() to shade column
    style = list(gt::cell_fill(color = "white")),
    locations = gt::cells_body(columns = c(estimate_3,ci_3,p.value_3))
  ) %>%
  gt::tab_style(
    style=cell_text(size=11),
    locations=gt::cells_body(columns=c(estimate_1,ci_1,p.value_1,
                                       estimate_2,ci_2,p.value_2,
                                       estimate_3,ci_3,p.value_3),
                             rows=everything())
  ) %>%
  gt::tab_style(
    style=list(cell_text(size=11,weight="bold",color="black"),
               cell_fill(color="white"),
               cell_borders(sides="bottom",weight=px(3))),
    locations=gt::cells_row_groups()
  ) %>%
  gt::tab_style(
    style=list(cell_text(size=11,weight="bold",color="black"),
               cell_fill(color="white"),
               cell_borders(sides="bottom",weight=px(3),color="white")),
    locations=list(gt::cells_column_labels(),gt::cells_column_spanners())
  ) %>%
  gt::tab_style(
    style=list(cell_text(size=11,weight="bold"),
               cell_borders(sides="bottom",weight=px(2))),
    locations=list(gt::cells_body(columns=everything(),rows=c(which(stk$table_body$row_type=="label"))))
  ) %>% gt::tab_style(
    style=list(cell_borders(sides="bottom",weight=px(1))),
    #locations=list(gt::cells_body(columns=everything(),rows=c(1,4,8,11,14,18,21,24,28,31,34,38)))
    locations=list(gt::cells_body(columns=everything(),rows=which(stk$table_body$row_type!="label")))) 

stk2 |> gtsave(file.path(outdir,"Parameter_Table","contact_rps_parms.html"))

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

#make keys (use for filtering in plot code)
keys=unique(preds$key)

# Plot results -----------------------------------------------------------------

#preds for dot/whisker
#pdiffs5 for charts showing differences
#trt.keys=keys[grep("area_aer_sex",keys)]
#trt.col=aer.hex
#preds_ylab_list=""
#xlabtitle=""

# Make function for Dot Whiskers -----------------------------------------------------------------

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


# Make barchart function -----------------------------------------------------------------
#keys=keys[grep("speed_aer_whole",keys)]
#trt.keys=keys[grep("ncon_aer_whole",keys)]
#trt.col=aer.hex
#preds_ylab_list=""
#xlabtitle=""
#xaxislabels=TRUE
#legend=TRUE
#facet_label=TRUE
#yaxislabels=TRUE
#pt_size=4
#l_size=1
#axtitlesize=14
#sexfilter="whole"
#sep_scales=TRUE
#MakeBarchart(keys[grep("ncon_aer_whole",keys)],aer.hex,"N. contacts","Removal Period",FALSE,FALSE,FALSE,TRUE,4,1,10,"whole",TRUE)


MakeBarchart=function(trt.keys,trt.col,preds_ylab_list,xlabtitle,xaxislabels,legend,facet_label,yaxislabels,pt_size,l_size,axtitlesize,sexfilter,sep_scales){
  if(missing(sexfilter)){
    sexfilter="whole"
  }
  
  if(missing(sep_scales)){
    sep_scales=FALSE
    resp=preds %>% filter(key==trt.keys) %>% select(response) %>% unique() %>% as.character
    #sex=preds2 %>% filter(key==trt.keys) %>% select(sm) %>% unique() %>% as.character
    ymax=preds %>% filter(response==resp,sex==sexfilter) %>% select(conf.high) %>% max()
    
  } else{
    resp=preds %>% filter(key==trt.keys) %>% select(response) %>% unique() %>% as.character
    remov=preds %>% filter(key==trt.keys) %>% select(rem) %>% unique() %>% as.character
    ymax=preds %>% filter(response==resp,sex==sexfilter,rem==remov) %>% select(conf.high) %>% max()
  }
  
  preds2=preds %>% filter(sex==sexfilter)


  if(xaxislabels==TRUE){
    p1=preds2 %>% filter(key==trt.keys) %>%
      ggplot(aes(x=per,y=predicted,group=trt,fill=trt))+
      #geom_segment(aes(y=conf.low),linewidth=l_size,position = position_dodge(width = 0.2))+
      geom_col(position = position_dodge(),
               width=0.5)+
      geom_errorbar(aes(ymin = conf.low, ymax = conf.high,group=trt),
                    position = position_dodge(width=0.5),
                    width=0.1)+
      theme_ipsum(axis_title_size=axtitlesize,
                  axis_text_size=8,
                  strip_text_size=axtitlesize,
                  base_size=15,
                  caption_margin=2,
                  plot_margin = margin(3, 3, 3, 3))+
      theme(legend.text=element_text(size=15))+
      scale_fill_manual(name="Removal\ntreatment",
                        values=c("#FFD065",trt.col))+
      ylab(preds_ylab_list)+
      xlab(xlabtitle)+
      facet_wrap(~sex)+
      ylim(0,ymax)
  }
  
  if(xaxislabels==FALSE){
    p1=preds2 %>% filter(key==trt.keys) %>%
      ggplot(aes(x=per,y=predicted,group=trt,fill=trt))+
      #geom_segment(aes(y=conf.low),linewidth=l_size,position = position_dodge(width = 0.2))+
      geom_col(position = position_dodge(),
               width=0.5)+
      geom_errorbar(aes(ymin = conf.low, ymax = conf.high,group=trt),
                    position = position_dodge(width=0.5),
                    width=0.1)+
      theme_ipsum(axis_title_size=axtitlesize,
                  axis_text_size=8,
                  strip_text_size =axtitlesize,
                  base_size=15,
                  caption_margin=2,
                  plot_margin = margin(3, 3, 3, 3))+
      theme(legend.text=element_text(size=15),
            axis.text.x = element_blank(),
            axis.title.x = element_blank())+
      scale_fill_manual(name="Removal\ntreatment",
                         values=c("#FFD065",trt.col))+
      ylab(preds_ylab_list)+
      xlab(xlabtitle)+
      facet_wrap(~sex)+
      ylim(0,ymax)
  }
  
  if(legend==FALSE){
    p1=p1 + theme(legend.position="none")
  }
  
  if(facet_label==FALSE){
    p1=p1 + theme(strip.text=element_blank())
  }
  
  if(yaxislabels==FALSE&sep_scales==FALSE){
    p1=p1+theme(axis.title.y=element_blank(),
                axis.text.y=element_blank())
  } else{
    if(yaxislabels==FALSE&sep_scales==TRUE){
      p1=p1+theme(axis.title.y=element_blank())
    }
  }
  
  return(p1)
}

#a_are_w=MakeBarchart(keys[grep("area_aer_whole",keys)],aer.hex,"area (km^2)","Removal Period",TRUE,FALSE,FALSE,TRUE,4,1,14)

#MakeBarchart(keys[grep("ncon_aer_whole",keys)],aer.hex,"N. contacts","Removal Period",FALSE,FALSE,FALSE,TRUE,4,1,10,"whole",TRUE)

# Set hex keys -------------------------
#Run function to make paired plots
tox.hex="#FF985A"
trap.hex="#FF57A9"
aer.hex="#822AFF"
ctrl.hex="#ffc400"

# Bar chart plots -----------------------------------------------------------------

## Movement Bar charts -----------------------------------------------------------------

# area
a_are_w=MakeBarchart(keys[grep("area_aer_whole",keys)],aer.hex,"area (km^2)","Removal Period",FALSE,FALSE,FALSE,TRUE,4,1,10)
t_are_w=MakeBarchart(keys[grep("area_trap_whole",keys)],trap.hex,"area (km^2)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,10)
x_are_w=MakeBarchart(keys[grep("area_tox_whole",keys)],tox.hex,"area (km^2)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,10)

# speed
a_sp_w=MakeBarchart(keys[grep("speed_aer_whole",keys)],aer.hex,"speed (km/hr)","Removal Period",FALSE,FALSE,FALSE,TRUE,4,1,10)
t_sp_w=MakeBarchart(keys[grep("speed_trap_whole",keys)],trap.hex,"speed (km/hr)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,10)
x_sp_w=MakeBarchart(keys[grep("speed_tox_whole",keys)],tox.hex,"speed (km/hr)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,10)

# NSD
a_nsd_w=MakeBarchart(keys[grep("nsd_aer_whole",keys)],aer.hex,"NSD (m^2)","Removal Period",FALSE,FALSE,FALSE,TRUE,4,1,10)
t_nsd_w=MakeBarchart(keys[grep("nsd_trap_whole",keys)],trap.hex,"NSD (m^2)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,10)
x_nsd_w=MakeBarchart(keys[grep("nsd_tox_whole",keys)],tox.hex,"NSD (m^2)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,10)

# distance
a_dis_w=MakeBarchart(keys[grep("distance_aer_whole",keys)],aer.hex,"distance (km)","Removal Period",TRUE,FALSE,FALSE,TRUE,4,1,10)
t_dis_w=MakeBarchart(keys[grep("distance_trap_whole",keys)],trap.hex,"distance (km)","Removal Period",TRUE,FALSE,FALSE,FALSE,4,1,10)
x_dis_w=MakeBarchart(keys[grep("distance_tox_whole",keys)],tox.hex,"distance (km)","Removal Period",TRUE,FALSE,FALSE,FALSE,4,1,10)

#combine removal types
nsd_pg=cowplot::plot_grid(a_nsd_w,t_nsd_w,x_nsd_w,nrow=1,rel_widths=c(0.65,0.65,0.65))
dis_pg=cowplot::plot_grid(a_dis_w,t_dis_w,x_dis_w,nrow=1,rel_widths=c(0.65,0.65,0.65))
sp_pg=cowplot::plot_grid(a_sp_w,t_sp_w,x_sp_w,nrow=1,rel_widths=c(0.65,0.65,0.65))
are_pg=cowplot::plot_grid(a_are_w,t_are_w,x_are_w,nrow=1,rel_widths=c(0.65,0.65,0.65))

#combine responses
allp=cowplot::plot_grid(are_pg,sp_pg,nsd_pg,dis_pg,ncol=1,scale=c(0.9,0.9,0.9,0.9))
ggsave("~/Downloads/barchart_mvmt.png",allp,width=5,height=5,units="in")

## Movement Bar charts (sex) -----------------------------------------------------------------

#Male
# area
a_are_m=MakeBarchart(keys[grep("area_aer_sex",keys)],aer.hex,"area (km^2)","Removal Period",FALSE,FALSE,FALSE,TRUE,4,1,10,"Male")
x_are_m=MakeBarchart(keys[grep("area_tox_sex",keys)],tox.hex,"area (km^2)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,10,"Male")

# speed
a_sp_m=MakeBarchart(keys[grep("speed_aer_sex",keys)],aer.hex,"speed (km/hr)","Removal Period",FALSE,FALSE,FALSE,TRUE,4,1,10,"Male")
t_sp_m=MakeBarchart(keys[grep("speed_trap_sex",keys)],trap.hex,"speed (km/hr)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,10,"Male")
x_sp_m=MakeBarchart(keys[grep("speed_tox_sex",keys)],tox.hex,"speed (km/hr)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,10,"Male")

# NSD
a_nsd_m=MakeBarchart(keys[grep("nsd_aer_sex",keys)],aer.hex,"NSD (m^2)","Removal Period",FALSE,FALSE,FALSE,TRUE,4,1,10,"Male")
t_nsd_m=MakeBarchart(keys[grep("nsd_trap_sex",keys)],trap.hex,"NSD (m^2)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,10,"Male")
x_nsd_m=MakeBarchart(keys[grep("nsd_tox_sex",keys)],tox.hex,"NSD (m^2)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,10,"Male")

# distance
a_dis_m=MakeBarchart(keys[grep("distance_aer_sex",keys)],aer.hex,"distance (km)","Removal Period",TRUE,FALSE,FALSE,TRUE,4,1,10,"Male")
t_dis_m=MakeBarchart(keys[grep("distance_trap_sex",keys)],trap.hex,"distance (km)","Removal Period",TRUE,FALSE,FALSE,FALSE,4,1,10,"Male")
x_dis_m=MakeBarchart(keys[grep("distance_tox_sex",keys)],tox.hex,"distance (km)","Removal Period",TRUE,FALSE,FALSE,FALSE,4,1,10,"Male")

#combine removal types, male
nsd_pgm=cowplot::plot_grid(a_nsd_m,t_nsd_m,x_nsd_m,nrow=1,rel_widths=c(0.65,0.65,0.65))
dis_pgm=cowplot::plot_grid(a_dis_m,t_dis_m,x_dis_m,nrow=1,rel_widths=c(0.65,0.65,0.65))
sp_pgm=cowplot::plot_grid(a_sp_m,t_sp_m,x_sp_m,nrow=1,rel_widths=c(0.65,0.65,0.65))
are_pgm=cowplot::plot_grid(a_are_m,NULL,x_are_m,nrow=1,rel_widths=c(0.65,0.65,0.65))

#Female
# area
a_are_f=MakeBarchart(keys[grep("area_aer_sex",keys)],aer.hex,"area (km^2)","Removal Period",FALSE,FALSE,FALSE,TRUE,4,1,10,"Female")
x_are_f=MakeBarchart(keys[grep("area_tox_sex",keys)],tox.hex,"area (km^2)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,10,"Female")

# speed
a_sp_f=MakeBarchart(keys[grep("speed_aer_sex",keys)],aer.hex,"speed (km/hr)","Removal Period",FALSE,FALSE,FALSE,TRUE,4,1,10,"Female")
t_sp_f=MakeBarchart(keys[grep("speed_trap_sex",keys)],trap.hex,"speed (km/hr)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,10,"Female")
x_sp_f=MakeBarchart(keys[grep("speed_tox_sex",keys)],tox.hex,"speed (km/hr)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,10,"Female")

# NSD
a_nsd_f=MakeBarchart(keys[grep("nsd_aer_sex",keys)],aer.hex,"NSD (m^2)","Removal Period",FALSE,FALSE,FALSE,TRUE,4,1,10,"Female")
t_nsd_f=MakeBarchart(keys[grep("nsd_trap_sex",keys)],trap.hex,"NSD (m^2)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,10,"Female")
x_nsd_f=MakeBarchart(keys[grep("nsd_tox_sex",keys)],tox.hex,"NSD (m^2)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,10,"Female")

# distance
a_dis_f=MakeBarchart(keys[grep("distance_aer_sex",keys)],aer.hex,"distance (km)","Removal Period",TRUE,FALSE,FALSE,TRUE,4,1,10,"Female")
t_dis_f=MakeBarchart(keys[grep("distance_trap_sex",keys)],trap.hex,"distance (km)","Removal Period",TRUE,FALSE,FALSE,FALSE,4,1,10,"Female")
x_dis_f=MakeBarchart(keys[grep("distance_tox_sex",keys)],tox.hex,"distance (km)","Removal Period",TRUE,FALSE,FALSE,FALSE,4,1,10,"Female")

#combine removal types, female
nsd_pgf=cowplot::plot_grid(a_nsd_f,t_nsd_f,x_nsd_f,nrow=1,rel_widths=c(0.65,0.65,0.65))
dis_pgf=cowplot::plot_grid(a_dis_f,t_dis_f,x_dis_f,nrow=1,rel_widths=c(0.65,0.65,0.65))
sp_pgf=cowplot::plot_grid(a_sp_f,t_sp_f,x_sp_f,nrow=1,rel_widths=c(0.65,0.65,0.65))
are_pgf=cowplot::plot_grid(a_are_f,NULL,x_are_f,nrow=1,rel_widths=c(0.65,0.65,0.65))

#combine responses
allpf=cowplot::plot_grid(are_pgf,sp_pgf,nsd_pgf,dis_pgf,ncol=1,scale=c(0.9,0.9,0.9,0.9))
#ggsave("~/Downloads/barchart_mvmt_female.png",allpf,width=4.5,height=5,units="in")

#combine sexes
mvmt_sex=cowplot::plot_grid(allpf,allpm,ncol=2,scale=c(0.92,0.92),labels=c("Female","Male"))
ggsave("~/Downloads/barchart_mvmt_sex.png",mvmt_sex,width=9,height=5,units="in")


## Contact Bar charts -----------------------------------------------------------------

#adjust preds for contact

#separate scales for contact ones

# area
a_nc_w=MakeBarchart(keys[grep("ncon_aer_whole",keys)],aer.hex,"N. contacts","Removal Period",FALSE,FALSE,FALSE,TRUE,4,1,10,"whole",TRUE)
t_nc_w=MakeBarchart(keys[grep("ncon_trap_whole",keys)],trap.hex,"N. contacts","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,10,"whole",TRUE)
x_nc_w=MakeBarchart(keys[grep("ncon_tox_whole",keys)],tox.hex,"N. contacts","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,10,"whole",TRUE)

# speed
a_ni_w=MakeBarchart(keys[grep("nind_aer_whole",keys)],aer.hex,"N. unique contact","Removal Period",TRUE,FALSE,FALSE,TRUE,4,1,10,"whole",TRUE)
t_ni_w=MakeBarchart(keys[grep("nind_trap_whole",keys)],trap.hex,"N. unique contact","Removal Period",TRUE,FALSE,FALSE,FALSE,4,1,10,"whole",TRUE)
x_ni_w=MakeBarchart(keys[grep("nind_tox_whole",keys)],tox.hex,"N. unique contact","Removal Period",TRUE,FALSE,FALSE,FALSE,4,1,10,"whole",TRUE)

#combine removal types
nc_pg=cowplot::plot_grid(a_nc_w,t_nc_w,x_nc_w,nrow=1,rel_widths=c(0.65,0.65,0.65))
ni_pg=cowplot::plot_grid(a_ni_w,t_ni_w,x_ni_w,nrow=1,rel_widths=c(0.65,0.65,0.65))

#combine responses
all_contact=cowplot::plot_grid(nc_pg,ni_pg,ncol=1,scale=c(0.9,0.9))
ggsave("~/Downloads/barchart_contact.png",all_contact,width=5,height=3,units="in")

## Contact Bar charts (sex) -----------------------------------------------------------------

#Male
# ncon
a_ncon_m=MakeBarchart(keys[grep("ncon_aer_sex",keys)],aer.hex,"N. contacts","Removal Period",FALSE,FALSE,FALSE,TRUE,4,1,10,"Male",TRUE)
t_ncon_m=MakeBarchart(keys[grep("ncon_trap_sex",keys)],trap.hex,"N. contacts","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,10,"Male",TRUE)
x_ncon_m=MakeBarchart(keys[grep("ncon_tox_sex",keys)],tox.hex,"N. contacts","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,10,"Male",TRUE)

# nind
a_nind_m=MakeBarchart(keys[grep("nind_aer_sex",keys)],aer.hex,"N. unique contacts","Removal Period",FALSE,FALSE,FALSE,TRUE,4,1,10,"Male",TRUE)
t_nind_m=MakeBarchart(keys[grep("nind_trap_sex",keys)],trap.hex,"N. unique contacts","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,10,"Male",TRUE)
x_nind_m=MakeBarchart(keys[grep("nind_tox_sex",keys)],tox.hex,"N. unique contacts","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,10,"Male",TRUE)

#combine removal types, male
ncon_pgm=cowplot::plot_grid(a_ncon_m,t_ncon_m,x_ncon_m,nrow=1,rel_widths=c(0.65,0.65,0.65))
nind_pgm=cowplot::plot_grid(a_nind_m,t_nind_m,x_nind_m,nrow=1,rel_widths=c(0.65,0.65,0.65))

allpm=cowplot::plot_grid(ncon_pgm,nind_pgm,ncol=1,scale=c(0.9,0.9))

#Female
# Ncon
a_ncon_f=MakeBarchart(keys[grep("ncon_aer_sex",keys)],aer.hex,"N. contacts","Removal Period",FALSE,FALSE,FALSE,TRUE,4,1,10,"Female",TRUE)
t_ncon_f=MakeBarchart(keys[grep("ncon_trap_sex",keys)],trap.hex,"N. contacts","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,10,"Female",TRUE)
x_ncon_f=MakeBarchart(keys[grep("ncon_tox_sex",keys)],tox.hex,"N. contacts","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,10,"Female",TRUE)

# Nind
a_nind_f=MakeBarchart(keys[grep("nind_aer_sex",keys)],aer.hex,"N. unique contacts","Removal Period",FALSE,FALSE,FALSE,TRUE,4,1,10,"Female",TRUE)
t_nind_f=MakeBarchart(keys[grep("nind_trap_sex",keys)],trap.hex,"N. unique contacts","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,10,"Female",TRUE)
x_nind_f=MakeBarchart(keys[grep("nind_tox_sex",keys)],tox.hex,"N. unique contacts","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,10,"Female",TRUE)

#combine removal types, female
ncon_pgf=cowplot::plot_grid(a_ncon_f,t_ncon_f,x_ncon_f,nrow=1,rel_widths=c(0.65,0.65,0.65))
nind_pgf=cowplot::plot_grid(a_nind_f,t_nind_f,x_nind_f,nrow=1,rel_widths=c(0.65,0.65,0.65))

#combine responses
allpf=cowplot::plot_grid(ncon_pgf,nind_pgf,ncol=1,scale=c(0.9,0.9))
#ggsave("~/Downloads/barchart_mvmt_female.png",allpf,width=4.5,height=5,units="in")

#combine sexes
con_sex=cowplot::plot_grid(allpf,allpm,ncol=2,scale=c(0.90,0.90),labels=c("Female","Male"))
ggsave("~/Downloads/barchart_contact_sex.png",con_sex,width=9,height=3,units="in")

# Dot Whisker plots -----------------------------------------------------------------

## Make movement dot whiskers -----------------------------------------------------------------
#trt.keys,trt.col,preds_ylab_list,xlabtitle,xaxislabels,legend,facet_label,yaxislabels,pt_size,l_size,axtitlesize

# NSD
a_nsd_w=MakeDotWhisker(keys[grep("nsd_aer_whole",keys)],aer.hex,"NSD (m^2)","Removal Period",FALSE,FALSE,FALSE,TRUE,4,1,14)
t_nsd_w=MakeDotWhisker(keys[grep("nsd_trap_whole",keys)],trap.hex,"NSD (m^2)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,14)
x_nsd_w=MakeDotWhisker(keys[grep("nsd_tox_whole",keys)],tox.hex,"NSD (m^2)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,14)

# distance
a_dis_w=MakeDotWhisker(keys[grep("distance_aer_whole",keys)],aer.hex,"distance (km)","Removal Period",FALSE,FALSE,FALSE,TRUE,4,1,14)
t_dis_w=MakeDotWhisker(keys[grep("distance_trap_whole",keys)],trap.hex,"distance (km)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,14)
x_dis_w=MakeDotWhisker(keys[grep("distance_tox_whole",keys)],tox.hex,"distance (km)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,14)

# speed
a_sp_w=MakeDotWhisker(keys[grep("speed_aer_whole",keys)],aer.hex,"speed (km/hr)","Removal Period",FALSE,FALSE,FALSE,TRUE,4,1,14)
t_sp_w=MakeDotWhisker(keys[grep("speed_trap_whole",keys)],trap.hex,"speed (km/hr)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,14)
x_sp_w=MakeDotWhisker(keys[grep("speed_tox_whole",keys)],tox.hex,"speed (km/hr)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,14)

# area
a_are_w=MakeDotWhisker(keys[grep("area_aer_whole",keys)],aer.hex,"area (km^2)","Removal Period",TRUE,FALSE,FALSE,TRUE,4,1,14)
t_are_w=MakeDotWhisker(keys[grep("area_trap_whole",keys)],trap.hex,"area (km^2)","Removal Period",TRUE,FALSE,FALSE,FALSE,4,1,14)
x_are_w=MakeDotWhisker(keys[grep("area_tox_whole",keys)],tox.hex,"area (km^2)","Removal Period",TRUE,FALSE,FALSE,FALSE,4,1,14)

#combine removal types
nsd_pg=cowplot::plot_grid(a_nsd_s,t_nsd_s,x_nsd_s,nrow=1,rel_widths=c(0.65,0.65,0.65))
dis_pg=cowplot::plot_grid(a_dis_s,t_dis_s,x_dis_s,nrow=1,rel_widths=c(0.65,0.65,0.65))
sp_pg=cowplot::plot_grid(a_sp_s,t_sp_s,x_sp_s,nrow=1,rel_widths=c(0.65,0.65,0.65))
are_pg=cowplot::plot_grid(a_are_s,t_are_s,x_are_s,nrow=1,rel_widths=c(0.65,0.65,0.65))

#combine responses
allp=cowplot::plot_grid(nsd_pg,dis_pg,sp_pg,are_pg,ncol=1)
ggsave("~/Downloads/allp.png",allp,width=9,height=6.5,units="in")

#put asterisk on significant interactions
#doing manually in preview
  #intxn=intxn[,c(8,4,5,6,7,1,2,3,9)]
  intxn[intxn$sex=="whole"&intxn$alpha==1,]

## Make movement dot whiskers  (sex) ------------------------------------------------------------

  # NSD
  #MakeDotWhisker=function(trt.keys,trt.col,preds_ylab_list,xlabtitle,xaxislabels,legend,facet_label,yaxislabels,pt_size,l_size){
  #trt.keys,trt.col,preds_ylab_list,xlabtitle,xaxislabels,legend,facet_label,yaxislabels,sex_model,pt_size,l_size
  a_nsd_w_s=MakeDotWhisker(keys[grep("nsd_aer_sex",keys)],aer.hex,"NSD (m^2)","Removal Period", FALSE,FALSE,TRUE,TRUE,4,1,14)
  t_nsd_w_s=MakeDotWhisker(keys[grep("nsd_trap_sex",keys)],trap.hex,"NSD (m^2)","Removal Period",FALSE,FALSE,TRUE,FALSE,4,1,14)
  x_nsd_w_s=MakeDotWhisker(keys[grep("nsd_tox_sex",keys)],tox.hex,"NSD (m^2)","Removal Period",FALSE,FALSE,TRUE,FALSE,4,1,14)
  
  # distance
  a_dis_w_s=MakeDotWhisker(keys[grep("distance_aer_sex",keys)],aer.hex,"distance (km)","Removal Period",FALSE,FALSE,FALSE,TRUE,4,1,14)
  t_dis_w_s=MakeDotWhisker(keys[grep("distance_trap_sex",keys)],trap.hex,"distance (km)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,14)
  x_dis_w_s=MakeDotWhisker(keys[grep("distance_tox_sex",keys)],tox.hex,"distance (km)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,14)
  
  # speed
  a_sp_w_s=MakeDotWhisker(keys[grep("speed_aer_sex",keys)],aer.hex,"speed (km/hr)","Removal Period",FALSE,FALSE,FALSE,TRUE,4,1,14)
  t_sp_w_s=MakeDotWhisker(keys[grep("speed_trap_sex",keys)],trap.hex,"speed (km/hr)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,14)
  x_sp_w_s=MakeDotWhisker(keys[grep("speed_tox_sex",keys)],tox.hex,"speed (km/hr)","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,14)
    
  # area
  a_are_w_s=MakeDotWhisker(keys[grep("area_aer_sex",keys)],aer.hex,"area (km^2)","Removal Period",TRUE,FALSE,FALSE,TRUE,4,1,14)
  x_are_w_s=MakeDotWhisker(keys[grep("area_tox_sex",keys)],tox.hex,"area (km^2)","Removal Period",TRUE,FALSE,FALSE,FALSE,4,1,14)
  
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

  colnames(pdiffs5)[2]="trt"
  pdiffs5$sex<-tolower(pdiffs5$sex)
  pdiff_int=left_join(pdiffs5,intxn,by=c("sex","trt","period","response"))
  pdiff_int=pdiff_int[,c(1:8,13)]
  
  #pdiff_int
  pdiff_int=pdiff_int[complete.cases(pdiff_int),]
  pdiff_int=pdiff_int[pdiff_int$alpha==1,]
  
## Make contact dot whiskers  ------------------------------------------------------------
  #trt.keys,trt.col,preds_ylab_list,xlabtitle
  #xaxislabels,legend,facet_label,yaxislabels,
  #pt_size,l_size,axtitlesize

  # num contacts
  a_nc_w=MakeDotWhisker(keys[grep("ncon_aer_whole",keys)],aer.hex,"Num. total contacts","Removal Period",FALSE,FALSE,FALSE,TRUE,4,1,14)
  t_nc_w=MakeDotWhisker(keys[grep("ncon_trap_whole",keys)],trap.hex,"Num. total contacts","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,14)
  x_nc_w=MakeDotWhisker(keys[grep("ncon_tox_whole",keys)],tox.hex,"Num. total contacts","Removal Period",FALSE,FALSE,FALSE,FALSE,4,1,14)
  
  #num uniq contacts
  a_ni_w=MakeDotWhisker(keys[grep("nind_aer_whole",keys)],aer.hex,"Num. unique contacts","Removal Period",TRUE,FALSE,FALSE,TRUE,4,1,14)
  t_ni_w=MakeDotWhisker(keys[grep("nind_trap_whole",keys)],trap.hex,"Num. unique contacts","Removal Period",TRUE,FALSE,FALSE,FALSE,4,1,14)
  x_ni_w=MakeDotWhisker(keys[grep("nind_tox_whole",keys)],tox.hex,"Num. unique contacts","Removal Period",TRUE,FALSE,FALSE,FALSE,4,1,14)

  #combine response plots
  nc_pg=cowplot::plot_grid(a_nc_w,t_nc_w,x_nc_w,nrow=1,rel_widths=c(0.65,0.65,0.65))
  ni_pg=cowplot::plot_grid(a_ni_w,t_ni_w,x_ni_w,nrow=1,rel_widths=c(0.65,0.65,0.65))
  
  #combine responses
  allp_c=cowplot::plot_grid(nc_pg,ni_pg,ncol=1,rel_heights=c(1,1))

  ggsave("~/Downloads/allcontact_whole.png",allp_c,width=9,height=3.25,units="in")
  
  #trt.keys,trt.col,preds_ylab_list,xlabtitle
  #xaxislabels,legend,facet_label,yaxislabels,
  #pt_size,l_size,axtitlesize
  
## Make contact dot whiskers (sex) ------------------------------------------------------------
  
  # num contacts
  a_nc_w_s=MakeDotWhisker(keys[grep("ncon_aer_sex",keys)],aer.hex,"Num. total contacts","Removal Period",FALSE,FALSE,TRUE,TRUE,4,1,14)
  t_nc_w_s=MakeDotWhisker(keys[grep("ncon_trap_sex",keys)],trap.hex,"Num. total contacts","Removal Period",FALSE,FALSE,TRUE,FALSE,4,1,14)
  x_nc_w_s=MakeDotWhisker(keys[grep("ncon_tox_sex",keys)],tox.hex,"Num. total contacts","Removal Period",FALSE,FALSE,TRUE,FALSE,4,1,14)
  
  #num uniq contacts
  a_ni_w_s=MakeDotWhisker(keys[grep("nind_aer_sex",keys)],aer.hex,"Num. unique contacts","Removal Period",TRUE,FALSE,FALSE,TRUE,4,1,14)
  t_ni_w_s=MakeDotWhisker(keys[grep("nind_trap_sex",keys)],trap.hex,"Num. unique contacts","Removal Period",TRUE,FALSE,FALSE,FALSE,4,1,14)
  x_ni_w_s=MakeDotWhisker(keys[grep("nind_tox_sex",keys)],tox.hex,"Num. unique contacts","Removal Period",TRUE,FALSE,FALSE,FALSE,4,1,14)
  
  #combine response plots
  nc_pg_s=cowplot::plot_grid(a_nc_w_s,t_nc_w_s,x_nc_w_s,nrow=1,rel_widths=c(0.65,0.65,0.65))
  ni_pg_s=cowplot::plot_grid(a_ni_w_s,t_ni_w_s,x_ni_w_s,nrow=1,rel_widths=c(0.65,0.65,0.65))
  
  #combine responses
  allp_c=cowplot::plot_grid(nc_pg_s,ni_pg_s,ncol=1,rel_heights=c(1,1))
  
  ggsave("~/Downloads/allcontact_sex.png",allp_c,width=9,height=3.25,units="in")
  

  
# Make heatmaps of predictions -----------------------------------
  
## Make coded pdiffs table -----------------
  
  #calc treatment mag diffs
  pdiffs3$during_trt_ctrl<-pdiffs3$before_during_trt/pdiffs3$before_during_ctrl
  pdiffs3$after_trt_ctrl<-pdiffs3$before_after_trt/pdiffs3$before_after_ctrl
  
  
  #trim unneeded cols, don't need diffs within control/treatment
  ptbl=pdiffs3[,c("sex",
                     "rem",
                     "response",
                     "before_during_ctrl",
                     "before_during_trt",
                     "before_after_ctrl",
                     "before_after_trt",
                     "during_trt_ctrl",
                     "after_trt_ctrl")]
  
  #pivot_logner
  ptbl1=tidyr::pivot_longer(ptbl,cols=4:9,names_to="per",values_to="pred")
  
  ptbl1$period=NA
  ptbl1$period[grep("during",ptbl1$per)]<-"during"
  ptbl1$period[grep("after",ptbl1$per)]<-"after"
  ptbl1$period<-forcats::fct_relevel(ptbl1$period,c("during","after"))
  
  ptbl1$trt=NA
  ptbl1$trt[grep("after_ctrl",ptbl1$per)]<-"ctrl"
  ptbl1$trt[grep("during_ctrl",ptbl1$per)]<-"ctrl"
  ptbl1$trt[grep("after_trt",ptbl1$per)]<-"trt"
  ptbl1$trt[grep("during_trt",ptbl1$per)]<-"trt"
  ptbl1$trt[grep("trt_ctrl",ptbl1$per)]<-"fulldiff"
  
  #remove NA rows (aerial, no during period)
  ptbl1=ptbl1[!is.na(ptbl1$pred),]
  
  #replace NA in sex with whole model
  ptbl1$sex[is.na(ptbl1$sex)]<-"whole"
  
  #make key to match to pdiffs later
  ptbl1$sex=tolower(ptbl1$sex)
  ptbl1$key=paste(ptbl1$response,ptbl1$rem,ptbl1$sex,ptbl1$period,sep="_")
  
  #get just the fulldiff ones
  ptbl1=ptbl1[ptbl1$trt=="fulldiff",]
  
  #rename pred to mag
  ptbl1$mag=ptbl1$pred
  
  #make just key and mag
  ptbl1=ptbl1[,c(8,9)]
  
  #intxn2=intxn
  str_sub(intxn[intxn$response=="nsd",]$model[grep("rps",intxn[intxn$response=="nsd",]$model)],4L,8L)<-"_nsd_rps_"
  str_sub(intxn[intxn$response=="area",]$model[grep("rps",intxn[intxn$response=="area",]$model)],4L,8L)<-"_area_rps_"
  
  str_sub(intxn[intxn$response=="nsd",]$model[grep("res.rp_",intxn[intxn$response=="nsd",]$model)],4L,7L)<-"_nsd_rp_"
  str_sub(intxn[intxn$response=="area",]$model[grep("res.rp_",intxn[intxn$response=="area",]$model)],4L,7L)<-"_area_rp_"
  
  #make key in intxn to match to pdiffs5
  intxn$key=paste(intxn$response,intxn$trt,intxn$sex,intxn$period,sep="_")
  #str_sub(intxn[grep("female",intxn$key),]$key,-6L)<-"whole"
  #str_sub(intxn[grep("male",intxn$key),]$key,-4L)<-"whole"
  
  #remake pdiffs5 key
  pdiffs5$sex=tolower(pdiffs5$sex)
  pdiffs5$key=paste(pdiffs5$response,pdiffs5$rem,pdiffs5$sex,pdiffs5$period,sep="_")
  
  
  intxn2=intxn[,c(9,10)] #trim to alpha and key only
  pd=left_join(pdiffs5,intxn2,by="key")
  
  #significant intxns for male ref group models:
  pdm=pd[pd$sex=="male",]
  pdm[pdm$key=="nsd_trap_male_during"|pdm$key=="distance_trap_male_during",]$alpha=1
  pdm[pdm$key!="nsd_trap_male_during"&pdm$key!="distance_trap_male_during",]$alpha=0.1
  pdm[pdm$key=="ncon_trap_male_during"|
        pdm$key=="ncon_trap_male_after"|
        pdm$key=="nind_trap_male_after"|
        pdm$key=="ncon_tox_male_after",]$alpha=1
  
  pdf=pd[pd$sex=="female",]
  pdw=pd[pd$sex=="whole",]
  pda=rbind(pdf,pdm)
  pda=rbind(pda,pdw)
  pda=pda[complete.cases(pda),]

  #join to get magnitude diffs
  pdj=left_join(pda,ptbl1,by="key")
  
  pdj$mag<-abs(pdj$mag)
  
  #change magnitude sign to be related to direction of effect, movement
  pdj[pdj$pred<0,]$mag=pdj[pdj$pred<0,]$mag*(-1)
  
  pdj$response[pdj$response=="distance"]<-"distance (km)"
  pdj$response[pdj$response=="speed"]<-"speed (km/hr)"
  pdj$response[pdj$response=="area"]<-"area (km^2)"
  pdj$response[pdj$response=="nsd"]<-"nsd (m^2)"
  
  #change 0.1 alpha to 0, gets confusing otherwise
  pdj[pdj$alpha==0.1,]$alpha<-0
  
  #clean digits
  pdj$pred<-round(pdj$pred,digits=2)
  
  #relevel
  pdj$rem<-factor(pdj$rem,levels=c("aer","trap","tox"))
  
  #remove contact ones, plot separately
  pdjmv=pdj[pdj$response!="nind"&pdj$response!="ncon",]
  pdjco=pdj[pdj$response=="nind"|pdj$response=="ncon",]
  
  #relevel contact responses
  pdjco$response<-factor(pdjco$response,levels=c("ncon","nind"))
  
  
  pdjmv_f=pdjmv[pdjmv$sex=="female",]
  pdjmv_m=pdjmv[pdjmv$sex=="male",]
  pdjmv_w=pdjmv[pdjmv$sex=="whole",]
  
  pdjco_f=pdjco[pdjco$sex=="female",]
  pdjco_m=pdjco[pdjco$sex=="male",]
  pdjco_w=pdjco[pdjco$sex=="whole",]
  
  
## Plot individual heatmaps -----------------------------------

  #make max mag 20+ to keep scale of others observable
  pdjmv_w$mag[pdjmv_w$mag==max(pdjmv_w$mag)]<-20 
  pdjmv_w$mag[pdjmv_w$mag==max(pdjmv_w$mag)]<-20 
  
### Female movement heatmap ------------------------  
  heatmap_mvmt_f <- 
    ggplot(pdjmv_f,aes(x = period, y = response, fill = mag, label = formatC(pred, format = "e",digits=2))) +
    geom_tile(aes(alpha=alpha),color = "white",show.legend=TRUE) + # Create heatmap
    geom_text(color = "black", size = 2.5) + # Add text labels
    scale_alpha(range=c(0,1))+
    scale_fill_continuous_divergingx(palette = 'RdBu', mid = 1,limits=c(-20,20),h1=6,c1=10,c2=10,c3=10) + 
    theme_ipsum() + # Set theme
    labs(x = "Removal Period", y = "Movement Response",fill="mag. diff.") + # Labels
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    facet_wrap(~rem)+
    theme(legend.position="none",
          axis.text.x=element_blank(),
          axis.title.x=element_blank())

### Male movement heatmap ------------------------  
  
  heatmap_mvmt_m <- 
    ggplot(pdjmv_m,aes(x = period, y = response, fill = mag, label = formatC(pred, format = "e",digits=2))) +
    geom_tile(aes(alpha=alpha),color = "white",show.legend=TRUE) + # Create heatmap
    geom_text(color = "black", size = 2.5) + # Add text labels
    scale_alpha(range=c(0,1),guide="none")+
    scale_fill_continuous_divergingx(palette = 'RdBu', mid = 1,limits=c(-20,20),h1=6,c1=10,c2=10,c3=10) + 
    theme_ipsum() + # Set theme
    labs(x = "Removal Period", y = "Movement Response",fill="mag. diff.") + # Labels
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    facet_wrap(~rem)+
    theme(legend.position="bottom")

### Whole movement heatmap ------------------------  

  heatmap_mvmt_w <- 
    ggplot(pdjmv_w,aes(x = period, y = response, fill = mag, label = formatC(pred, format = "e",digits=2))) +
    geom_tile(aes(alpha=alpha),color = "white",show.legend=TRUE) + # Create heatmap
    geom_text(color = "black", size = 2.5) + # Add text labels
    scale_alpha(range=c(0,1))+
    scale_fill_continuous_divergingx(palette = 'RdBu', mid = 1,limits=c(-20,20),h1=6,c1=10,c2=10,c3=10) + 
    theme_ipsum() + # Set theme
    labs(x = "Removal Period", y = "Movement Response",fill="mag. diff.") + # Labels
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    facet_wrap(~rem)+
    theme(legend.position="none",
          axis.text.x=element_blank(),
          axis.title.x=element_blank())
  
### Female contact heatmap ------------------------  
  
  heatmap_con_f <- 
    ggplot(pdjco_f,aes(x = period, y = response, fill = mag, label = formatC(pred, format = "e",digits=2))) +
    geom_tile(aes(alpha=alpha),color = "white",show.legend=TRUE) + # Create heatmap
    geom_text(color = "black", size = 2.5) + # Add text labels
    scale_alpha(range=c(0,1))+
    scale_fill_continuous_divergingx(palette = 'RdBu', mid = 1,limits=c(-20,20),h1=6,c1=10,c2=10,c3=10) + 
    theme_ipsum() + # Set theme
    labs(x = "Removal Period", y = "Movement Response",fill="mag. diff.") + # Labels
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    facet_wrap(~rem)+
    theme(legend.position="none",
          axis.text.x=element_blank(),
          axis.title.x=element_blank())
  
  
### Male contact heatmap ------------------------
  pdjco_m$mag[pdjco_m$mag==min(pdjco_m$mag)]=-20
  pdjco_m$mag[pdjco_m$mag==min(pdjco_m$mag)]=-20
  
  heatmap_con_m <- 
    ggplot(pdjco_m,aes(x = period, y = response, fill = mag, label = formatC(pred, format = "e",digits=2))) +
    geom_tile(aes(alpha=alpha),color = "white",show.legend=TRUE) + # Create heatmap
    geom_text(color = "black", size = 2.5) + # Add text labels
    scale_alpha(range=c(0,1),guide="none")+ #none significant here
    scale_fill_continuous_divergingx(palette = 'RdBu', mid = 1,limits=c(-20,20),h1=6,c1=10,c2=10,c3=10) + 
    theme_ipsum() + # Set theme
    labs(x = "Removal Period", y = "Movement Response",fill="mag. diff.") + # Labels
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    facet_wrap(~rem)+
    theme(legend.position="bottom")
  
  
### Whole contact heatmap ------------------------  
  
  heatmap_con_w <- 
    ggplot(pdjco_w,aes(x = period, y = response, fill = mag, label = formatC(pred, format = "e",digits=2))) +
    geom_tile(aes(alpha=alpha),color = "white",show.legend=TRUE) + # Create heatmap
    geom_text(color = "black", size = 2.5) + # Add text labels
    scale_alpha(range=c(0,1))+
    scale_fill_continuous_divergingx(palette = 'RdBu', mid = 1,limits=c(-20,20),h1=6,c1=10,c2=10,c3=10) + 
    theme_ipsum() + # Set theme
    labs(x = "Removal Period", y = "Movement Response",fill="mag. diff.") + # Labels
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    facet_wrap(~rem)+
    theme(legend.position="none",
          axis.text.x=element_blank(),
          axis.title.x=element_blank())
  
  
## Group heatmaps with cowplot and save to file -----------------------------------
  
### Movement heatmaps ------------------------------------  
  heatmaps_mvmt=cowplot::plot_grid(heatmap_mvmt_w,NULL,heatmap_mvmt_f,NULL,heatmap_mvmt_m, ncol=1, rel_heights=c(0.5,-0.1,0.5,-0.1,0.75))
  ggsave(file.path(outdir,"HeatMaps","heatmap_mvmt.png"),heatmaps_mvmt,width=6.5,height=9,units="in")


### contact heatmaps ------------------------------------  
  
  heatmaps_con=cowplot::plot_grid(heatmap_con_w,NULL,heatmap_con_f,NULL,heatmap_con_m, ncol=1, rel_heights=c(0.5,-0.1,0.5,-0.1,0.75))
  ggsave(file.path(outdir,"HeatMaps","heatmap_contact.png"),heatmaps_con,width=6.5,height=9,units="in")
  
  
  
  
  