#The purpose of this script is to do some diagnostics on the akde areas
#Some of the areas output had very large est bias, 
#with correspondingly large CI values around the estimates,
#and unnaturally large akde area estimates

#Based on recommendations from ctmm developer, going to do the following
#1. Use a cutoff of estbias (3%) by looking at forest plots/effective Ns
#2. For those with estbias above cutoff, look at paths/variograms to determine if there are large displacements affecting results
#3. For those tracks, segment paths using seg2clust, and run akde separately on each cluster
#4.Do the following: 
  #a. Remove clusters within track that have >3% est bias (implies displacement behavior, but use tracks to visualize this as well)
  #b. Keep clusters with <3% bias
  #c. Average remaining clusters to get hr area size
#5. Then, merge pigs with averaged clustered paths back to original akde data frame output


#This script is broken up into the following sections:
#1. Script setup and data formatting
#2. Paths/Variogram/Forest plot diagnostic outputs
#3. Do segmentation and check akde output of segments to make decisions
#4. Loop through pigs and use AKDE function to get hr areas

###########################################
#### Script setup and data formatting #####
###########################################
#git clone --progress --verbose . /Volumes/Projects/MUDD/ASF_NIFA/Pipelines/Removals_Movement
{
#Load libraries
library(ctmm)
library(segclust2d)
#library(NIFApackagev1.1)
library(stringr)

#Load and format geo data
home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/1_Data/"
geo=read.csv(paste0(home,"Objects/geo_remtyp.csv"))
geo=geo[,-1]
geo=geo[!is.na(geo$Removal.Type),]
geo$date_only<-Neat.Dates.POSIXct(geo$date_only,tz="UTC")

#Read input objects
input="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/1_Data/Input/"
activities.raw=readRDS(paste0(input,"activities.rds"))

#source needed functions
funcdir<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/2_Scripts/Functions"
func.list=list.files(funcdir,full.names=TRUE)
for(f in 1:length(func.list)){
  source(func.list[f])
}

#Load df outputs from Get_Period_AKDEs
akaer=readRDS("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/outdf_akde_aerial.rds")
aktrap=readRDS("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/outdf_akde_trap.rds")
aktox=readRDS("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/outdf_akde_tox.rds")

#Fix area units-- sometimes area rec'd as hectares, meters
Fix.Area.Units<-function(akrem){
  akrem[grep("square meters",akrem$akde.hrarea.units),][3:5]=akrem[grep("square meters",akrem$akde.hrarea.units),][3:5]/1e6
  akrem[grep("hectares",akrem$akde.hrarea.units),][3:5]=akrem[grep("hectares",akrem$akde.hrarea.units),][3:5]/100
  return(akrem)
}
aktox=Fix.Area.Units(aktox)
aktrap=Fix.Area.Units(aktrap)
akaer=Fix.Area.Units(akaer)

#Get est bias statistic. Want to keep below 2%
#originally made cutoff 1%, but this seemed too stringent
#with 1%, had a lot of pigs that were truly range resident but with larger hr areas
akaer$estbias=(1/(akaer$Neff)^2)*100
aktrap$estbias=(1/(aktrap$Neff)^2)*100
aktox$estbias=(1/(aktox$Neff)^2)*100


#Get removal dates
{
  #trap dates
  trap.act=activities.raw[activities.raw$activity=="trap",]
  trap.act$trap_datetime=Neat.Dates.POSIXct(trap.act$trap_datetime,tz="UTC")
  trap.start.date=as.Date(min(trap.act$trap_datetime))
  trap.end.date=as.Date(max(trap.act$trap_datetime))
  
  #tox dates
  tox.act=activities.raw[activities.raw$activity=="toxic",]
  tox.act$tox_datetime=Neat.Dates.POSIXct(tox.act$tox_datetime,tz="UTC")
  tox.start.date=as.Date(min(tox.act$tox_datetime))
  tox.end.date=as.Date("2023-03-09")
  
  #aerial dates
  aer.start.date=as.Date("2023-03-27")
  aer.end.date=as.Date("2023-03-29")
  
}

#Get cutoffs for each removal type
{
#Aerial 
#determined 4 weeks from sensitivity analysis
#geo.aer=geo2[geo2$Removal.Type=="aer"|geo2$Removal.Type=="ctrl",]
aer.cutoff1=aer.start.date-37
aer.cutoff2=aer.end.date+37

#Trap
#use length of time in trap after period (smaller than during period)
trap.len=as.Date(max(geo[geo$Removal.Type=="trap",]$date_only))-as.Date(trap.end.date)-1
#use as trap start date, 56 days before trap end date
#think this is when trapping really started more intensively, anyways
trap.start.date.akde=trap.end.date-trap.len-1
trap.cutoff1=trap.start.date-trap.len-1
trap.cutoff2=trap.end.date+trap.len+1

#Tox
#use length of time in tox period
#tox.len=tox.end.date-tox.start.date
tox.len=36
tox.cutoff1=tox.start.date-tox.len-1
tox.cutoff2=tox.end.date+tox.len+1

#removal="tox"
#cutoff1=tox.cutoff1
#cutoff2=tox.cutoff2
#startdate=tox.start.date
#enddate=tox.end.date
Do.Date.Cutoffs<-function(geo,removal,cutoff1,cutoff2,startdate,enddate){
  geo.rem=geo[geo$Removal.Type==removal|geo$Removal.Type=="ctrl",]
  geo.rem.before=geo.rem[geo.rem$date_only>=cutoff1&geo.rem$date_only<startdate,]
  geo.rem.after=geo.rem[geo.rem$date_only>enddate&geo.rem$date_only<=cutoff2,]
  geo.rem.during=geo.rem[geo.rem$date_only>=startdate&geo.rem$date_only<=enddate,]
  geo.rem.before$removal.period.akdecalc="before"
  geo.rem.after$removal.period.akdecalc="after"
  geo.rem.during$removal.period.akdecalc="during"
  geo.rem=rbind(geo.rem.before,geo.rem.during,geo.rem.after)
  return(geo.rem)
}

geo.aer=Do.Date.Cutoffs(geo,"aer",aer.cutoff1,aer.cutoff2,aer.start.date,aer.end.date)
#geo.trap=Do.Date.Cutoffs(geo,"trap",trap.cutoff1,trap.cutoff2,trap.start.date.akde,trap.end.date)
geo.tox=Do.Date.Cutoffs(geo,"tox",tox.cutoff1,tox.cutoff2,tox.start.date,tox.end.date)

#do trap ones manually
geo.rem=geo[geo$Removal.Type=="trap"|geo$Removal.Type=="ctrl",]
geo.rem.before=geo.rem[geo.rem$date_only>=trap.cutoff1&geo.rem$date_only<trap.start.date,]
geo.rem.after=geo.rem[geo.rem$date_only>trap.end.date&geo.rem$date_only<=trap.cutoff2,]
geo.rem.during=geo.rem[geo.rem$date_only>trap.start.date.akde&geo.rem$date_only<=trap.end.date,]
geo.rem.before$removal.period.akdecalc="before"
geo.rem.after$removal.period.akdecalc="after"
geo.rem.during$removal.period.akdecalc="during"
geo.rem=rbind(geo.rem.before,geo.rem.during,geo.rem.after)
geo.trap=geo.rem

#verification
geo.tox.sums=
  geo.tox %>% 
  dplyr::group_by(animalid,removal.period.akdecalc) %>% 
  dplyr::summarise(num.locs=n(),
                   data_from=first(data_from),
                   strt.date=min(date_only),
                   end.date=max(date_only))
geo.tox.sums$difftime=as.Date(geo.tox.sums$end.date)-as.Date(geo.tox.sums$strt.date)

geo.trap.sums=
  geo.trap %>% 
  dplyr::group_by(animalid,removal.period.akdecalc) %>% 
  dplyr::summarise(num.locs=n(),
                   data_from=first(data_from),
                   strt.date=min(date_only),
                   end.date=max(date_only))
geo.trap.sums$difftime=as.Date(geo.trap.sums$end.date)-as.Date(geo.trap.sums$strt.date)

geo.aer.sums=
  geo.aer %>% 
  dplyr::group_by(animalid,removal.period.akdecalc) %>% 
  dplyr::summarise(num.locs=n(),
                   data_from=first(data_from),
                   strt.date=min(date_only),
                   end.date=max(date_only))
geo.aer.sums$difftime=as.Date(geo.aer.sums$end.date)-as.Date(geo.aer.sums$strt.date)

######Trim based on above summaries
#remove any pigs that were tracked less than half the number of days of other pigs in any period
datestooshort.trap=geo.trap.sums[geo.trap.sums$difftime<(trap.len/2),]$animalid
geo.trap=geo.trap[!geo.trap$animalid%in%datestooshort.trap,]

datestooshort.tox=c("86070_H2_H2")
geo.tox=geo.tox[!geo.tox$animalid%in%datestooshort.tox,]

}


###### set needed functions
Convert.Telemetry<-function(geolocs){
  id.n=which(colnames(geolocs)=="animalid")
  id.t=which(colnames(geolocs)=="datetime")
  id.lon=which(colnames(geolocs)=="latitude")
  id.lat=which(colnames(geolocs)=="longitude")
  tk=geolocs[,c(id.n,id.t,id.lon,id.lat)]
  colnames(tk)<-c("ID","timestamp","longitude","latitude")
  crs_str="+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  te1 <- ctmm::as.telemetry(tk,projection=crs_str)
  return(te1)
}

}

###########################################
#### Paths/Variogram Diagnostic plots #####
###########################################
{
#Make function to plot each pig/period: paths with ggplot and variogram
#geo.rem=geo.aer
#Paths.Variogram.Diagnostics(geo.rem,"aer")
#outdir="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Exploratory_Outputs/Paths_Variograms/"
#eo.rem=geo.aerb
Paths.Variogram.Diagnostics<-function(geo.rem,remtype,outdir){
  IDs=unique(geo.rem$animalid)
  #periods=unique(geo.rem$removal.period.akdecalc)
  for(i in 1:length(IDs)){
    geo.pig=geo.rem[geo.rem$animalid==IDs[i],]
    pers=unique(geo.pig$removal.period.akdecalc)
    if(remtype=="aer"&"during"%in%pers){pers=pers[-which(pers=="during")]}
    for(p in 1:length(pers)){
    geo.pig.per=geo.pig[geo.pig$removal.period.akdecalc==pers[p],]
    paths=ggplot(geo.pig.per,aes(x=X,y=Y))+
      geom_path()
    te1=Convert.Telemetry(geo.pig.per)
    SVF <- variogram(te1,CI="Gauss")
    if(!dir.exists(paste0(outdir,remtype))){dir.create(paste0(outdir,remtype))}
    ggsave(paste0(outdir,remtype,"/",IDs[i],"_",pers[p],"_ggpaths",".png"),plot=paths)
    png(paste0(outdir,remtype,"/",IDs[i],"_",pers[p],"_var",".png"), width=600, height=500, res=120) # start export
    plot(SVF)
    dev.off()
    }
    
  }
  
  
  #return(geo.rem)
}


#outdir.pv="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Exploratory_Outputs/Paths_Variograms/"
#Paths.Variogram.Diagnostics(geo.aer,"aer",outdir.pv)
#Paths.Variogram.Diagnostics(geo.trap,"trap",outdir.pv)
#Paths.Variogram.Diagnostics(geo.tox,"tox",outdir.pv)

#Gex expected order of bias-- cutoff of 1% for where we will try to improve akde area estimates
aktrap$expbias=(1/aktrap$Neff^2)*100 #convert to percentage
akaer$expbias=(1/akaer$Neff^2)*100
aktox$expbias=(1/aktox$Neff^2)*100

#Look at those with large order of bias-- >2%
  #Look at variograms and paths to ID range residency

#make subset of IDs/periods with larger bias
aktrap.bias=aktrap[aktrap$expbias>2&!is.na(aktrap$expbias),]
akaer.bias=akaer[akaer$expbias>2&!is.na(akaer$expbias),,]
aktox.bias=aktox[aktox$expbias>2&!is.na(aktox$expbias),,]

#make joinkeys
aktrap.bias$joinkey=paste(aktrap.bias$animalid,aktrap.bias$period,sep="_")
akaer.bias$joinkey=paste(akaer.bias$animalid,akaer.bias$period,sep="_")
aktox.bias$joinkey=paste(aktox.bias$animalid,aktox.bias$period,sep="_")

#remove during for tox
aktox.bias=aktox.bias[aktox.bias$period!="during",]

#make joinkeys
geo.trap$joinkey=paste(geo.trap$animalid,geo.trap$removal.period.akdecalc,sep="_")
geo.aer$joinkey=paste(geo.aer$animalid,geo.aer$removal.period.akdecalc,sep="_")
geo.tox$joinkey=paste(geo.tox$animalid,geo.tox$removal.period.akdecalc,sep="_")

#do subset
geo.trapb=geo.trap[geo.trap$joinkey%in%aktrap.bias$joinkey,]
geo.aerb=geo.aer[geo.aer$joinkey%in%akaer.bias$joinkey,]
geo.toxb=geo.tox[geo.tox$joinkey%in%aktox.bias$joinkey,]

#Get above plots for diagnostics for only high errors
#outdir="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Exploratory_Outputs/Paths_Variograms_Biased/"
#Paths.Variogram.Diagnostics(geo.aerb,"aer",outdir)
#Paths.Variogram.Diagnostics(geo.trapb,"trap",outdir)
#Paths.Variogram.Diagnostics(geo.toxb,"tox",outdir)

#Another diagnostic plot:
#recommended in google groups discussion
#help to find cutoff
#make sure cutoff of 1% also corresponds to ones with high error
aktrap$is.bias=NA
aktrap[aktrap$expbias>2,]$is.bias<-"yes"
aktrap[aktrap$expbias<=2,]$is.bias<-"no"
ggplot(aktrap,
       aes(x=area.CI.low,xend=area.CI.high,y=animalid,color=is.bias))+
  geom_segment()+geom_point(aes(x=area.est,color=is.bias))+
  facet_wrap(~period)

#this might be a good plot for the supps
#aktrap_onepercent
akaer$is.bias=NA
akaer[akaer$expbias>2&!is.na(akaer$expbias),]$is.bias<-"yes"
akaer[akaer$expbias<=2&!is.na(akaer$expbias),]$is.bias<-"no"
ggplot(akaer,
       aes(x=area.CI.low,xend=area.CI.high,y=animalid,color=is.bias))+
  geom_segment()+geom_point(aes(x=area.est,color=is.bias))+
  facet_wrap(~period)

aktox$is.bias=NA
aktox[aktox$expbias>2&!is.na(aktox$expbias),]$is.bias<-"yes"
aktox[aktox$expbias<=2&!is.na(aktox$expbias),]$is.bias<-"no"
ggplot(aktox,
       aes(x=area.CI.low,xend=area.CI.high,y=animalid,color=is.bias))+
  geom_segment()+geom_point(aes(x=area.est,color=is.bias))+
  facet_wrap(~period)

}

##################################
#### Do segmentation process #####
##################################

#Following advice from ctmm akde developer:
#https://groups.google.com/g/ctmm-user/c/gfS4JvwGefA
  #and other posts.. searched 'home range shift', 'displacement' 'range residency' and 'forays' to get common advice from developer on these topics

{
#Function to automate segmetnation
Get.SegNum<-function(pig,kmax,lmin){
  pig$indice=1:nrow(pig)
  pig2=pig[,c(17,9,10,5)]
  colnames(pig2)=c("indice","x","y","dateTime")
  shift_seg <- segmentation(pig2,
                            seg.var = c("x","y"),
                            lmin = lmin, Kmax = kmax)
  plot_likelihood(shift_seg)
  
}

Do.Seg<-function(pig,kmax,nclust,lmin){
pig$indice=1:nrow(pig)
pig2=pig[,c(17,9,10,5)]
colnames(pig2)=c("indice","x","y","dateTime")
#shift_seg <- segmentation(pig2,
#                          seg.var = c("x","y"),
#                          lmin = 240, Kmax = 15,
#                          subsample_by = 60)
#plot_likelihood(shift_seg)
mode_segclust <- segclust(pig2,
                          Kmax = kmax, lmin=lmin, 
                          ncluster = nclust,
                          seg.var = c("x","y"),
                          scale.variable = FALSE)
#plot_BIC(mode_segclust)
#segmap(mode_segclust)
#plot(mode_segclust)
pig3=augment(mode_segclust)
pig3=pig3[,c(1,12)]
pig=left_join(pig,pig3,by="indice")

return(pig)
}
#test
#Get.SegNum(pig,40,5)
#pigseg=Do.Seg(pig,10,2:3,5)

#See how many for each
nrow(akaer.bias) #8
nrow(aktrap.bias) #4
nrow(aktox.bias) #3

#Make function to output original area ests, cluster ests, cluster combined ests

DoClustAKDE<-function(pigseg){
  clusts=unique(pigseg$state_ordered)
  outdf=data.frame(matrix(nrow=(length(clusts)+2),ncol=8))
  colnames(outdf)=c("animalid","period","vers","Neff","estbias","area.low","area.est","area.high")
  #pigseg=pigseg[pigseg$state_ordered==clusts[c],]
  outdf$animalid=pigseg$animalid[1]
  outdf$period=pigseg$removal.period.akdecalc[1]
  #First do as before to see issue
  te1=Convert.Telemetry(pigseg)
  GUESS1 <- ctmm.guess(te1, interactive = FALSE)
  print("fitting ctmm model for full path")
  FIT1_pHREML <- ctmm.select(te1, GUESS1, method = 'pHREML')
  Neff=summary(FIT1_pHREML)$DOF["area"]
  UD1_pHREML <- akde(te1, FIT1_pHREML)
  area.est.full=summary(UD1_pHREML)$CI
  
  #assign to outdf
  outdf$vers[1]="full"
  outdf$Neff[1]=Neff
  outdf$estbias[1]=1/Neff^2
  outdf[1,6:8]=area.est.full
  
  #Then do in loop for each cluster
  
  for(c in 1:length(clusts)){
    print(clusts[c])
    pigseg.clust=pigseg[pigseg$state_ordered==clusts[c],]
    print(head(pigseg.clust))
    #First do as before to see issue
    #if(nrow(pigseg.clust))
    print("converting telemetry")
    te1=Convert.Telemetry(pigseg.clust)
    print(head(te1))
    print(paste0("fitting ctmm model for cluster ",c," of ",length(clusts)))
    GUESS1 <- ctmm.guess(te1, interactive = FALSE)
    print("did Guess")
    FIT1_pHREML <- ctmm.select(te1, GUESS1, method = 'pHREML')
    print("fit phreml")
    Neff=summary(FIT1_pHREML)$DOF["area"]
    UD1_pHREML <- akde(te1, FIT1_pHREML)
    area.est.clust=summary(UD1_pHREML)$CI
    
    #assign to outdf
    outdf$vers[c+1]=clusts[c]
    outdf$Neff[c+1]=Neff
    outdf$estbias[c+1]=1/Neff^2
    outdf[c+1,6:8]=area.est.clust
    
  }
  
  #At end of loop, get summed values for area.est
  outdf$vers[nrow(outdf)]="clustsums"
  outdf[nrow(outdf),6:8]=colSums(outdf[2:(length(clusts)+1),6:8])
  #colSums(rbind(area.est.p3,area.est.p2,area.est.p1))
  beep()
  return(outdf)
  
}

#outdir="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/BiasPigsSeg/trap/"
#vers="full"
#outdf=outdf2
#append=FALSE
#remtype="trap"
#Keep.Clust(outdf,vers,outdir,append,remtype)
#Function that saves or appends list of decisions for whether to use clustered or full akde
Keep.Clust<-function(outdf,vers,outdir,append,remtype){
  outdf.keep=outdf[outdf$vers==vers,]
  if(file.exists(paste0(outdir,"biasdecisions_",remtype,".rds"))&append){
    biasdecisions.outdf=readRDS(paste0(outdir,"biasdecisions_",remtype,".rds"))
    biasdecisions.outdf=rbind(biasdecisions.outdf,outdf.keep)
    saveRDS(biasdecisions.outdf,paste0(outdir,"biasdecisions_",remtype,".rds"))
  } else{
    saveRDS(outdf.keep,paste0(outdir,"biasdecisions_",remtype,".rds"))
  }
}

###################
#Do ops here
#cant do automated loop because need to manually adjust vals
i=3
pig=geo.toxb[geo.toxb$joinkey==aktox.bias$joinkey[i],]
ggplot(pig,aes(x=X,y=Y))+
  geom_path()
Get.SegNum(pig,50,5)

kmax=9
nclust=3
lmin=5
pigseg=Do.Seg(pig,kmax,nclust,lmin)

#vizualize
ggplot(pigseg,aes(x=X,y=Y,color=as.factor(state_ordered)))+
  geom_path()
ggplot(pigseg,aes(x=X,y=Y,color=indice))+
  geom_path()+facet_wrap(~state_ordered)

#saveRDS(pigseg,paste0("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/BiasPigsSeg/tox/",aktox.bias$joinkey[i],".rds"))
###################

####Notes

i=3

#Redo akde with segs
outdir="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/BiasPigsSeg/tox/"
pigseg=readRDS(paste0("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/BiasPigsSeg/tox/",aktox.bias$joinkey[i],".rds"))
ggplot(pigseg,aes(x=X,y=Y,color=as.factor(state_ordered)))+
  geom_path()
#ggplot(pigseg[pigseg$state_ordered==1,],aes(x=X,y=Y,color=as.factor(state_ordered)))+
#  geom_path()

outdf2=DoClustAKDE(pigseg)

#vers=2
#Keep.Clust(outdf2,vers,outdir,TRUE,"tox")

}

#finished aerial
#finished trap
#finished tox


###########################################################
#### Re-calculate AKDE and merge with original output #####
###########################################################

#set outdir for akde output df with corrected hr areas
outdir.corr.akde="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/"

#function to do averages
Average.outdfclust=function(out.df.clust){
  out.df.clust2=Fix.Area.Units(out.df.clust)
  out.df.clust2$akde.hrarea.units="area (square kilometers)"
  out.df.clust.avg=out.df.clust2 %>% group_by(animalid,period) %>% dplyr::summarise(area.CI.low=mean(area.CI.low),
                                                                                    area.est=mean(area.est),
                                                                                    area.CI.high=mean(area.CI.high),
                                                                                    ctr.x=mean(ctr.x),
                                                                                    ctr.y=mean(ctr.y),
                                                                                    rem.overlaps.CI.low=NA,
                                                                                    rem.overlaps.est=NA,
                                                                                    rem.overlaps.CI.high=NA,
                                                                                    Neff=mean(Neff),
                                                                                    akde.hrarea.units=first(akde.hrarea.units),
                                                                                    expbias=mean(expbias),
                                                                                    clust="full_average") %>% as.data.frame()
  return(out.df.clust.avg)
}


#get list of each pig/cluster
#biasdecisionsaer.rds

removal="tox" #toggle this for each type
file=paste0("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/BiasPigsSeg/",removal,"/biasdecisions_",removal,".rds")
clusts=readRDS(file)

#Remove any that say full, these will not change
clusts=clusts[clusts$vers!="full",]
clusts$joinkey=paste(clusts$animalid,clusts$period,sep="_")
joinkeys=unique(clusts$joinkey)

for(i in 1:length(joinkeys)){
#for(i in 1:2){ #testing
  akde.outdir=paste0("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/AKDE_clusters/",removal,"/")
  print(paste0("starting ",joinkeys[i]," (",i," of ",length(joinkeys),")"))
  pigseg=readRDS(paste0("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline/Data/BiasPigsSeg/",removal,"/",joinkeys[i],".rds"))
  pigclust=clusts[clusts$joinkey==joinkeys[i],]
  for(j in 1:nrow(pigclust)){
    pigID=pigclust$animalid[j]
    print(paste0("starting segment ",j," of ",nrow(pigclust)))
    pigseg.clust=pigseg[pigseg$state_ordered==pigclust$vers[j],]
    period=pigseg.clust$removal.period.akdecalc[1]
    te1=Convert.Telemetry(pigseg.clust)
    print("Getting akde hr")
    out.df.clust.j=GetAKDE_AC(te1,akde.outdir,pigID,period,removal)
    out.df.clust.j$clust=pigclust$vers[j]
    
    if(i==1&j==1){
      out.df.clust=out.df.clust.j
    } else{
      out.df.clust=rbind(out.df.clust,out.df.clust.j)
    }
    
    } 
  
  }

####hr/disp exploratory/descriptive plot ideas
#1-do paths/variograms with akde polygons AND centerpoints (show displ) for each pig, with facet for before/during/after combined
#2-forest plots with area ests

#######TRAP
{
######Now read original akde output, format each, and combine
aktrap$clust="full"
aktrap=aktrap[,c(1:12,14,16)]

out.df.clust$expbias=(1/(out.df.clust$Neff^2))*100
out.df.clust=out.df.clust[,c(1:12,14,13)]
all(colnames(out.df.clust)==colnames(aktrap))

#Do averages for out.df.clust before combining with main
out.df.clust.avg=Average.outdfclust(out.df.clust)
#remove IDs in akaer that are in out.df.clust
aktrap$joinkey=paste(aktrap$animalid,aktrap$period,sep="_")
out.df.clust.avg$joinkey=paste(out.df.clust.avg$animalid,out.df.clust.avg$period,sep="_")
aktrap.remc=aktrap[-which(aktrap$joinkey%in%out.df.clust.avg$joinkey),]
aktrap.c=rbind(out.df.clust.avg,aktrap.remc)
aktrap.c=left_join(aktrap.c,geo.id,by="animalid")
#saveRDS(aktrap.c,paste0(outdir.corr.akde,"outdf_akde_trap_corrected.rds")) #Need save it

#plot forest plot
ggplot(aktrap.c,
       aes(x=area.CI.low,xend=area.CI.high,y=animalid,color=period))+
  geom_segment()+geom_point(aes(x=area.est,color=period))+
  facet_wrap(~Removal.Type)+theme_ipsum()
#ggsave("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Exploratory_Outputs/Forest_Plots/trap_forestplot.png",height=2400,width=2000,units="px")

}

#######AER
{
######Now read original akde output, format each, and combine
akaer$clust="full"
out.df.clust$expbias=(1/(out.df.clust$Neff^2))*100
akaer=akaer[,c(1:12,14,16)]
out.df.clust=out.df.clust[,c(1:12,14,13)]
all(colnames(out.df.clust)==colnames(akaer))

#Do averages for out.df.clust before combining with main
out.df.clust=Average.outdfclust(out.df.clust)
#remove IDs in akaer that are in out.df.clust
akaer$joinkey=paste(akaer$animalid,akaer$period,sep="_")
out.df.clust.avg$joinkey=paste(out.df.clust.avg$animalid,out.df.clust.avg$period,sep="_")
akaer.remc=akaer[-which(akaer$joinkey%in%out.df.clust.avg$joinkey),]
akaer.c=rbind(out.df.clust.avg,akaer.remc)
akaer.c=left_join(akaer.c,geo.id,by="animalid")

#Need save it, change path later
#saveRDS(akaer.c,paste0(outdir.corr.akde,"outdf_akde_aer_corrected.rds")) 

#plot forest plot
ggplot(akaer.c,
       aes(x=area.CI.low,xend=area.CI.high,y=animalid,color=period))+
  geom_segment()+geom_point(aes(x=area.est,color=period))+
  facet_wrap(~Removal.Type)+theme_ipsum()
#ggsave("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Exploratory_Outputs/Forest_Plots/aer_forestplot.png",height=2400,width=2000,units="px")

}

#######TOX
{
  ######Now read original akde output, format each, and combine
  aktox$clust="full"
  out.df.clust$expbias=(1/(out.df.clust$Neff^2))*100
  aktox=aktox[,c(1:12,14,16)]
  out.df.clust=out.df.clust[,c(1:12,14,13)]
  all(colnames(out.df.clust)==colnames(aktox))
  
  #Do averages for out.df.clust before combining with main
  out.df.clust.avg=Average.outdfclust(out.df.clust)
  #remove IDs in akaer that are in out.df.clust
  aktox$joinkey=paste(aktox$animalid,aktox$period,sep="_")
  out.df.clust.avg$joinkey=paste(out.df.clust.avg$animalid,out.df.clust.avg$period,sep="_")
  aktox.remc=aktox[-which(aktox$joinkey%in%out.df.clust.avg$joinkey),]
  aktox.c=rbind(out.df.clust.avg,aktox.remc)
  aktox.c=left_join(aktox.c,geo.id,by="animalid")
  
  #Need save it, change path later
  #saveRDS(aktox.c,"~/Downloads/aktox_corr.rds") 
  
  #plot forest plot
  ggplot(aktox.c3,
         aes(x=area.CI.low,xend=area.CI.high,y=animalid,color=period))+
    geom_segment()+geom_point(aes(x=area.est,color=period))+
    facet_wrap(~Removal.Type)+theme_ipsum()
#ggsave("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Exploratory_Outputs/Forest_Plots/tox_forestplot.png",height=2400,width=2000,units="px")

###looks pretty good, but need to remove 85411_C3_C3 
  #and 86075_Y5_Y5
#could not subset to range resident behaviors

aktox.c2=aktox.c[aktox.c$animalid!="85411_C3_C3"&aktox.c$animalid!="86075_Y5_Y5",]
aktox.c3=aktox.c[aktox.c$period!="during",]


#saveRDS(aktox.c3,paste0(outdir.corr.akde,"outdf_akde_tox_corrected.rds")) 


}





