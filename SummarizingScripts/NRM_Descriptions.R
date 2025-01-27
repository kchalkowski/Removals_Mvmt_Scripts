library(ggplot2)
library(vistime)

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

#Get data start and ends, to trim data
#determined in AKDE_HR_Diagnostics:
aer.start.date=as.Date("2023-03-27")
aer.end.date=as.Date("2023-03-29")
aer.data.start=aer.start.date-37
aer.data.end=aer.end.date+37

trap.start.date=as.Date("2023-02-17")
trap.end.date=as.Date("2023-05-04")
trap.data.start=as.Date("2022-12-23")
trap.data.end=as.Date("2023-06-30")

tox.start.date=as.Date("2023-02-17")
tox.end.date=as.Date("2023-03-09")
tox.data.start=as.Date("2023-01-11")
tox.data.end=as.Date("2023-04-15")

start.vec=c(aer.start.date,trap.start.date,tox.start.date)
end.vec=c(aer.end.date,trap.end.date,tox.end.date)
data.start.vec=c(aer.data.start,trap.data.start,tox.data.start)
data.end.vec=c(aer.data.end,trap.data.end,tox.data.end)

timeline.df=data.frame("Treatment"=c("aer","trap","tox"),
           "Tracking.Start"=data.start.vec,
           "Tracking.End"=data.end.vec,
           "Treatment.Start"=start.vec,
           "Treatment.End"=end.vec)

timeline.df.before=data.frame("Treatment"=c("aer","trap","tox"),
                              "Period"="before",
                              "Start"=data.start.vec,
                              "End"=start.vec)
timeline.df.during=data.frame("Treatment"=c("aer","trap","tox"),
                              "Period"="during",
                              "Start"=start.vec,
                              "End"=end.vec)
timeline.df.after=data.frame("Treatment"=c("aer","trap","tox"),
                              "Period"="after",
                              "Start"=end.vec,
                              "End"=data.end.vec)

timeline.df=rbind(timeline.df.before,timeline.df.during,timeline.df.after)
timeline.df$color= NA
timeline.df$color[timeline.df$Period=="before"&timeline.df$Treatment=="aer"]<-"#cdabff"
timeline.df$color[timeline.df$Period=="during"&timeline.df$Treatment=="aer"]<-"#822AFF"
timeline.df$color[timeline.df$Period=="after"&timeline.df$Treatment=="aer"]<-"#311854"
timeline.df$color[timeline.df$Period=="before"&timeline.df$Treatment=="trap"]<-"#fcaed4"
timeline.df$color[timeline.df$Period=="during"&timeline.df$Treatment=="trap"]<-"#FF57A9"
timeline.df$color[timeline.df$Period=="after"&timeline.df$Treatment=="trap"]<-"#632242"
timeline.df$color[timeline.df$Period=="before"&timeline.df$Treatment=="tox"]<-"#ffbf99"
timeline.df$color[timeline.df$Period=="during"&timeline.df$Treatment=="tox"]<-"#FF985A"
timeline.df$color[timeline.df$Period=="after"&timeline.df$Treatment=="tox"]<-"#8a4820"

#library(vistime)
timeline=gg_vistime(timeline.df, col.event = "Period", col.color="color", col.group = "Treatment", 
        col.start = "Start", col.end = "End",optimize_y=TRUE, show_labels=FALSE)+
  theme_ipsum(axis_text_size=40)
#ggsave("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Descriptive_Outputs/Methods_Descriptive_Fig/timeline.png",bg="white",height=3,width=11,units="in")
