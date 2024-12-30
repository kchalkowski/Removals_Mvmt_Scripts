#test:
#td.input=geo$datetime
#td.input=temp$datetime
#td.input<-geolocs$datetime
#Input a vector of datetimes with different formats, returns vector of datetimes in character format, with same structure

#assumptions of function:
#organized date space time (ie 2023-02-14 05:55:02)
#separators in dates are all same (ie 2023-02-13 OR 2023/02/13, NOT 2023/02-13)
#dates are month-day-year

#assumptions of function:
#organized as date space time (ie 2023-02-14 05:55:02)
#separators in dates are all same (ie 2023-02-13 OR 2023/02/13, NOT 2023/02-13)
#dates are month-day-year

#Can handle:
#if some in vector are 24 hr, some 12 hr (with AM PM)
#year in front or back of date
#mixed double and single digits for month/day

NeatDates<-function(td.input){
  if(all(is.na(td.input))){return(td.input)
  } else{
    td<-td.input[!is.na(td.input)]
    jdt=rep(NA,length(td))
    jt=rep(NA,length(td))
    jm=rep(NA,length(td))
    jday=rep(NA,length(td))
    jy=rep(NA,length(td))
    jtx=rep(NA,length(td))
    #ask if any spaces
    has.space=grep("\\s",td,perl=TRUE)
    
    #ask if no space in string
    no.space=which(!(1:length(td)%in%grep("\\s",td,perl=TRUE)))
    
    #extract first word before space
    jdt[has.space]=str_sub(str_match(td,regex("^\\S+\\s")),1L,-2L)[has.space]
    #jdt[1142348:1142415]
    
    #paste ones without space, date only
    jdt[no.space]=td[no.space]
    
    #extract second word before space, put in just time col jt
    jt=str_sub(str_match(td,regex("\\s\\S+$")),2L,-1L)
    #jt[1142348:1142415]
    
    #check if there are any with two separate spaces (e.g. UTC or PM)
    two.spaces=grep("\\s\\S+\\s",td,perl=TRUE)
    jt[two.spaces]=str_sub(str_match(td[two.spaces],regex("\\s\\S+\\s")),2L,-2L)
    jtx[two.spaces]=str_sub(str_match(td[two.spaces],regex("\\s\\S+$")),2L,-1L)
    #jt[1142348:1142415]
    
    #ask if there is a four digit sequence of numbers
    
    #for single digits, replace with digits with 0 in front
    #jdt=str_replace_all(jdt,regex("\\b\\d{1}\\b"),paste0(0,str_match(jdt,regex("\\b\\d{1}\\b"))))
    #jdt[1142348:1142415]
    #^^^^ this is problem
    jdt=str_replace(jdt,regex("\\b\\d{1}\\D"),paste0(0,str_sub(str_match(jdt,regex("\\b\\d{1}\\D")),1L,-2L),str_match(jdt,regex("\\D"))))
    jdt=str_replace(jdt,regex("\\b\\d{1}\\D"),paste0(0,str_sub(str_match(jdt,regex("\\b\\d{1}\\D")),1L,-2L),str_match(jdt,regex("\\D"))))
    
    #ask where there is a four digit sequence of numbers
    is.four=grep("\\d{4}",jdt,perl=TRUE)
    no.four=which(!(1:length(td)%in%is.four))
    #return four digit sequence
    
    #return all two digit numbers in each match, for those without four digits
    all.two=str_match_all(jdt[no.four],regex("\\b\\d{2}\\b"))
    jm[no.four]=sapply(all.two,"[[",1)
    jday[no.four]=sapply(all.two,"[[",2)
    jy[no.four]=sapply(all.two,"[[",3)
    
    #turn all two digit jy's to four digit jy's
    #don't have anything <2000, so just paste 20 to all two digit years
    jy[!is.na(jy)]=paste0("20",jy[!is.na(jy)])
    
    #return four digit years to year vector
    jy[is.four]=str_match(jdt[is.four],regex("\\d{4}"))
    
    #occasionally have dates represented like: 0012
    #get all years formatted this way, starting with two 0's
    jy[grep("^[0]{2}",jy,perl=TRUE)]<-str_replace(jy[grep("^[0]{2}",jy,perl=TRUE)],"00","20")
    
    #return all two digit numbers
    all.two.isfour=str_match_all(jdt[is.four],regex("\\b\\d{2}\\b"))
    jm[is.four]=sapply(all.two.isfour,"[[",1)
    jday[is.four]=sapply(all.two.isfour,"[[",2)
    
    #do some jt formatting
    #if only 00:00 (hour and minute), add second to make consistent format
    jt[grep("^[0-9]*:[0-9]*$",jt,perl=TRUE)]<-paste0(jt[grep("^[0-9]*:[0-9]*$",jt,perl=TRUE)],":00")
    
    #also want hour to be double digits, for consistency
    jt[grep("^[0-9]{1}:",jt,perl=TRUE)]<-paste0(0,jt[grep("^[0-9]{1}:",jt,perl=TRUE)])
    
    #if any jt==00:00:00, add one second
    #this is because R will OBNOXIOUSLY drop times that are equal to 00:00 or 00:00:00 in POSIXct format
    #sometimes, this results in R dropping all times from the vector.... :) SO fun
    #adding one second is easy way to avoid this without changing data much
    #(for animal movement data, one second is inconsequential)
    jt[jt=="00:00:00"]<-"00:00:01"
    
    #paste together as same format date strings
    jdt.new=paste(jm,jday,jy,sep="/") 
    
    neat.dates=jdt.new
    neat.dates[has.space]=paste(neat.dates[has.space],jt[has.space],sep=" ")
    neat.dates[two.spaces]=paste(neat.dates[two.spaces],jtx[two.spaces],sep=" ")
    
    #####Times, remove am/pm and convert to 24 hr clock
    ##only acts on neat.dates if there is 'am' or 'pm' in string
    #remove AM
    am.index=grep("AM",neat.dates,perl=TRUE)
    front=str_sub(neat.dates[am.index],1L,11L)
    back=str_sub(neat.dates[am.index],14L)
    hours=as.numeric(str_sub(str_match(neat.dates[am.index],regex("\\s\\d{2}")),2L))
    new.hours=hours
    new.hours[which(hours==12)]=0
    neat.dates[am.index]=paste0(front,new.hours,back)
    neat.dates[am.index]=str_sub(neat.dates[am.index],1L,-4L)
    
    #remove pm, but add 12 to hours
    pm.index=grep("PM",neat.dates,perl=TRUE)
    front=str_sub(neat.dates[pm.index],1L,11L)
    back=str_sub(neat.dates[pm.index],14L)
    hours=as.numeric(str_sub(str_match(neat.dates[pm.index],regex("\\s\\d{2}")),2L))
    new.hours=hours+12
    new.hours[which(hours==12)]=hours[which(hours==12)]
    neat.dates[pm.index]=paste0(front,new.hours,back)
    neat.dates[pm.index]=str_sub(neat.dates[pm.index],1L,-4L)
    
    #td<-td.input[!is.na(td.input)]
    td.input[!is.na(td.input)]<-neat.dates
    
    return(td.input)
  }
}


Neat.Dates.POSIXct<-function(td.input,tz){
  if(all(is.na(td.input))){return(td.input)
  } else{
    neat.dates<-NeatDates(td.input)
    neat.dates.posixct=as.POSIXct(neat.dates,tz=tz,tryFormats=c("%m/%d/%Y %H:%M:%OS",
                                                                "%m/%d/%Y %H:%M:%S",
                                                                "%m/%d/%Y %H:%M",
                                                                "%m/%d/%Y"))
    
    return(neat.dates.posixct)
  }
}


#td.input<-geolocs.t5$datetime

#Input a vector of datetimes with different formats, returns vector of datetimes in POSIXct format, with same structure
Neat.Dates.POSIXct<-function(td.input,tz){
  if(all(is.na(td.input))){return(td.input)
  } else{
  neat.dates<-NeatDates(td.input)
  neat.dates.posixct=as.POSIXct(neat.dates,tz=tz,tryFormats=c("%m/%d/%Y %H:%M:%OS",
                                                              "%m/%d/%Y %H:%M:%S",
                                                              "%m/%d/%Y %H:%M",
                                                              "%m/%d/%Y"))

  return(neat.dates.posixct)
  }
}

#length(geolocs$datetime)
#length(neat.dates)

#neat.dates=NeatDates(geolocs$datetime)
#as.POSIXct(geolocs$datetime[236086:236087],tz=tz,tryFormats=c("%m/%d/%Y %H:%M:%OS",
#                                                            "%m/%d/%Y %H:%M:%S",
#                                                            "%m/%d/%Y %H:%M",
#                                                            "%m/%d/%Y"))
#as.POSIXct(geolocs$datetime[236088:236089],tz=tz,tryFormats=c("%m/%d/%Y %H:%M:%OS",
#                                                        "%m/%d/%Y %H:%M:%S",
#                                                        "%m/%d/%Y %H:%M",
#                                                        "%m/%d/%Y"))
#236087-236088
#geolocs[236088,]
#neat.dates[236088,]
