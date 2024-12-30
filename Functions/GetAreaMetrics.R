GetAreaMetrics<-function(trknested,out,mcp.iso,IDs){
  
  kde_list=vector(mode="list",length=length(trknested))
  area_list=vector(mode="list",length=length(trknested))
  mcp_list=vector(mode="list",length=length(trknested))
  date_list=vector(mode="list",length=length(trknested))
  for(i in 1:length(trknested)){
    #print(i)
    if(!(
      (max(trknested[[i]]$x_)-min(trknested[[i]]$x_)<0.01)| #if didn't move at all
      (max(trknested[[i]]$y_)-min(trknested[[i]]$y_)<0.01)| #or less than 10 points, skip
      length(trknested[[i]]$x_)<10)){
      if("kde"%in%out|"area"%in%out){
        kde_list[[i]]=hr_kde(trknested[[i]])
        area_list[[i]]=as.numeric(hr_area(kde_list[[i]])$area)
      }
      if("mcp"%in%out){
        mcp_list[[i]]=hr_mcp(trknested[[i]],levels=c(mcp.iso))
        date_list[[i]]=as.Date(mcp_list[[i]]$data$t_[1])
      }
    } else{
      kde_list[[i]]=0
      area_list[[i]]=0
      mcp_list[[i]]=0
    }
  }
  
  #get clean list of just the mcps, areas, whatever
  if("mcp"%in%out){
    mcp.list.2=vector(mode="list",length=length(mcp_list$mcp))
    for(i in 1:length(mcp_list)){
      if(length(mcp_list[[i]])>1){
        mcp.list.2[[i]]=mcp_list[[i]]$mcp
      } else{mcp.list.2[[i]]=NA}
    }
    
    #na.IDs=IDs[which(is.na(mcp.list.2))]
    notna.IDs=IDs[which(!(is.na(mcp.list.2)))]
    mcp.list.2=mcp.list.2[!is.na(mcp.list.2)]
    
    trpmcps=dplyr::bind_rows(mcp.list.2)
    dates=as.Date(unlist(date_list))
    
    trpmcps$animalid=notna.IDs
    trpmcps$date=dates
    mcp_list=trpmcps
  }
  
  out.list=list("kde"=kde_list,"area"=area_list,"mcp"=mcp_list)
  #which(names(out.list)%in%out)
  
  return(out.list[which(names(out.list)%in%out)])
  
}