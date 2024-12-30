
GetRemovalOverlaps<-function(mcps,IDs,removal.object){
  
  overlap.props=c()
  overlap.areas=c()
  for(x in 1:nrow(mcps)){
    print(x)
    t=st_intersection(mcps[x,], removal.object) %>% 
      mutate(intersect_area = st_area(.)) #%>%   # create new column with shape area
    if(nrow(t)>0){
      p_overlap=st_area(t)/st_area(mcps[x,])
      area_overlap=st_area(t)
    } else{
      p_overlap=0
      area_overlap=0
    }
    
    overlap.props=c(overlap.props,p_overlap)
    overlap.areas=c(overlap.areas,area_overlap)
  }
  
  overlap.df=data.frame("animalid"=IDs,"prop_overlap"=overlap.props,"area_overlap"=overlap.areas)
  
  return(overlap.df)
  
}