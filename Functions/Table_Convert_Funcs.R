# Function to pull needed info for table from sdmTMB object
#res is sdmTMB model object
get_tidy_coefs<-function(res){
  tib=tidy(res)
  tib$statistic=tib$estimate/tib$std.error
  tib$p_value=exp(-0.717*tib$statistic-0.416*tib$statistic^2)
  tib$cl_low=tib$estimate-2*tib$std.error
  tib$cl_high=tib$estimate+2*tib$std.error
  tib$estimate=exp(tib$estimate)
  tib$cl_low=round(exp(tib$cl_low),digits=2)
  tib$cl_high=round(exp(tib$cl_high),digits=2)
  tib$ci=glue::glue("{tib$cl_low},{tib$cl_high}")
  tib=tib[,c(1,2,8,5)]
  #colnames(tib)[c(1,2,4)]<-c("Characteristic","exp(Beta)","p-value")
  
  #remove intercept to match other gt
  tib=tib[-1,]
  
  #format numbers
  tib[,2]=round(tib[,2],digits=2)
  tib[,4]=formatC(as.numeric(tib[,4][[1]]), format = "e",digits=2)
  return(tib)
}

#res=sdmTMB model object
#dat=data used to run sdmTMB model object
# Function to format sdmTMB model info into gt table object
make_sdmTMB_gt<-function(res,dat){
  tib=get_tidy_coefs(res)
  call=eval(parse(text=as.character(res$call)[2]))
  res=glmmTMB(call,data=dat)
  tbl_gt <- tbl_regression(res, exponentiate = TRUE)
  #tbl_gt$table_body %>% as.data.frame()
  tbl_gt$table_body$estimate[!is.na(tbl_gt$table_body$estimate)]=tib$estimate
  tbl_gt$table_body$ci[!is.na(tbl_gt$table_body$ci)]=tib$ci
  tbl_gt$table_body$p.value[!is.na(tbl_gt$table_body$p.value)]=as.numeric(tib$p_value)
  return(tbl_gt)
}
