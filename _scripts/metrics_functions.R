######
#####
##### functions used during training process

### read and resample radar data 

optical_metrics = function(tile,init_date,end_date,
                           optical_bands = c("blue","green","red","vrededg1","vrededg2","vrededg3","nir","narrownir","swir1","swir2"),
                           bad_pixels_limit = 15){
  
  ## get sentinel 2 images
  
  inventory = get_inventory(tile)
  optical_images = get_optical_imagery(inventory , init_date = DateStart, end_date = DateEnd, 
                                       bands = optical_bands,bad_pixels_limit = bad_pixels_limit, optical_mission = c("L2A","S2"))
  
  ########## ----- > Create features for optical images 
  
  dataModel_Param_optical=stack(
    lapply (optical_bands,function(band){
      
      band_info = stack(lapply(optical_images , function(x) x[[band]]))
      
      band_ifo = metrics_(band_info)
      
      return(band_ifo)
    })
  )
  ## set a satellite image as reference
  satimage = optical_images[[1]][[1]]
  rm(optical_images)
  ## transform optical data, extract data and storage in a ff data frame
  
  dataImages_optical= ff(vmode="double",dim=c(ncell(dataModel_Param_optical), nlayers(dataModel_Param_optical)))
  
  for(raster_layer in 1:nlayers(dataModel_Param_optical)){
    dataImages_optical[,raster_layer]=dataModel_Param_optical[[raster_layer]][]
  }
  rm(dataModel_Param_optical)
  optical_names_database = do.call(c,lapply(optical_bands, function(band) paste0(band ,c("_SD","_MEAN","_MIN","_MAX"))))
  
  return(list(features_data = dataImages_optical, feature_names = optical_names_database, satimage = satimage))
  
}

radar_metrics = function(tile,init_date,end_date,
                         radar_Bands = c("db"),
                         polarization = "VV"){
  
  ## set band names 
  
  if(length(polarization) == 1){
    radar_Bands = paste0("Sigma0_",polarization, "_",radar_Bands)
  }else{
    cat("TODO: module is not implemented yet ..")
  }
  
  ## read images 
  inventory = get_inventory(tile)
  
  radar_images = get_rada_images(inventory , init_date = DateStart, end_date = DateEnd, 
                                 bands = radar_Bands)
  
  ## set a satellite image as reference
  satimage = radar_images[[1]][[1]]
  
  
  ########## ----- > Create Parametrics for optical images
  quantile_imags = calculate_raster_quantiles(radar_images,radar_Bands = radar_Bands, division = 25, ncores = 5)
  
  ## database is transformed to ff format, which is ideal to optimize space
  quantile_imags = stack(quantile_imags)
  
  dataImages_radar= ff(vmode="double",dim=c(ncell(quantile_imags), nlayers(quantile_imags)))
  
  for(raster_layer in 1:nlayers(quantile_imags)){
    dataImages_radar[,raster_layer]=quantile_imags[[raster_layer]][]
  }
  rm(quantile_imags)
  radar_names_database = do.call(c,lapply(radar_Bands, function(band) paste0(band,"_pa_",1:length(c(0.05,.25, .50, .75,0.95)))))
  
  return(list(features_data = dataImages_radar, feature_names = radar_names_database, satimage = satimage))
}


calculate_raster_quantiles = function(radar_images,radar_Bands = "db", division = 20, ncores = 2){
  

  SentinelimageReference = radar_images[[1]][[1]]
  ###split rows
  levelsDiv=cut_number(1:ncell(SentinelimageReference), division)
  
  quantile_imags =(lapply(radar_Bands, function(band) {
    
    extract_band_info = stack(lapply(radar_images, 
                                     function(x) x[[band]]))
    
    dataImages= ff(vmode="double",dim=c(ncell(extract_band_info), nlayers(extract_band_info)))
    
    for(raster_layer in 1:nlayers(extract_band_info)){
      dataImages[,raster_layer]=extract_band_info[[raster_layer]][]
    }
    
    
    cl =  makeCluster(ncores)
    registerDoParallel(cl)
    Sys.time()->start
    quantile_imags = do.call(rbind,foreach(step_ = 1:length(unique(levelsDiv)), 
                                           .export = "levelsDiv") %dopar% {
                                             
                                             library(ff)
                                             rowstoSelect=which(levelsDiv%in%unique(levelsDiv)[step_])
                                             
                                             Q10 = t(apply(dataImages[rowstoSelect,], 1,quantile, probs = c(0.05,.25, .50, .75,0.95), na.rm=TRUE))
                                             cat(Sys.time()-start)
                                             Q10
                                           })
    quantile_imags = raster::stack(lapply(1:ncol(quantile_imags), function(rastval){
      SentinelimageReference[] = NA
      SentinelimageReference[] = quantile_imags[,rastval]
      return(SentinelimageReference)
    }))
    stopCluster(cl)
    cat("\n",band," processed time:")
    cat(Sys.time()-start)
    
    return(quantile_imags)}))
  return(quantile_imags)
}

std_raster = function(rimages){
  meanValues=mean(rimages,na.rm=T)
  SumbandDiff=list()
  for (j in 1:nlayers(rimages)){
    SumbandDiff[[j]]=(rimages[[j]]-meanValues)*(rimages[[j]]-meanValues)
  }
  SumbandDiff=sum(stack(SumbandDiff),na.rm=T)
  
  stdValues=((SumbandDiff*(1/(nlayers(rimages)-1)))^(1/2))
  
}

metrics_ = function(images, metrics =c("std","mean", "min", "max")){
  
  ## define options
  list_metrics = list("mean" = "mean(images,na.rm=T)",
                      "std" = "std_raster(images)",
                      "min" = "min(images,na.rm=T)",
                      "max" = "max(images,na.rm=T)")
  
  ## export raster stacked 
  return(stack(lapply(list_metrics[metrics],function(metric) eval(parse(text = metric)))))
}


calculate_rastermetrics = function(DateStart, DateEnd,radar_bands, SentinelimageReference,
                                   fpath = "process/rdata/radar_info/", metrics = c("mean", "min", "max")){
  ## set capturing period
  Date_interval = DateStart : DateEnd
  # filter filenames
  fnames = list.files(fpath, pattern = paste0(Locality , ".RData$"))
  
  ## Filter files by dates
  
  fnames = fnames[as.Date(IdentiyImageDate(fnames), format = "%Y%m%d")%in% Date_interval]
  fnames = fnames[order(as.Date(IdentiyImageDate(fnames), format = "%Y%m%d"))]
  
  ## read imagery
  radar_images = read_radar_Rdata (fnames, radar_bands,
                                   SentinelimageReference, paste0(Main_Folder , fpath))
  
  ## calculate metrics
  metrics_features = lapply(radar_bands, function(band){
    band_info = stack(lapply(radar_images , function(x) x[[band]]))

    band_ifo = metrics_(band_info,metrics)
  })
  
  stack(metrics_features)
  
}

