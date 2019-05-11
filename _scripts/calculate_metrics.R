######
#####
##### functions used during training process

### read and resample radar data 

read_radar_Rdata = function(list_radar_names,radar_Bands, SentinelimageReference, fpath =  "process/rdata/radar_info/"){

  radar_images = lapply(1:length(list_radar_names), function(imag_index){
    file_name = list_radar_names[imag_index]
    
    load(file = paste0(fpath,file_name))
    
    bands_toSelect = names(radar_info)[names(radar_info)%in%radar_Bands]
    radar_info = radar_info[[bands_toSelect]]
    ##project raster
    
    radar_info_projected = projectRaster(radar_info, crs = crs(SentinelimageReference))
    ## resample
    radar_info_resample = resample(radar_info_projected, SentinelimageReference, method = "ngb")
    
    
    cat(str_sub(file_name ,1, -7) ," loaded\n")
    return(radar_info_resample)
  })
  names(radar_images) = list_radar_names
  return(radar_images)
}




calculate_raster_quantiles = function(radar_images,radar_Bands = "db",polarization = "VV", division = 20, ncores = 5){
  
  ## set band names 
  if(length(polarization) == 1){
    bands = paste0("Sigma0_",polarization, "_",radar_Bands)
  }else{
    cat("TODO: module is not implemented yet ..")
  }
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

