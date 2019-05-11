################ 
### 
### This is a code which contains all the function created for the project Rice Crop Identification
### 
###                   
###   Author:  Andr?s Aguilar
###     CIAT - DAPA - AEPS
#########


### Function used for plotting Images and constructing Quality pixels table
setPathStyle=function(pathOrig){ifelse(str_sub(pathOrig,-1)=="/",pathOrig,paste0(pathOrig,"/"))}

### 
initconf = 
  read.table(paste0(dirname(sys.frame(1)$ofile),"//init_configuration.txt"), header = FALSE, sep = "=")



#### 

### Set stallite parameters
get_sat_list_parameters = function(mission){

  Sentinel2_Folder=setPathStyle(as.character(initconf[as.character(initconf[,1])%in%'sentinel2_images',2]))
  Sentinel1_Folder=setPathStyle(as.character(initconf[as.character(initconf[,1])%in%'sentinel1_images',2]))
  Landsat_Folder=setPathStyle(as.character(initconf[as.character(initconf[,1])%in%'landsat_images',2]))
  
  ### bands names for each mission
  
  # original names
  S1_VV_Bands=c("Sigma0_VV_db",'Sigma0_VV_Contrast', 'Sigma0_VV_Dissimilarity',
                'Sigma0_VV_Energy', 'Sigma0_VV_Entropy',
                'Sigma0_VV_GLCMCorrelation',
                'Sigma0_VV_GLCMMean', 'Sigma0_VV_GLCMVariance',
                'Sigma0_VV_Homogeneity',
                'Sigma0_VV_MAX')
  
  
  S2_Bands = c(  "B2",    "B3",    "B4","B8", "B5" ,"B6","B7", "B8A", "B11",   "B12","quality_scene_classification","quality_cloud_confidence")
  L8_Bands = c( "sr_band2",    "sr_band3",    "sr_band4",    "sr_band5", "sr_band6", "sr_band7",    "pixel_qa",    "sr_aerosol")
  L7_Bands = c("sr_band1",   "sr_band2",    "sr_band3",    "sr_band4",    "sr_band5", "sr_band7",    "pixel_qa","sr_cloud_qa")
  
  # Standar names
  
  Standar_bandNames_Landsat = c("blue","green","red","nir","swir1","swir2","QL","SLC")
  Standar_bandNames_Sentinel2 = c("blue","green","red","nir","vrededg1","vrededg2","vrededg3","narrownir","swir1","swir2","SLC","QL")
  
  
  ##### Quality threshoolds
  
  S1QualityPixels=c(1000)
  S2QualityPixels=list("QL" =9,
                       "SLC" = c(0, 1, 3, 8, 9, 10, 11))
  
  L8QualityPixels=list("QL" =c(323),
                       "SLC" = c(8, 16,80,208, 144, 100))
  L7QualityPixels=list("QL" = c(67),
                       "SLC" = c(2, 4,8,12, 20, 34,36,52))
  
  
  #### Class information
  
  S1_VV_Param=SatelliteParameters("S1_VV",Sentinel1_Folder,S1_VV_Bands,"tif",S1QualityPixels)
  S2_Param=SatelliteParameters("L2A",Sentinel2_Folder,S2_Bands,"tif",S2QualityPixels , Standar_bandNames_Sentinel2)
  L8_Param=SatelliteParameters("LC08",Landsat_Folder,L8_Bands,"tif",L8QualityPixels , Standar_bandNames_Landsat)
  L7_Param=SatelliteParameters("LE07",Landsat_Folder,L7_Bands,"tif",L7QualityPixels , Standar_bandNames_Landsat)
  
  sat_Params=switch(mission,
                    "L2A"=S2_Param,"S2"=S2_Param,
                    "LC08"=L8_Param,"L8" =L8_Param ,
                    "LE07"=L7_Param, "L7" = L7_Param,
                    "S1_VV" = S1_VV_Param)
  
  return(sat_Params)
}

### Class Definition
SatelliteParameters = function(SatIndex,DirFolder,BandNames,ExtImg,QualityLimits , standar_band_names = ""){
  SatParam = list(
    Sat_Index = SatIndex,
    Band_Names = BandNames,
    Dir_Path = DirFolder,
    ExtImg=ExtImg,
    QualityLimits=QualityLimits,
    Standar_bandNames = standar_band_names
  )
  class(SatParam) = append(class(SatParam),"SatelliteParemetes")
  return(SatParam)
}


### create quichlooks


create_quicklooks = function(inventory,
                             start_date, end_date,
                             quicklooks_path, remove_clouds =T, 
                             band_combination = c("red", "green","blue"),
                             exportasrdata = F){
  
  satellite_image = get_optical_imagery(inventory,start_date, end_date,band_combination,bad_pixels_limit =100,remove_clouds,
                                 imagesto_delete = NA)
  if(exportasrdata){
    sat_imag = satellite_image[[1]]
    tif_tordata(satellite_image[[1]],names(satellite_image),tile,
                optical_folder)
  }
  lapply(1:length(satellite_image), function(i){
    PlotRGBImage(satellite_image[[i]] , names(satellite_image)[i], OutputDirectoty = quicklooks_path, bands = band_combination)
    
  }
    
  )
    
}

###
PackageReading=function(libs){
  lapply(libs,function(package){
    if(require(package,character.only = T)){
      require(package,character.only = T)
    }else{
      install.packages(package)
      require(package,character.only = T)
    }
  })
}


UnzipLandsatFiles = function(dirFolPrincipal,nameFolder,outputLandsat7,outputLandsat8,L8Pattern="LC08",extentRef){
  tempna=tempfile()
  cat("UnTar: " , i ," \n")
  untar(paste0(dirFolPrincipal,nameFolder), exdir = tempna)
  
  ### Create Directory to export the crop rasters
  NameFile=list.files(path = tempna, pattern = "*.tif$")[1]
  rx = gregexpr(".*.T1",NameFile, perl = TRUE)
  nameImage=str_sub(unlist(regmatches(NameFile, rx)),1,-1)
  if(length(nameImage) ==0){
    rx = gregexpr(".*.RT",NameFile, perl = TRUE)
    nameImage=str_sub(unlist(regmatches(NameFile, rx)),1,-1)
  }
  if(grepl(NameFile, pattern = L8Pattern)){
    DirLandsat=setPathStyle(outputLandsat8)
    NamesBands=L8_Bands
  }else{
    DirLandsat=setPathStyle(outputLandsat7)
    NamesBands=L7_Bands
  } 
  
  ### Crop Correct Bands
  cat("Crop : " , nameImage ," \n")
  Cropraster=MaskLandsatImages(tempna,extentRef,NamesBands)
  
  ### Export reduced rasters
  cat("Export : " , nameImage ," \n")
  writeRaster(Cropraster,filename=paste0(DirLandsat,nameImage,".tif") , format="GTiff", overwrite=T)
  
  
  unlink(tempna, recursive=TRUE)
  cat("Done !!!!\n\n")
}

calculate_vi = function(sat_imag , names_idexes, equations){
  
  NewLayers=lapply(equations,function(x){
    eval(parse(text = replace_equation(x)))})
  
  imageI=stack(NewLayers)
  names(imageI) = names_idexes
  return(imageI)
}

identify_index_bands= function(equation){
  # this are the bands commonly used in this study
  standar_bandnames =  c("blue","green","red","nir","vrededg1","vrededg2","vrededg3",
                         "narrownir","swir1","swir2")
  
  return(standar_bandnames[sapply(standar_bandnames,function(band){grepl(equation,pattern=band)})])
}


## get vegetation indexes layers
get_vi_layers = function(inventory , veg_indexes, startdate, enddate, 
                         satelliteimages_path,PercMaxBadPixels, add_new_vi = NA){
  
  # vegetations indexes available
  ListIndex=list(NDVI="(nir-red)/(nir+red)",
                 LSWI="(nir-swir1)/(nir+swir1)",
                 GNDVI="(green-nir)/(green+nir)",
                 BSI="(((2*swir1/(swir1+nir))-((nir/(nir+red))+(green)/(green + swir1) ))/((2*swir1/(swir1+nir))+((nir/(nir+red))+(green)/(green + swir1) )))")
  
  ## selcet bands to use
  bands_touse = do.call(c, lapply(veg_indexes, function(index) identify_index_bands(ListIndex[[index]])))
  
  ## read optical data
  
  optical_images = get_optical_imagery(inventory , init_date = startdate, end_date = enddate, bands = bands_touse,bad_pixels_limit =  PercMaxBadPixels)
  
  ## calculate vegetation indexes
  equations = do.call(c, lapply(veg_indexes, function(index) ListIndex[[index]]))
  
  vi_images = lapply(optical_images, function(sat_imag){
    calculate_vi(sat_imag,veg_indexes,equations)
  })
  return(vi_images)
}





## save rdata
tif_tordata = function(sat_imag,file_name,locality,
                       ouput_folder){
  # identify mission
  
  sub_imag = switch(str_sub(file_name , 1, 4) , 
                    "L2A_" = "L2A", 
                    "LC08" = "LC8",
                    "LE07" = "LE7",
                    "S1_2" = "S1_VV",
                    "S1_V" = "S1_VV")
  # create new file name
  file_name = paste0(sub_imag , "_" ,IdentiyImageDate(file_name) , "_" , tolower(locality) , ".RData")
  
  #export as RData
  save(sat_imag, file = paste0(setPathStyle(ouput_folder),file_name))
  cat(file_name, " was saved\n")
}

### read rada images

get_rada_images = function(inventory,init_date, end_date, bands = c("Sigma_VV_db"),
                           imagesto_delete = NA){
  
  
  # filter filenames
  init_date=CheckDateFormat(init_date)
  end_date=CheckDateFormat(end_date)
  inventory = inventory[CheckDateFormat(inventory$Date)%in%init_date:end_date,]
  ## identify polarization
  if(grepl(bands[1], pattern = "VV")){
    polarization = "VV"
  }else{
    polarization = "VH"
  }
  ##  s1 inventory
  inventory = inventory[inventory[,1] %in% paste0("S1_" , polarization),]
  
  # filter 
  inventory = inventory[!inventory$Delete %in% TRUE,]
  
  if(!is.na(imagesto_delete)){
    inventory = inventory[inventory$NameImage %in%imagesto_delete]
  }
  
  
  satellite_images = do.call(c,lapply(1:nrow(inventory), function(i){
    # configure satellite parameters
    sat_Params = get_sat_list_parameters(as.character(inventory$Mission_Type[i]))
    
    fname = do.call(c,lapply(inventory$NameImage[i], function(fname) paste0(fname,".", sat_Params$ExtImg)))
    
    ## set band names 

    satellite_image = Read_SatellitalImagery(fname, sat_Params, 
                                             filter_bands = bands)
    
    return(satellite_image)
  }))
  satellite_images = satellite_images[order(as.Date(IdentiyImageDate(names(satellite_images)), format = "%Y%m%d"))]
  
  return(satellite_images)
}

# read and filter optical images

get_optical_imagery = function(inventory,init_date, end_date, bands = c("red","nir"),bad_pixels_limit =100,remove_clouds = TRUE,
                               imagesto_delete = NA, optical_mission = c( "L2A","LC8","LE7","S2","L7","L8","LC08","LE07")){
  

  # filter filenames
  init_date=CheckDateFormat(init_date)
  end_date=CheckDateFormat(end_date)
  inventory = inventory[CheckDateFormat(inventory$Date)%in%init_date:end_date,]
  ## remove s1 inventory
  inventory = inventory[!inventory[,1] %in% c("S1_VV" , "S1_VH"),]
  ## filter by optical mission
  inventory = inventory[inventory[,1] %in% optical_mission,]
  
  # filter by cloud pixels and manually criteriums
  inventory = inventory[inventory$BadPixels<bad_pixels_limit,]
  inventory = inventory[!inventory$Delete %in% TRUE,]
  if(!is.na(imagesto_delete)){
    inventory = inventory[inventory$NameImage %in%imagesto_delete]
  }
  
  cloud_bands = c("QL",   "SLC")
  satellite_images = do.call(c,lapply(1:nrow(inventory), function(i){
    # configure satellite parameters
    sat_Params = get_sat_list_parameters(as.character(inventory$Mission_Type[i]))
    
    fname = do.call(c,lapply(inventory$NameImage[i], function(fname) paste0(fname,".", sat_Params$ExtImg)))
    
    
    
    # mask image
    if(remove_clouds){
      # get image with mask layers
      
      satellite_image = Read_SatellitalImagery(fname, sat_Params, 
                                               filter_bands = c(bands,cloud_bands))
      
      mask_sentence = get_mask_condition(sat_Params)
      ## remove clouds and shadows
      satellite_image = RemoveClouds(imag = satellite_image[[1]], condition = mask_sentence, numbands = length(bands))
      ## convert in list
      satellite_image = list(satellite_image)
      names(satellite_image) = fname
    }else{
      # get image without mask layers
      satellite_image = Read_SatellitalImagery(fname, sat_Params, 
                                               filter_bands = c(bands))
    }
    
    return(satellite_image)
  }))
  satellite_images = satellite_images[order(as.Date(IdentiyImageDate(names(satellite_images)), format = "%Y%m%d"))]
  
  return(satellite_images)
}

## 

get_mask_condition = function(sat_Params){
  if(sat_Params$Sat_Index=="LE07"){
    Limit=paste0("imag[['QL']][] > ",sat_Params$QualityLimits$QL)
    
  }else if(sat_Params$Sat_Index=="LC08"){
    Limit=paste0("imag[['QL']][] > ",sat_Params$QualityLimits$QL , 
                 " | imag[['SLC']][] %in% c(",paste0(sat_Params$QualityLimits$SLC,collapse = ","), ")")
  }else{
    Limit=paste0("imag[['QL']][] > ",sat_Params$QualityLimits$QL , 
                 " | imag[['SLC']][] %in% c(",paste0(sat_Params$QualityLimits$SLC,collapse = ","), ")")
  }
  return(Limit)
}

update_inventory = function(sat_Params, tile, start_date = "1990-01-01", end_date = Sys.Date(),
                                    inventory_path = setPathStyle(as.character(initconf[as.character(initconf[,1])%in%'inventory_folder',2]))){

  ## check dates format
  start_date = CheckDateFormat(start_date)
  end_date = CheckDateFormat(end_date)
  
  ## filter list files that correspond to the mission

  
  list_files = list.files(path = sat_Params$Dir_Path , pattern = paste0(sat_Params$ExtImg, "$"))
  list_missionfiles = list_files[grepl(list_files, pattern = sat_Params$Sat_Index)]
  list_missionfiles = list_missionfiles[CheckDateFormat(IdentiyImageDate(list_missionfiles))%in%start_date:end_date]
  
  # select new images to extract cloud percentage
  
  sat_images_names = sapply(list_missionfiles, function(fname)
    stringr::str_sub(fname,1,-(stringr::str_length(sat_Params$ExtImg)+2)))
  
  newimages = check_newimages(names_images = unname(sat_images_names), locality = tile )
  if(length(newimages)>0){
  
    ## this part is only for optical data: read the image and calculate the cloud percentage
    
    if(sat_Params$Sat_Index %in% c("L2A","LC08","LE07","S2","L8","L7")){
      # organice file names
      newimages = do.call(c,lapply(newimages, function(fname) paste0(fname,".", sat_Params$ExtImg)))
      ## read satellite data
      cloud_bands = c("QL",   "SLC")
      satellite_data = Read_SatellitalImagery(newimages, sat_Params, 
                                              filter_bands = cloud_bands)
      
      # get mask condition
      mask_condition = get_mask_condition(sat_Params)
      
      # get cloud pixels percentage
      cloud_percentages = lapply(satellite_data, function(imag){CountingPixelsWithClouds(imag,mask_condition)})
      
    }else{
      cloud_percentages = lapply(newimages, function(name)NA)
      names(cloud_percentages) = newimages
    }

    # addnew information to the inventory
    image_name = names(cloud_percentages)[1]
    lapply(names(cloud_percentages),function(image_name){
      quality_table = data.frame("Mission_Type" = sat_Params$Sat_Index, 
                                 "Date" = as.Date(IdentiyImageDate(image_name), format = "%Y%m%d"),
                                 "BadPixels" = cloud_percentages[[image_name]], NameImage = image_name, "Delete" = "")
      add_newdata(quality_table,tile)
      cat(image_name, "was added in the inventory","\n")
    })
  }else{
    cat("\nthere is no new image to add")
  }
  
 
 }

###### 3. Rice Classification

###-> Split DataSet


TransformAndBuffer=function(shapetoTransform,
                            crs_system = "+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0",
                            buffer_dist = -20){
  ## Buffer polygon
  planarSystem = spTransform(shapetoTransform,crs(crs_system )) 
  planarSystem_buffer=raster::buffer(planarSystem,buffer_dist,dissolve=F)
  
  return(planarSystem_buffer)
}



### --> Metrics


Extract_raster_data = function(raster_images , spatial_points , spectral_bands = ""){
  require(raster)
  if(!is.list(raster_images)){
    raster_images = list(raster_images)
  }
  return(lapply(1:length(raster_images), function(count){
    image_i = raster_images[[count]]
    if("" %in%spectral_bands){
      spectral_bands=names(image_i)
    }
    image_i=image_i[[spectral_bands]]
    
    ## compare coordinates system
    
    if(as.character(crs(spatial_points)) != as.character(image_i@crs)){
      spatial_points = spTransform(spatial_points, image_i@crs)
    }
    
    dataExtracted=data.frame(raster::extract(image_i,spatial_points))
    
    if(length(spectral_bands) == 1){
      names(dataExtracted) = spectral_bands
    }
    dataExtracted=cbind(spatial_points@data,dataExtracted)
    
    
    cat("The information in", names(raster_images)[count],"was extracted\n")
    return(dataExtracted)
  }))
}



## Divide the paramaters into windows


CalculateNumDivisions=function(dates_radar,numDiv=2){
  datapixel=data.frame(dates=dates_radar,
                       month=as.numeric(as.character(as.Date(dates_radar,format="%Y%m%d"),format="%m")))
  
  ## define groups
  count=1
  datanew=array()
  
  for(i in 2:(length(datapixel$month))){
    if(i-1 == 1) datanew[i-1]=count 
    if(datapixel$month[i]!=datapixel$month[i-1]){
      count=count+1
      datanew[i]=count
    }else{
      datanew[i]=count
    }
  }
  datapixel$groups=round((datanew+(datanew%%(12/numDiv)))/(12/numDiv))
  if(numDiv==1){datapixel$groups=1}
  return(datapixel)
}

### Models


SplitDataSet=function(mainDataSet,n.ite,split_number=0.7){
  output <- ncol(mainDataSet)
  setseed <- .Random.seed[1:n.ite]
  
  inTrain <- createDataPartition(y=mainDataSet[,output], p=0.7, list=F)
  training <- mainDataSet[inTrain,]
  testing <- mainDataSet[-inTrain,]
  return(list(training,testing))
  
}

ClassificationMetrics=function(conf_matrix){
  n = sum(conf_matrix) # number of instances
  nc = nrow(conf_matrix) # number of classes
  diag = diag(conf_matrix) # number of correctly classified instances per class 
  rowsums = apply(conf_matrix, 1, sum) # number of instances per class
  colsums = apply(conf_matrix, 2, sum) # number of predictions per class
  p = rowsums / n # distribution of instances over the actual classes
  q = colsums / n # distribution of instances over the predicted classes
  expAccuracy = sum(p*q)
  accuracy = sum(diag) / n
  kappa = (accuracy - expAccuracy) / (1 - expAccuracy)
  precision = diag / colsums 
  return(list(kappa,accuracy,precision))
}

SVMRadialForestModel=function(training , testing){
  
  
  ctrl <- trainControl(method="repeatedcv",   # 10fold cross validation
                       repeats=5,  	    # do 5 repititions of cv
                       # summaryFunction=twoClassSummary,	# Use AUC to pick the best model
                       classProbs=F)
  
  
  model <- train(x=training[,-ncol(training)],
                 y= training[,ncol(training)],
                 method = "svmRadial",   # Radial kernel
                 tuneLength = 9,					# 9 values of the cost function
                 preProc = c("center","scale"),  # Center and scale data
                 trControl=ctrl,
                 importance=T)
  varImportance=varImp(model)
  
  VarImportance=row.names(varImportance$importance[order(varImportance$importance[,1],decreasing = T),])[1:20]
  training[,ncol(training)]=factor(training[,ncol(training)])
  predictedValues = factor(as.character(predict(model, testing[,-ncol(testing)])))
  CM=table(testing[,ncol(testing)],predictedValues)
  
  Metrics=ClassificationMetrics(CM)
  
  return(list(model,Metrics[[1]],Metrics[[2]],Metrics[[3]],CM,training,testing,VarImportance))
}


RandomForestModel=function(training , testing){
  
  grid <- expand.grid(mtry=round((ncol(training)-1)/3))
  
  
  model <- train(training[,-ncol(training)], training[,ncol(training)],
                 method="rf", tuneGrid=grid,importance = TRUE,ntree = 3000)
  
  
  varImportance=varImp(model)
  
  VarImportance=row.names(varImportance$importance[order(varImportance$importance[,1],decreasing = T),])[1:20]
  training[,ncol(training)]=factor(training[,ncol(training)])
  predictedValues = factor(as.character(predict(model, testing[,-ncol(testing)])))
  CM=table(testing[,ncol(testing)],predictedValues)
  
  Metrics=ClassificationMetrics(CM)
  
  return(list(model,Metrics[[1]],Metrics[[2]],Metrics[[3]],CM,training,testing,VarImportance))
}


XgBoostingModel=function(training , testing){
  
  xgb_grid_1 = expand.grid(
    nrounds = 1000,
    eta = c(.1,0.01, 0.001),
    max_depth = c(2, 4, 6, 8),
    gamma = 1,
    colsample_bytree = c(.7,1),    #default=1
    min_child_weight = 1,
    subsample=c(.8, 1)
  )
  
  xgb_trcontrol_1 = trainControl(
    method = "cv",
    number = 5,
    verboseIter = TRUE,
    returnData = FALSE,
    returnResamp = "all",                                                        # save losses across all models
    classProbs = TRUE,                                                           # set to TRUE for AUC to be computed
    #summaryFunction = twoClassSummary,
    allowParallel = TRUE
  )
  
  model <- train(training[,-ncol(training)], training[,ncol(training)],
                 trControl = xgb_trcontrol_1,
                 tuneGrid = xgb_grid_1,
                 method = "xgbTree",importance = TRUE)
  
  
  varImportance=varImp(model)
  
  VarImportance=row.names(varImportance$importance[order(varImportance$importance[,1],decreasing = T),])[1:20]
  training[,ncol(training)]=factor(training[,ncol(training)])
  predictedValues = factor(as.character(predict(model, testing[,-ncol(testing)])))
  CM=table(testing[,ncol(testing)],predictedValues)
  
  Metrics=ClassificationMetrics(CM)
  
  return(list(model,Metrics[[1]],Metrics[[2]],Metrics[[3]],CM,training,testing,VarImportance))
}


#############

###### check inventory

get_inventory=  function(locality,
                         inventory_path = setPathStyle(as.character(initconf[as.character(initconf[,1])%in%'inventory_folder',2]))){
  
  ## create the inventory file name
  qtablename= paste0(setPathStyle(inventory_path),"Table_QualityPixels_",locality,".csv")
  # check the file existance
  if(file.exists(qtablename)){
    inventory = read.csv(qtablename)
    datadates = CheckDateFormat(inventory$Date)
    inventory$Date= datadates
  }
  ## create a inventory
  else{
    dir.create(inventory_path)
    # create columns and structure
    inventory = data.frame("Mission_Type" = "", 
                           "Date" = as.Date("1990/01/01", format = "%Y/%m/%d"),
                           "BadPixels" = 100 , NameImage = "first_element", "Delete" = TRUE)
    write.csv(inventory, qtablename, row.names = F)
    ## 
    cat("the inventory for ", locality,"was created")
    
  }
  return(inventory)
}

check_newimages = function(names_images,locality,
                           inventory_path = setPathStyle(as.character(initconf[as.character(initconf[,1])%in%'inventory_folder',2]))){
  inventory = get_inventory(locality)
  names_images = names_images[grepl(names_images , pattern = locality)]
  new_images = names_images
  ## check that images are not repeated
  img_dup = names_images%in% as.character(inventory$NameImage)
  if (length(which(img_dup))>0){
    cat("Following images are repeated:\n",names_images[img_dup],"\n")
    new_images = names_images[!img_dup]
  }
  return(new_images)
}

add_newdata = function(new_data,locality,
                       inventory_path = setPathStyle(as.character(initconf[as.character(initconf[,1])%in%'inventory_folder',2]))){
  
  inventory = get_inventory(locality)
  qtablename= paste0(setPathStyle(inventory_path),"Table_QualityPixels_",locality,".csv")
  # update table
  inventory = rbind(inventory,new_data)
  write.csv(inventory, qtablename, row.names = F)
}


cropRasters=function(rast_tomod, rasterRef,sps){

  newrast=crop(rast_tomod,extent(sps))
  newrast@extent=extent(sps)
  
  newrast=resample(newrast, rasterRef)
  return(newrast)
}


IdentiyImageDate=function(nameImage,patternDates="201",numPosition=7){
  posPattern=regexpr(patternDates, nameImage)
  return(str_sub(nameImage,posPattern,(numPosition+posPattern)))
}

replace_equation=function(equation, rast_varname = "sat_imag"){
  BandsToReplace = identify_index_bands(equation)
  functionIndex=equation
  for(coun in 1 : length(BandsToReplace)){
    functionIndex=gsub(BandsToReplace[coun], paste0(rast_varname,"$",BandsToReplace[coun]), functionIndex)
  }
  return(functionIndex)
  
}

PlotRGBImage=function(imageI,nameImage,OutputDirectoty,bands = c("red","green","blue")){
  PackageReading("RStoolbox")
  OutputDirectoty = setPathStyle(OutputDirectoty)
  nameOuput=paste0(OutputDirectoty,"RGB_",IdentiyImageDate(nameImage),"_",nameImage,".png")
  RStoolbox::ggRGB(imageI,r=bands[1], g=bands[2], b=bands[3],stretch = "lin", quantiles = c(0.2, 0.97))
  ggsave(nameOuput,width = 16,height = 8,units="cm")
}

Plot_Indexes= function(NamesIndexes,imageI,nameImage,OutputDirectoty){
  breaks = seq(0, 1, by=0.1)
  nb = length(breaks)-1 
  cols = rev(terrain.colors(nb))
  lapply(NamesIndexes,function(NameIndex){
    ImaIndex=imageI[[NameIndex]]
    nameOuput=paste0(OutputDirectoty,NameIndex,"_",IdentiyImageDate(nameImage),"_",nameImage,".png")
    png(filename=nameOuput,width = 1100,height = 950)
    plot(ImaIndex,breaks=breaks, zlim=c(0,1), col=cols)
    dev.off()
  })
}



Read_SatellitalImagery=function(filenames , sat_Params , init_date = "1990-08-20",
                                end_date = "2020-01-01", filter_bands = "all"){
  
  cat(sat_Params$Sat_Index, ": \n")
  patternExt=paste0(".",sat_Params$ExtImg,"$")
  
  ## Filter Images by dates
  ImagesDates=as.Date(IdentiyImageDate(filenames) , format = "%Y%m%d")
  
  init_date=CheckDateFormat(init_date)
  end_date=CheckDateFormat(end_date)
  
  filenames = filenames[ImagesDates %in% init_date:end_date]
  
  ## Read Images
  RasterImages=lapply(filenames,function(RastName){
    ImagesStack=stack(paste0(setPathStyle(sat_Params$Dir_Path),RastName))
    
    cat("Image Read: ",RastName[1],"\n")
    return(ImagesStack)
  })
  if(!sat_Params$Standar_bandNames[1] == ""){
    RasterImages=lapply(RasterImages,function(imag){names(imag)=sat_Params$Standar_bandNames
    return(imag)})
  }else{
    RasterImages=lapply(RasterImages,function(imag){names(imag)=sat_Params$Band_Names
    return(imag)})
  }
  
  # select only bands of interest
  if(filter_bands[1] != "all"){
    RasterImages=lapply(RasterImages,function(imag){imag[[which(names(imag)%in%filter_bands)]]})
    
    # this is only for sentinel -2 images, some images before of 2016 swifted their 
    # cloud probability  position
    if(length(which(filter_bands %in% c("QL", "SLC")))>1 & (sat_Params$Sat_Index == "L2A")){
      RasterImages = lapply(RasterImages,function(imag){
        if(max(imag[['QL']][],na.rm=T)<12 & max(imag[['SLC']][],na.rm=T)>10){
          pos_slc = grep(names(imag), pattern = "SLC")
          pos_ql =grep(names(imag), pattern = "QL")
          names(imag)[pos_slc] = ""
          names(imag)[pos_ql] = "SLC"
          names(imag)[pos_slc] = "QL"
        }
        return(imag)
        
      })
    }

  }
     
  # rename images
  posPattern=regexpr("[.]", filenames)
  names(RasterImages)=str_sub(filenames,1,(posPattern-1))
  return(RasterImages)
}

#### Create features Process
CountingPixelsWithClouds= function(imag,condition){
  CloudsPixels=which(eval(parse(text=paste0(condition))))
  percentBadPixels=length(CloudsPixels)/ncell(imag)*100
  return(percentBadPixels)
}

RemoveClouds = function(imag,condition=">9",numbands=NA){
  if (is.na(numbands)){
    numbands = length(names(imag))
  }
  
  ## mask images using condition
  
  CloudsPixels=which(eval(parse(text=paste0(condition))))

  return(stack(lapply(1:numbands,function(band){
    imgBand=imag[[band]]
    if(length(CloudsPixels)>0)imgBand[CloudsPixels]=NA
    cat("Band: ",names(imag)[band]," Was processed \n")
    return(imgBand)
  })))
}

rice_layerclassification = function(step_, levelsDiv, model_ML,dataImages, features_names,
                        model_variables){
  
  cat(step_ ,"\n")
  rowstoSelect=which(levelsDiv%in%unique(levelsDiv)[step_])
  
  # subset and join 
  dataToClass = dataImages[rowstoSelect, ]
  
  # remove na rows
  LevelsNA=which(rowSums(is.na(dataToClass))>0)
  if(length(LevelsNA)>0){
    dataToClass = dataToClass[-LevelsNA,]
    
  }
  
  # assign names
  dataToClass =data.frame(dataToClass)
  names(dataToClass) = features_names
  
  # Classify the database using a machine learning model
  
  Classify_ML=predict(model_ML,dataToClass[,model_variables[-length(model_variables)]])
  
  # export output
  final_class = data.frame(pixel = 1:length(rowstoSelect), class = 1)
  final_class[!final_class$pixel%in%LevelsNA,] = as.numeric(Classify_ML)
  
  return(final_class$class)}


######## Remove Noise

RemoveNoisefromaLayer=function(RiceLayer,MatrixSize,NumberMinumofPixels){
  ### Transform the raster information to a matrix
  
  dataLayerMatrix=as.matrix(RiceLayer)
  
  ### 
  
  if(!grepl(tolower(MatrixSize), pattern ="x")){
    stop("The matrix size doesn't have the correct format, please insert as follows '11 X 11'")
  }
  formatMatrix=unlist(strsplit(tolower(gsub(" ", "" ,formatMatrix)),split = "x"))
  N_rows=as.numeric(formatMatrix[1])
  N_Columns=as.numeric(formatMatrix[length(N_rows)])
  
  if(!N_rows%%2==0) N_rowsHalf=(N_rows-1)/2 else N_rowsHalf=(N_rows)/2
  if(!N_Columns%%2==0) N_ColumnsHalf=(N_Columns-1)/2 else N_ColumnsHalf=(N_Columns)/2
  
  MatrixModified=dataLayerMatrix
  
  RowLevelswithRice=apply(dataLayerMatrix,1,function(x){which(x==2)})
  
  for(i in (N_rowsHalf+1):(nrow(dataLayerMatrix)-N_rowsHalf)){
    for(j in unname(RowLevelswithRice[[i]])){
      if(j>(N_ColumnsHalf)&j<(ncol(dataLayerMatrix)-(N_ColumnsHalf-1))){
        if(sum((dataLayerMatrix[((i-N_rowsHalf):(i+N_rowsHalf)),((j-N_ColumnsHalf):(j+N_ColumnsHalf))])==2)<NumberMinumofPixels & dataLayerMatrix[i,j]==2){
          MatrixModified[i,j]=1
        }
      }
    }
    cat(i,"\n")
  }
  ImageClass=RiceLayer
  ImageClass[]=NA
  ImageClass[]=as.vector(t(MatrixModified))
  return(ImageClass)
  
}


### fill Empties



FillEmptiesLayer=function(RiceLayer,MatrixSize,NumberMinumofPixels){
  ### Transform the raster information to a matrix
  
  dataLayerMatrix=as.matrix(RiceLayer)
  
  ### Validate format
  
  if(!grepl(tolower(MatrixSize), pattern ="x")){
    stop("The matrix size doesn't have the correct format, please insert as follows '11 X 11'")
  }
  formatMatrix=unlist(strsplit(tolower(gsub(" ", "" ,formatMatrix)),split = "x"))
  N_rows=as.numeric(formatMatrix[1])
  N_Columns=as.numeric(formatMatrix[length(N_rows)])
  
  if(!N_rows%%2==0) N_rowsHalf=(N_rows-1)/2 else N_rowsHalf=(N_rows)/2
  if(!N_Columns%%2==0) N_ColumnsHalf=(N_Columns-1)/2 else N_ColumnsHalf=(N_Columns)/2
  
  MatrixModified=dataLayerMatrix
  
  RowLevelswithOutRice=apply(dataLayerMatrix,1,function(x){which(x==1)})
  
  for(i in (N_rowsHalf+1):(nrow(dataLayerMatrix)-N_rowsHalf)){
    for(j in unname(RowLevelswithOutRice[[i]])){
      if(j>(N_ColumnsHalf)&j<(ncol(dataLayerMatrix)-(N_ColumnsHalf-1))){
        if(sum((dataLayerMatrix[((i-N_rowsHalf):(i+N_rowsHalf)),((j-N_ColumnsHalf):(j+N_ColumnsHalf))])==2)>numberMaxPixels & dataLayerMatrix[i,j]==1){
          MatrixModified[i,j]=2
        }
      }
    }
    cat(i,"\n")
  }
  
  
  ImageClass=RiceLayer
  ImageClass[]=NA
  ImageClass[]=as.vector(t(MatrixModified))
  return(ImageClass)
  
}


########## Rice Grow Map


CheckDateFormat=function(DateToEvaluate){
  if(length(DateToEvaluate)>1){
    arrayaux=as.Date("1990-10-10" )
    count=1
    while(count<=length(DateToEvaluate)){
      arrayaux=c(arrayaux,CheckDateFormat(as.character(DateToEvaluate[count])))
      count=count+1
    }
    
    DateFormat=arrayaux[-1]
  }else{
    Possiblestrings=c("/","-"," ","[.]")
    formatDates=unlist(sapply(Possiblestrings,function(x){
      grep(DateToEvaluate,pattern = x)
    }))
    names(formatDates)
    ## evaluate the number of possible strings that are present in the date
    if(length(formatDates)== 0){
      DateFormat=as.Date(DateToEvaluate,format = "%Y%m%d")
      if(is.na(DateFormat)){
        DateFormat=as.Date(DateToEvaluate,format = "%Y%d%m")
        if(is.na(DateFormat))
          stop("Check the date format. it must be as follows: Year-month-day")
      }
      
    }else{
      splitDate=unlist(strsplit(as.character(DateToEvaluate),names(formatDates)))
      YearCond=which(as.numeric(splitDate)>1000)
      year=as.numeric(splitDate[YearCond])
      splitDate=as.numeric(splitDate[-ifelse(length(YearCond)>0,YearCond,3)])
      if(splitDate[1]>12){
        month=splitDate[2]
        day=splitDate[1]
      }else if(splitDate[2]>12){
        month=splitDate[1]
        day=splitDate[2]
      }else{
        month=splitDate[1]
        day=splitDate[2]
      }
      if(month<10)month=paste0('0',month)
      DateFormat=as.Date(paste0(year,month,day),format = "%Y%m%d")
      if(is.na(DateFormat)){
        stop("Check the date format. it must be as follows: Year-month-day")
      }
    }
   
  }
  
  return(DateFormat)
}


selectImages=function(TableQualityPixels,DateStart,DateEnd,PercMaxBadPixels){
  DateStart=CheckDateFormat(DateStart)
  DateEnd=CheckDateFormat(DateEnd)
  TableQualityPixels$Date=as.Date(CheckDateFormat(TableQualityPixels$Date))
  filter1=TableQualityPixels[TableQualityPixels$Date%in%DateStart:DateEnd & TableQualityPixels$BadPixels <PercMaxBadPixels,]
  return(as.character(filter1$NameImage)[!is.na(as.character(filter1$NameImage))])
}


#### get spatial poitns that will be used for extracting images information

get_spatialpoints = function(image_reference){
  
  #if there is no rice layer, the entire image will be processed
  
  LevelsWNA_FirstFilter=which(is.na(image_reference[]))
  
  # 
  points_info = SpatialPointsDataFrame(xyFromCell(image_reference,1:ncell(image_reference))[-LevelsWNA_FirstFilter,],
                                       data = data.frame(ID = c(1:ncell(image_reference))[-LevelsWNA_FirstFilter]),
                                       proj4string = crs(image_reference))
  return(points_info)
  
}


calculate_interpolation_series = function(images_data, date_int, veg_indexes,
                                          image_reference_path = NA, parallel_process = TRUE, 
                                          num_process = 5, summaryts = T){
  
  ## Set range time
  
  rangeDate=c(as.Date(date_int,format="%Y%m%d")-96,(as.Date(date_int,format="%Y%m%d")+0))
  nDaysPoints=length(rangeDate[1]:rangeDate[2])/16
  
  #### reorganice data
  
  if(is.na(image_reference_path)){
    image_reference =  images_data[[1]][[1]]
  }else{
    # if the file is geven it means that this is the rice layer
    image_reference =  raster(image_reference_path)
    image_reference[image_reference[] == 1]=NA
  }
  
  ## get spatial grid for extracting satellite data
  spatial_grid =  get_spatialpoints(image_reference)
  
  training_points = Extract_raster_data(images_data, spatial_grid, veg_indexes)
  
  ## merge data
  training_points = do.call(rbind,lapply(1:length(names(images_data)), function(image_index){
    data_perdate = training_points[[image_index]]
    #assign date
    data_perdate$Image_date = IdentiyImageDate(names(images_data)[image_index])
    return(data_perdate)
  }))
  
  ## change date format
  training_points$Image_date = as.Date(training_points$Image_date, format = "%Y%m%d")
  # reshape data
  data_reshape = lapply(veg_indexes , function(VI){
    
    ## filter by interest date
    data_WithOutDuplicated_date=
      training_points[training_points$Image_date %in% (as.Date(date_int,format="%Y%m%d")-130):as.Date(date_int,format=" %Y%m%d"),]
    
    
    data_WithOutDuplicated_date$Image_date = paste0(VI,"_", as.character(data_WithOutDuplicated_date$Image_date, format = "%Y%m%d"))
    ## reshape dataset
    
    data_reshape = reshape2::dcast( data_WithOutDuplicated_date , ID ~ Image_date,value.var = VI, fun= mean , na.rm =T)
    ## assign ID code
   
    return(data_reshape)
  })  
  
  ## merge
  data_reshape = plyr::join_all(data_reshape, by = "ID")
  
  ### remove those pixels with a high percentage of NA in the time series
  
  pixelsPercentageNA=data.frame(Pixel=data_reshape$ID,
                                Percentage=apply(data_reshape[,-which(names(data_reshape)%in% c("ID"))],
                                                 1,function(x){(sum(is.na(x))/length(x)*100)}))
  
  pixelsHighPercentage=pixelsPercentageNA[which(pixelsPercentageNA$Percentage>55),"Pixel"]
  
  if(length(pixelsHighPercentage)>0){
    data_reshape=data_reshape[!data_reshape$ID%in%as.character(pixelsHighPercentage),]
    
  }
  
  
  #### Smooth process
  
  rm(pixelsPercentageNA)
  rm(training_points)
  
  row.names(data_reshape) = data_reshape$ID
  
  KernelSmootValues=SmoothAllPixels(DateInt = date_int,
                                    ImagesInfo = data_reshape,
                                    Veg_Indexes = veg_indexes,
                                    rangeDate = rangeDate,nDaysPoints = nDaysPoints,
                                    parallelProcess = parallel_process,ncores = num_process)
  
  ######## Organize the information in a table
  
  TableValuesToClassify=do.call(cbind,lapply(1:length(veg_indexes),function(index){
    ValTableperIndex=data.frame(do.call(rbind,lapply(lapply(KernelSmootValues,function(x){x[[index]]}),
                                                     function(x){x[[2]]})))
    names(ValTableperIndex)=paste0(veg_indexes[index],"_Date_",1:length(names(ValTableperIndex)))
    return(ValTableperIndex)
  }))
  row.names(TableValuesToClassify)=row.names(data_reshape)
  
  
  ## get derivative values
  
  derivativeValues=do.call(cbind,lapply(1:1,function(index){
    ValTableperIndex=data.frame(do.call(rbind,lapply(lapply(KernelSmootValues,function(x){x[[index]]}),
                                                     function(x){x[[1]]})))
    names(ValTableperIndex)=paste0("NDVI_derivative_",1:length(names(ValTableperIndex)))
    return(ValTableperIndex)
  }))
  row.names(derivativeValues)=row.names(data_reshape)
  TableValuesToClassify = cbind(TableValuesToClassify,derivativeValues)
  
  ## get summaryvariables
  if (summaryts ){
    summaryValues=do.call(cbind,lapply(1:1,function(index){
      ValTableperIndex=data.frame(do.call(rbind,lapply(lapply(KernelSmootValues,function(x){x[[index]]}),
                                                       function(x){x$summary_ts})))
      names(ValTableperIndex)=paste0(veg_indexes[index],"_",names(ValTableperIndex))
      return(ValTableperIndex)
    }))
    row.names(summaryValues)=row.names(data_reshape)
    
    TableValuesToClassify = cbind(TableValuesToClassify,derivativeValues,summaryValues)
    
  }

  
  return(TableValuesToClassify)
}

#################
####-------------------------> Smooth process


timeseriesmetrics = function(data_PerPixel_VI,DateInt,rangeDate,nDaysPoints, crosscorr = T ){
  ## check dates order
  pixel_time = as.Date(IdentiyImageDate(names(data_PerPixel_VI)),format = "%Y%m%d")
  data_PerPixel_VI = data_PerPixel_VI[order(pixel_time)]
  
  ## remove NA pixels 
  namesColumns_withoutNA=names(data_PerPixel_VI)[!is.na(data_PerPixel_VI)]
  data_PerPixel_WithOutNA=data_PerPixel_VI[!is.na(data_PerPixel_VI)]
  
  ## get vi dates
  sat_dates=as.Date(IdentiyImageDate(namesColumns_withoutNA),format="%Y%m%d")
  ## identify last date
  last_date = sat_dates[length(sat_dates)]
  
  ## get date of interest intersection
  posDate=which(as.character(sat_dates)%in%as.character(as.Date(DateInt,format="%Y%m%d")))
  
  ## find in which position the date of interest is in the array
  if(length(posDate)==0){
    for(z in 2 : length(sat_dates)){
      if(sat_dates[z-1]<as.Date(DateInt,format="%Y%m%d") & sat_dates[z]>as.Date(DateInt,format="%Y%m%d")) 
        posDate=z
    }
    if(length(posDate)==0){
      posDate=length(sat_dates) 
    }
  }
  
  width_values = seq(20,70,1)
  ## set pixel value in 0 if the difference betweenm the reference date and the last info date is more than 20 days
  if(!((as.Date(DateInt,format="%Y%m%d") - last_date)>10 & max(diff(sat_dates))<=55)){
    
    data_tosmooth = data_PerPixel_WithOutNA[!is.na(data_PerPixel_WithOutNA[1:length(data_PerPixel_WithOutNA)])][1:posDate]
    bd = 36
    ## apply a 5 values window if there is more than 4 values
    if(length(data_tosmooth)>11){
      as_2 = signal::sgolayfilt(data_tosmooth, 5, 7)
      bd = 28
      
    }else if(length(data_tosmooth)>7){
      as_2 = signal::sgolayfilt(data_tosmooth, 4, 7)
      
    }else if(length(data_tosmooth)>4){
      as_2 = signal::sgolayfilt(data_tosmooth, 3, 5)
      
    }else{
      as_2 = signal::sgolayfilt(data_tosmooth, 2, 3)
    }
    df_to_fit = data.frame(id = row.names(data_PerPixel_VI),
                           y = as_2, x = as.numeric(as.Date(sat_dates)[1:posDate]))
    
    #36
    kernel_regres=ksmooth(df_to_fit$x, df_to_fit$y, "normal",
                          bandwidth = bd,n.points = nDaysPoints,range.x = rangeDate)
    # test
    regression_vals = kernel_regres$y
    # 
    plot(kernel_regres$x,kernel_regres$y, col = 'blue',ylim=c(0, 0.95),
         xlim= c(min(df_to_fit$x),max(kernel_regres$x)), pch = 16)

    points(df_to_fit$x,data_tosmooth, col="red")

    points(df_to_fit$x,df_to_fit$y, col = 'gray30')
    # points(kernel_regres_cross$x, kernel_regres_cross$y)
    # points(kernel_regres_cross$x, kernel_regres_cross$y, col ='green')

    
    ## ## cross correlation
    pos_torandom = (2:(length(df_to_fit$x)-1))
    newvals = sample(pos_torandom, round((length(df_to_fit$x)-2)*0.6))
    pos_training = c(1,newvals[order(newvals)],length(df_to_fit$x))
    pos_evaluation = c(1,pos_torandom[!pos_torandom%in%pos_training],length(df_to_fit$x))
    
    if(crosscorr){
      rsquare_cross = do.call(c,lapply(width_values,function(width){
        
        kernel_regres_cross=ksmooth(df_to_fit$x[pos_training], df_to_fit$y[pos_training], "normal",
                                    bandwidth = width,range.x = c(min(df_to_fit$x),max(df_to_fit$x)))
        predvalues = kernel_regres_cross$y[which(round(as.numeric(kernel_regres_cross$x))%in%
                                                   df_to_fit$x[pos_evaluation])]
        realvalues = df_to_fit$y[pos_evaluation]
        
        if(!T%in%is.na(predvalues) & !grepl(row.names(data_PerPixel_VI), pattern = "other")){
          rmse = caret::RMSE(predvalues,realvalues )
        }else{
          rmse = NA
        }
        return(rmse)
        
      }))
      rsquare_cross = (rsquare_cross-min(rsquare_cross))/(max(rsquare_cross)-min(rsquare_cross))
      
      
    }else{
      rsquare_cross = NA
    }
    
    ### step derivatives
    
    doi = as.Date(DateInt, "%Y%m%d")
    interval_90days = df_to_fit$x%in%doi:(doi-100)
    interval_40days = df_to_fit$x%in%doi:(doi-60)
    if((length(which(interval_90days))>3)){

      df_2condpoly = df_to_fit[interval_90days,]
      poly2grade_second = lm(y~ 1 + x + I(x^2)  , data = df_2condpoly)
      ## get coefficients : y = a + bx + cx^2
      coefficients_second = unname(poly2grade_second$coefficients)

      #df_toprodict = data.frame(x = as.numeric(kernel_regres$x[as.numeric(kernel_regres$x)>=min(df_2condpoly$x)]))
      df_toprodict = data.frame(x = as.numeric(kernel_regres$x))
      ## calculate second derivative last section: y' = b + 2c*x
      derivative_end =coefficients_second[2] + 2*coefficients_second[3]*max(df_toprodict$x)

      derivative_first =  coefficients_second[2] + 2*coefficients_second[3]*(df_toprodict$x[2])
      
      modelderi =c(derivative_first,derivative_end)
      ##test
      # modelderi <- sapply(df_toprodict$x, function(x)
      #   coefficients_second[2] + 2*coefficients_second[3]*x)
      # plot(df_toprodict$x, modelderi, col ='black')
      
      #### mean dif , sd and day of maximum
      
      #mean_val = mean(df_to_fit[interval_90days,"y"], na.rm = T)
      
      sd_val = sd(df_to_fit[interval_90days,"y"], na.rm = T)
      
      max_val = max(df_to_fit[interval_40days,"y"], na.rm = T)
      
      min_val = min(df_to_fit[interval_40days,"y"], na.rm = T)
      
      #dif_val = max_val - min_val
      
      df_to_fit40 = df_to_fit[interval_40days,]
      pos_max40daysbefore = df_to_fit[interval_40days,"y"] %in% max(df_to_fit[interval_40days,"y"], na.rm = T)
      
      pos_min40daysbefore = df_to_fit[interval_40days,"y"] %in% min(df_to_fit[interval_40days,"y"], na.rm = T)
      
      day_ofmaximum = abs(df_to_fit40[pos_max40daysbefore,"x"] - max(df_toprodict$x))[1]
      
      day_ofminimum = abs(df_to_fit40[pos_min40daysbefore,"x"] - max(df_toprodict$x))[1]
      
      summary_ts = data.frame(min_ts = min_val, 
                              max_ts = max_val,
                 sd_ts = sd_val,
                 maximumvi_day = day_ofmaximum,
                 minimumvi_day = day_ofminimum)
      
    }else{
      modelderi = rep(NA, 2)
      regression_vals= rep(NA, 7)
      summary_ts =  data.frame(min_ts = NA, 
                               max_ts = NA,
                               sd_ts = NA,
                               maximumvi_day = NA,
                               minimumvi_day = NA)
    }
    
    #interval_60days = df_to_fit$x%in%doi:(doi-75)
    # 
    # interval_70days = df_to_fit$x%in%(doi-120):(doi-35)
    # regression_vals = kernel_regres$y
    # 
    # if((length(which(interval_60days))>2)&(length(which(interval_70days))>2)){
    #   
    #   df_2condpoly = df_to_fit[interval_60days,]
    #   poly2grade_second = lm(y~ 1 + x + I(x^2), data = df_2condpoly) 
    #   ## get coefficients : y = a + bx + cx^2
    #   coefficients_second = unname(poly2grade_second$coefficients)
    #   
    #   #df_toprodict = data.frame(x = as.numeric(kernel_regres$x[as.numeric(kernel_regres$x)>=min(df_2condpoly$x)]))
    #   df_toprodict = data.frame(x = as.numeric(kernel_regres$x))
    #   ## calculate second derivative last section: y' = b + 2c*x
    #   derivative_end =coefficients_second[2] + 2*coefficients_second[3]*max(df_toprodict$x)
    #   
    #   ## middle derivative
    #   
    #   df_2condpoly = df_to_fit[interval_70days,]
    #   poly2grade_first = lm(y~ 1 + x + I(x^2), data = df_2condpoly) 
    #   coefficients_first = unname(poly2grade_first$coefficients)
    #   
    #   # ## calculate second derivative last section: y'' = 2c
    #   derivative_first =  coefficients_first[2] + 2*coefficients_first[3]*(df_toprodict$x[1])
    #   
    #   derivative_middlefirst =  coefficients_first[2] + 2*coefficients_first[3]*(df_toprodict$x[4])
    #   derivative_middlesecond = coefficients_second[2] + 2*coefficients_second[3]*(df_toprodict$x[4])
    #   
    #   modelderi =c(derivative_first,mean(c(derivative_middlefirst,derivative_middlesecond)),derivative_end)
    #   ##test 
    #   first_portion = data.frame(x= df_toprodict[1:4,])
    #   new_first = predict(poly2grade_first,first_portion)
    #   points(first_portion$x, new_first, col ='forestgreen', pch = 16)
    # 
    #   second_portion = data.frame(x= df_toprodict[4:7,])
    #   new_second = predict(poly2grade_second,second_portion)
    #   points(second_portion$x, new_second, col ='green', pch = 16)
    # 
    #   regression_vals = c(new_first[1:3],mean(c(new_second[1],new_first[4])),new_second[2:4])
    #   points(df_toprodict$x, regression_vals, col ='black', pch = 16)
    #   modelderi <- sapply(df_toprodict$x, function(x)
    #     coefficients_second[2] + 2*coefficients_second[3]*x)
    #   plot(df_toprodict$x, modelderi, col ='black')
    #   
    # }else{
    #   modelderi = rep(NA, 3)
    #   regression_vals= rep(NA, 7)
    # }
    
    ## results data frame
    
  }else{
    modelderi = rep(NA, 2)
    regression_vals= rep(NA, 7)
    rsquare_cross = rep(NA, length(width_values))
    df_to_fit = data.frame(id = NA,y = NA, x = NA)
    summary_ts = data.frame(min_ts = NA, 
                            max_ts = NA,
                            sd_ts = NA,
                            maximumvi_day = NA,
                            minimumvi_day = NA)
    
  }
  return(list(
    derivative_values = modelderi,
    regression_values = regression_vals,
    sg_smooth = df_to_fit,
    cross_corr_kernel = rsquare_cross,
    summary_ts = summary_ts))
}

SmoothAllPixels=function(DateInt,ImagesInfo,Veg_Indexes,rangeDate,nDaysPoints,parallelProcess=F,ncores=2,crosscorr = F){
  ## Stimate time
  if(parallelProcess){
    sfInit(parallel=T,cpus=ncores)
    sfLibrary("snowfall", character.only=TRUE)
    sfLibrary(stringr)
    sfLibrary(ff)

    
    
    sfExport("DateInt")
    sfExport("ImagesInfo")
    
    sfExport("timeseriesmetrics")
    sfExport("Veg_Indexes")
    sfExport("IdentiyImageDate")
    sfExport("rangeDate") 
    sfExport("nDaysPoints")
    sfExport("crosscorr")
    
    Sys.time()->start
    resultsKernel=sfLapply(1:nrow(ImagesInfo),function(j){
      
      ImagesInfo_PerPixel=ImagesInfo[j,]
      
      return(lapply(Veg_Indexes,function(VI){
        data_PerPixel_VI  = ImagesInfo_PerPixel[grepl(names(ImagesInfo_PerPixel), pattern =VI)]
        Smooth_VI = timeseriesmetrics(data_PerPixel_VI,DateInt,rangeDate,nDaysPoints)
        # graph_ = data.frame(x = Smooth_VI$x , y = Smooth_VI$y)
        # ggplot(graph_, aes(x , y ))+geom_point()+geom_line(alpha = 0.8)
        # 
        
        return(Smooth_VI)
      }))
    })
    
    
    print(Sys.time()-start)
    sfStop()
    
  }else{
    resultsKernel=lapply(1:nrow(ImagesInfo),function(j){
      #4805978
      which(row.names(ImagesInfo)%in% "3429023")
      ImagesInfo_PerPixel=ImagesInfo[j,]
      cat(paste("interation ", j))
      return(lapply(Veg_Indexes,function(VI){
        data_PerPixel_VI  = ImagesInfo_PerPixel[grepl(names(ImagesInfo_PerPixel), pattern =VI)]
        timeseriesmetrics(data_PerPixel_VI,DateInt,rangeDate,nDaysPoints,crosscorr= F)
      }))
    })
    
  }
  return(resultsKernel) 
}


######## Extract information e,ploying coordinates

ExtractInfo=function(SentinelimageReference,ImagesRiceZones,coordLevels,LevelsWNA_FirstFilter=0){
  if(as.character(is(ImagesRiceZones)[1])=="RasterBrick"){limitCount=nlayers(ImagesRiceZones)
  }else{limitCount=length(ImagesRiceZones)}
  return(lapply(1:length(ImagesRiceZones), function(count){
    coordRef=coordLevels
    if(length(LevelsWNA_FirstFilter)>1) coordRef=coordLevels[-LevelsWNA_FirstFilter,]
    nameImage=names(ImagesRiceZones)[count]
    imag=ImagesRiceZones[[count]]
    if(as.character(crs(imag))!=as.character(crs(SentinelimageReference))){
      spdf=SpatialPoints(coordRef,proj4string = crs(SentinelimageReference))
      coordRef = spTransform(spdf, crs(imag))
      coordRef=coordRef@coords
      cat("a change in the coordinates reference was done \n")
    }
    
    dataExtract=data.frame(raster::extract(imag,coordRef))
    names(dataExtract)=paste0(names(dataExtract),"_",IdentiyImageDate(nameImage))
    cat("The information in", nameImage,"was extracted\n")
    return(dataExtract)
  }))
}
