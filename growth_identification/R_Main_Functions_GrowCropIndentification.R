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


#### 1. Pre - Processing - Landsat

MaskLandsatImages=function(dirRaster,extentRef,NamesBands){
  tifFiles=list.files(path = dirRaster, pattern = "*.tif$")
  tifFiles=tifFiles[unlist(sapply(NamesBands,function(i){grep(tifFiles,pattern = i)}))]
  dirRaster=setPathStyle(dirRaster)
  ReadFiles=stack(paste0(dirRaster,tifFiles))
  cropRaster=crop(ReadFiles,extentRef)
  names(cropRaster)=NamesBands
  return(cropRaster)
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

calculate_vi = function(imageI , list_indices){
  NewLayers=lapply(list_indices,function(x){
    eval(parse(text = Identiffy_Band(names(imageI),x)))})
  imageI=stack(c(imageI,NewLayers))
  return(imageI)
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


get_optical_imagery = function(optical_path, inventory, DateStart, DateEnd, PercMaxBadPixels,band_names, select_s2 = FALSE){
  optical_images = selectImages(TableQualityPixels = inventory,DateStart,DateEnd,PercMaxBadPixels)
  
  ### delete images with bad quality
  ImagestoDelete=as.character(inventory[inventory$Delete%in%T,4])
  if(length(which(optical_images%in%ImagestoDelete))>0) optical_images=optical_images[-which(optical_images%in%ImagestoDelete)]
  
  ### read optical imagery information
  
  ## list of files into the folder
  
  list_optical_names = list.files(paste0( 
    optical_path) , pattern = paste0(tile , ".RData$"))
  
  if(select_s2){
    ## Filter only Sentinel 2 imagery
    optical_images = optical_images[grepl(optical_images, pattern ="L2A")]
    
  }
  
  # change names
  sub_imag = sapply(optical_images , function(imag){
    switch(str_sub(imag , 1, 4) , 
           "L2A_" = "L2A", 
           "LC08" = "LC8",
           "LE07" = "LE7")})
  
  images_names = paste0(sub_imag , "_" ,IdentiyImageDate(optical_images) , "_" , tolower(tile) , ".RData")
  
  # filter file names
  list_optical_names = list_optical_names[list_optical_names %in% images_names]
  
  # delete manually
  #list_optical_names = list_optical_names[!list_optical_names%in%"S2_20151211_saldana.RData"]
  
  # load RData
  optical_images = lapply(1:length(list_optical_names), function(imag_index){
    file_name = list_optical_names[imag_index]
    
    load(file = paste0( optical_path,file_name))
    cat(str_sub(file_name ,1, -7) ," loaded\n")
    bands_toSelect = names(sat_imag)[names(sat_imag)%in%band_names]
    sat_imag = sat_imag[[bands_toSelect]]
    return(sat_imag)
  })
  names(optical_images) = list_optical_names
  # order optical imagery
  optical_images = optical_images[order(as.Date(IdentiyImageDate(names(optical_images)), format = "%Y%m%d"))]
  return(optical_images)
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

get_inventory=  function(inventory_path, locality){
  
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

check_newimages = function(names_images,inventory_path,locality){
  inventory = get_inventory(inventory_path,locality)
  new_images = names_images
  ## check that images are not repeated
  img_dup = names_images%in% as.character(inventory$NameImage)
  if (length(which(img_dup))>0){
    cat("Following images are repeated:\n",names_images[img_dup],"\n")
    new_images = names_images[!img_dup]
  }
  return(new_images)
}

add_newdata = function(new_data,inventory_path,locality){
  
  inventory = get_inventory(inventory_path, locality)
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

export_rdata_images = function(Satellital_imagery,index,output_folder,sub_imag){}

replaceMultypleChar=function(BandToReplace,Func){
  return(gsub(BandToReplace, paste0("imageI$",BandToReplace), Func))
}

IdentiyImageDate=function(nameImage,patternDates="201",numPosition=7){
  posPattern=regexpr(patternDates, nameImage)
  return(str_sub(nameImage,posPattern,(numPosition+posPattern)))
}

CountingPixelsWithClouds= function(imag,NameBandCloud="CloudQA1",Limit=">9"){
  CloudsPixels=which(eval(parse(text=paste0("imag[[NameBandCloud]][]",Limit))))
  percentBadPixels=length(CloudsPixels)/ncell(imag)*100
  return(percentBadPixels)
}

Identiffy_Band=function(Standar_bandNames,Equation){
  BandsToReplace=Standar_bandNames[sapply(Standar_bandNames,function(band){grepl(Equation,pattern=band)})]
  functionIndex=Equation
  for(coun in 1 : length(BandsToReplace)){
    functionIndex=replaceMultypleChar(BandsToReplace[coun],functionIndex)
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


Read_SatellitalImagery=function(sat_Params , init_date = "1990-08-20",
                                end_date = "2020-01-01", tile = NA){
  cat(sat_Params$Sat_Index, ": \n")
  patternExt=paste0(".",sat_Params$ExtImg,"$")
  
  ## List file names
  names(ListRasterNames)=ListRasterNames=list.files(path=sat_Params$Dir_Path,pattern=patternExt)
  ## second filter by mission pattern
  ListRasterNames = ListRasterNames[grep(ListRasterNames, pattern = sat_Params$Sat_Index)]
  ## consider a third fileter when the tile is given
  if(!is.na(tile)){
    ListRasterNames = ListRasterNames[grep(ListRasterNames, pattern = tile)]
  }
  
  ## There are special conditions for S2 and S1, due to the data are into folders
  
  #if(sat_Params$Sat_Index == "S2"){ListRasterNames=SentinelL2A_Process(sat_Params)}
  #if(sat_Params$Sat_Index == "S1"){ListRasterNames=SentinelL1A_Process(sat_Params)}
  
  
  ## Filter Images by dates
  ImagesDates=as.Date(IdentiyImageDate(ListRasterNames) , format = "%Y%m%d")
  
  init_date=CheckDateFormat(init_date)
  end_date=CheckDateFormat(end_date)
  
  ListRasterNames = ListRasterNames[ImagesDates %in% init_date:end_date]
  
  ## Read Images
  RasterImages=lapply(ListRasterNames,function(RastName){
    ImagesStack=stack(paste0(setPathStyle(sat_Params$Dir_Path),RastName))
    
    cat("Image Readed: ",RastName[1],"\n")
    return(ImagesStack)
  })
  if(!sat_Params$Standar_bandNames[1] == ""){
    RasterImages=lapply(RasterImages,function(imag){names(imag)=sat_Params$Standar_bandNames
    return(imag)})
  }else{
    RasterImages=lapply(RasterImages,function(imag){names(imag)=sat_Params$Band_Names
    return(imag)})
    }
  posPattern=regexpr("[.]", names(ListRasterNames))
  names(ListRasterNames)=str_sub(names(ListRasterNames),1,(posPattern-1))
  names(RasterImages)=names(ListRasterNames)
  return(RasterImages)
}
#### Create features Process

RemoveClouds = function(imag,condition=">9",numbands=10){
  #CloudsPixels=which(eval(parse(text=paste0("imag[[NameBandCloud]][]",Limit))))
  ## mask images using condition
  
  CloudsPixels=which(eval(parse(text=paste0(condition))))
  percentBadPixels=length(CloudsPixels)/ncell(imag[[1]])*100
  return(list(percentBadPixels,stack(lapply(1:numbands,function(band){
    imgBand=imag[[band]]
    if(length(CloudsPixels)>0)imgBand[CloudsPixels]=NA
    cat("Band: ",names(imag)[band]," Was processed \n")
    return(imgBand)
  }))))
}

PredictPixelsClass=function(step_){
  rowstoSelect=which(levelsDiv%in%unique(levelsDiv)[step_])
  ImageProcess=data.frame(dataImages[rowstoSelect,])
  names(ImageProcess)=ColumnNames
  DataDef=ImageProcess[,names(ImageProcess)%in%VariablesUsed[-length(VariablesUsed)]]
  
  levelsData=1:nrow(DataDef)
  LevelsNA=which(rowSums(is.na(DataDef))>0)
  LevelsWNA=levelsData[!levelsData%in%LevelsNA]
  DataDef=DataDef[,VariablesUsed[-length(VariablesUsed)]]
  Classify_ML=predict(model_ML,DataDef[LevelsWNA,])
  levels(Classify_ML)=1:length(levels(Classify_ML))
  Classify_ML=as.numeric(as.character(Classify_ML))
  
  dfPredict=data.frame(PixelN=levelsData,value=1)
  
  dfPredict[LevelsWNA,2]=Classify_ML
  
  return(dfPredict$value)
  cat(paste0(step_, " to ",division,"\n"))
}

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
      stop("Check the date format. it must have as follows: Year-month-day")
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



#################
####-------------------------> Smooth process



smoothPerIndex=function(data_PerPixel_VI,DateInt,rangeDate,nDaysPoints ){
  
  namesColumns=names(data_PerPixel_VI)
  namesColumns_withoutNA=namesColumns[!is.na(data_PerPixel_VI)]
  data_PerPixel_WithOutNA=data_PerPixel_VI[!is.na(data_PerPixel_VI)]
  
  difDays=as.Date(IdentiyImageDate(namesColumns_withoutNA),format="%Y%m%d")
  infoPerPixel_Filled=data.frame(t(data.frame(row.names=c(namesColumns_withoutNA),vals=c(data_PerPixel_WithOutNA)) ))
  infoPerPixel_Filled=infoPerPixel_Filled[,order(names(infoPerPixel_Filled))]
  
  # order the names by dates
  namesColumns=names(infoPerPixel_Filled)
  namesColumns=namesColumns[order(as.Date(IdentiyImageDate(namesColumns),format = "%Y%m%d"))]
  infoIm_perIndex=infoPerPixel_Filled[,namesColumns]
  
  datesIm=as.Date(IdentiyImageDate(names(infoIm_perIndex)),format="%Y%m%d")
  datesImAux=datesIm[!is.na(infoIm_perIndex)]
  
  datesImAux=datesImAux[order(datesImAux)]
  posDate=which(as.character(datesImAux)%in%as.character(as.Date(DateInt,format="%Y%m%d")))
  
  if(length(posDate)==0){
    for(z in 2 : length(datesImAux)){
      if(datesImAux[z-1]<as.Date(DateInt,format="%Y%m%d") & datesImAux[z]>as.Date(DateInt,format="%Y%m%d")) 
        posDate=z
    }
    if(length(posDate)==0){
      posDate=length(datesImAux) 
      #infoIm_perIndex[!is.na(infoIm_perIndex[1:length(infoIm_perIndex)])]=0
    }
  }
  
  #sgolayfilt(x, p = 2, n = p + 3 )
  #33
  data_tosmooth = infoIm_perIndex[!is.na(infoIm_perIndex[1:length(infoIm_perIndex)])][1:posDate]
  if(length(data_tosmooth)>4){
    as_2 = signal::sgolayfilt(data_tosmooth, 3, 5)
    
  }else{
    as_2 = signal::sgolayfilt(data_tosmooth, 2, 3)
  }
  df_to_fit = data.frame(y = as_2, x = as.numeric(as.Date(datesImAux)[1:posDate]))
  
  Smooth_VI=ksmooth(df_to_fit$x, df_to_fit$y, "normal",
                    bandwidth = 25,n.points = nDaysPoints,range.x = rangeDate)

  #splpredict = spline(x = df_to_fit$x, y=df_to_fit$y, method="natural", xout = Smooth_VI$x)
  model <- lm(y ~ x + I(x^2)  + I(x^3),data = df_to_fit)
  df_to_predict = data.frame(x = as.numeric(Smooth_VI$x))
  new_s = predict(model, newdata = df_to_predict)
  modelderi <- sapply(as.numeric(Smooth_VI$x), function(x)model$coefficients[2] + 2*model$coefficients[3]*x)
  
  
  #ant 33
  #
  return(list(c(modelderi[1] , modelderi[length(modelderi)]),Smooth_VI$y))
}

SmoothAllPixels=function(DateInt,ImagesInfo,Veg_Indexes,rangeDate,nDaysPoints,parallelProcess=F,ncores=2){
  ## Stimate time
  if(parallelProcess){
    sfInit(parallel=T,cpus=ncores)
    sfLibrary("snowfall", character.only=TRUE)
    sfLibrary(stringr)
    sfLibrary(ff)

    
    
    sfExport("DateInt")
    sfExport("ImagesInfo")
    
    sfExport("smoothPerIndex")
    sfExport("Veg_Indexes")
    sfExport("IdentiyImageDate")
    sfExport("rangeDate") 
    sfExport("nDaysPoints")
    
    Sys.time()->start
    resultsKernel=sfLapply(1:nrow(ImagesInfo),function(j){
      
      ImagesInfo_PerPixel=ImagesInfo[j,]
      
      return(lapply(Veg_Indexes,function(VI){
        data_PerPixel_VI  = ImagesInfo_PerPixel[grepl(names(ImagesInfo_PerPixel), pattern =VI)]
        Smooth_VI = smoothPerIndex(data_PerPixel_VI,DateInt,rangeDate,nDaysPoints)
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
      ImagesInfo_PerPixel=ImagesInfo[j,]
      cat(paste("interation ", j))
      return(lapply(Veg_Indexes,function(VI){
        data_PerPixel_VI  = ImagesInfo_PerPixel[grepl(names(ImagesInfo_PerPixel), pattern =VI)]
        smoothPerIndex(data_PerPixel_VI,DateInt,rangeDate,nDaysPoints)
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
