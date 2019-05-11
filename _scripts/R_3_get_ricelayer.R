################ 
### 
### This is a code, that was designed to Map Rice Zones for an spicific region
### To do this part, you had to select four Sentinel 2 images
###                   
###                   
###   Author:  Andres Aguilar
###     CIAT - DAPA - AEPS
#########


rm(list=ls())

##### Set Locality Name

tile = "col_t1"


## signal to use c("radar", "optical", "optical-radar")
signal_touse = "optical-radar"

######### Load functions


source(paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/_scripts/R_Main_Functions_GrowCropIndentification.R"))
source("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/_scripts/metrics_functions.R")

############ Load Libraries

libs=c("snowfall","caret","nnet","SDMTools","stringr","raster","ff", "rgeos","rgdal","dplyr" , "tidyverse","doParallel")
PackageReading(libs)

############## MAIN CODE

##### Set Dates

DateEnd="2018-09-12"
DateStart=CheckDateFormat(DateEnd) - 190




### Get features

### define signal to use
if(signal_touse == "optical-radar"){
  
  
  ## get metrics
  dataImages_optical =optical_metrics(tile, DateStart, DateEnd,bad_pixels_limit = 12)
  dataImages_radar = radar_metrics(tile, DateStart, DateEnd, radar_Bands = "db")

  dataImages= ff(vmode="double",dim=c(nrow(dataImages_optical$features_data), 
                                      ncol(dataImages_optical$features_data) + ncol(dataImages_radar$features_data)))
  
  for(col_i in 1:ncol(dataImages)){
    if(col_i <= ncol(dataImages_optical$features_data)){
      dataImages[,col_i]=dataImages_optical$features_data[,col_i]
    }else{
      dataImages[,col_i]=dataImages_radar$features_data[,col_i-ncol(dataImages_optical$features_data)]
    }
    
  } 

  features_names = c(dataImages_optical$feature_names,dataImages_radar$feature_names)
  
  satimage = dataImages_optical$satimage
  
}
if(signal_touse == "optical"){
  # get optical feautres
  dataImages_optical =optical_metrics(tile, DateStart, DateEnd)
  dataImages= dataImages_optical$features_data
  ##model file path
  model_fname = paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/sup_class_models/rice_identification/",ml_model,"_optical_metrics_4im_6months.RData")
  ## feature names
  features_names = dataImages_optical$feature_names
  ## set a sentinel 2 image as reference
  satimage = dataImages_optical$satimage
  
}
if(signal_touse == "radar"){
  
  dataImages_radar = radar_metrics(tile, DateStart, DateEnd, radar_Bands = "db")
  
  dataImages= dataImages_radar$features_data
  model_fname = paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/sup_class_models/rice_identification/",ml_model,"_radar_optical_6months.RData")
  ## feature names
  features_names = dataImages_radar$feature_names
  ## set a sentinel 1 image as reference
  satimage = dataImages_radar$satimage
  
}



####### Classify layer

### -- > Load Models
ml_model = "svmRadial"
ml_model = "rf"
ml_model = "xgbTree"
model_fname = paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/sup_class_models/rice_identification/",ml_model,"_radar_optical_onlydb_6months.RData")

load(file = paste0(model_fname))
allModels   <- lapply(ClassificationModels,function(x){x[[1]]})    
model_variables = names(ClassificationModels[[1]][[6]])
rm(ClassificationModels)

satimage[] =NA 
division = 20
levelsDiv=cut_number(1:ncell(satimage), division)

final_class = lapply(1:length(allModels), function(model_i){
  
  model_ML = allModels[[model_i]]
  
  ncores = 5
  cl =  makeCluster(ncores)
  registerDoParallel(cl)
  Sys.time()->start
  ML_class = do.call(c,foreach(step_ = 1:length(unique(levelsDiv)), 
                               .export = c("levelsDiv","rice_layerclassification","model_ML","dataImages",
                                           "features_names","model_variables")) %dopar% {
                                             
                                             library(ff)
                                             library(caret)
                                             
                                             ## select rows 
                                             rowstoSelect=which(levelsDiv%in%unique(levelsDiv)[step_])
                                             
                                             ## call function to classify the scene
                                             subset_class = rice_layerclassification(step_, levelsDiv, model_ML,
                                                                         dataImages, features_names,model_variables)
                                             # export
                                             subset_class
                                             
                                           })
  stopCluster(cl)
  doParallel::stopImplicitCluster()
  cat("\nmodel " , model_i , "finsihed time:", Sys.time()-start)
  
  satimage[] = ML_class
  plot(satimage)
  return(satimage)
})



final_class = stack(final_class)
## the final classification is stacked and changed its units to percentaje 
stack_maps = (((sum(final_class)/nlayers(final_class))*100)-100)

## export final layer

plot(stack_maps)
datingPeriod=paste0(gsub(x = as.character(DateStart),pattern = "-",replacement = ""),"_",
                    gsub(x = as.character(DateEnd),pattern = "-",replacement = ""))


filename=paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/model_ouputs/rice_layers/",tolower(tile),allModels[[1]]$method,"_",datingPeriod,"_",signal_touse,".tif")

writeRaster(stack_maps, filename=filename,format="GTiff",overwrite=T)

 
