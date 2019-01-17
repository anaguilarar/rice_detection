################ 
### 
### This is a code, that was designed to Map Rice Zones for an spicific region
### To do this part, you had to select four Sentinel 2 images
###                   
###                   
###   Author:  Andr?s Aguilar
###     CIAT - DAPA - AEPS
#########


rm(list=ls())

######## Functions
### Load functions
source("D:/phen_identification/_scripts/r/novembrer_2018/R_Main_Functions_GrowCropIndentification.R")
source("D:/phen_identification/_scripts/r/novembrer_2018/process_main_function.R")

############## MAIN CODE

##### Load Libraries
libs=c("snowfall","caret","nnet","SDMTools","stringr","raster","ff", "rgeos","rgdal","dplyr" , "tidyverse","parallel", "foreach","doParallel")
PackageReading(libs)
##### Set Locality Name

tile = "north_tolima"
## Period
DateStart="2017-04-04"
DateEnd="2017-10-21"
as.Date(DateStart) - as.Date(DateEnd)
Date_interval = as.Date(DateStart) : as.Date(DateEnd)
### Set folders path

Main_Folder = paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/")
optical_data = "D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/satellite_imagery/optical_data/"
radar_data = "D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/satellite_imagery/radar_data/"

################################
################# CLASIFICATION PROCESS
################# 
# use a raster image reference 



################################
################# READ OPTICAL DATA
#################

##### Read quality table

qpath = "D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/quality_info/"
# get inventory
inventory = get_inventory(qpath, tile )

### set vegetation indexes names

standarBand_names = c("blue","green","red","vrededg1","vrededg2","vrededg3","nir","narrownir","swir1","swir2")

######## Select Images for the study: there are two critriums 1) dates 2) quality and 3) visual 


### Set Criteriums

PercMaxBadPixels=38.3


###


########## ----- > Create Parametrics for optical images 


dataModel_Param_optical=stack(
    lapply (standarBand_names,function(band){
  
    band_info = stack(lapply(optical_images , function(x) x[[band]]))
    
    band_ifo = metrics_(band_info)
    
    return(band_ifo)
  })
)
rm(optical_images)

## transform optical data, extract data and storage in a ff data frame

dataImages_optical= ff(vmode="double",dim=c(ncell(dataModel_Param_optical), nlayers(dataModel_Param_optical)))

for(raster_layer in 1:nlayers(dataModel_Param_optical)){
  dataImages_optical[,raster_layer]=dataModel_Param_optical[[raster_layer]][]
}
rm(dataModel_Param_optical)
optical_names_database = do.call(c,lapply(standarBand_names, function(band) paste0(band ,c("_SD","_MEAN","_MIN","_MAX"))))

################################
################# READ RADAR DATA
#################

radar_Bands=c("Sigma0_VV_db","Sigma0_VV_Energy","Sigma0_VV_GLCMCorrelation","Sigma0_VV_GLCMVariance","Sigma0_VV_Contrast")
radar_Bands="Sigma0_VV_db"
#radar_Bands=c("Sigma0_VV_db","Sigma0_VV_Energy","Sigma0_VV_Entropy","Sigma0_VV_GLCMVariance","Sigma0_VV_Dissimilarity","Sigma0_VV_GLCMMean")

## list of files into the folder
list_radar_names = list.files(paste0(radar_data) , pattern = paste0(tile , ".RData$"))

## Filter files by dates

list_radar_names = list_radar_names[as.Date(IdentiyImageDate(list_radar_names), format = "%Y%m%d")%in% Date_interval]


## Load radar files

radar_images = lapply(1:length(list_radar_names), function(imag_index){
  file_name = list_radar_names[imag_index]
  
  load(file = paste0(  radar_data,file_name))
  
  bands_toSelect = names(sat_imag)[names(sat_imag)%in%radar_Bands]
  sat_imag = sat_imag[[bands_toSelect]]
  cat(str_sub(file_name ,1, -7) ," loaded\n")
  return(sat_imag)
})
names(radar_images) = list_radar_names

radar_images = radar_images[order(as.Date(IdentiyImageDate(names(radar_images)), format = "%Y%m%d"))]
sat_imageref = radar_images[[1]][[1]]


## Change data frame to a ff variable, to optimize space

ncores = 5
calculate_raster_quantiles = function(radar_images,radar_Bands,SentinelimageReference , division = 20, ncores = 5){
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

quantile_imags = calculate_raster_quantiles(radar_images,radar_Bands,sat_imageref, division = 25, ncores = 5)
rm(radar_images)
## database is transformed to ff format, which is ideal to optimize space
quantile_imags = stack(quantile_imags)

dataImages_radar= ff(vmode="double",dim=c(ncell(quantile_imags), nlayers(quantile_imags)))

for(raster_layer in 1:nlayers(quantile_imags)){
  dataImages_radar[,raster_layer]=quantile_imags[[raster_layer]][]
}
rm(quantile_imags)
radar_names_database = do.call(c,lapply(radar_Bands, function(band) paste0(band,"_pa_",1:length(c(0.05,.25, .50, .75,0.95)))))


## a parallel process was implemented, using a division 


dataImages= ff(vmode="double",dim=c(nrow(dataImages_optical), ncol(dataImages_optical) + ncol(dataImages_radar)))

for(col_i in 1:ncol(dataImages)){
  if(col_i <= ncol(dataImages_optical)){
    dataImages[,col_i]=dataImages_optical[,col_i]
  }else{
    dataImages[,col_i]=dataImages_radar[,col_i-ncol(dataImages_optical)]
  }
  
}


predict_tile = function(step_, levelsDiv, model_ML,dataImages, features_names,
                        model_variables){
  
  cat(step_ ,"\n")
  rowstoSelect=which(levelsDiv%in%unique(levelsDiv)[step_])
 
  # subset and join 
  dataToClass = dataImages[rowstoSelect, ]
  
  # remove na rows
  LevelsNA=which(rowSums(is.na(dataToClass))>0)
  dataToClass = dataToClass[-LevelsNA,]
  
  # assign names
  dataToClass =data.frame(dataToClass)
  names(dataToClass) = features_names
  
  # Classify the database using a machine learning model
  
  Classify_ML=predict(model_ML,dataToClass[,model_variables[-length(model_variables)]])
  
  # export output
  final_class = data.frame(pixel = 1:length(rowstoSelect), class = 1)
  final_class[!final_class$pixel%in%LevelsNA,] = as.numeric(Classify_ML)
  
  return(final_class$class)}


#### the data will be joined. the final dabase is classified by the best model.
signal_used = "optical"
library(doParallel)
library(foreach)

if(signal_used == "radar"){
  model_fname = "results/rice_mapping_models/rice_identification/svmRadial_radar_conf2_quantile_GLCM2_conf2_6months_V2.RData"
  dataImages_ = dataImages_radar
  sat_imageref = radar_images[[1]][[1]]
  features_names = radar_names_database
  
}else if(signal_used == "optical"){
  #model_fname = "results/rice_mapping_models/rice_identification/svmRadial_optical_metrics_4im_conf2_6months_V2.RData"
  model_fname = "D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/sup_class_models/rice_identification/xgbTree_optical_metrics_4im_6months.RData"
  dataImages_ = dataImages_optical
  load(paste0(optical_data, list_optical_names[1]))
  sat_imageref = sat_imag[[1]]
  rm(sat_imag)
  sat_imageref[] = NA
  features_names = optical_names_database
}else{
  model_fname = "results/rice_mapping_models/rice_identification/svmRadial_optical_metrics_conf2_6months_V2.RDataradar_optical_2conf2_6months_V2.RData"
  model_fname = "D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/sup_class_models/rice_identification/svmRadial_radar_optical_onlydb_6months.RData"
  #model_fname = "D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/sup_class_models/rice_identification/xgbTree_radar_optical_onlydb_6months.RData"
  
  dataImages_ = dataImages
  features_names = c(optical_names_database,radar_names_database)
}

### -- > Load Models
setwd("D:/phen_identification/")
load(file = paste0(model_fname))
allModels   <- lapply(ClassificationModels,function(x){x[[1]]})    
model_variables = names(ClassificationModels[[1]][[6]])
rm(ClassificationModels)
division = 20
levelsDiv=cut_number(1:ncell(sat_imageref), division)

final_class = lapply(1:length(allModels), function(model_i){
  
  model_ML = allModels[[model_i]]
  
  ncores = 5
  cl =  makeCluster(ncores)
  registerDoParallel(cl)
  Sys.time()->start
  ML_class = do.call(c,foreach(step_ = 1:length(unique(levelsDiv)), 
    .export = c("levelsDiv","predict_tile","model_ML","dataImages_",
                "features_names","model_variables")) %dopar% {
     
    library(ff)
    library(caret)
    
    ## select rows 
    rowstoSelect=which(levelsDiv%in%unique(levelsDiv)[step_])
    
    ## call function to classify the scene
    subset_class = predict_tile(step_, levelsDiv, model_ML,
                                dataImages_, features_names,model_variables)
    # export
    subset_class
    
    })
  stopCluster(cl)
  cat("\nmodel " , model_i , "finsihed time:", Sys.time()-start)

  sat_imageref[] = ML_class
  plot(sat_imageref)
  return(sat_imageref)
})

final_class = stack(final_class)
## the final classification is stacked and changed its units to percentaje 
stack_maps = (((sum(final_class)/nlayers(final_class))*100)-100)

## export final layer


datingPeriod=paste0(gsub(x = as.character(DateStart),pattern = "-",replacement = ""),"_",
                    gsub(x = as.character(DateEnd),pattern = "-",replacement = ""))


filename=paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/model_ouputs/rice_layers/",tolower(tile),allModels[[1]]$method,"_conf2_",datingPeriod,"_optical_db.tif")

writeRaster(stack_maps, filename=filename,format="GTiff",overwrite=T)
