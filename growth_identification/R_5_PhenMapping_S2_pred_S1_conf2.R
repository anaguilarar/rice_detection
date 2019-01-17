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
server=Sys.info()[1]
serverPath=switch(server, "Linux"="/mnt","Windows"="//dapadfs")
source("D:/phen_identification/_scripts/r/novembrer_2018/R_Main_Functions_GrowCropIndentification.R")



##### Load Libraries

libs=c("snowfall","caret","nnet","SDMTools","stringr","raster","ff", "rgeos","rgdal","dplyr" , "tidyverse")
PackageReading(libs)
##### Set Locality Name

tile = "north_tolima"

optical_data = "D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/satellite_imagery/optical_data/"
### Set folders path

Main_Folder = paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification")

setwd(Main_Folder)


#### Set the date where will be identified the grow stages



#### Load Rice Layer



rice_layer=raster::raster("model_ouputs/rice_layers/Filled_north_tolimaxgbTree_conf2_20170404_20171021_optical_db_15 x 15_100.tif")

#SentinelimageReference = imageI[[1]]
#rm(imageI)
plot(rice_layer)
LevelsWNA_FirstFilter=which(rice_layer[]==1)
#LevelsWNA_FirstFilter=which(is.na(SentinelimageReference[]))

# 
points_info = SpatialPointsDataFrame(xyFromCell(rice_layer,1:ncell(rice_layer))[-LevelsWNA_FirstFilter,],
                                     data = data.frame(ID = c(1:ncell(rice_layer))[-LevelsWNA_FirstFilter]),
                                     proj4string = crs(rice_layer))


#################
#################
#### ---------->> Create Parameters for optical images

##### get inventory
inventory = get_inventory("quality_info/", tile)


### set vegetation indexes names

indexes = c("NDVI")

######## Select Images for the study: there are two critriums 1) dates 2) quality and 3) visual 


DateInt="20171008"
## Period
DateStart=as.Date(DateInt, format = "%Y%m%d")-120
DateEnd=as.Date(DateInt, format = "%Y%m%d")

PercMaxBadPixels=38.4
#PercMaxBadPixels=17.3


optical_images = get_optical_imagery(optical_data, inventory, DateStart, DateEnd, PercMaxBadPixels,indexes)

######### ---> Extract data information using the Rice Layer

training_points = Extract_raster_data(raster_images = optical_images , spatial_points = points_info , 
                                      spectral_bands = indexes)

names(training_points) = sapply(names(optical_images ), function(name_imag) IdentiyImageDate(name_imag))


### Join all points 

training_points = do.call(rbind , lapply(1:length(training_points) , function(i){
  dataper_=training_points[[i]]
  dataper_$Image_date = names(training_points)[i]
  
  return(dataper_)
}))


training_points$Image_date = as.Date(training_points$Image_date , format = "%Y%m%d")
rm(optical_images)
######### ---> Kernel Smooth process

## Set parameters

rangeDate=c(as.Date(DateInt,format="%Y%m%d")-96,(as.Date(DateInt,format="%Y%m%d")+0))
nDaysPoints=length(rangeDate[1]:rangeDate[2])/16

############

VI = "NDVI"
ListIndex = c("NDVI")
data_reshape = lapply(ListIndex , function(VI){
  
  ## filter by interest date
  data_WithOutDuplicated_date=
    training_points[training_points$Image_date %in% (as.Date(DateInt,format=" %Y%m%d")-120):as.Date(DateInt,format=" %Y%m%d"),]
  
  
  data_WithOutDuplicated_date$Image_date = paste0(VI,"_", as.character(data_WithOutDuplicated_date$Image_date, format = "%Y%m%d"))
  ## reshape dataset
  
  data_reshape = reshape2::dcast( data_WithOutDuplicated_date , ID ~ Image_date, value.var = VI)
  ## assign ID code
  
  return(data_reshape)
})  


data_reshape = plyr::join_all(data_reshape, by = "ID")

### remove those pixels with high percentage of NA in the time series

pixelsPercentageNA=data.frame(Pixel=data_reshape$ID,
                              Percentage=apply(data_reshape[,-which(names(data_reshape)%in% c("ID"))],
                                               1,function(x){(sum(is.na(x))/length(x)*100)}))

pixelsHighPercentage=pixelsPercentageNA[which(pixelsPercentageNA$Percentage>50),"Pixel"]

if(length(pixelsHighPercentage)>0){
  imagesToProcess=data_reshape[!data_reshape$ID%in%as.character(pixelsHighPercentage),]
  
}else{
  imagesToProcess=data_reshape
}

###############

#### Smooth process

rm(pixelsPercentageNA)
rm(training_points)

row.names(imagesToProcess) = imagesToProcess$ID

KernelSmootValues=SmoothAllPixels(DateInt = DateInt,ImagesInfo = imagesToProcess,indexes,
                                  rangeDate = rangeDate,nDaysPoints = nDaysPoints,parallelProcess = T,ncores = 4)

#-- End smooth process



######## Organize the information in a table

TableValuesToClassify=do.call(cbind,lapply(1:length(ListIndex),function(index){
  ValTableperIndex=data.frame(do.call(rbind,lapply(lapply(KernelSmootValues,function(x){x[[index]]}),
                                                   function(x){x[[2]]})))
  names(ValTableperIndex)=paste0(ListIndex[index],"_Date_",1:length(names(ValTableperIndex)))
  return(ValTableperIndex)
}))
row.names(TableValuesToClassify)=row.names(imagesToProcess)


## get derivative values

derivativeValues=do.call(cbind,lapply(1:1,function(index){
  ValTableperIndex=data.frame(do.call(rbind,lapply(lapply(KernelSmootValues,function(x){x[[index]]}),
                                                   function(x){x[[1]]})))
  names(ValTableperIndex)=paste0(ListIndex[index],"_derivative_",1:length(names(ValTableperIndex)))
  return(ValTableperIndex)
}))
row.names(derivativeValues)=row.names(imagesToProcess)
TableValuesToClassify = cbind(TableValuesToClassify,derivativeValues)

rm(KernelSmootValues)

## check point
dir.create(paste0(setPathStyle(Main_Folder),"temp"))
write.csv(TableValuesToClassify, paste0(setPathStyle(Main_Folder),"temp/timeseriessmoothed_",tile,DateInt,"_v1.csv"))
#rm(TableValuesToClassify)
########
#----------------------------------
######## 

####### -------> Calculate radar metrics

### Set Radar Bands

radar_Bands=c("Sigma0_VV_db",
              "Sigma0_VV_Energy","Sigma0_VV_GLCMCorrelation","Sigma0_VV_GLCMVariance",
              "Sigma0_VV_MAX","Sigma0_VV_GLCMMean")


###filter those bands that were not used in the model

radar_Bands = radar_Bands[sapply(radar_Bands , function(band){
  T%in% (str_sub(namesVar ,1,str_length(band) ) %in%band )
})]

## list of files into the folder
list_radar_names = list.files(paste0(Main_Folder , "process/rdata/radar_info/") , pattern = paste0(Locality , ".RData$"))

## Filter files by dates

list_radar_names = list_radar_names[order(as.Date(IdentiyImageDate(list_radar_names), format = "%Y%m%d"))]

list_radar_names_used=list_radar_names[as.Date(IdentiyImageDate(list_radar_names), format = "%Y%m%d")%in%(as.Date(DateInt, format = "%Y%m%d"):(as.Date(DateInt, format = "%Y%m%d")-90))]

## create groups
Groups_radar = data.frame(dates =  IdentiyImageDate(list_radar_names_used) , groups = 1)
dates_radar = as.Date(Groups_radar$dates, format = "%Y%m%d")
limitsDates = (dates_radar[length(dates_radar)]-dates_radar[1])/3

for(i in 1:3){
  Groups_radar[dates_radar%in%(dates_radar[1]+(limitsDates*i)):(dates_radar[1]+limitsDates*(i-1)),"groups"]=i
}


## calculate parameters per periods
group_ref = 3
band =radar_Bands[1]

### Load data
radar_images_used = lapply(1:length(list_radar_names_used), function(imag_index){
  file_name = list_radar_names_used[imag_index]
  
  load(file = paste0(Main_Folder ,  "process/rdata/radar_info/",file_name))
  cat(str_sub(file_name ,1, -7) ," loaded\n")
  bands_toSelect = names(radar_info)[names(radar_info)%in%radar_Bands]
  radar_info = radar_info[[bands_toSelect]]
  return(radar_info)
})


names(radar_images_used) = list_radar_names_used
radar_images_used = radar_images_used[order(as.Date(IdentiyImageDate(names(radar_images_used)), format = "%Y%m%d"))]


### calculate parameter second configuration

dataModel_Param_radar=stack(lapply (radar_Bands,function(band){
  
  infoRaster = stack(lapply(radar_images_used , function(x) x[[band]]))
  
  meanValues=mean(infoRaster,na.rm=T)
  
  SumbandDiff=list()
  meanValues=mean(infoRaster,na.rm=T)
  for (j in 1:nlayers(infoRaster)){
    SumbandDiff[[j]]=(infoRaster[[j]]-meanValues)*(infoRaster[[j]]-meanValues)
  }
  SumbandDiff=sum(stack(SumbandDiff),na.rm=T)
  
  stdValues=((SumbandDiff*(1/(nlayers(infoRaster)-1)))^(1/2))
  minValues=min(infoRaster,na.rm=T)
  maxValues=max(infoRaster,na.rm=T)
  
  dataModel_Param_radar = stack(lapply(3:1, function(group_ref){
    radar_toCalc = infoRaster[[names(infoRaster)[which(Groups_radar$groups%in%group_ref)]]]
    meanValues=mean(radar_toCalc,na.rm=T)
    return(meanValues)
  }))
  cat(band, " processed \n")
  return(stack(dataModel_Param_radar, stdValues, meanValues, minValues, maxValues))
  
}))

rm(radar_images_used)

names(dataModel_Param_radar) =  unlist(lapply(radar_Bands, function(band){
  paste0(band , c("_MEAN_3","_MEAN_2","_MEAN_1","_SD","_MEAN","_MIN","_MAX"))
}))


radar_data_extracted = Extract_raster_data(raster_images = dataModel_Param_radar, spatial_points = points_info, 
                                           spectral_bands = names(dataModel_Param_radar) )[[1]]

rm(dataModel_Param_radar)

#write.csv(radar_data_extracted, paste0("process/temp/radarallglc_",DateInt,"_conf1.csv"))
#dataModel_Param_radar$diff_sgima_mean = dataModel_Param_radar$Sigma0_VV_db_MEAN_2 - dataModel_Param_radar$Sigma0_VV_db_MEAN_1
######################
#-- start classification process

TableValuesToClassify = read.csv(paste0("temp/timeseriessmoothed_",tile,DateInt,"_v1.csv"), row.names = 1)

######### ---> Join Data Base
TableValuesToClassify$ID = row.names(TableValuesToClassify)


#TableValuesToClassify = plyr::join_all(list(TableValuesToClassify, radar_data_extracted), by = "ID",match = "first")
TableValuesToClassify = TableValuesToClassify[,-which(names(TableValuesToClassify)%in%"ID")]
######### ---> Predict Values

LevelsWNA_SecondFilter=unique(unlist(sapply(1:ncol(TableValuesToClassify),function(columnVal){which(is.na(TableValuesToClassify[,columnVal]))})))
LevelsWNA_SecondFilter=LevelsWNA_SecondFilter[order(LevelsWNA_SecondFilter)]

if(length(LevelsWNA_SecondFilter)!=0){
  TableValuesToClassify_WithoutNA=TableValuesToClassify[-c(LevelsWNA_SecondFilter),]
}else{
  TableValuesToClassify_WithoutNA=TableValuesToClassify

}

## Organice database columns
TableValuesToClassify_WithoutNA=TableValuesToClassify_WithoutNA[,namesVar]
#######
#TableValuesToClassify_WithoutNA$diff_sgima_mean = TableValuesToClassify_WithoutNA$Sigma0_VV_db_MEAN_2 - TableValuesToClassify_WithoutNA$Sigma0_VV_db_MEAN_1

TableValuesToClassify_WithoutNA_ff=ff(vmode="double",dim=dim(TableValuesToClassify_WithoutNA))

for(Band in 1:ncol(TableValuesToClassify_WithoutNA)){
  TableValuesToClassify_WithoutNA_ff[,Band]=TableValuesToClassify_WithoutNA[,Band]
  cat("Done \n")
}


#######
########## Load rice Grow Identification Model


load(file="sup_class_models/phen_identification/rf_conf_8optical_ndvi_deri.RData")

allModels=lapply(ClassificationModels,function(x){x[[1]]})

namesVar=names(ClassificationModels[[1]][[6]])
namesVar=namesVar[-length(namesVar)]
####


datatoPredict=data.frame(TableValuesToClassify_WithoutNA_ff[1:nrow(TableValuesToClassify_WithoutNA_ff),])
names(datatoPredict)=names(TableValuesToClassify)
rm(TableValuesToClassify)

lapply(ClassificationModels, function(x)x[[2]])
model_ML=allModels[[3]]

PheStagePre=predict(model_ML,datatoPredict[,namesVar])

levels(PheStagePre)=1:length(levels(PheStagePre))
PheStagePre=as.numeric(as.character(PheStagePre))

testImage=rice_layer
testImage[]=NA

testImage[as.numeric(row.names(TableValuesToClassify_WithoutNA))]=PheStagePre


plot(testImage)

dir.create("model_ouputs/growth_stages")
modelMethod=allModels[[1]]$method

writeRaster(testImage, filename=paste0("model_ouputs/growth_stages/",tile,"_",modelMethod,"_conf8_",DateInt,"ndvi.tif"), format="GTiff",overwrite=T)





