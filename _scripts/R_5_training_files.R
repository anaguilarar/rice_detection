################ 
### 
###   Create data training and validatinfg for clasification modeling
### 
###                   
###                   
###   Author:  Andr?s Aguilar
###     CIAT - DAPA - AEPS
#########


rm(list=ls())
count = 2
######## Functions
### Load functions
server=Sys.info()[1]

Main_Folder = paste0("D:/phen_identification/")

setwd(Main_Folder)

serverPath=switch(server, "Linux"="/mnt","Windows"="//dapadfs")

source(paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/_scripts/R_Main_Functions_GrowCropIndentification.R"))


source("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/_scripts/process_main_function.R")


##### Load Libraries

libs=c("snowfall","caret","nnet","SDMTools","stringr","raster","ff", "rgeos","rgdal","dplyr" , "foreach","doParallel")
PackageReading(libs)
##### Set Locality Name

tile = "col_t3"
### Set Criteriums

## Period
DateStart="2015-07-30"
DateEnd="2016-01-16"

inventory = get_inventory( tile)

vi_images = get_vi_layers(inventory = inventory , veg_indexes = c("NDVI","LSWI"), DateStart, DateEnd,PercMaxBadPixels = 34)

plot(vi_images$LE07_20150925_col_t3_10m.tif)

## set radar band names
# 

### Read Spatial data for one iteration
setwd("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/")
crossval_iter = 3

typeData = "training"

typeData = "validation"
data_p= c("training", "validation")

training_points_ = do.call(rbind,lapply(data_p, function(typeData){
  training_points_ = do.call(rbind,lapply(1:3 , function(crossval_iter){
    training_points = readOGR(paste0("model_inputs/phen_identification/partitions/",typeData,"_phenIdenti/v8",typeData,"_iterationCV_", crossval_iter ,".shp"))
    
    training_points@data$iteration = crossval_iter
    
    ### extract raster data
    # training_points = Extract_raster_data(raster_images = optical_images , spatial_points = training_points , 
    #                                       spectral_bands = c("NDVI" , "LSWI" , "BSI"))
    return(training_points@data)
    
  }))
  training_points_$data_type = typeData
  return(training_points_)
}))



clean_ts_function = function(data_WithOutDuplicated, Typ_Stg=NULL, crossval_iter=NULL){
  ############### Clean Process
  ###### Graphics
  data_WithOutDuplicated$ID = with(data_WithOutDuplicated,paste(x , y , plygn_n,  Date , Typ_Stg ,data_type , "i", iteration, sep = "_"))
  
  Date_Int = "20151228"

  
  ##### Vegetative Ending
  
  training_points_veg = data_WithOutDuplicated[ data_WithOutDuplicated$Typ_Stg %in% c("late_vegetative"), ]
  
  training_points_veg_date = training_points_veg[ training_points_veg$Date %in% Date_Int, ]
  
  #ggplot(training_points_veg_date , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2015-12-11") & training_points_veg_date$NDVI <0.45)
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = training_points_veg_date$ID[wrongLevels]
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2016-01-10") & training_points_veg_date$NDVI <0.55)
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = c(wronglevels_ID,training_points_veg_date$ID[wrongLevels])
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2015-10-19") & training_points_veg_date$NDVI >0.45)
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = c(wronglevels_ID,training_points_veg_date$ID[wrongLevels])
  
  wrongLevels = which(training_points_veg_date$ID%in%wronglevels_ID)
  #m=ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  #ggsave(plot = m, filename = paste0("process/exploratory_images/time_Series/raw_data/plot_",Date_Int,"Ending",crossval_iter,"_", typeData,".png"))
  
  #write.csv(training_points_veg_date[-wrongLevels,], paste0("process/exploratory_images/time_Series/raw_data/data_",Date_Int,"Ending",crossval_iter,"_", typeData,".csv"))
  
  ## remove bad registers
  data_WithOutDuplicated=data_WithOutDuplicated[!data_WithOutDuplicated$ID %in% wronglevels_ID,]
  
  
  ###### early_vegetative
  training_points_veg = data_WithOutDuplicated[ data_WithOutDuplicated$Typ_Stg %in% c("early_vegetative"), ]
  
  training_points_veg_date = training_points_veg[ training_points_veg$Date %in% Date_Int, ]
  #ggplot(training_points_veg_date , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2015-10-22") & (training_points_veg_date$NDVI >0.65 |
                                                                                      training_points_veg_date$NDVI <0.4))
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = training_points_veg_date$ID[wrongLevels]
  
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2015-10-19") & (training_points_veg_date$NDVI >0.65 |
                                                                                      training_points_veg_date$NDVI <0.4))
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = c(wronglevels_ID,training_points_veg_date$ID[wrongLevels])
  
  wrongLevels = which(training_points_veg_date$ID%in%wronglevels_ID)
  m=ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  #ggsave(plot = m, filename = paste0("process/exploratory_images/time_Series/raw_data/plot_",Date_Int,"_early_vegetative_",crossval_iter,"_", typeData,".png"))
  
  #write.csv(training_points_veg_date[-wrongLevels,], paste0("process/exploratory_images/time_Series/raw_data/data_",Date_Int,"_early_vegetative_",crossval_iter,"_", typeData,".csv"))
  ## remove bad registers
  data_WithOutDuplicated=data_WithOutDuplicated[!data_WithOutDuplicated$ID %in% wronglevels_ID,]
  
  ###### harvested
  training_points_veg = data_WithOutDuplicated[ data_WithOutDuplicated$Typ_Stg %in% c("harvested"), ]
  
  training_points_veg_date = training_points_veg[ training_points_veg$Date %in% Date_Int, ]
  #ggplot(training_points_veg_date , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2015-10-22") & (training_points_veg_date$NDVI <0.5))
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = training_points_veg_date$ID[wrongLevels]
  
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2015-11-20") & (training_points_veg_date$NDVI <0.55))
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = c(wronglevels_ID,training_points_veg_date$ID[wrongLevels])
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2015-09-01") & (training_points_veg_date$LSWI >0.51))
  #ggplot(training_points_veg_date[!training_points_veg_date$ID%in%wronglevels_ID,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = c(wronglevels_ID, training_points_veg_date$ID[wrongLevels])
  
  wrongLevels = which(training_points_veg_date$ID%in%wronglevels_ID)
  m = ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  #ggsave(plot = m, filename = paste0("process/exploratory_images/time_Series/raw_data/plot_",Date_Int,"_harvested_",crossval_iter,"_", typeData,".png"))
  
  #write.csv(training_points_veg_date[-wrongLevels,], paste0("process/exploratory_images/time_Series/raw_data/data_",Date_Int,"_harvested_",crossval_iter,"_", typeData,".csv"))
  
  ## remove bad registers
  
  data_WithOutDuplicated=data_WithOutDuplicated[!data_WithOutDuplicated$ID %in% wronglevels_ID,]
  
  ###### reproductive
  training_points_veg = data_WithOutDuplicated[ data_WithOutDuplicated$Typ_Stg %in% c("reproductive"), ]
  
  training_points_veg_date = training_points_veg[ training_points_veg$Date %in% Date_Int, ]
  #ggplot(training_points_veg_date , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2016-01-10") & (training_points_veg_date$NDVI <0.4))
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = training_points_veg_date$ID[wrongLevels]
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2015-10-19") & (training_points_veg_date$NDVI >0.55))
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = c(wronglevels_ID,training_points_veg_date$ID[wrongLevels])
  wrongLevels = which(training_points_veg_date$ID%in%wronglevels_ID)
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2015-10-19") & (training_points_veg_date$LSWI >0.51))
  #ggplot(training_points_veg_date[!training_points_veg_date$ID%in%wronglevels_ID,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = c(wronglevels_ID, training_points_veg_date$ID[wrongLevels])
  
  m=ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  #ggsave(plot = m, filename = paste0("process/exploratory_images/time_Series/raw_data/plot_",Date_Int,"_reproductive_",crossval_iter,"_", typeData,".png"))
  #write.csv(training_points_veg_date[-wrongLevels,], paste0("process/exploratory_images/time_Series/raw_data/data_",Date_Int,"_reproductive_",crossval_iter,"_", typeData,".csv"))
  
  ## remove bad registers
  data_WithOutDuplicated=data_WithOutDuplicated[!data_WithOutDuplicated$ID %in% wronglevels_ID,]
  
  
  ###### ripening
  training_points_veg = data_WithOutDuplicated[ data_WithOutDuplicated$Typ_Stg %in% c("ripening"), ]
  
  training_points_veg_date = training_points_veg[ training_points_veg$Date %in% Date_Int, ]
  #ggplot(training_points_veg_date , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2015-10-22") & (training_points_veg_date$NDVI <0.5))
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = c(training_points_veg_date$ID[wrongLevels])
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2015-11-20") & (training_points_veg_date$NDVI <0.55))
  #ggplot(training_points_veg_date[!training_points_veg_date$ID%in%wronglevels_ID,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = c(wronglevels_ID, training_points_veg_date$ID[wrongLevels])
  

  
  
  
  wrongLevels = which(training_points_veg_date$ID%in%wronglevels_ID)
  m = ggplot(training_points_veg_date[!training_points_veg_date$ID%in%wronglevels_ID,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  #  ggsave(plot = m, filename = paste0("process/exploratory_images/time_Series/raw_data/plot_",Date_Int,"_ripening_",crossval_iter,"_", typeData,".png"))
  #write.csv(training_points_veg_date[!training_points_veg_date$ID%in%wronglevels_ID,], paste0("process/exploratory_images/time_Series/raw_data/data_",Date_Int,"_ripening_",crossval_iter,"_", typeData,".csv"))
  
  ## remove bad registers
  data_WithOutDuplicated=data_WithOutDuplicated[!data_WithOutDuplicated$ID %in% wronglevels_ID,]
  
  ###### soil
  training_points_veg = data_WithOutDuplicated[ data_WithOutDuplicated$Typ_Stg %in% c("soil"), ]
  training_points_veg_date = training_points_veg[ training_points_veg$Date %in% Date_Int, ]
  #ggplot(training_points_veg_date , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2015-12-21") & (training_points_veg_date$NDVI >0.4))
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = training_points_veg_date$ID[wrongLevels]
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2015-12-11") & (training_points_veg_date$NDVI >0.5))
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = c(wronglevels_ID,training_points_veg_date$ID[wrongLevels])
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2015-11-20") & (training_points_veg_date$NDVI >0.6))
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = c(wronglevels_ID,training_points_veg_date$ID[wrongLevels])
  
  
  wrongLevels = which(training_points_veg_date$ID%in%wronglevels_ID)
  m = ggplot(training_points_veg_date[!training_points_veg_date$ID%in%wronglevels_ID,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  #ggsave(plot = m, filename = paste0("process/exploratory_images/time_Series/raw_data/plot_",Date_Int,"_soil_",crossval_iter,"_", typeData,".png"))
  #write.csv(training_points_veg_date[-wrongLevels,], paste0("process/exploratory_images/time_Series/raw_data/data_",Date_Int,"_soil_",crossval_iter,"_", typeData,".csv"))
  
  ## remove bad registers
  data_WithOutDuplicated=data_WithOutDuplicated[!data_WithOutDuplicated$ID %in% wronglevels_ID,]
  
  
  ###### Vegetation
  training_points_veg = data_WithOutDuplicated[ data_WithOutDuplicated$Typ_Stg %in% c("other"), ]
  training_points_veg_date = training_points_veg[ training_points_veg$Date %in% Date_Int, ]
  #ggplot(training_points_veg_date , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  
  ###### Urban_Zones
  training_points_veg = data_WithOutDuplicated[ data_WithOutDuplicated$Typ_Stg %in% c("Urban_Zones"), ]
  training_points_veg_date = training_points_veg[ training_points_veg$Date %in% Date_Int, ]
  #ggplot(training_points_veg_date , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  
  
  
  ######## Date Change
  Date_Int = "20151127"
  
  ##### early_vegetative
  sdsad = data_WithOutDuplicated[data_WithOutDuplicated$Date %in% Date_Int,]
  unique(sdsad$Typ_Stg)
  training_points_veg = data_WithOutDuplicated[ data_WithOutDuplicated$Typ_Stg %in% c("early_vegetative"), ]
  training_points_veg_date = training_points_veg[ training_points_veg$Date %in% Date_Int, ]
  #ggplot(training_points_veg_date , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2015-12-11") & training_points_veg_date$NDVI <0.45)
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = training_points_veg_date$ID[wrongLevels]
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2016-01-10") & training_points_veg_date$NDVI <0.55)
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = c(wronglevels_ID,training_points_veg_date$ID[wrongLevels])
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2015-10-19") & training_points_veg_date$NDVI >0.50)
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = c(wronglevels_ID,training_points_veg_date$ID[wrongLevels])
  
  wrongLevels = which(training_points_veg_date$ID%in%wronglevels_ID)
  m = ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  #ggsave(plot = m, filename = paste0("process/exploratory_images/time_Series/raw_data/plot_",Date_Int,"early_vegetative",crossval_iter,"_", typeData,".png"))
  #write.csv(training_points_veg_date[-wrongLevels,], paste0("process/exploratory_images/time_Series/raw_data/data_",Date_Int,"early_vegetative",crossval_iter,"_", typeData,".csv"))
  
  ## remove bad registers
  data_WithOutDuplicated=data_WithOutDuplicated[!data_WithOutDuplicated$ID %in% wronglevels_ID,]
  
  ##### late_vegetative
  sdsad = data_WithOutDuplicated[data_WithOutDuplicated$Date %in% Date_Int,]
  unique(sdsad$Typ_Stg)
  training_points_veg = data_WithOutDuplicated[ data_WithOutDuplicated$Typ_Stg %in% c("late_vegetative"), ]
  training_points_veg_date = training_points_veg[ training_points_veg$Date %in% Date_Int, ]
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2015-11-20") & training_points_veg_date$NDVI <0.6)
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = training_points_veg_date$ID[wrongLevels]
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2015-10-22") & training_points_veg_date$NDVI >0.6)
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = c(wronglevels_ID,training_points_veg_date$ID[wrongLevels])
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2015-10-19") & training_points_veg_date$LSWI >0.51)
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = c(wronglevels_ID,training_points_veg_date$ID[wrongLevels])
  
  
  wrongLevels = which(training_points_veg_date$ID%in%wronglevels_ID)
  m = ggplot(training_points_veg_date[!training_points_veg_date$ID%in%wronglevels_ID,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  #ggsave(plot = m, filename = paste0("process/exploratory_images/time_Series/raw_data/plot_",Date_Int,"late_vegetative",crossval_iter,"_", typeData,".png"))
  #write.csv(training_points_veg_date[-wrongLevels,], paste0("process/exploratory_images/time_Series/raw_data/data_",Date_Int,"late_vegetative",crossval_iter,"_", typeData,".csv"))
  
  
  ## remove bad registers
  data_WithOutDuplicated=data_WithOutDuplicated[!data_WithOutDuplicated$ID %in% wronglevels_ID,]
  
  
  ##### ripening
  sdsad = data_WithOutDuplicated[data_WithOutDuplicated$Date %in% Date_Int,]
  unique(sdsad$Typ_Stg)
  training_points_veg = data_WithOutDuplicated[ data_WithOutDuplicated$Typ_Stg %in% c("ripening"), ]
  training_points_veg_date = training_points_veg[ training_points_veg$Date %in% Date_Int, ]
  m =ggplot(training_points_veg_date , aes(Image_date ,NDVI ,colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  
  #ggsave(plot = m, filename = paste0("process/exploratory_images/time_Series/raw_data/plot_",Date_Int,"ripening",crossval_iter,"_", typeData,".png"))
  #write.csv(training_points_veg_date[-wrongLevels,], paste0("process/exploratory_images/time_Series/raw_data/data_",Date_Int,"ripening",crossval_iter,"_", typeData,".csv"))
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2015-09-01") & training_points_veg_date$LSWI >0.55)
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = c(training_points_veg_date$ID[wrongLevels])
  
  ## remove bad registers
  data_WithOutDuplicated=data_WithOutDuplicated[!data_WithOutDuplicated$ID %in% wronglevels_ID,]
  
  ##### reproductive
  sdsad = data_WithOutDuplicated[data_WithOutDuplicated$Date %in% Date_Int,]
  unique(sdsad$Typ_Stg)
  training_points_veg = data_WithOutDuplicated[ data_WithOutDuplicated$Typ_Stg %in% c("reproductive"), ]
  training_points_veg_date = training_points_veg[ training_points_veg$Date %in% Date_Int, ]
  m =ggplot(training_points_veg_date , aes(Image_date ,NDVI ,colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2015-10-22") & training_points_veg_date$NDVI <0.55)
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = training_points_veg_date$ID[wrongLevels]
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2015-11-20") & training_points_veg_date$NDVI <0.55)
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = c(wronglevels_ID,training_points_veg_date$ID[wrongLevels])
  
  wrongLevels = which(training_points_veg_date$ID%in%wronglevels_ID)
  m =ggplot(training_points_veg_date[!training_points_veg_date$ID%in%wronglevels_ID,] , aes(Image_date ,NDVI ,colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  
  #ggsave(plot = m, filename = paste0("process/exploratory_images/time_Series/raw_data/plot_",Date_Int,"reproductive",crossval_iter,"_", typeData,".png"))
  #write.csv(training_points_veg_date[-wrongLevels,], paste0("process/exploratory_images/time_Series/raw_data/data_",Date_Int,"reproductive",crossval_iter,"_", typeData,".csv"))
  
  ## remove bad registers
  data_WithOutDuplicated=data_WithOutDuplicated[!data_WithOutDuplicated$ID %in% wronglevels_ID,]
  
  ##### soil
  sdsad = data_WithOutDuplicated[data_WithOutDuplicated$Date %in% Date_Int,]
  unique(sdsad$Typ_Stg)
  training_points_veg = data_WithOutDuplicated[ data_WithOutDuplicated$Typ_Stg %in% c("soil"), ]
  training_points_veg_date = training_points_veg[ training_points_veg$Date %in% Date_Int, ]
  m =ggplot(training_points_veg_date , aes(Image_date ,NDVI ,colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2015-10-22") & training_points_veg_date$NDVI <0.3)
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = training_points_veg_date$ID[wrongLevels]
  
  wrongLevels = which(training_points_veg_date$ID%in%wronglevels_ID)
  m =ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI ,colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  
  #ggsave(plot = m, filename = paste0("process/exploratory_images/time_Series/raw_data/plot_",Date_Int,"soil",crossval_iter,"_", typeData,".png"))
  #write.csv(training_points_veg_date[-wrongLevels,], paste0("process/exploratory_images/time_Series/raw_data/data_",Date_Int,"soil",crossval_iter,"_", typeData,".csv"))
  
  ## remove bad registers
  data_WithOutDuplicated=data_WithOutDuplicated[!data_WithOutDuplicated$ID %in% wronglevels_ID,]
  
  
  
  ######## Date Change
  Date_Int = "20160114"
  sdsad = data_WithOutDuplicated[data_WithOutDuplicated$Date %in% Date_Int,]
  unique(sdsad$Typ_Stg)
  
  ##### early_vegetative
  
  training_points_veg = data_WithOutDuplicated[ data_WithOutDuplicated$Typ_Stg %in% c("early_vegetative"), ]
  training_points_veg_date = training_points_veg[ training_points_veg$Date %in% Date_Int, ]
  #ggplot(training_points_veg_date , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2015-10-22") & (training_points_veg_date$NDVI >0.6 |
                                                                                      training_points_veg_date$NDVI <0.2))
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = training_points_veg_date$ID[wrongLevels]
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2016-01-10") & training_points_veg_date$NDVI >0.7)
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = c(wronglevels_ID,training_points_veg_date$ID[wrongLevels])
  
  wrongLevels = which(training_points_veg_date$ID%in%wronglevels_ID)
  m = ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  #ggsave(plot = m, filename = paste0("process/exploratory_images/time_Series/raw_data/plot_",Date_Int,"early_vegetative",crossval_iter,"_", typeData,".png"))
  #write.csv(training_points_veg_date[-wrongLevels,], paste0("process/exploratory_images/time_Series/raw_data/data_",Date_Int,"early_vegetative",crossval_iter,"_", typeData,".csv"))
  
  ## remove bad registers
  data_WithOutDuplicated=data_WithOutDuplicated[!data_WithOutDuplicated$ID %in% wronglevels_ID,]
  
  ##### ripening
  training_points_veg = data_WithOutDuplicated[ data_WithOutDuplicated$Typ_Stg %in% c("ripening"), ]
  training_points_veg_date = training_points_veg[ training_points_veg$Date %in% Date_Int, ]
  #ggplot(training_points_veg_date , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2015-10-22") & training_points_veg_date$NDVI >0.5)
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = training_points_veg_date$ID[wrongLevels]
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2016-01-10") & training_points_veg_date$NDVI <0.30)
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = c(wronglevels_ID,training_points_veg_date$ID[wrongLevels])
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2015-10-19") & training_points_veg_date$NDVI >0.5)
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = c(wronglevels_ID,training_points_veg_date$ID[wrongLevels])
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2015-10-19") & (training_points_veg_date$LSWI >0.51))
  #ggplot(training_points_veg_date[!training_points_veg_date$ID%in%wronglevels_ID,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = c(wronglevels_ID, training_points_veg_date$ID[wrongLevels])
  
  wrongLevels = which(training_points_veg_date$ID%in%wronglevels_ID)
  m = ggplot(training_points_veg_date[!training_points_veg_date$ID%in%wronglevels_ID,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  #ggsave(plot = m, filename = paste0("process/exploratory_images/time_Series/raw_data/plot_",Date_Int,"ripening",crossval_iter,"_", typeData,".png"))
  #write.csv(training_points_veg_date[-wrongLevels,], paste0("process/exploratory_images/time_Series/raw_data/data_",Date_Int,"ripening",crossval_iter,"_", typeData,".csv"))
  
  ## remove bad registers
  data_WithOutDuplicated=data_WithOutDuplicated[!data_WithOutDuplicated$ID %in% wronglevels_ID,]
  
  ##### soil
  training_points_veg = data_WithOutDuplicated[ data_WithOutDuplicated$Typ_Stg %in% c("soil"), ]
  training_points_veg_date = training_points_veg[ training_points_veg$Date %in% Date_Int, ]
  #ggplot(training_points_veg_date , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  
  
  wrongLevels = (training_points_veg_date$Image_date==as.Date("2016-01-10") & training_points_veg_date$NDVI>0.30)
  #ggplot(training_points_veg_date[!wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID =training_points_veg_date$ID[wrongLevels]
  
  wrongLevels = (training_points_veg_date$ID%in%wronglevels_ID)
  m = ggplot(training_points_veg_date[!wrongLevels,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  #ggsave(plot = m, filename = paste0("process/exploratory_images/time_Series/raw_data/plot_",Date_Int,"soil",crossval_iter,"_", typeData,".png"))
  #write.csv(training_points_veg_date[!wrongLevels,], paste0("process/exploratory_images/time_Series/raw_data/data_",Date_Int,"soil",crossval_iter,"_", typeData,".csv"))
  
  ## remove bad registers
  data_WithOutDuplicated=data_WithOutDuplicated[!data_WithOutDuplicated$ID %in% wronglevels_ID,]
  
  ##### harvested
  training_points_veg = data_WithOutDuplicated[ data_WithOutDuplicated$Typ_Stg %in% c("harvested"), ]
  training_points_veg_date = training_points_veg[ training_points_veg$Date %in% Date_Int, ]
  #ggplot(training_points_veg_date , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  
  
  wrongLevels = (training_points_veg_date$Image_date==as.Date("2015-10-22") & training_points_veg_date$NDVI<0.50)
  #ggplot(training_points_veg_date[!wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID =training_points_veg_date$ID[wrongLevels]
  
  wrongLevels = which(training_points_veg_date$Image_date==as.Date("2015-11-20") & training_points_veg_date$NDVI <0.7)
  #ggplot(training_points_veg_date[-wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = c(wronglevels_ID,training_points_veg_date$ID[wrongLevels])
  
  wrongLevels = (training_points_veg_date$ID%in%wronglevels_ID)
  m = ggplot(training_points_veg_date[!wrongLevels,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  #ggsave(plot = m, filename = paste0("process/exploratory_images/time_Series/raw_data/plot_",Date_Int,"harvested",crossval_iter,"_", typeData,".png"))
  #write.csv(training_points_veg_date[!wrongLevels,], paste0("process/exploratory_images/time_Series/raw_data/data_",Date_Int,"harvested",crossval_iter,"_", typeData,".csv"))
  
  ## remove bad registers
  data_WithOutDuplicated=data_WithOutDuplicated[!data_WithOutDuplicated$ID %in% wronglevels_ID,]
  
  
  ######## Date Change
  Date_Int = "20151211"
  sdsad = data_WithOutDuplicated[data_WithOutDuplicated$Date %in% Date_Int,]
  unique(sdsad$Typ_Stg)
  
  ##### early_vegetative
  training_points_veg = data_WithOutDuplicated[ data_WithOutDuplicated$Typ_Stg %in% c("early_vegetative"), ]
  training_points_veg_date = training_points_veg[ training_points_veg$Date %in% Date_Int, ]
  #ggplot(training_points_veg_date , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  
  
  wrongLevels = (training_points_veg_date$Image_date==as.Date( "2015-12-11") & training_points_veg_date$NDVI <0.5)
  #ggplot(training_points_veg_date[!wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = training_points_veg_date$ID[wrongLevels]
  
  wrongLevels = (training_points_veg_date$ID%in%wronglevels_ID)
  m = ggplot(training_points_veg_date[!wrongLevels,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  #ggsave(plot = m, filename = paste0("process/exploratory_images/time_Series/raw_data/plot_",Date_Int,"early_vegetative",crossval_iter,"_", typeData,".png"))
  #write.csv(training_points_veg_date[!wrongLevels,], paste0("process/exploratory_images/time_Series/raw_data/data_",Date_Int,"early_vegetative",crossval_iter,"_", typeData,".csv"))
  
  ## remove bad registers
  data_WithOutDuplicated=data_WithOutDuplicated[!data_WithOutDuplicated$ID %in% wronglevels_ID,]
  
  ##### ripening
  training_points_veg = data_WithOutDuplicated[ data_WithOutDuplicated$Typ_Stg %in% c("ripening"), ]
  training_points_veg_date = training_points_veg[ training_points_veg$Date %in% Date_Int, ]
  #ggplot(training_points_veg_date , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  
  
  wrongLevels = (training_points_veg_date$Image_date==as.Date( "2015-12-11") & training_points_veg_date$NDVI <0.35)
  #ggplot(training_points_veg_date[!wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = training_points_veg_date$ID[wrongLevels]
  
  wrongLevels = (training_points_veg_date$Image_date==as.Date( "2015-11-20") & training_points_veg_date$NDVI <0.7)
  #ggplot(training_points_veg_date[!wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = c(wronglevels_ID,training_points_veg_date$ID[wrongLevels])
 
  wrongLevels = (training_points_veg_date$Image_date==as.Date( "2015-09-01") & training_points_veg_date$LSWI >0.5)
  #ggplot(training_points_veg_date[!wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = c(wronglevels_ID,training_points_veg_date$ID[wrongLevels])
  
  
  wrongLevels = (training_points_veg_date$ID%in%wronglevels_ID)
  m = ggplot(training_points_veg_date[!wrongLevels,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  #ggsave(plot = m, filename = paste0("process/exploratory_images/time_Series/raw_data/plot_",Date_Int,"ripening",crossval_iter,"_", typeData,".png"))
  #write.csv(training_points_veg_date[!wrongLevels,], paste0("process/exploratory_images/time_Series/raw_data/data_",Date_Int,"ripening",crossval_iter,"_", typeData,".csv"))
  
  ## remove bad registers
  data_WithOutDuplicated=data_WithOutDuplicated[!data_WithOutDuplicated$ID %in% wronglevels_ID,]
  
  ##### reproductive
  training_points_veg = data_WithOutDuplicated[ data_WithOutDuplicated$Typ_Stg %in% c("reproductive"), ]
  training_points_veg_date = training_points_veg[ training_points_veg$Date %in% Date_Int, ]
  ggplot(training_points_veg_date , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  
  
  wrongLevels = (training_points_veg_date$Image_date==as.Date( "2015-11-20") & training_points_veg_date$NDVI <0.65)
  #ggplot(training_points_veg_date[!wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = training_points_veg_date$ID[wrongLevels]
  
  wrongLevels = (training_points_veg_date$Image_date==as.Date( "2015-10-22") & training_points_veg_date$NDVI <0.35)
  #ggplot(training_points_veg_date[!wrongLevels,] , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)
  wronglevels_ID = c(wronglevels_ID,training_points_veg_date$ID[wrongLevels])
  
  wrongLevels = (training_points_veg_date$ID%in%wronglevels_ID)
  m = ggplot(training_points_veg_date[!wrongLevels,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+geom_point()+geom_line(alpha = 0.8)
  
  ## remove bad registers
  data_WithOutDuplicated=data_WithOutDuplicated[!data_WithOutDuplicated$ID %in% wronglevels_ID,]
  
  ###############End cleaning process
  data_WithOutDuplicated = data_WithOutDuplicated[,-ncol(data_WithOutDuplicated)]
  
  return(data_WithOutDuplicated)
}

names_optical = names(vi_images)
SentinelimageReference = vi_images[[1]][[1]]
SentinelimageReference[] = NA

#training_points =readOGR(paste0("data/",typeData,"_phenIdenti/v5",typeData,"_iterationCV_", crossval_iter ,".shp"))
training_points = SpatialPointsDataFrame(training_points_[,1:2],data = training_points_,
                                         proj4string = crs(SentinelimageReference))
training_points = Extract_raster_data(raster_images = vi_images , spatial_points = training_points , 
                                      spectral_bands = names(vi_images[[1]]))
names(training_points) = sapply(names_optical, function(name_imag) IdentiyImageDate(name_imag))


rm(vi_images)

### Join all points
i = 1
training_points = do.call(rbind , lapply(1:length(training_points) , function(i){
  dataper_=training_points[[i]]
  dataper_$Image_date = names(training_points)[i]
  dataper_$plygn_n = factor(as.character(dataper_$plygn_n))
  return(dataper_)
}))


head(training_points)
data_WithOutDuplicated = data.frame(training_points %>%
                                      group_by(x,  y , data_type,iteration, plygn_n,  Date , Typ_Stg , Image_date) %>% filter(row_number( y) == 1))



data_WithOutDuplicated$Image_date = as.Date(data_WithOutDuplicated$Image_date , format = "%Y%m%d")


training_points_veg = data_WithOutDuplicated[ data_WithOutDuplicated$Typ_Stg %in%  c("late_vegetative"  ,"early_vegetative","harvested","reproductive","ripening"), ]

training_points_veg$ID = with(training_points_veg,paste(x , y , plygn_n,  Date , Typ_Stg ,data_type , "i", iteration, sep = "_"))

m=ggplot(training_points_veg[training_points_veg$iteration %in% "1", ], aes(Image_date ,NDVI , colour=plygn_n, group = ID))+
  geom_point()+geom_line(alpha = 0.8)+ facet_grid(rows = vars(Typ_Stg), cols = vars(Date))

ggsave(plot = m , filename = paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/temp/testuncleanNDVI.jpg"),
       width = 34, height = 25, units = "cm")
#ggplot(training_points_veg_date , aes(Image_date ,NDVI , colour=Typ_Stg, group = ID))+geom_point()+geom_line(alpha = 0.8)


data_WithOutDuplicated = clean_ts_function(data_WithOutDuplicated, crossval_iter)
training_points_veg = data_WithOutDuplicated[ data_WithOutDuplicated$Typ_Stg %in%  c("late_vegetative"  ,"early_vegetative","harvested","reproductive","ripening"), ]

training_points_veg$ID = with(training_points_veg,paste(x , y , plygn_n,  Date , Typ_Stg ,data_type , "i", iteration, sep = "_"))
# 
# data_date = training_points_veg[training_points_veg$Date %in% "20151127" &
#                                   training_points_veg$Image_date<as.Date("2015-11-27"),]
data_date = do.call(rbind,lapply(c("20151127","20151228","20160114"), function(date_campaign){
  training_points_veg[training_points_veg$Date %in% date_campaign & 
                        training_points_veg$Image_date<CheckDateFormat(date_campaign),]
}))
lines_plot = data.frame(group_by(data_date, Typ_Stg,Date,Image_date) %>%summarise(ndvi = median(NDVI, na.rm = T)))
data_date$temp = paste0(data_date$Image_date, data_date$Typ_Stg)
m = ggplot(data_date, aes(as.factor(Image_date), NDVI ,group = temp)) + 
  geom_boxplot(aes(fill = Typ_Stg), alpha = .60)+lims(y = c(0,1))+ theme_bw() + labs( fill = "Growth Stage", x = "Satellite Date")+ 
  stat_summary(fun.y=median, geom="line", aes(colour = Typ_Stg,group = Typ_Stg), size = 1.4)+ 
  stat_summary(fun.y=median, geom="point", aes(colour = Typ_Stg,group = Typ_Stg), size = 1.8)+
  facet_grid(Date ~.)


m+ guides(colour=FALSE)+scale_colour_manual(values = c ('harvested'='firebrick3', "ripening"="goldenrod", 
                                                        "early_vegetative"="cyan4", "reproductive" = "darkgreen","late_vegetative" = "chartreuse3"
                                                        ))+
  scale_fill_manual(values = c ('harvested'='firebrick3', "ripening"="goldenrod", 
                                    "early_vegetative"="cyan4", "reproductive" = "darkgreen","late_vegetative" = "chartreuse3"
  ))


###### lines 
CheckDateFormat("20151127")
data_date = do.call(rbind,lapply(c("20151127","20151228","20160114"), function(date_campaign){
  training_points_veg[training_points_veg$Date %in% date_campaign & 
                        training_points_veg$Image_date<CheckDateFormat(date_campaign),]
}))


lines_plot = data.frame(group_by(data_date, Typ_Stg,Date,Image_date) %>%summarise(ndvi = median(NDVI, na.rm = T)))

lines_plot$temp = paste0(lines_plot$Image_date, lines_plot$Typ_Stg)


m = ggplot(lines_plot, aes(as.factor(Image_date), ndvi ,group = Typ_Stg)) + 
  geom_boxplot( alpha = .60)+
  geom_line(aes(color = Typ_Stg), alpha = .9)+
  geom_point(aes(color = Typ_Stg), alpha = .9)+
  lims(y = c(0,1))+ theme_bw() + labs( color = "Etapa de Crecimiento",
                                       y = "promedio de NDVI",
                                       x = "Fecha de captura de las imágenes satelitales")+ 
  facet_grid(Date ~.)

m+scale_colour_manual(values = c ('harvested'='firebrick3', "ripening"="goldenrod", 
                                                        "early_vegetative"="cyan4", "reproductive" = "darkgreen","late_vegetative" = "chartreuse3"
))+
  scale_fill_manual(values = c ('harvested'='firebrick3', "ripening"="goldenrod", 
                                "early_vegetative"="cyan4", "reproductive" = "darkgreen","late_vegetative" = "chartreuse3"
  ))



m=ggplot(training_points_veg[training_points_veg$iteration%in%1,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+
  geom_point()+geom_line(alpha = 0.8)+ facet_grid(rows = vars(Typ_Stg), cols = vars(Date))



ggsave(plot = m , filename = paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/temp/testcleanNDVI.jpg"),
       width = 34, height = 25, units = "cm")


head(data_WithOutDuplicated)
############################### end cleaning dataset
#------------------------------------------------------------------------

#-- Start smooth process

######### ---> Kernel Smooth process

## Date of Interest
ListIndex = c("NDVI")
## Set date period
datestoEvaluate = c("20151127","20151211","20151228", "20160114")
DateInt="20151127"
dataSmoothed = lapply(datestoEvaluate, function(DateInt){
  as.Date(DateInt,format=" %Y%m%d")-96
  rangeDate=c(as.Date(DateInt,format="%Y%m%d")-96,(as.Date(DateInt,format="%Y%m%d")+0))
  nDaysPoints=length(rangeDate[1]:rangeDate[2])/16
  
  VI = "NDVI"
  data_reshape = lapply(ListIndex , function(VI){
    ## Select data for the interested date
    
    data_WithOutDuplicated_date = data_WithOutDuplicated[data_WithOutDuplicated$Date %in% DateInt,]
    
    ## filter by interest date
    data_WithOutDuplicated_date=
      data_WithOutDuplicated_date[data_WithOutDuplicated_date$Image_date %in% (as.Date(DateInt,format=" %Y%m%d")-120):as.Date(DateInt,format=" %Y%m%d"),]
    
    if(dim(data_WithOutDuplicated_date)[1]>0){
      data_WithOutDuplicated_date$Image_date = paste0(VI,"_", as.character(data_WithOutDuplicated_date$Image_date, format = "%Y%m%d"))
      ## reshape dataset
      data_WithOutDuplicated_date$ID =
        with(data_WithOutDuplicated_date,paste0(as.character(x), as.character(y),"_p_" ,as.character(plygn_n),"_s_",as.character(type),"_type_",as.character(Typ_Stg),"_d" ,as.character(Date), "_i_", as.character(iteration), "_t_",data_type))
      
      pos_names = sapply(c("x","y", "plygn_n", "type","Typ_Stg","Date","iteration","data_type"), function(x)
        which(names(data_WithOutDuplicated_date)%in%x))
      data_WithOutDuplicated_date[,-pos_names]
      data_reshape = reshape2::dcast( data_WithOutDuplicated_date[,c(VI,"Image_date","ID")] , ID ~ Image_date, value.var = VI ,fun.aggregate=mean)
      ## assign ID code
      return(data_reshape)
      
    }else{
      return(NULL)
    }
  })
  
  if(!is.null(data_reshape[[1]])){
    ###### join
    
    data_reshape = plyr::join_all(data_reshape, by = "ID")
    
    ### remove those pixels with high percentage of NA in the time series
    
    pixelsPercentageNA=data.frame(Pixel=data_reshape$ID,
                                  Percentage=apply(data_reshape[,-which(names(data_reshape)%in% c("ID"))],
                                                   1,function(x){(sum(is.na(x))/length(x)*100)}))
    
    pixelsHighPercentage=which(pixelsPercentageNA$Percentage>50)
    if(length(pixelsHighPercentage)>0){
      imagesToProcess=data_reshape[-pixelsHighPercentage,]
      
    }else{
      imagesToProcess=data_reshape
    }
    
    
    ###############
    
    #### Smooth process
    
    rm(pixelsPercentageNA)
    row.names(imagesToProcess) = imagesToProcess$ID
    
    imagesToProcess[150:200,]
    
    KernelSmootValues=SmoothAllPixels(DateInt = DateInt,ImagesInfo = imagesToProcess,
                                      Veg_Indexes = ListIndex,
                                      rangeDate = rangeDate,
                                      nDaysPoints = nDaysPoints,parallelProcess =T,ncores = 5)
    
    #-- End smooth process
    ########
    #----------------------------------
    ########
    #-- Start export step
    
    ####Organize the information in a table
    
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
    
    ####### Assign class label to data frame
    
    num_pos = regexpr(pattern = "_type_",row.names(TableValuesToClassify))
    num_pos_end = regexpr(pattern = "_d2",row.names(TableValuesToClassify))
    TableValuesToClassify$Class =str_sub(row.names(TableValuesToClassify),num_pos+6, num_pos_end-1)
    
    ### export training data
    ## remove NA
    TableValuesToClassify = TableValuesToClassify[!apply(TableValuesToClassify , 1 , function(column)T%in%is.na(column)),]
    

    return(TableValuesToClassify)
  }else{
    return(NULL)
  }
  
})

TableValuesToClassify = do.call(rbind,dataSmoothed)


graphic_smoothed = function(test, vi = "NDVI", iter = "1", final = FALSE){
  test$ID = row.names(test)
  test = (reshape2::melt(test))
  num_pos = regexpr(pattern = "_p_",test$ID)
  num_pos_end = regexpr(pattern = "_s_",test$ID)
  test$plygn_n = str_sub(test$ID,num_pos+3, num_pos_end-1)
  
  num_pos = regexpr(pattern = "_type_",test$ID)
  num_pos_end = regexpr(pattern = "_d2",test$ID)
  test$Typ_Stg= str_sub(test$ID,num_pos+6, num_pos_end-1)
  
  num_pos = regexpr(pattern = "_d2",test$ID)
  num_pos_end = regexpr(pattern = "_i_",test$ID)
  test$Date= str_sub(test$ID,num_pos+2, num_pos_end-1)
  
  num_pos = regexpr(pattern = "_i_",test$ID)
  num_pos_end = regexpr(pattern = "_t_",test$ID)
  test$iteration= str_sub(test$ID,num_pos+3, num_pos_end-1)
  
  num_pos = regexpr(pattern = "_i_",test$ID)
  num_pos_end = regexpr(pattern = "_t_",test$ID)
  test$training= str_sub(test$ID,num_pos+3, num_pos_end-1)
  if (final){
    test = test[test$Typ_Stg%in%unique(test$Typ_Stg)[c(1:5, 7)],]
  }else{
    test = test[test$Typ_Stg%in%unique(test$Typ_Stg)[c(1:4, 6)],]
  }
  
  m=ggplot(test[test$iteration %in%iter, ] , aes(variable ,value , colour=plygn_n, group = ID))+
    geom_point()+geom_line(alpha = 0.8)+ facet_grid(rows = vars(Typ_Stg), cols = vars(Date))
  m
}


test = TableValuesToClassify[,grepl(names(TableValuesToClassify), pattern = "NDVI_Da")]
m = graphic_smoothed (test, vi = "NDVI", iter = "3")

ggsave(plot = m , filename = paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/temp/llData_v8_1_NDVI_da.png"),
       width = 34, height = 25, units = "cm")



## second Clean Proccess.
cleandata = TableValuesToClassify[TableValuesToClassify$Class %in%  "late_vegetative" & (!TableValuesToClassify$NDVI_Date_7<0.45),]
TableValuesToClassify = TableValuesToClassify[!TableValuesToClassify$Class %in%  "late_vegetative",]
TableValuesToClassify = rbind(cleandata,TableValuesToClassify)

cleandata = TableValuesToClassify[TableValuesToClassify$Class %in%  "late_vegetative" & (!TableValuesToClassify$NDVI_Date_2>0.5),]
TableValuesToClassify = TableValuesToClassify[!TableValuesToClassify$Class %in%  "late_vegetative",]
TableValuesToClassify = rbind(cleandata,TableValuesToClassify)


cleandata = TableValuesToClassify[TableValuesToClassify$Class %in%  "late_vegetative" & (!TableValuesToClassify$NDVI_Date_5>0.6),]
TableValuesToClassify = TableValuesToClassify[!TableValuesToClassify$Class %in%  "late_vegetative",]
TableValuesToClassify = rbind(cleandata,TableValuesToClassify)


cleandata = TableValuesToClassify[TableValuesToClassify$Class %in%  "early_vegetative" & (!TableValuesToClassify$NDVI_Date_7<0.35),]
TableValuesToClassify = TableValuesToClassify[!TableValuesToClassify$Class %in%  "early_vegetative",]
TableValuesToClassify = rbind(cleandata,TableValuesToClassify)

cleandata = TableValuesToClassify[TableValuesToClassify$Class %in%  "early_vegetative" & (!TableValuesToClassify$NDVI_Date_4>0.5),]
TableValuesToClassify = TableValuesToClassify[!TableValuesToClassify$Class %in%  "early_vegetative",]
TableValuesToClassify = rbind(cleandata,TableValuesToClassify)



cleandata = TableValuesToClassify[TableValuesToClassify$Class %in%  "ripening" & (!TableValuesToClassify$NDVI_Date_1>0.6),]
TableValuesToClassify = TableValuesToClassify[!TableValuesToClassify$Class %in%  "ripening",]
TableValuesToClassify = rbind(cleandata,TableValuesToClassify)


cleandata = TableValuesToClassify[TableValuesToClassify$Class %in%  "late_vegetative" & (!TableValuesToClassify$NDVI_Date_3>0.55),]
TableValuesToClassify = TableValuesToClassify[!TableValuesToClassify$Class %in%  "late_vegetative",]
TableValuesToClassify = rbind(cleandata,TableValuesToClassify)

cleandata = TableValuesToClassify[TableValuesToClassify$Class %in%  "reproductive" & (!TableValuesToClassify$NDVI_Date_5<0.30),]
TableValuesToClassify = TableValuesToClassify[!TableValuesToClassify$Class %in%  "reproductive",]
TableValuesToClassify = rbind(cleandata,TableValuesToClassify)

cleandata = TableValuesToClassify[TableValuesToClassify$Class %in%  "harvested" & (!TableValuesToClassify$NDVI_Date_7>0.55),]
TableValuesToClassify = TableValuesToClassify[!TableValuesToClassify$Class %in%  "harvested",]
TableValuesToClassify = rbind(cleandata,TableValuesToClassify)

test = TableValuesToClassify[,grepl(names(TableValuesToClassify), pattern = "NDVI_Da")]


graphic_s = graphic_smoothed (test, vi =  "NDVI", iter = "2", final = TRUE)


head(graphic_s$data)

m = ggplot(graphic_s$data, aes(as.factor(variable), value )) + 
  geom_boxplot(aes(fill = Typ_Stg), alpha = 0.4)+
  lims(y = c(-0.035,1))+ theme_bw() + 
  labs( fill = "Growth Stage", y = "NDVI", x = "")+ 
  stat_summary(fun.y=median, geom="line", aes(colour = Typ_Stg,group = Typ_Stg), size = 1.4)+ 
  stat_summary(fun.y=median, geom="point", aes(colour = Typ_Stg,group = Typ_Stg), size = 1.8)      



m+ guides(colour=FALSE)+scale_colour_manual(values = c ('harvested'='firebrick3', "ripening"="goldenrod", "other" = "gray", 
                                                        "early_vegetative"="cyan4", "reproductive" = "darkgreen","late_vegetative" = "chartreuse3"
))+
  scale_fill_manual(values =  c ('harvested'='firebrick3', "ripening"="goldenrod","other" = "gray", 
                                 "early_vegetative"="cyan4", "reproductive" = "darkgreen","late_vegetative" = "chartreuse3"
  ))+ stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                      geom="crossbar", width=0.5)+
   stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="pointrange", color="red")

test = TableValuesToClassify[,grepl(names(TableValuesToClassify), pattern = "NDVI_Da")]
m = graphic_smoothed (test, vi = "NDVI", iter = "3")


ggsave(plot = m , filename = paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/temp/allData_v8_smoothed_NDVI_date.png"),
       width = 34, height = 25, units = "cm")


##### dot plot

my_mean=aggregate(graphic_s$data$value , by=list(graphic_s$data$Typ_Stg,graphic_s$data$variable) , mean)
colnames(my_mean)=c("type","vari" , "mean")
my_sd=aggregate(graphic_s$data$value , by=list(graphic_s$data$Typ_Stg,graphic_s$data$variable) , sd) 
colnames(my_sd)=c("type","vari" , "sd")

my_info=merge(my_mean , my_sd , by.x=c(1,2) , by.y=c(1,2))

ggplot(graphic_s$data) + 
  geom_point(aes(x = variable, y = value, colour = Typ_Stg) , size=0.1) + 
  geom_point(data = my_info, aes(x=vari , y = mean, colour = type) , size = 2) +
  geom_errorbar(data = my_info, aes(x = vari, y = sd, ymin = mean - sd, ymax = mean + sd,  colour = type) , width = 0.7 , size=1)


# m = ggplot(graph_ , aes(variable , value , colour = Class,group = ID))+geom_point()+geom_line(alpha = 0.8)
# ggsave(plot = m , filename = paste0("process/exploratory_images/time_Series/phen_identification_",typeData,"_iterationCV_", crossval_iter ,"_class_todas_predicted_V2.png"),width = 30, height = 14, units = "cm")

print(dim(TableValuesToClassify))
head(TableValuesToClassify)
table(TableValuesToClassify$Class)
write.csv(TableValuesToClassify ,paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/model_inputs/phen_identification/optical_data_conf8_ndvi_derstepndvi.csv"))

#############################################3
crossval_iter = 2



typeData = "Validation"
training_points_radar = SpatialPointsDataFrame(training_points_[,1:2],data = training_points_,
                                               proj4string = crs(SentinelimageReference))

########## ----- > Create Parametrics for radar images 
training_points = training_points_radar

datestoEvaluate = c("20151127","20151211","20151228", "20160114")
DateInt="20151127"
radar_Bands=c("Sigma0_VV_db", 'Sigma0_VV_GLCMMean', 'Sigma0_VV_GLCMVariance')


#### metrics per month

library(doParallel)
library(foreach)
ncores = 2
cl =  makeCluster(ncores)
registerDoParallel(cl)
Sys.time()->start
radar_path = "D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/satellite_imagery/radar_data/"
calculate_rastermetrics = function(DateStart, DateEnd,radar_bands,tile,
         fpath = "D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/satellite_imagery/radar_data/",
         metrics = c("mean", "min", "max")){
  ## set capturing period
  Date_interval = DateStart : DateEnd
  # filter filenames
  fnames = list.files(fpath, pattern = paste0(tile , ".RData$"))
  
  ## Filter files by dates
  
  fnames = fnames[as.Date(IdentiyImageDate(fnames), format = "%Y%m%d")%in% Date_interval]
  fnames = fnames[order(as.Date(IdentiyImageDate(fnames), format = "%Y%m%d"))]
  
  ## read imagery
  radar_images = lapply(1:length(fnames), function(imag_index){
    file_name = fnames[imag_index]
    
    load(file = paste0(  fpath,file_name))
    
    bands_toSelect = names(sat_imag)[names(sat_imag)%in%radar_Bands]
    sat_imag = sat_imag[[bands_toSelect]]
    cat(str_sub(file_name ,1, -7) ," loaded\n")
    return(sat_imag)
  })
  
  ## calculate metrics
  metrics_features = lapply(radar_bands, function(band){
    band_info = stack(lapply(radar_images , function(x) x[[band]]))
    
    band_ifo = metrics_(band_info,metrics)
  })
  
  stack(metrics_features)
  
}


data_radar = foreach(DateInt = datestoEvaluate, 
                     .export = c("calculate_rastermetrics","radar_Bands", "training_points")) %dopar% {
                       
  library(raster)
  library(stringr)

  ## 30 days metrics
  DateEnd=as.Date(DateInt, format = "%Y%m%d")
  DateStart= DateEnd - 30

  metrics_30days = calculate_rastermetrics ( DateEnd - 30,DateEnd ,radar_bands , tile)

  ## 60  days metrics
  DateStart= (DateEnd ) - ( 60)

  metrics_60days = calculate_rastermetrics (DateStart, ( DateEnd - 31),radar_Bands, tile)
  
  ## 90  days metrics
  DateStart= (DateEnd ) - ( 90)

  metrics_90days = calculate_rastermetrics (DateStart, ( DateEnd - 61),radar_Bands, tile)
  
  ## stack metrics
  metrics= c("mean","min","max")
  radar_data = stack(metrics_30days, metrics_60days, metrics_90days)
  names(radar_data) = do.call(c,lapply(c("30days","60days","90days"),function(y)
    paste0(do.call(c,lapply(radar_Bands,function(x) 
      (paste0(x, "_", metrics)))),"_",y)))

  radar_data_extracted = Extract_raster_data(raster_images = radar_data, spatial_points = training_points [training_points$Date %in% DateInt,] , 
                                             spectral_bands = names(radar_data))
  
  # graph_values_limits = reshape2::melt(radar_data_extracted[[1]])
  # 
  # ggplot(graph_values_limits[grepl(graph_values_limits$variable, pattern = "Sigma0_VV_db_min"),] , aes(variable, value, fill = Typ_Stg)) + geom_boxplot()
  

  radar_data_extracted[[1]]
  
  
}
print(Sys.time()-start)
stopCluster(cl)

#### metrics quantile


Sys.time()->start

data_radar = foreach(DateInt = datestoEvaluate, 
                     .export = c("calculate_raster_quantiles","radar_Bands", "training_points")) %do% {
                       
    library(raster)
    library(stringr)
                       
                       
  DateEnd=as.Date(DateInt, format = "%Y%m%d")
  DateStart= DateEnd - 60
  Date_interval = DateStart : DateEnd
  # filter filenames
  fnames = list.files(radar_path, pattern = paste0(tile , ".RData$"))
  
  ## Filter files by dates
  
  fnames = fnames[as.Date(IdentiyImageDate(fnames), format = "%Y%m%d")%in% Date_interval]
  fnames = fnames[order(as.Date(IdentiyImageDate(fnames), format = "%Y%m%d"))]
  
  ## read imagery
  radar_images = lapply(1:length(fnames), function(imag_index){
    file_name = fnames[imag_index]
    
    load(file = paste0(  radar_path,file_name))
    
    bands_toSelect = names(sat_imag)[names(sat_imag)%in%radar_Bands]
    sat_imag = sat_imag[[bands_toSelect]]
    cat(str_sub(file_name ,1, -7) ," loaded\n")
    return(sat_imag)
  })
  
  ncores = 4
  
  quantile_imags = calculate_raster_quantiles(radar_images,radar_Bands,SentinelimageReference, division = 20, ncores = 4)
  
  quantile_imags = stack(quantile_imags)
  
  names(quantile_imags) = do.call(c,lapply(radar_Bands,function(x) 
      (paste0(x, c("_p05","_p25","_p50","_p75","_p95")))))
  
  
  
  radar_data_extracted = Extract_raster_data(raster_images = quantile_imags, 
                                             spatial_points = training_points [training_points$Date %in% DateInt,] , 
                                             spectral_bands = names(quantile_imags))[[1]]
  
  radar_data_extracted
  #radar_data_extracted$diff = radar_data_extracted$Sigma0_VV_db_p95 - radar_data_extracted$Sigma0_VV_db_p05
  
}

graph_values_limits = reshape2::melt(radar_data_extracted[[1]])
# 
ggplot(graph_values_limits[grepl(graph_values_limits$variable, pattern = "saldana_rasterref"),] , aes(variable, value, fill = Typ_Stg)) + geom_boxplot()

data_ras = data_radar[[1]]
data_ras = do.call(rbind,lapply(data_radar, function(data_ras){
  data_ras$Class = data_ras$Typ_Stg
  
  data_ras$ID = 
    with(data_ras,paste0(as.character(x), as.character(y),"_p_" ,as.character(plygn_n),"_s_",as.character(type),"_type_",as.character(Typ_Stg),"_d" ,as.character(Date), "_i_", as.character(iteration), "_t_",data_type))
  data_ras = data_ras[,-which(names(data_ras)%in% c("x", "y", "plygn_n", "type","Typ_Stg", "Date" , "iteration","data_type"))]
  duplicated_info = which(duplicated(data_ras$ID ))
  if(length(duplicated_info)>0)
    data_ras = data_ras[-duplicated_info,]
  row.names(data_ras) = data_ras$ID
  data_ras = data_ras[,-ncol(data_ras)]
  print(dim(data_ras))
  return(data_ras)
}))

head(data_ras) 
dataras_back = data_ras
data_ras = dataras_back
for(band in radar_Bands){
  data_ras= data_ras[,-ncol(data_ras)]
  data_ras[, paste0(band, "_diffmin_30_60days")] = data_ras[, paste0(band, "_min_30days")] - data_ras[, paste0(band, "_min_60days")]
  data_ras[, paste0(band, "_diffmax_30_60days")] =data_ras[, paste0(band, "_max_30days")] - data_ras[, paste0(band, "_max_60days")]
  
  data_ras[, paste0(band, "_diffmin_60_90days")] =data_ras[, paste0(band, "_min_60days")] - data_ras[, paste0(band, "_min_90days")]
  data_ras[, paste0(band, "_diffmax_60_90days")] =data_ras[, paste0(band, "_max_60days")] - data_ras[, paste0(band, "_max_90days")]
  
  data_ras[, paste0(band, "_diffmin_30_90days")] =data_ras[, paste0(band, "_min_30days")] - data_ras[, paste0(band, "_min_90days")]
  data_ras[, paste0(band, "_diffmax_30_90days")] =data_ras[, paste0(band, "_max_30days")] - data_ras[, paste0(band, "_max_90days")]
  data_ras$Class = dataras_back$Class
}

band = radar_Bands[1]
lapply(radar_Bands ,function(band){
  graph_ = data_ras[,c(grep(names(data_ras) , pattern = paste0(band,"_m")), ncol(data_ras))]
  
  graph_ = reshape2::melt(graph_ )
  t_print = unique(graph_$Class)[1]
  m=ggplot(graph_, aes(variable , value , fill = Class))+geom_boxplot()+
    theme(text=element_text(size=10),
          axis.title.y=element_text(size = rel(1.2),colour = "#999999"),
          axis.title.x=element_text(size = rel(1.2),colour = "#888888"),
          axis.text.x  = element_text(angle=60, hjust=1),
          panel.background=element_rect(fill="white",colour = "black"),
          panel.grid.major = element_line(colour = "gray"))
  ggsave(plot = m , filename = paste0("process/exploratory_images/radar_images/phen_identification_",typeData,"_iterationCV_", crossval_iter ,"_",band,".png"),width = 30, height = 14, units = "cm")
  
} )

str_length(row.names(TableValuesToClassify))
write.csv(data_ras ,paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/model_inputs/phen_identification/RADAR_conf8_radar.csv"))

table(str_sub(row.names(data_ras),str_length(row.names(data_ras))-2,str_length(row.names(data_ras))))
windows()



# 
# DateStart= (DateEnd ) - ( 90)
# Date_interval = DateStart :( DateEnd - 61)
# 
# list_radar_names = list.files(paste0(Main_Folder , "process/rdata/radar_info/") , pattern = paste0(Locality , ".RData$"))
# 
# ## Filter files by dates
# 
# list_radar_names = list_radar_names[as.Date(IdentiyImageDate(list_radar_names), format = "%Y%m%d")%in% Date_interval]
# list_radar_names = list_radar_names[order(as.Date(IdentiyImageDate(list_radar_names), format = "%Y%m%d"))]
# 
# 
# radar_images = read_radar_Rdata (list_radar_names, radar_Bands,
#                                  SentinelimageReference, paste0(Main_Folder , "process/rdata/radar_info/"))
# 
# 
# max_60 = max(stack(radar_images))
# min_60 = min(stack(radar_images))
# mean_60 = mean(stack(radar_images))
# plot(max_30 - max_60)
# radar_data_extracted = Extract_raster_data(raster_images = stack(abs(min_30- max_60),abs(min_60- max_90)), spatial_points = training_points , 
#                                            spectral_bands = names( stack((min_30- min_60),(min_60- min_90))))
# 
# graph_values_limits = reshape2::melt(radar_data_extracted[[1]])
# 
# ggplot(graph_values_limits[grepl(graph_values_limits$variable, pattern = "layer"),] , aes(variable, value, fill = Typ_Stg)) + geom_boxplot()
# 
# 
# dates_images = as.Date(IdentiyImageDate(list_radar_names), format = "%Y%m%d")
# 
# wavelengths =  sapply(dates_images , function(x) x - dates_images[1])
# names(radar_images) = IdentiyImageDate(list_radar_names)
# 
# Sys.time()->start
# data_to_integrated = do.call(cbind,lapply(radar_data_extracted, function(x) x[,"Sigma0_VV_db"] ))
# data_integral = apply(data_to_integrated, 1, function(valraster)calc_integral(wavelengths,valraster))
# 
# head()[,1:6]
# datatoplot = data.frame(value = data_integral, typecov = radar_data_extracted[[1]]$Typ_Stg)
# ggplot(datatoplot, aes(typecov, value))+geom_boxplot()
# 
# dataModel_inegral=lapply (radar_data_extracted,function(image_info){
#   image_info = image_info[,names(central_wavelengths)]/2^12
# })
# print(Sys.time()-start)
# image_info = optical_data_extracted[[1]]
# 
# dataModel_inegral=lapply (optical_data_extracted,function(image_info){
#   image_info = image_info[,names(central_wavelengths)]/2^12
#   data_integral = apply(image_info, 1, function(valraster)calc_integral(wavelengths,valraster))
# })
# 
# rm(quantile_imags)
# 
# rm(radar_images)
# quantile_imags = raster::stack(quantile_imags)
# radar_names_database = do.call(c,lapply(radar_Bands, function(band) paste0(band,"_pa_",1:length(c(0.05,.25, .50, .75,0.95)))))
# names(quantile_imags) = radar_names_database
# 
# radar_min =  list(min_30, min_60,min_90)
# meanValues = mean(stack(radar_min))
# SumbandDiff=list()
# for (j in 1:length(radar_min)){
#   SumbandDiff[[j]]=(radar_min[[j]]-meanValues)*(radar_min[[j]]-meanValues)
# }
# SumbandDiff=sum(stack(SumbandDiff),na.rm=T)
# 
# stdValues=((SumbandDiff*(1/(length(radar_min)-1)))^(1/2))
# 
# 
# rm(raster_temp)
# 
# dataModel_Param_radar = cbind(dataModel_Param_radar, dataModel_Param_radar_conf2)

