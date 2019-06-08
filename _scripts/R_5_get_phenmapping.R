################ 
### 
### This is a code, that was designed to Map Rice Zones for an spicific region
### To do this part, you had to select four Sentinel 2 images
###                   
###                   
###   Author:  Andr?s Aguilar
###     CIAT - DAPA - AEPS
#########


rm(list = ls())
##### Set Locality Name



get_featuresmin_maxvals = function(dataf){
  ## number of features
  num_features = ncol(dataf)
  
  #calculate minimun and maximun values
  minmax_list =lapply(1:num_features, function(x){
    
    return(c(min(dataf[,x]),max(dataf[,x])))
  })
  
  ### rename list
  names(minmax_list) = names(dataf)
  return(minmax_list)
}
# 
feature_scaling =  function(feature, limits = c(0,1), type = "min_max"){
  
  ## calculate min max feature scaling
  (feature  - limits[1])/(limits[2]-limits[1])
  
}



######### Load functions



source(paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/_scripts/R_Main_Functions_GrowCropIndentification.R"))

##### Load Libraries
setwd("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/")
libs=c("snowfall","caret","nnet","SDMTools","stringr","raster","ff", "rgeos","rgdal","dplyr" , "tidyverse")
PackageReading(libs)



### set vegetation indexes names

indexes = c("NDVI")
tile = "north_tolima"
tiles = c("col_t2","col_t1")
#tiles = c("north_tolima")
######## Select Images for the study: there are two critriums 1) dates 2) quality and 3) visual 

DateInt="2016-09-06"
summaryts = T

start_date = "2016-06-10"
end_date = "2016-11-11"
dates_percent = 10
dates_percent = 37


### interpolate multiple dates 
lapply(tiles, function(tile){
  PercMaxBadPixels=50
  if(tile == "north_tolima"){
    PercMaxBadPixels=38.7
    dates_percent = 38

    
  }
  inventory = get_inventory(tile)
  inventory$BadPixels[is.na( inventory$BadPixels)] = 100
  dates_topredict = inventory$Date[CheckDateFormat(inventory$Date)>as.Date(start_date) &
                                     inventory$BadPixels < dates_percent &
                                     is.na(inventory$Delete)  &
                                     CheckDateFormat(inventory$Date)<as.Date(end_date)]
  

  print(PercMaxBadPixels)
  # 
  dates_topredict = dates_topredict[!is.na(dates_topredict)]
  dates_topredict = c(dates_topredict[1],dates_topredict[2:length(dates_topredict)][diff(as.Date(dates_topredict))>4])
  print(dates_topredict)
  DateInt = as.character(dates_topredict)[4]
  lapply(as.character(dates_topredict), function(DateInt){
    cat(DateInt,' \n')
    ##### get inventory
    inventory = get_inventory( tile)
    
    
    
    ## Period
    startdate=CheckDateFormat(DateInt)-130
    enddate=CheckDateFormat(DateInt)
    

    #
    vi_images = get_vi_layers(inventory = inventory , veg_indexes = indexes, startdate, enddate,
                              PercMaxBadPixels = PercMaxBadPixels)
    
    ## the smooth process is done with at least 5 images
    if(length(vi_images)>4){
      ### get interpolation time series per pixel
      if(tile == "north_tolima"){
        sat_ref_path =  'D:/OneDrive - Universidad Nacional de Colombia/MScPhil/disease_identification/spatial_data/field_ref_north_tolima.tif'
        TableValuesToClassify = calculate_interpolation_series (images_data = vi_images, date_int = enddate, veg_indexes = indexes,summaryts = T,
                                                                image_reference_path =sat_ref_path,num_process = 6)
      }
      if(tile == "col_t1"){
        sat_ref_path = "D:/OneDrive - Universidad Nacional de Colombia/MScPhil/disease_identification/spatial_data/field_disesase_ref_col_t1_fields.tif"
        #sat_ref_path = "D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/model_ouputs/rice_layers/Filled_col_t1xgbTree_20180306_20180912_optical-radar_15 x 15_100.tif"
        rasteref = raster(sat_ref_path)
        #raster::xyFromCell(rasteref,3264996)
        raster::cellFromXY(rasteref,c(695694.76,1158995.33))
        TableValuesToClassify = calculate_interpolation_series (images_data = vi_images, 
                                                                date_int = enddate, 
                                                                veg_indexes = indexes,
                                                                image_reference_path = sat_ref_path,
                                                                num_process = 6)
      }
      if(tile == "col_t3"){
        #sat_ref_path = "D:/OneDrive - Universidad Nacional de Colombia/MScPhil/disease_identification/spatial_data/field_disesase_ref_col_t3.tif"
        sat_ref_path = "D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/model_ouputs/rice_layers/Filled_col_t3xgbTree_conf2_20151016_20160420_optical_db_15 x 15_100.tif"
        # rasteref = raster(sat_ref_path)
        # raster::xyFromCell(rasteref,3264996)
        # raster::cellFromXY(rasteref,c(500457.6503,427736.7741))
        TableValuesToClassify = calculate_interpolation_series (images_data = vi_images, date_int = enddate, veg_indexes = indexes,
                                                                image_reference_path = sat_ref_path,num_process = 6)
      }
      if(tile == "col_t7"){
        sat_ref_path = "D:/OneDrive - Universidad Nacional de Colombia/MScPhil/disease_identification/spatial_data/field_disesase_ref_col_t7.tif"
        TableValuesToClassify = calculate_interpolation_series (images_data = vi_images, date_int = enddate,
                                                                veg_indexes = indexes,
                                                                image_reference_path = sat_ref_path,num_process = 6)
      }
      if(tile == "col_t2"){
        sat_ref_path = "D:/OneDrive - Universidad Nacional de Colombia/MScPhil/disease_identification/spatial_data/field_disesase_ref_col_t2.tif"
        
        rasteref = raster(sat_ref_path)
        raster::xyFromCell(rasteref,3264996)
        raster::cellFromXY(rasteref,c(622017.736,1089370.003))
        
        TableValuesToClassify = calculate_interpolation_series (images_data = vi_images, date_int = enddate, veg_indexes = indexes,
                                                                image_reference_path = sat_ref_path,num_process = 6)
      }
      
      if(nrow(TableValuesToClassify) !=  length(which(is.na(TableValuesToClassify$NDVI_Date_1)))){
        #dir.create(paste0(setPathStyle(mainfolder),"temp"))
        ## create backup file
        
        file_name = paste0("temp/timeseriessmoothed_",tile,CheckDateFormat(DateInt),"_filled.csv")
        write.csv(TableValuesToClassify, file_name)
        
        TableValuesToClassify = read.csv(file_name,row.names = 1)
        
        ## create backup file
        #TableValuesToClassify = read.csv( file_name,row.names = 1)
        
        
        ##### remove rows with na
        
        LevelsWNA_SecondFilter=unique(unlist(sapply(1:ncol(TableValuesToClassify),function(columnVal){which(is.na(TableValuesToClassify[,columnVal]))})))
        LevelsWNA_SecondFilter=LevelsWNA_SecondFilter[order(LevelsWNA_SecondFilter)]
        
        if(length(LevelsWNA_SecondFilter)!=0){
          TableValuesToClassify_WithoutNA=TableValuesToClassify[-c(LevelsWNA_SecondFilter),]
        }else{
          TableValuesToClassify_WithoutNA=TableValuesToClassify
        }
        testImage = raster(sat_ref_path)
        # 
        rm(vi_images)
        # # predict
        load( file = "D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/classification_models/xgboost_phen_identification_veg_newfeatures.RData")
        library(xgboost)
        model_ML=model
        namesVar=model_ML$feature_names
        data_toclassxgboost = xgboost::xgb.DMatrix(as.matrix(TableValuesToClassify_WithoutNA[,namesVar]),  missing = NA)
        
        
        PheStagePre=predict(model_ML,data_toclassxgboost)
        PheStagePre <- matrix(PheStagePre, ncol = 6, byrow = TRUE)
        PheStagePre <- apply(PheStagePre, 1, which.max)
        modelMethod="xgboost"
        testImage[]=NA
        testImage[as.numeric(row.names(TableValuesToClassify_WithoutNA))]=PheStagePre
        plot(testImage)
        dir.create(paste0('D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/model_ouputs/growth_stages/',tile,'/xgboost_veg/'),showWarnings = F)
        
        writeRaster(testImage, filename=paste0('D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/model_ouputs/growth_stages/',tile,'/xgboost_veg/',tile,"_",modelMethod,"_ndvits_summary_",as.character(CheckDateFormat(DateInt),format("%Y%m%d")),".tif"), format="GTiff",overwrite=T)
        
        ## random Forest
        load( file = "D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/classification_models/rf_phen_identification_veg_newfeatures.RData")
        library(randomForest)
        namesVar = row.names(model$importance)
        data_toclassrf = TableValuesToClassify_WithoutNA[,namesVar]
        
        
        model_ML=model
        
        
        PheStagePre=as.numeric(as.character(predict(model_ML,data_toclassrf)))
        testImage[]=NA
        
        testImage[as.numeric(row.names(TableValuesToClassify_WithoutNA))]=(PheStagePre+1)
        plot(testImage)
        modelMethod="rf"
        dir.create(paste0('D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/model_ouputs/growth_stages/',tile,'/rf_veg/'),showWarnings = F)
        writeRaster(testImage, filename=paste0('D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/model_ouputs/growth_stages/',tile,'/',modelMethod,'_veg/',tile,"_",modelMethod,"_ndvits_summary_",as.character(CheckDateFormat(DateInt),format("%Y%m%d")),".tif"), format="GTiff",overwrite=T)
        
        ## svm poly

        library(e1071)
        load( file = "D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/classification_models/svm_polynomial_phen_identification_veg_newfeatures.RData")

        
        numfeatures = length(namesVar)
        dates_pos = grepl("Date", names(model$min_maxvalues))
        

        #
        # ###
        #
        min_dates = min(unlist(model$min_maxvalues[dates_pos]))
        max_dates = max(unlist(model$min_maxvalues[dates_pos]))
        # 
        # ## calculate min max derivatives
        # 
        deriva_pos = grepl("deriva",  names(model$min_maxvalues))
        
        min_der = min(unlist(model$min_maxvalue[deriva_pos]))
        max_der = max(unlist(model$min_maxvalue[deriva_pos]))

        
        dates_scaled = do.call(cbind,lapply(which(dates_pos),function(x){
          feature_scaling(TableValuesToClassify_WithoutNA[,namesVar][,x], c(min_dates, max_dates))
        }))
        
        derivatives_scaled = do.call(cbind,lapply(which(deriva_pos),function(x){
          feature_scaling(TableValuesToClassify_WithoutNA[,namesVar][,x], c(min_der, max_der))
        }))
        
        x_variables_scaled = data.frame(cbind(dates_scaled,
                                              derivatives_scaled))
        
        if(summaryts){
          
          summarize_scaled = do.call(cbind,lapply((1:numfeatures)[-c(which(dates_pos),which(deriva_pos))],function(x){
            feature_scaling(TableValuesToClassify_WithoutNA[,namesVar][,x], model$min_maxvalue[[x]])
          }))
          
          #
          # ##
          x_variables_scaled = data.frame(cbind(dates_scaled,derivatives_scaled,summarize_scaled
                                                ))
          
        }
        names(x_variables_scaled) = namesVar
        
        
        predictedValues = as.numeric(unname(predict(model$model,x_variables_scaled)))
        
        
        testImage[]=NA
        
        testImage[as.numeric(row.names(TableValuesToClassify_WithoutNA))]=predictedValues
        plot(testImage)
        modelMethod="svm_pol"
        dir.create(paste0('D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/model_ouputs/growth_stages/',tile,'/',modelMethod,'_veg/'),showWarnings = F)
        
        writeRaster(testImage, filename=paste0('D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/model_ouputs/growth_stages/',tile,'/',modelMethod,'_veg/',tile,"_",modelMethod,"_ndvits_summary_",as.character(CheckDateFormat(DateInt),format("%Y%m%d")),".tif"), format="GTiff",overwrite=T)
        
        ## svm radial
        load( file = "D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/classification_models/svm_radial_phen_identification_veg_newfeatures.RData")
        predictedValues = as.numeric(unname(predict(model$model,x_variables_scaled)))
        
        
        testImage[]=NA
        
        testImage[as.numeric(row.names(TableValuesToClassify_WithoutNA))]=predictedValues
        plot(testImage)
        modelMethod="svm_radial"
        dir.create(paste0('D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/model_ouputs/growth_stages/',tile,'/',modelMethod,'_veg/'),showWarnings = F)
        
        writeRaster(testImage, filename=paste0('D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/model_ouputs/growth_stages/',tile,'/',modelMethod,'_veg/',tile,"_",modelMethod,"_ndvits_summary_",as.character(CheckDateFormat(DateInt),format("%Y%m%d")),".tif"), format="GTiff",overwrite=T)
        
        
      }
      
         #print(paste0('D:/temp/',tile,"_",modelMethod,"_conf8_",as.character(CheckDateFormat(DateInt),format("%Y%m%d")),"ndvi_test.tif"))
    }
  })
    
    
  
})


#----------------------------------
######## 

paste0('D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/model_ouputs/growth_stages/',tile,'/',tile,"filled_",modelMethod,"_conf8_",as.character(CheckDateFormat(DateInt),format("%Y%m%d")),"ndvideri_field.tif")

data_rast = raster(paste0('D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/model_ouputs/growth_stages/',tile,'/',tile,"filled_",modelMethod,"_conf8_",as.character(CheckDateFormat(DateInt),format("%Y%m%d")),"ndvi_test.tif"))

data_count = table(data_rast[])

(data_count*100)/10000

row.names(model$importance)[order(model$importance)]
