################ 
### 
###   Proccess Satellite data 
### 
###                   
###                   
###   Author:  Andres Aguilar
###     CIAT - DAPA - AEPS - UoM
#########


rm(list=ls())

##### Set Locality Name

tile = "col_t3"

######### Load functions


source(paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/_scripts/R_Main_Functions_GrowCropIndentification.R"))


source(paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/_scripts/R_Main_Functions_GrowCropIndentification.R"))
source("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/_scripts/metrics_functions.R")

############ Load Libraries





libs=c("snowfall","caret","nnet","SDMTools","stringr","raster","ff", "rgeos","rgdal","dplyr" , "tidyverse")
PackageReading(libs)


############## MAIN CODE

DateEnd="2016-04-20"
DateStart=CheckDateFormat(DateEnd) - 190

inventory = get_inventory(tile)

optical_bands = c("blue","green","red","vrededg1","vrededg2","vrededg3","nir","narrownir","swir1","swir2")

bad_pixels_limit = 15
optical_images = get_optical_imagery(inventory , init_date = DateStart, 
                                     end_date = DateEnd, 
                                     bands = optical_bands,bad_pixels_limit = bad_pixels_limit, optical_mission = c("L2A","S2"))


optical_images = optical_images[order(as.Date(IdentiyImageDate(names(optical_images)), format = "%Y%m%d"))]


SentinelimageReference = optical_images[[1]][[1]]

SentinelimageReference[] = NA

#### unsupervised classification

data_ndvi = do.call(cbind,lapply(1:length(optical_images), function(i){
  ndvi = (optical_images[[i]]$nir-optical_images[[i]]$red) / (optical_images[[i]]$nir + optical_images[[i]]$red)
  ndvi[]
}))

head(data_ndvi)

napixels = unique(do.call(c,apply(data_ndvi, 2,function(x)which(is.na(x)))))
if(length(napixels)>0){
  data_ndvi = data_ndvi[-napixels,]
}
datakmeans = kmeans(data_ndvi ,10)
SentinelimageReference[!(1:ncell(SentinelimageReference))%in%napixels] = datakmeans$cluster
plot(SentinelimageReference)


writeRaster(SentinelimageReference, filename=paste0("D:/temp/",tile,"_clusterndvi2.tif"), format="GTiff",overwrite=T)



optical_images
dataImages_optical =optical_metrics(tile, DateStart, DateEnd)



################################
################# READ RADAR DATA
#################

## set radar band names

# radar_Bands=c("Sigma0_VV_db","Sigma0_VV_ASM","Sigma0_VV_Dissimilarity",
#            "Sigma0_VV_Energy","Sigma0_VV_Contrast",
#                    "Sigma0_VV_Entropy","Sigma0_VV_GLCMCorrelation","Sigma0_VV_GLCMVariance",
#                    "Sigma0_VV_Homogeneity","Sigma0_VV_MAX","Sigma0_VV_GLCMMean")

radar_Bands=c("Sigma0_VV_db")
inventory = get_inventory(tile)

radar_images = get_rada_images(inventory , init_date = DateStart, end_date = DateEnd, 
                               bands = "Sigma0_VV_db")

plot(radar_images$S1_VV_20151017_col_t3_10m)
radar_images = radar_images[order(as.Date(IdentiyImageDate(names(radar_images)), format = "%Y%m%d"))]

#### unsupervised classification

data_radar = do.call(cbind,lapply(1:length(radar_images), function(i){
  radar_da = radar_images[[i]]
  radar_da[]
}))

head(data_radar)

napixels = unique(do.call(c,apply(data_radar, 2,function(x)which(is.na(x)))))
if(length(napixels)>0){
  data_radar = data_radar[-napixels,]
}
datakmeans = kmeans(data_radar ,10)
SentinelimageReference[!(1:ncell(SentinelimageReference))%in%napixels] = datakmeans$cluster
plot(SentinelimageReference)


writeRaster(SentinelimageReference, filename=paste0("D:/temp/",tile,"_clusterrada.tif"), format="GTiff",overwrite=T)



##################################### 

# use a raster image reference 

rep = 1
setwd(Main_Folder)

set_types = c("training","validation")

### read data training
mergedata = do.call(rbind,lapply(set_types, function(set_type){
  do.call(rbind,lapply(1:10 , function(rep){
    data_model = read.csv(paste0("data/",set_type,"_riceidentification/" ,set_type,"_iterationCV_",rep,"_addpoints.csv")) 
    
    data_model$iteration = rep
    
    data_model$feed = set_type
    return(data_model)}))
}))

print(unique(mergedata$feed))
sp_data= SpatialPointsDataFrame(mergedata[,1:2],data = mergedata,
                                proj4string = crs(SentinelimageReference))


### Extract information
optical_data_extracted = Extract_raster_data(optical_images , sp_data , optical_bands)
names(optical_data_extracted) = names(optical_images)

radar_data_extracted = Extract_raster_data(radar_images , sp_data , radar_Bands)
names(radar_data_extracted) = names(radar_images)


########## ----- > Create Parametrics for optical images 

band = standarBand_names [1]

dataModel_Param_optical=do.call(cbind , lapply (optical_bands,function(band){
  
  band_info = do.call(cbind,lapply(optical_data_extracted , function(x) x[,band]))
  
  band_info = data.frame(do.call(rbind,lapply(1:nrow(band_info) ,  function(j){
    
    meanBand = stdBand = minBand = maxBand = NA
    row_i = band_info[j,]
    count_na = sum(is.na(row_i))
    if(!count_na > 1){
      
      meanBand=mean(row_i,na.rm = T)
      stdBand=sd(row_i,na.rm = T)
      minBand=min(row_i,na.rm = T)
      maxBand=max(row_i,na.rm = T)
      
    }
    
    return(c(stdBand,meanBand,minBand,maxBand ))
  })))
  names(band_info) = paste0(band , c("_SD","_MEAN","_MIN","_MAX"))
  return(band_info)
}))

dataModel_Param_optical= cbind(optical_data_extracted[[1]][,1:7],dataModel_Param_optical)

rm(optical_images)

############# parameters 2 integral


###


library(pracma)
rm(radar_images)


central_wavelengths = list('blue'=496.6,
                           'green'=560.0,
                           'red'=664.5,
                           'vrededg1'=703.9,
                           'vrededg2'=740.2,
                           'vrededg3'=782.5,
                           'nir'=835.1,
                           'narrownir'=864.8,
                           'swir1'=1613.7,
                           'swir2'=2202.4)

calc_integral = function(wavelengths,valraster){
  
  if(T%in%is.na(valraster)){
    return(NA)
  }else{
    return(trapz(wavelengths,valraster))
  }
  
}

wavelengths = do.call(c, central_wavelengths)
image_info = optical_data_extracted[[1]]
Sys.time()->start
dataModel_inegral=lapply (optical_data_extracted,function(image_info){
  image_info = image_info[,names(central_wavelengths)]/2^12
  data_integral = apply(image_info, 1, function(valraster)calc_integral(wavelengths,valraster))
})

data_integral = data.frame(do.call(rbind,lapply(1:nrow(do.call(cbind,dataModel_inegral)), function(i){
  x = do.call(cbind,dataModel_inegral)[i,]
  partial_data =c(mean(x, na.rm = TRUE ),
                  sd(x, na.rm = TRUE ),
                  max(x, na.rm = TRUE ),
                  min(x, na.rm = TRUE ))
  
  
  return(partial_data)
  
})))
names(data_integral) = c("ref_int_mean", "ref_int_sd", "ref_int_max", "ref_int_min")

print(Sys.time()-start)

dataModel_Param_optical = cbind(dataModel_Param_optical,data_integral)

########## ----- > Create Parametrics for radar images 

## number of splits
names(radar_data_extracted) = IdentiyImageDate(names(radar_data_extracted))

dates_radar = (names(radar_data_extracted))

Groups_radar = CalculateNumDivisions(dates_radar, numDiv = 6)

datapixel=data.frame(dates=dates_radar,
                     date=as.Date(dates_radar,format="%Y%m%d"))

x
(lapply(c(60, 120), function(x)as.Date(datapixel$date[length(datapixel$date)]-x)))

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
## calculate parameters per periods 1 
## summarisse by mean sd max and min for data gather each 2 months


dataModel_Param_radar = do.call(cbind, lapply(unique(Groups_radar$groups), function(group_ref){
  radar_toCalc = radar_data_extracted[Groups_radar[Groups_radar$groups%in%group_ref,1]]
  
  
  dataModel_Param_radar=do.call(cbind , lapply (radar_Bands,function(band){
    
    band_info = do.call(cbind,lapply(radar_toCalc , function(x) x[,band]))
    
    band_info = data.frame(do.call(rbind,lapply(1:nrow(band_info) ,  function(j){
      
      meanBand = stdBand = minBand = maxBand = NA
      row_i = band_info[j,]
      count_na = sum(is.na(row_i))
      if(!count_na > 1){
        
        meanBand=mean(row_i,na.rm = T)
        stdBand=sd(row_i,na.rm = T)
        minBand=min(row_i,na.rm = T)
        maxBand=max(row_i,na.rm = T)
        
      }
      
      return(c(stdBand,meanBand,minBand,maxBand ))
    })))
    names(band_info) = paste0(band ,"_", group_ref , c("_SD","_MEAN","_MIN","_MAX"))
    return(band_info)
  }))
  
  return(dataModel_Param_radar)
  
}))

dataModel_Param_radar= cbind(radar_data_extracted[[1]][,1:7],dataModel_Param_radar)

dataModel_all= cbind(dataModel_Param_optical,dataModel_Param_radar)

## calculate parameters per periods 2
## order that from less to highest after calculate summatory under the curve

sigma_value= radar_Bands[1]
dataModel_Param_radar_conf2 = do.call(cbind, lapply(radar_Bands, function(sigma_value){
  
  extract_band_info = sapply(radar_data_extracted, 
         function(x) x[[sigma_value]])

  
  dataOrdered  = lapply(1:nrow(extract_band_info), function(j){
    pixel_val = extract_band_info[j,]
    if(sum(is.na(pixel_val))==length(pixel_val)){
      sorted = unname(pixel_val)
      data_toreturn = t(data.frame(rep(NA,5), row.names = paste0(sigma_value,"_pa_",1:5)))
    }else{

      limits = quantile(pixel_val, c(0.05,.25, .50, .75,0.95), na.rm = T)
      data_toreturn = t(data.frame(limits, row.names = paste0(sigma_value,"_pa_",1:length(limits))))
      #extremes = quantile(sorted, c(.05, .95) )
      #data_toreturn = cbind(data_toreturn , data.frame(dife = abs(extremes[1] - extremes[2])))
    }
    return(data_toreturn)
  })
  
  test_ = do.call(rbind,dataOrdered)
  return(test_)
}))
library(doParallel)
quantile_imags = calculate_raster_quantiles(radar_images,radar_Bands,SentinelimageReference, division = 20, ncores = 5)

names(quantile_imags) = radar_Bands
lapply(radar_Bands, function(radar_band){
  
  data_extracted = Extract_raster_data(quantile_imags[[radar_band]] , sp_data , names(quantile_imags[[radar_band]]))

  })


dataModel_Param_radar_conf2 = cbind(radar_data_extracted[[1]][,1:7],dataModel_Param_radar_conf2)
graph_values_limits = reshape2::melt(cbind(radar_data_extracted[[1]][,3:5],dataModel_Param_radar_conf2[,-(1:7)]))
dataModel_Param_radar=dataModel_Param_radar_conf2
ggplot(graph_values_limits[grepl(graph_values_limits$variable, pattern = "Sigma0_VV_Dissimilarity"),] , aes(variable, value, fill = type)) + geom_boxplot()
#dataModel_Param_radar = cbind(dataModel_Param_radar, dataModel_Param_radar_conf2)
### export data
#export_tables = function(df, columns_to_delete, )
tye_confopt = "optical_metrics_4im"
tye_confrad = "radar_conf2_quantile_db_GLCM"
iter = 1
feed = "training"
lapply(c(1:10), function(iter){
  lapply(set_types, function(feed){
    df_export = dataModel_Param_optical[(dataModel_Param_optical$iteration %in% iter) & (as.character(dataModel_Param_optical$feed) %in% feed),]

    df_export = df_export[,-which(names(df_export) %in% c("feed", "iteration"))]
    file_name = paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/model_inputs/rice_identification/",tye_confopt,"_",feed,"_iterationCV_",iter,"_dataV2.csv")
    write.csv(df_export , file_name)

    df_export = dataModel_Param_radar[(dataModel_Param_radar$iteration %in% iter) & (as.character(dataModel_Param_radar$feed) %in% feed),]

    df_export = df_export[,-which(names(df_export) %in% c("feed", "iteration"))]
    file_name = paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/model_inputs/rice_identification/",tye_confrad,"_",feed,"_iterationCV_",iter,".csv")
    write.csv(df_export , file_name)
  })
})

  

#### graphics
optical_data = dataModel_Param_optical 
bands_order = optical_bands_order
graphic_data = function(optical_data,bands_order, metric="_MEAN"){
  data_red = optical_data[,grepl(names(optical_data), pattern = metric)]
  data_red = data_red[,sapply(names(data_red),function(name_col)T%in%grepl(paste0(bands_order,metric), pattern = name_col))]
  
  data_red$output_model = optical_data$type2
  data_red$output_model_2 = optical_data$type
  data_red$output_model = as.character(data_red$output_model)
  bool_pos = do.call(c, lapply(c("Ripening", "Vegetative", "Reproductive","Harvested"),function(patt){
    which(grepl(unique(data_red$output_model), pattern = patt))}))
  data_red$output_model[data_red$output_model%in%unique(data_red$output_model)[unique(bool_pos)]] = "Rice" 
  data_reshape = reshape2::melt( data_red )
  
  data_reshape$variable = factor(data_reshape$variable , levels = paste0(bands_order, metric))
  
  return(data_reshape)
  
}

optical_bands_order = c("blue","green","red","vrededg1","vrededg2","vrededg3",
                        "nir","narrownir","swir1","swir2" )
cover_types = c("RipeningHarvestedHarvestedSoil",
                "Vegetative_EndingReproductiveRipening",
                "Forest","Mountain","Bare_Soil",
                "Vegetation","Vegetative_EndingRipening","SoilVegetative_Starting" ,"Water",'Roads',"Urban_Zones")

graphdata_ = graphic_data(dataModel_Param_optical[dataModel_Param_optical$type2 %in% cover_types,], optical_bands_order, "_MAX")

graphdata_$aux = paste0(graphdata_$output_model, graphdata_$output_model_2)

ggplot(graphdata_ , aes(variable, value/(2^12)))+geom_boxplot(aes(fill = output_model)) + labs(x = "Spectral Bands Average", y = "Reflectance")



radar_bands_order = c("Sigma0_VV_GLCMCorrelation_1","Sigma0_VV_GLCMCorrelation_2","Sigma0_VV_GLCMCorrelation_3")
graphdata_radar = graphic_data(dataModel_Param_radar[dataModel_Param_radar$type2 %in% cover_types,], radar_bands_order, "_SD")

ggplot(graphdata_radar , aes(variable, value))+geom_boxplot(aes(fill = output_model_2)) + labs(x = "Spectral Bands Average", y = "Reflectance")



head(dataModel_Param_radar)

dataModel_Param_radar$dif_1_MIN = with(dataModel_Param_radar,(Sigma0_VV_db_2_MIN - Sigma0_VV_db_1_MIN))
dataModel_Param_radar$dif_2_MIN = with(dataModel_Param_radar,(Sigma0_VV_db_3_MIN - Sigma0_VV_db_2_MIN))

dataModel_Param_radar$dif_1_MAX = with(dataModel_Param_radar,(Sigma0_VV_db_2_MAX - Sigma0_VV_db_1_MAX))
dataModel_Param_radar$dif_2_MAX = with(dataModel_Param_radar,(Sigma0_VV_db_3_MAX - Sigma0_VV_db_2_MAX))

dataModel_Param_radar$dif_1_MEAN = with(dataModel_Param_radar,(Sigma0_VV_db_2_MEAN - Sigma0_VV_db_1_MEAN))
dataModel_Param_radar$dif_2_MEAN = with(dataModel_Param_radar,(Sigma0_VV_db_3_MEAN - Sigma0_VV_db_2_MEAN))

dataModel_Param_radar$dif_1_SD = with(dataModel_Param_radar,(Sigma0_VV_db_2_SD - Sigma0_VV_db_1_SD))
dataModel_Param_radar$dif_2_SD = with(dataModel_Param_radar,(Sigma0_VV_db_3_SD - Sigma0_VV_db_2_SD))

radar_bands_order_2 = c("dif_1","dif_2")

graphdata_radar = graphic_data(dataModel_Param_radar[dataModel_Param_radar$type2 %in% cover_types,], radar_bands_order_2
                               , "_MEAN")

ggplot(graphdata_radar , aes(variable, value))+geom_boxplot(aes(fill = output_model_2)) + labs(x = "Spectral Bands Average", y = "Reflectance")



###


library(pracma)
rm(radar_images)

central_wavelengths = list('blue'=496.6,
                           'green'=560.0,
                           'red'=664.5,
                           'vrededg1'=703.9,
                           'vrededg2'=740.2,
                           'vrededg3'=782.5,
                           'nir'=835.1,
                           'narrownir'=864.8,
                           'swir1'=1613.7,
                           'swir2'=2202.4)

calc_integral = function(wavelengths,valraster){
  
  if(T%in%is.na(valraster)){
    return(NA)
  }else{
    return(trapz(wavelengths,valraster))
  }
  
}
wavelengths = do.call(c, central_wavelengths)
image_info = optical_data_extracted[[1]]
Sys.time()->start
dataModel_inegral=lapply (optical_data_extracted,function(image_info){
  image_info = image_info[,names(central_wavelengths)]/2^12
  data_integral = apply(image_info, 1, function(valraster)calc_integral(wavelengths,valraster))
})

data_integral = data.frame(do.call(rbind,lapply(1:nrow(do.call(cbind,dataModel_inegral)), function(i){
  x = do.call(cbind,dataModel_inegral)[i,]
  partial_data =c(mean(x, na.rm = TRUE ),
                sd(x, na.rm = TRUE ),
                max(x, na.rm = TRUE ),
                min(x, na.rm = TRUE ))
  
  
  return(partial_data)
  
})))
names(data_integral) = c("ref_int_mean", "ref_int_sd", "ref_int_max", "ref_int_min")

print(Sys.time()-start)


data_joined = cbind(optical_data_extracted[[1]][1:7], data_integral)


ggplot(data_joined, aes(type,ref_int_min))+geom_boxplot()

length(data_integral)

image_info$integral =  data_integral

cover_types = c("RipeningHarvestedHarvestedSoil",
                "Vegetative_EndingReproductiveRipening",
                "Bare_Soil",
                "Vegetation","Vegetative_EndingRipening","SoilVegetative_Starting" )


data_joined_ = do.call(rbind,lapply(1:length(dataModel_inegral), function(index_list){
  data_joined = cbind(optical_data_extracted[[index_list]][3:7], dataModel_inegral[[index_list]])
  data_joined$date = IdentiyImageDate(names(optical_data_extracted)[index_list])
  return(data_joined)
}))

ggplot(data_joined_[data_joined_$type2 %in%cover_types,], aes(type, data_joined_[data_joined_$type2 %in%cover_types,][,6]))+geom_boxplot() +
  facet_wrap(vars(date), ncol = 1)
