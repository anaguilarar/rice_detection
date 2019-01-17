################ 
### 
###   Proccess SAtellte Iamgery 
### 
###                   
###                   
###   Author:  Andr?s Aguilar
###     CIAT - DAPA - AEPS
#########


rm(list=ls())

##### Set Locality Name

tile = "north_tolima"

######### Load functions
server=Sys.info()[1]

serverPath=switch(server, "Linux"="/mnt","Windows"="//dapadfs")

#source(paste0(serverPath,"/workspace_cluster_6/TRANSVERSAL_PROJECTS/MADR/COMPONENTE_2/FEDEARROZ/REMOTE_SENSING/Sentinel_Saldanna/_scripts/r/R_Main_Functions_GrowCropIndentification.R"))

source("D:/phen_identification/_scripts/r/novembrer_2018/R_Main_Functions_GrowCropIndentification.R")

############ Load Libraries

libs=c("snowfall","caret","nnet","SDMTools","stringr","raster","ff", "rgeos")
PackageReading(libs)

## Use a sentinel image as reference to stack both mission radar and optical,
# regardless the images were previously resampled, the sezes are diferent. so a vector file is used 
# to reduce the borders

shpvector = readOGR(paste0("F:/satellite_imagery/refdata/",tile, ".shp"))
sentinelref = raster(paste0("F:/satellite_imagery/refdata/", tile, ".tif"))
shpvector =  spTransform(shpvector, sentinelref@crs)

## check borders


### Set folders path

Sentinel2_Folder=paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/satellite_imagery/s2_processed/")
Sentinel1_Folder=paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/satellite_imagery/s1_processed/")
#Sentinel1_FolderTextures=paste0("satellite_imagery/sentinel_1/",Locality,"_GRDH/9.TextureAnalysis/")
Landsat_Folder=paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/satellite_imagery/landsat_processed/")

optical_folder = paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/satellite_imagery/optical_data/")
radar_folder = paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/satellite_imagery/radar_data/")

### bands names for each mission

# original names
S1_Bands=c("Sigma0_VV_db",'Sigma0_VV_Contrast', 'Sigma0_VV_Dissimilarity',
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

S1QualityPixels=c(100)
S2QualityPixels=list("QL" =9,
  "SLC" = c(0, 1, 3, 8, 9, 10, 11))

L8QualityPixels=list("QL" =c(323),
                     "SLC" = c(8, 208, 144, 100))
L7QualityPixels=list("QL" = c(67))

#### Class information

S1_Param=SatelliteParameters("S1",Sentinel1_Folder,S1_Bands,"tif",S1QualityPixels)
#S1_texture_Param=SatelliteParameters("S1",Sentinel1_FolderTextures,S1Textures_Bands,"img",1000)
S2_Param=SatelliteParameters("L2A",Sentinel2_Folder,S2_Bands,"tif",S2QualityPixels , Standar_bandNames_Sentinel2)
L8_Param=SatelliteParameters("LC08",Landsat_Folder,L8_Bands,"tif",L8QualityPixels , Standar_bandNames_Landsat)
L7_Param=SatelliteParameters("LE07",Landsat_Folder,L7_Bands,"tif",L7QualityPixels , Standar_bandNames_Landsat)


################
##
##get inventory
qpath = "D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/quality_info/"

### Filter images List by Dates

DateStart="2016-01-01"
DateEnd="2017-12-30"


### Process Radar Information

##############
S1_Param


radar_imagery=Read_SatellitalImagery(sat_Params = S1_Param,init_date = DateStart,
                                      end_date = DateEnd, tile = tile)
names_images = names(radar_imagery)
newimages = check_newimages(names_images, qpath, tile )

# update images to process
if(length(newimages) != length(names_images)){
  radar_imagery = radar_imagery[newimages]
  names_images = names(radar_imagery)
}

cropinidicator = "None"

if(length(radar_imagery)>0){
  
  
  ### Crop images
  if(length(unique(sapply(radar_imagery,ncell))) != ncell(sentinelref)){
    cat("The images differs in their extent, clip process is necessary to do \n")
    
    if (ncell(radar_imagery[[1]])<ncell(sentinelref)){
      refextent = extent(radar_imagery[[1]])
      cropinidicator = "optical"
    }else{
      refextent = extent(sentinelref)
      cropinidicator = "radar"
    }
    
    # coerce to a SpatialPolygons object
    ref_pol = as(refextent, 'SpatialPolygons') 
    crs(ref_pol) = sentinelref@crs
  }
  
  ### Delete Manually any image 
  
  ### Export radar information
  radar_path = paste0( "D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/satellite_imagery/radar_data/")
  dir.create(paste0(radar_path,"quicklooks"))
  
  lapply(1:length(radar_imagery),function(count){
    imag=radar_imagery[[count]]
    nameImage=names(radar_imagery)[count]
    if(cropinidicator == "radar"){
      imag = crop(imag, ref_pol)
    }
    
    tif_tordata(imag,nameImage,tile,
                radar_path)
    
    ## export quick true color view
    PlotRGBImage(imag , nameImage , OutputDirectoty = paste0(setPathStyle(radar_path),"quicklooks"),
                 bands = c('Sigma0_VV_db','Sigma0_VV_Energy','Sigma0_VV_Entropy'))
    
    ## add new information 
    quality_table = data.frame("Mission_Type" = "S1_VV", 
                               "Date" = as.Date(IdentiyImageDate(nameImage), format = "%Y%m%d"),
                               "BadPixels" = NA , NameImage = nameImage, "Delete" = "")
    
    add_newdata(quality_table,qpath,tile)
  })
}



##############

### Process Optical Information
### Set Vegetation Index list

ListIndex=list(NDVI="(nir-red)/(nir+red)",
               LSWI="(nir-swir1)/(nir+swir1)",
               GNDVI="(green-nir)/(green+nir)",
               BSI="(((2*swir1/(swir1+nir))-((nir/(nir+red))+(green)/(green + swir1) ))/((2*swir1/(swir1+nir))+((nir/(nir+red))+(green)/(green + swir1) )))")



######### Read Satellital imagery

missions=c("S2","L8","L7")
Mission = "L8"

## read images
ImagesInfo=unlist(lapply(missions,function(Mission){

  sat_Params=switch(Mission,"S2"=S2_Param,"L8"=L8_Param,
                    "L7"=L7_Param)
  cat(Mission,"\n")
  return(Read_SatellitalImagery(sat_Params, DateStart , DateEnd, tile = tile))
  
}))

names_images = names(ImagesInfo)
##compare new images with the inventory

newimages = check_newimages(names_images, qpath, tile )
# update images to process
if(length(newimages) != length(names_images)){
  ImagesInfo = ImagesInfo[newimages]
  names_images = names(ImagesInfo)
}

### Remove Clouds and Shadows from images and export as RData

lapply(1:length(ImagesInfo),function(count){
  imag=ImagesInfo[[count]]
  nameImage=names(ImagesInfo)[count]
  cat("Remove Clouds Image:",nameImage," \n")
  ## match conditions with the satellite mission
  
  if(str_sub(nameImage,1,4)=="LE07"){
    mission = "L7"
    Limit=paste0("imag[['QL']][] > ",L7_Param$QualityLimits$QL)

  }else if(str_sub(nameImage,1,4)=="LC08"){
    mission = "L8"
    Limit=paste0("imag[['QL']][] > ",L8_Param$QualityLimits$QL , 
                 " | imag[['SLC']][] %in% c(",paste0(L8_Param$QualityLimits$SLC,collapse = ","), ")")
  }else{
    mission = "S2"
    if(max(imag[['QL']][],na.rm=T)<12 & max(imag[['SLC']][],na.rm=T)>10){
      pos_slc = grep(names(imag), pattern = "SLC")
      pos_ql =grep(names(imag), pattern = "QL")
      names(imag)[pos_slc] = ""
      names(imag)[pos_ql] = "SLC"
      names(imag)[pos_slc] = "QL"
    }
    
    Limit=paste0("imag[['QL']][] > ",S2_Param$QualityLimits$QL , 
                 " | imag[['SLC']][] %in% c(",paste0(S2_Param$QualityLimits$SLC,collapse = ","), ")")
  }
  ImageWitOutClouds = RemoveClouds(imag,Limit,length(names(imag)))
  names(ImageWitOutClouds[[2]]) = names(imag)
  # calculate indices
  imag_indices = calculate_vi(ImageWitOutClouds[[2]], ListIndex)
  cat("Vegetation Indexes Calculated for: ",nameImage,"\n")
  
  if(cropinidicator == "optical"){
    imag_indices = crop(imag_indices, ref_pol)
    cat("Image extension is greater than radar, a clip extent is applied")
  }
  
  ### Export optical information
  
  tif_tordata(imag_indices,nameImage,tile,
              optical_folder)
  
  ## export quick true color view
  PlotRGBImage(imag_indices , nameImage , OutputDirectoty = paste0(setPathStyle(optical_folder),"true_color"))
  
  
  ## add new information 
  quality_table = data.frame("Mission_Type" = mission, 
             "Date" = as.Date(IdentiyImageDate(nameImage), format = "%Y%m%d"),
             "BadPixels" = ImageWitOutClouds[[1]] , NameImage = nameImage, "Delete" = "")

  add_newdata(quality_table,qpath,tile)
  

})

##############

