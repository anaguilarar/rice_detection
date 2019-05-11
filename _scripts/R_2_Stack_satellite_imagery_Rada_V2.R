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

Locality = "saldana"

######### Load functions
server=Sys.info()[1]

serverPath=switch(server, "Linux"="/mnt","Windows"="//dapadfs")

#source(paste0(serverPath,"/workspace_cluster_6/TRANSVERSAL_PROJECTS/MADR/COMPONENTE_2/FEDEARROZ/REMOTE_SENSING/Sentinel_Saldanna/_scripts/r/R_Main_Functions_GrowCropIndentification.R"))

source("_scripts/r/R_Main_Functions_GrowCropIndentification.R")

############ Load Libraries

libs=c("snowfall","caret","nnet","SDMTools","stringr","raster","ff", "rgeos")
PackageReading(libs)


### Set folders path

Main_Folder = paste0("")
Sentinel_Folder=paste0("satellite_imagery/sentinel_2/L2A_Process_",Locality,"/")
Sentinel1_Folder=paste0("satellite_imagery/sentinel_1/",Locality,"_GRDH/8.ToDB/")
Sentinel1_FolderTextures=paste0("satellite_imagery/sentinel_1/",Locality,"_GRDH/9.TextureAnalysis/")
Landsat7_Folder=paste0("satellite_imagery/landsat_7/",Locality,"_corrected/")
Landsat8_Folder=paste0("satellite_imagery/landsat_8/",Locality,"_corrected/")


### bands names for each mission

# original names
S1_Bands=c("Sigma0_VV_db")

S1Textures_Bands=c("Sigma0_VV_ASM","Sigma0_VV_Dissimilarity","Sigma0_VV_Energy","Sigma0_VV_Contrast",
                   "Sigma0_VV_Entropy","Sigma0_VV_GLCMCorrelation","Sigma0_VV_GLCMVariance",
                   "Sigma0_VV_Homogeneity","Sigma0_VV_MAX","Sigma0_VV_GLCMMean")

S2_Bands = c(  "B2",    "B3",    "B4","B5","B6","B7",  "B8", "B8A", "B11",   "B12","quality_cloud_confidence","quality_scene_classification")
L8_Bands = c( "sr_band2",    "sr_band3",    "sr_band4",    "sr_band5", "sr_band6", "sr_band7",    "pixel_qa",    "sr_aerosol")
L7_Bands = c("sr_band1",   "sr_band2",    "sr_band3",    "sr_band4",    "sr_band5", "sr_band7",    "pixel_qa","sr_cloud_qa")

# Standar names

Standar_bandNames_Landsat = c("blue","green","red","nir","swir1","swir2","CloudQA1","CloudQA2")
Standar_bandNames_Sentinel2 = c("blue","green","red","vrededg1","vrededg2","vrededg3","nir","narrownir","swir1","swir2","CloudQA1","CloudQA2")


##### Quality threshoolds

S1QualityPixels=c(100)
S2QualityPixels=c(9,3,11)
L8QualityPixels=c(323)
L7QualityPixels=c(67)

#### Class information

S1_Param=SatelliteParameters("S1",Sentinel1_Folder,S1_Bands,"img",S1QualityPixels)
S1_texture_Param=SatelliteParameters("S1",Sentinel1_FolderTextures,S1Textures_Bands,"img",1000)
S2_Param=SatelliteParameters("S2",Sentinel_Folder,S2_Bands,"img",S2QualityPixels , Standar_bandNames_Sentinel2)
L8_Param=SatelliteParameters("L8",Landsat8_Folder,L8_Bands,"tif",L8QualityPixels , Standar_bandNames_Landsat)
L7_Param=SatelliteParameters("L7",Landsat7_Folder,L7_Bands,"tif",L7QualityPixels , Standar_bandNames_Landsat)


### Set Vegetation Index list

ListIndex=list(NDVI="(nir-red)/(nir+red)",
               LSWI="(nir-swir1)/(nir+swir1)",
               GNDVI="(nir-green)/(nir+green)",
               BSI="(((2*swir1/(swir1+nir))-((nir/(nir+red))+(green)/(green + swir1) ))/((2*swir1/(swir1+nir))+((nir/(nir+red))+(green)/(green + swir1) )))")


### Filter images List by Dates

DateStart="2015-06-01"
DateEnd="2016-05-30"

##############

### Process Optical Information

##############

######### Read Satellital imagery

missions=c("S2","L8","L7")
ImagesInfo=unlist(lapply(missions,function(Mission){

  sat_Params=switch(Mission,"S2"=S2_Param,"L8"=L8_Param,
                    "L7"=L7_Param)
  cat(Mission,"\n")
  return(Read_SatellitalImagery(sat_Params, DateStart , DateEnd))
  
}))



### Remove Clouds and Shadows from images
NameBandCloud="CloudQA1"
NameBandCloud2='CloudQA2'

MultipleImages_WithOutClouds=lapply(1:length(ImagesInfo),function(count){
  imag=ImagesInfo[[count]]
  nameImage=names(ImagesInfo)[count]
  cat("Remove Clouds Image:",nameImage," \n")
  if(str_sub(nameImage,1,4)=="LE07"){
    Limit=paste0(">",L7_Param$QualityLimits[1])
  }else if(str_sub(nameImage,1,4)=="LC08"){
    Limit=paste0(">",L8_Param$QualityLimits[1])
  }else{
    NameBandCloud2 = which(S2_Param$Standar_bandNames %in% NameBandCloud2 )
    Limit=paste0(">",S2_Param$QualityLimits[1] ," | imag[[",NameBandCloud2,"]][] %in% c(",paste0(S2_Param$QualityLimits[2:3],collapse = ","), ")")
  }
  ImageWitOutClouds = RemoveClouds(imag,NameBandCloud,Limit,length(names(imag)))
  return(ImageWitOutClouds)
})
names(MultipleImages_WithOutClouds)=names(ImagesInfo)




### Calculate Vegetation Indexes

MultipleImages_WithOutClouds_Indexes=lapply(1:length(MultipleImages_WithOutClouds),function(count){
  imageI=MultipleImages_WithOutClouds[[count]]
  
  NewLayers=lapply(ListIndex,function(x){eval(parse(text = Identiffy_Band(names(imageI),x)))})
  
  imageI=stack(c(imageI,NewLayers))
  
  cat("Vegetation Indexes Calculated for: ",names(MultipleImages_WithOutClouds)[count],"\n")
  return(imageI)
})

names(MultipleImages_WithOutClouds_Indexes)=names(MultipleImages_WithOutClouds)

### Export optical information
dir.create(paste0(Main_Folder , "process/rdata/optical_info/"))

lapply(1:length(MultipleImages_WithOutClouds_Indexes), function(imag_index){
  file_name = names(MultipleImages_WithOutClouds_Indexes)[imag_index]
  sub_imag = switch(str_sub(file_name , 1, 4) , 
         "Subs" = "S2", 
         "LC08" = "LC8",
         "LE07" = "LE7")
  file_name = paste0(sub_imag , "_" ,IdentiyImageDate(file_name) , "_" , tolower(Locality) , ".RData")
  imag = MultipleImages_WithOutClouds_Indexes [[imag_index]] 
  save(imag, file = paste0(Main_Folder , "process/rdata/optical_info/",file_name))
})

### 
dir.create(paste0(Main_Folder , "process/exploratory_images/optical_saldana/true_color"))
lapply(names(MultipleImages_WithOutClouds_Indexes) , function(nameImag){
  
  PlotRGBImage(MultipleImages_WithOutClouds_Indexes[[nameImag]] , nameImage = nameImag , OutputDirectoty = "process/exploratory_images/optical_saldana/true_color/")
})


## Remove images list
rm(MultipleImages_WithOutClouds)
rm(MultipleImages_WithOutClouds_Indexes)

##############

### Process Radar Information

##############
sat_Params

SentinelImages=Read_SatellitalImagery(sat_Params = S1_Param,init_date = "2015-07-01",
                                      end_date = "2016-05-20")
SentinelImages_text=Read_SatellitalImagery(sat_Params =S1_texture_Param,init_date = "2015-07-01",
                                              end_date = "2016-05-20")
plot(SentinelImages$Subset_S1A_IW_GRDH_1SDV_20151127_subset_VV)
### Join radar information

ImagesSentinel_Ready=lapply(names(SentinelImages_text),function(nameImage){
  cat(nameImage,"\n")
  stack1=SentinelImages_text[grepl(names(SentinelImages_text),pattern=nameImage)][[1]]
  stack2=SentinelImages[grepl(names(SentinelImages),pattern=nameImage)][[1]]
  return(stack(stack1,stack2))
})
rasterref = NULL
### Reescale Radar images
Satellital_imagery = ImagesSentinel_Ready
names(ImagesSentinel_Ready)=names(SentinelImages_text)
if(length(unique(sapply(ImagesSentinel_Ready,ncell)))>1){
  cat("The images differs in their extent, Cropping process is necessary to do \n")
  rasterref = create_raster_ref(ImagesSentinel_Ready)
}

### Delete Manually any image 

ImagesSentinel_Ready = ImagesSentinel_Ready[!IdentiyImageDate(names(ImagesSentinel_Ready)) %in%"20151110"]


### Export radar information
out_path = paste0(Main_Folder , "process/rdata/radar_info/")
dir.create(out_path)

export_radar_data(ImagesSentinel_Ready, out_path, Locality, rasterref = rasterref,multitask = TRUE, branches = 3)

load("D:/phen_identification/process/rdata/radar_info/S1_VV_Subset_S1A_IW_GRDH_1SSV_20151103_subset_VV_saldana.RData")

