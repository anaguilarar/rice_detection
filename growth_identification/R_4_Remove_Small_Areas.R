################ 
### 
### This is a code, that is designed to remove small areas identied as rice
### the requirement to use this code is h
###                   
###                   
###   Author:  Andr?s Aguilar
###     CIAT - DAPA - AEPS
#########

rm(list=ls())

######### Load functions
server=Sys.info()[1]
serverPath=switch(server, "Linux"="/mnt","Windows"="//dapadfs")
source(file = paste0("D:/phen_identification/_scripts/r/novembrer_2018/R_Main_Functions_GrowCropIndentification.R"))


############ Load Libraries
library(raster)
library(stringr)
########### set workspace directories

DirFolPrincipal=paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/")

############## Load rice Layer

setwd(DirFolPrincipal)
RiceLayerName = "north_tolimaxgbTree_conf2_20160310_20161007_optical_db.tif"
RiceLayer=raster(paste0("model_ouputs/rice_layers/",RiceLayerName))
plot(RiceLayer)

############## Select those pixels which have more probability to be rice, more than 30 percent
crop_threshold = 40

RiceLayer[RiceLayer[]<crop_threshold]=1
RiceLayer[RiceLayer[]>=crop_threshold]=2
RiceLayer[is.na(RiceLayer)] = 1
plot(RiceLayer)

############# Set size parameters

### rows per columns, Format "rows x columns"

formatMatrix= "15 x 15"
#NumberMinumofPixels = 90
NumberMinumofPixels = 70

Locality="north_tolima"

############ Remove pixels with a small area

CleanImage=RemoveNoisefromaLayer(RiceLayer = RiceLayer,MatrixSize = formatMatrix,NumberMinumofPixels)

plot(CleanImage)
############ Export the layer

writeRaster(CleanImage, filename=paste0("model_ouputs/rice_layers/Cleaned_",str_sub(RiceLayerName,1,-5),"_",formatMatrix,"_",NumberMinumofPixels,".tif"), format="GTiff",overwrite=T)


############ Fill empties between the areas

######Set size parameters

formatMatrix= "15 x 15"
numberMaxPixels=100

fillLayer=FillEmptiesLayer(CleanImage,formatMatrix,numberMaxPixels)

############ Export layer

plot(fillLayer)
writeRaster(fillLayer, filename=paste0("model_ouputs/rice_layers/Filled_",str_sub(RiceLayerName,1,-5),"_",formatMatrix,"_",numberMaxPixels,".tif"), format="GTiff",overwrite=T)

