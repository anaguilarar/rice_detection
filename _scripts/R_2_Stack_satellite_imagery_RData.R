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


source(paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/_scripts/R_Main_Functions_GrowCropIndentification.R"))

############ Load Libraries

libs=c("snowfall","stringr","raster", "rgdal")
PackageReading(libs)


################
## set folders
##get inventory


### Filter images List by Dates

StartDate="2016-01-01"
EndDate="2017-12-30"


## Update Inventory images
missions=c("L2A","LC08","LE07" , "S1_VV")
Mission = "S1_VV"

for(Mission in missions){
  
  sat_Params=get_sat_list_parameters( Mission)
  cat(Mission,"\n")
  update_inventory(sat_Params,tile,StartDate,EndDate)
  
}

### todo implement update sentinel 1 iventory
# create quicklooks

quicklooks_path = "D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/satellite_imagery/quicklooks/"

create_quicklooks(inventory = get_inventory(tile),start_date = StartDate, end_date = EndDate, quicklooks_path,
                             band_combination = c("red", "green","blue"))




