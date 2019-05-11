################ 
### 
###   Transform polygons to spatial grid:
###     
###                   
###                   
###   Author:  Andr?s Aguilar
###     CIAT - DAPA - AEPS
#########


rm(list=ls())

######## Functions

### Load functions
source("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/_scripts/R_Main_Functions_GrowCropIndentification.R")


##### Load Libraries

libs=c("stringr","raster","rgeos","rgdal")
PackageReading(libs)
##### Set Locality Name


### Set folders path

Main_Folder = paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/")

SHPFol=paste0(Main_Folder,"/SHP_Fields/")

################# ---> Read Spatial data information

setwd(SHPFol)

patternExt=".shp"
Filesin=list.files(pattern = paste0("*",patternExt, "$"))

DataShape=lapply(Filesin, readOGR )
names(DataShape)=str_sub(Filesin, 1,(-1*(str_length(patternExt)+1)))

## read real data names

cover_type = read.csv("tipodeCoberturas_2.csv" , stringsAsFactors = F)
cover_type = cover_type[!cover_type$Cover_Type %in% c("Forest"),]

### Extract plots Classification



# Select a raster reference. it must be a Sentinel 2 image
setwd(Main_Folder)
SentinelimageReference = raster("satellite_imagery/refdata/col_t3.tif")
SentinelimageReference[] = NA

### tranform all spatial files to spatial points

DataShape=DataShape[names(DataShape)%in%cover_type$File_Name]
SpatialPoints = lapply(1: length(DataShape) , function(z){
  ## get spatial information
  shpfile = DataShape[[z]]
  ## select the cover type
  cat(names(DataShape)[z],"\n")
  cover = as.character(cover_type[ cover_type[,1] %in%names(DataShape)[z],2])
  ## identify if the spatial inforamtion is polygons or points
  if(is(shpfile)[1]!="SpatialPointsDataFrame"){
    
    shpfile = TransformAndBuffer(shapetoTransform=shpfile,buffer_dist =  -5)
    if(!as.character(crs(shpfile)) == as.character(SentinelimageReference@crs)){
      shpfile = spTransform(shpfile, SentinelimageReference@crs)
    }
    intercepted_points = raster::extract(SentinelimageReference, shpfile,cellnumbers=T)
    names(intercepted_points) = as.character(shpfile@data[,ncol(shpfile@data)])
    dataper_pol = intercepted_points[[1]]
    intercepted_points_coord = do.call( rbind , lapply(names(intercepted_points),function(name_dataper_pol){
      ### 
      dataper_pol = intercepted_points[[name_dataper_pol]]
      ### get coordinates from a raster reference
      intercepted_points_coord  = data.frame(xyFromCell(SentinelimageReference,dataper_pol [,1]))
      
      ## define data frame name
      intercepted_points_coord$polygon_name =  name_dataper_pol
      
      return(intercepted_points_coord)
    }))
  }else{
    if(!as.character(crs(shpfile)) == as.character(SentinelimageReference@crs)){
      shpfile = spTransform(shpfile, SentinelimageReference@crs)
    }
    intercepted_points = raster::extract(SentinelimageReference, shpfile,cellnumbers=T)
    ### get coordinates from a raster reference
    intercepted_points_coord  = data.frame(xyFromCell(SentinelimageReference,intercepted_points[,1]))
    intercepted_points_coord$polygon_name =  NA
  }
  if(length(cover) == 0)stop("there is no landcover assignation for the file name ")
  intercepted_points_coord$type = cover
  return(intercepted_points_coord)
  
})

SpatialPoints = do.call(rbind,SpatialPoints)

### Join per Category
SpatialPoints = split(SpatialPoints , SpatialPoints$type)


#### Extract information and split in training and validation
phen_ident_data=read.csv("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/field_data/training_fields/growthphases_campaigns_veg_v2.csv")
phen_ident_data$Field = as.factor(as.character(phen_ident_data$Field))

table(phen_ident_data$state)
####

### organice phen data training

cover_type = "Rice"
Rice_SpatialPoints = SpatialPoints[[cover_type]]

## number of cross validation sub sets
phen_cat = as.character(unique(phen_ident_data$state))[2]

  
rice_classes = do.call(rbind,lapply(as.character(unique(phen_ident_data$state)), function(phen_cat){
  ## Select category
  cat(phen_cat , "\n")
  ## select those fields which belong to the category
  fieldNames = phen_ident_data[phen_ident_data$state%in% phen_cat,]
  fieldNames$Field = factor(as.character(fieldNames$Field))
  
  ### assign coordinates 
  
  field_name = as.character(fieldNames$Field)[6]
  Rice_SpatialPoints[Rice_SpatialPoints$polygon_name%in%field_name,'x']
  phen_fields = do.call(rbind,lapply(as.character(fieldNames$Field), function(field_name){
    cat("\n\t",field_name)
    ## check date duplicity
    capaign_dates = as.character(fieldNames$Date[fieldNames$Field %in% field_name])
    do.call(rbind,lapply(capaign_dates, function(capaign_date){
      data.frame(
        x = Rice_SpatialPoints[Rice_SpatialPoints$polygon_name%in%field_name,'x'],
        y = Rice_SpatialPoints[Rice_SpatialPoints$polygon_name%in%field_name,'y'],
        field_name = field_name,
        class = phen_cat,
        date = capaign_date)
    }))
  }))
  cat("\ttotal points :",nrow(phen_fields),'\n')
  return(phen_fields)
}))

#### Organice information of Soil and Other Vegetation 

cover_types = c("Water" ,  "Urban_Zones","other","soil")

cove = "Water"

  
otherclass = do.call(rbind,lapply(cover_types, function(cove){
  cat(cove , "\n")
  
  subset = SpatialPoints[[cove]]
  dateInt = "20151228"
  
  ## rename column
  names(subset)[names(subset)%in%"type"] = "class"
  names(subset)[names(subset)%in%"polygon_name"] = "field_name"
  ## assign values
  subset$class = cove
  subset$date = dateInt
  subset$field_name = "none"
  cat("\ttotal points :",nrow(subset),'\n')
  return(subset)
}))


table(rbind(rice_classes,otherclass)$class)
### Export Data

write.csv(rbind(rice_classes,otherclass),"model_inputs/phen_identification/data_coordinates.csv")

