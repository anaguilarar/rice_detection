################ 
### 
###   Split Spatial information in training and Validations subsets for phenological classification
### 
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
source("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/_scripts/R_Main_Functions_GrowCropIndentification.R")


##### Load Libraries

libs=c("snowfall","caret","nnet","SDMTools","stringr","raster","ff", "rgeos","rgdal")
PackageReading(libs)
##### Set Locality Name

Locality = "saldana"


### Set folders path

Main_Folder = paste0("D:/phen_identification/")

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
setwd(Main_Folder)


field_toDelete = ""



# Select a raster reference. it must be a Sentinel 2 image

SentinelimageReference = raster("process/sentine_image_ref/saldana_rasterref.tif")
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
    
    shpfile = TransformAndBuffer(shapetoTransform=shpfile,buffer_dist =  -30)
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
phen_ident_data=read.csv("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/field_data/training_fields/growthphases_campaigns.csv")
phen_ident_data$Field = as.factor(as.character(phen_ident_data$Field))

levels(phen_ident_data$Field )


sum(round(table(phen_ident_data$state)*0.7))
table(phen_ident_data$state)
####

training_per = 70 ## set percentage for training data

### organice phen data training

cover_type = "Rice"
Rice_SpatialPoints = SpatialPoints[[cover_type]]


### 
setwd(paste0(Main_Folder,"data"))

## number of cross validation sub sets
number_random = 3
run = 1
phen_cat = as.character(unique(phen_ident_data$state))[6]

dataRice = lapply(1:number_random , function(run){
  
  CrossValInfo = lapply(as.character(unique(phen_ident_data$state)), function(phen_cat){
    ## Select category
    cat(phen_cat , "\n")
    ## select those fields which belong to the category
    fieldNames = phen_ident_data[phen_ident_data$state%in% phen_cat,]
    fieldNames$Field = factor(as.character(fieldNames$Field))
    
    fieldNames$Date = factor(as.character(fieldNames$Date))
    feildDate = paste(fieldNames$Field,fieldNames$Date)
    ## create subsets for cross validations
    if(nrow(fieldNames)>7){
      field_perepoch = rep("field",8)
    }else{
      field_perepoch =rep("field",nrow(fieldNames))
    }
    set.seed(123)
    partitions = createDataPartition(y= field_perepoch, p =(training_per/100) , times = number_random, list=F)
    
    if(dim(fieldNames)[1]==5){
      set.seed(60)
      partitions = createDataPartition(y= field_perepoch, p =(0.7 ) , times = number_random, list=F)
      
      n_field = unname(sapply(c("62D835A"), function(x) grep (fieldNames$Field , pattern =x)))
      
      x= partitions[,2]
      fields_number = length(fieldNames$Field)
      partitions = apply(partitions,2, function(x) if(!T%in%(x%in%n_field)){
          validate = (1:fields_number)[!(1:fields_number)%in%x]
          validate = validate[!validate%in%n_field]
          set.seed(1990)
          x = x[!x%in%sample(x,1)]
          x = c(x,n_field)
        }else{
          return(x)
        })
      
    }
    if(dim(fieldNames)[1]==6){
      
      set.seed(1992)
      partitions = createDataPartition(y= field_perepoch, p =(60/100 ) , times = number_random, list=F)
      
      n_field = unname(sapply(c("11B011A"), function(x) grep (fieldNames$Field , pattern =x)))
      if(!is.list(n_field)){
        x= partitions[,1]
        fields_number = length(fieldNames$Field)
        partitions = apply(partitions,2, function(x) if(!T%in%(x%in%n_field)){
          validate = (1:fields_number)[!(1:fields_number)%in%x]
          validate = validate[!validate%in%n_field]
          set.seed(1992)
          x = x[!x%in%sample(x,1)]
          x = c(x,n_field)
        }else{
          return(x)
        })
      }
  
    }
    
    fieldNames$Field = as.character(fieldNames$Field)
    ## export coordinates
    print(partitions)
    dataTraining_per_Iteration = Rice_SpatialPoints[Rice_SpatialPoints$polygon_name %in%fieldNames$Field [partitions[,run]],]

    infostage = fieldNames [partitions[,run],]
    Iterations_training = do.call(rbind,lapply(split(dataTraining_per_Iteration, dataTraining_per_Iteration$polygon_name) , function(datapoints){
      cat(unique(datapoints$polygon_name), "\n")
      datapoints$Type_Stage  = as.character(phen_cat)
      
      data_state = as.character(infostage[infostage$Field %in% unique(datapoints$polygon_name),"Date"])
      if(!length(data_state) > 1){
        datapoints$Date  = data_state
      }else{
        datapoints = do.call(rbind,lapply(data_state, function(date){
          datapoints$Date = date
          return(datapoints)
        }))
      }
      if(((phen_cat %in% c( "harvested","late_vegetative","soil")) & nrow(datapoints)>600)){
        datapoints = datapoints[sample(1:nrow(datapoints),350),]
        print(nrow(datapoints))
      }
      if(((phen_cat %in% c("reproductive","ripening")) & nrow(datapoints)>500)){
        datapoints = datapoints[sample(1:nrow(datapoints),300),]
        print(nrow(datapoints))
      }
      
      if(nrow(datapoints)>800){
        datapoints = datapoints[sample(1:nrow(datapoints),650),]
        print(nrow(datapoints))
      } 
      return(datapoints)
    }))
    
    
        ### Export TRaining data base as spatial points
       
    ### create Validation data set
    dataVal_per_Iteration = Rice_SpatialPoints[Rice_SpatialPoints$polygon_name %in%fieldNames$Field [-partitions[,run]],]
    infostage = fieldNames [-partitions[,run],]
    Iterations_validation = do.call(rbind,lapply(split(dataVal_per_Iteration, dataVal_per_Iteration$polygon_name) , function(datapoints){
      datapoints$Type_Stage  = as.character(phen_cat)
      data_state = as.character(infostage[infostage$Field %in% unique(datapoints$polygon_name),"Date"])
      
      if(!length(data_state) > 1){
        datapoints$Date  = data_state
      }else{
        datapoints = do.call(rbind,lapply(data_state, function(date){
          datapoints$Date = date
          return(datapoints)
        }))
      }
      if(((phen_cat %in% c( "harvested","late_vegetative","soil")) & nrow(datapoints)>600)){
        datapoints = datapoints[sample(1:nrow(datapoints),350),]
        print(nrow(datapoints))
      }
      if(((phen_cat %in% c("reproductive")) & nrow(datapoints)>500)){
        datapoints = datapoints[sample(1:nrow(datapoints),350),]
        print(nrow(datapoints))
      }
      if(((phen_cat %in% c("ripening")) & nrow(datapoints)>400)){
        datapoints = datapoints[sample(1:nrow(datapoints),200),]
        print(nrow(datapoints))
      }
      
      if(nrow(datapoints)>800){
        datapoints = datapoints[sample(1:nrow(datapoints),650),]
        print(nrow(datapoints))
      } 
      
      return(datapoints)
    }))

    cat("\t training fields:",length(fieldNames$Field [partitions[,run]]))
    cat("training size:",nrow(Iterations_training))
    cat("\tvalidation fields:",length(fieldNames$Field [-partitions[,run]]),"\n")

    
    return(list(Iterations_training , Iterations_validation ))
  })

  Iterations_training = do.call(rbind,lapply(CrossValInfo, function(sp_info)return(sp_info[[1]])))
  
  Iterations_validation = do.call(rbind,lapply(CrossValInfo, function(sp_info)return(sp_info[[2]])))

  return(list(Iterations_training , Iterations_validation))
 
})
lapply(dataRice, function(x) lapply(x, function(y) table(y[,5])))

lapply(dataRice , function(CrossValInfo){
  lapply(CrossValInfo, function(sp_info)
    sp_info[is.na(sp_info$x),])})
#### Organice information of Soil and Other Vegetation 

cover_types = c("Water" ,  "Urban_Zones", "Roads","Bare_Soil","other")

c_type = "Water"
data_otherTypes = lapply(1:number_random , function(run){
  
  CrossValInfo = lapply(cover_types, function(c_type){
    cat(c_type , "\n")
    
    subset = SpatialPoints[[c_type]]
    dateInt = "20151228"
    if(c_type %in% c("Bare_Soil")){
      cov = "soil"
    }else{
      cov = c_type
    }
    
    ## create subsets for cross validations
    set.seed(1992)
    partitions = createDataPartition(y= subset[,4], p =(training_per/100) , times = number_random, list=F)
    
    points_training = subset[partitions[,run],]
    points_validation = subset[-partitions[,run],]
    
    if(c_type %in% c("other")){
      points_training = points_training[sample(1:nrow(points_training),800),]
      points_validation = points_validation[sample(1:nrow(points_validation),400),]
    }
    
    points_training$Type_Stage = cov
    points_training$Date = dateInt

    points_validation$Type_Stage = cov
    points_validation$Date = dateInt
    return(list(points_training , points_validation))

    
  })
  Iterations_training = do.call(rbind,lapply(CrossValInfo, function(sp_info)return(sp_info[[1]])))
  Iterations_validation = do.call(rbind,lapply(CrossValInfo, function(sp_info)return(sp_info[[2]])))
  return(list(Iterations_training , Iterations_validation))
})

lapply(data_otherTypes , function(CrossValInfo){
  lapply(CrossValInfo, function(sp_info)
    sp_info[is.na(sp_info$x),])})

lapply(data_otherTypes, function(x) lapply(x, function(y) table(y[,5])))

### Export Data

run= 1

setwd("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/model_inputs/phen_identification/parititions")
lapply(1:number_random , function(run){
  
  Iterations_training_Rice = dataRice[[run]][[1]]
  Iterations_training_NotRice = data_otherTypes[[run]][[1]]

  Iterations_training = rbind(Iterations_training_Rice, Iterations_training_NotRice)
  
  Iterations_training = SpatialPointsDataFrame(Iterations_training[,1:2],data = Iterations_training,
                                             proj4string = crs(SentinelimageReference))
  
  Iterations_training@data[is.na(Iterations_training@data$y),]
  nameFile = paste0("v8training_iterationCV_" , run)
  writeOGR (Iterations_training , dsn = "training_phenIdenti" , layer = nameFile , driver = "ESRI Shapefile" , overwrite_layer = T)
  
  ####
  Iterations_validation_Rice = dataRice[[run]][[2]]
  Iterations_validation_NotRice = data_otherTypes[[run]][[2]]
  
  Iterations_validation = rbind(Iterations_validation_Rice, Iterations_validation_NotRice)
  
  Iterations_validation = SpatialPointsDataFrame(Iterations_validation[,1:2],data = Iterations_validation,
                                               proj4string = crs(SentinelimageReference))

  nameFile = paste0("v8validation_iterationCV_" , run)
  writeOGR (Iterations_validation , dsn = "validation_phenIdenti" , layer = nameFile , driver = "ESRI Shapefile" , overwrite_layer = T)

})

