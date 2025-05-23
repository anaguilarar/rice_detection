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
count = 2
######## Functions
setwd("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/")

source(paste0("_scripts/R_Main_Functions_GrowCropIndentification.R"))


source("_scripts/process_main_function.R")


##### Load Libraries

libs=c("stringr","dplyr","caret")
PackageReading(libs)


################# ---> Read Spatial data information
training_data = read.csv(paste0("model_inputs/phen_identification/optical_data_ndvi_date_deri_veg.csv"),row.names =1)

training_data[grepl(row.names(training_data) , pattern =  "p_52B002_s_harvested_d20151228"),'class'] = 'soil'

table(training_data$class)
head(training_data)
## add columns

num_pos = regexpr(pattern = "_p_",row.names(training_data))
num_pos_end = regexpr(pattern = "_s_",row.names(training_data))
training_data$field_name = str_sub(row.names(training_data),num_pos+3, num_pos_end-1)


num_pos = regexpr(pattern = "_d2",row.names(training_data))

training_data$date= str_sub(row.names(training_data),num_pos+2, num_pos+10)


####
training_per = 70 ## set percentage for training data


phen_ident_data = training_data[training_data$class%in%
                                  c("vegetative","reproductive","soil","ripening","harvested"),]

## number of cross validation sub sets
kfolds = 6

sample_phendata = function(sample_data, validation = FALSE){
  lapply(sample_data, 
         function(datapoints){
           cat(unique(datapoints$field_name), "\n")
           
           data_state = unique(datapoints$date)
           if(length(data_state) > 1){
             print("%%%%%%%%%%%error############$$$$$$$$$$$")
           }
           if(((unique(datapoints$class) %in% c( "harvested")) & nrow(datapoints)>600)){
             datapoints = datapoints[sample(1:nrow(datapoints),350),]
             print(nrow(datapoints))
           }
           if(((unique(datapoints$class) %in% c("reproductive")) & nrow(datapoints)>600)){
             datapoints = datapoints[sample(1:nrow(datapoints),350),]
             print(nrow(datapoints))
           }
           if(((unique(datapoints$class) %in% c("ripening")) & nrow(datapoints)>450)){
             datapoints = datapoints[sample(1:nrow(datapoints),250),]
             print(nrow(datapoints))
           }
           
           if(nrow(datapoints)>800){
             datapoints = datapoints[sample(1:nrow(datapoints),700),]
             print(nrow(datapoints))
           } 
           if(validation){
             if(((unique(datapoints$class) %in% c("ripening")))& nrow(datapoints)>400){
               datapoints = datapoints[sample(1:nrow(datapoints),220),]
             }
             if(((unique(datapoints$class) %in% c("vegetative")))& nrow(datapoints)>600){
               datapoints = datapoints[sample(1:nrow(datapoints),350),]
             }
           }
           datapoints$id = row.names(datapoints)
           return(datapoints)
           
         })
}
run =1
phen_cat = as.character(unique(phen_ident_data$class))[5]

table(phen_ident_data$field_name,phen_ident_data$class)
resample_partitions = function(partitions, n_field,fields_number){
  for(i in 1:length(n_field)){
    n_field_i = n_field[i]
    partitions = apply(partitions,2, function(x) if(!T%in%(x%in%n_field_i)){
      validate = (1:fields_number)[!(1:fields_number)%in%x]
      validate = validate[!validate%in%n_field_i]
      set.seed(123)
      x = x[!x%in%sample(x,1)]
      x = c(x,n_field_i)
    }else{
      return(x)
    })
    
  }
  return(partitions)
}
kfolds= 6
dataRice = lapply(1:kfolds , function(run){
  
  CrossValInfo = lapply(as.character(unique(phen_ident_data$class)), function(phen_cat){
    ## Select category
    cat(phen_cat , "\n")
    ## select those fields which belong to the category
    fieldNames = phen_ident_data[phen_ident_data$class%in% phen_cat,]
    fieldNames$id = (paste0("s_",fieldNames$class , "_f_",fieldNames$field_name, 
                            "_d",fieldNames$date))
    
    subset = unique(fieldNames$id )
    
    ## create subsets for cross validations

    field_perepoch =rep("field",length(subset))

    set.seed(123)
    partitions = createDataPartition(y= field_perepoch, p =(training_per/100) , times = kfolds, list=F)
    n_field = unlist(unname(sapply(c("62D835A",'11B011A'), 
                            function(x) grep (subset, pattern =x))))
    
    fields_number = length(subset)

    if(length(n_field)>0){
      resample_partitions(partitions,n_field, fields_number )
    }

    
    if((phen_cat)=="soil"){

      n_field = unlist(unname(sapply(c("52D618","21D823"), 
                              function(x) grep (subset, pattern =x))))
      
      resample_partitions(partitions,n_field, fields_number )
     
    }
    # if(length(subset)==6 | (phen_cat == "early_vegetative")){
    #   set.seed(123)
    #   partitions = createDataPartition(y= field_perepoch, p =(0.7 ) , times = kfolds, list=F)
    #   
    #   n_field = unname(sapply(c("62D835A"), 
    #                           function(x) grep (subset, pattern =x)))
    #   
    #   x= partitions[,2]
    #   fields_number = length(subset)
    #   partitions = apply(partitions,2, function(x) if(!T%in%(x%in%n_field)){
    #     validate = (1:fields_number)[!(1:fields_number)%in%x]
    #     validate = validate[!validate%in%n_field]
    #     set.seed(60)
    #     x = x[!x%in%sample(x,1)]
    #     x = c(x,n_field)
    #   }else{
    #     return(x)
    #   })
    #   
    # }
    # if(phen_cat == "harvested"){
    # 
    #   n_field = unname(sapply(c("52B404A"), 
    #                           function(x) grep (subset, pattern =x)))
    #   
    #   x= partitions[,2]
    #   fields_number = length(subset)
    #   partitions = apply(partitions,2, function(x) if(!T%in%(x%in%n_field)){
    #     validate = (1:fields_number)[!(1:fields_number)%in%x]
    #     validate = validate[!validate%in%n_field]
    #     set.seed(60)
    #     x = x[!x%in%sample(x,1)]
    #     x = c(x,n_field)
    #   }else{
    #     return(x)
    #   })
    # }
    # if(length(subset)==6 & !(phen_cat == "early_vegetative")){
    #   
    #   set.seed(1992)
    #   partitions = createDataPartition(y= field_perepoch, p =(60/100 ) , times = kfolds, list=F)
    #   
    #   n_field = unname(sapply(c("11B011A"), function(x) 
    #     grep (subset, pattern =x)))
    #   
    #   if(!is.list(n_field)){
    #     x= partitions[,1]
    #     fields_number = length(subset)
    #     partitions = apply(partitions,2, function(x) if(!T%in%(x%in%n_field)){
    #       validate = (1:fields_number)[!(1:fields_number)%in%x]
    #       validate = validate[!validate%in%n_field]
    #       set.seed(1992)
    #       x = x[!x%in%sample(x,1)]
    #       x = c(x,n_field)
    #     }else{
    #       return(x)
    #     })
    #   }
    #   
    # }
    
    fieldNames$field_name = as.character(fieldNames$field_name)
    ## export coordinates
    print(partitions)
    for (x in 1:ncol(partitions)){
      for(j in (1:ncol(partitions))[-x]){
        if((unique(partitions[,x] == partitions[,j])[1] == T) & (length(unique(partitions[,x] == partitions[,j]))==1)){

          partitions[,x] = sample(1:fields_number, nrow(partitions))
        }
        
      }
      
    }
    for (x in 1:ncol(partitions)){
      for(j in (1:ncol(partitions))[-x]){
        if((unique(partitions[,x] == partitions[,j])[1] == T) & (length(unique(partitions[,x] == partitions[,j]))==1)){
          
          partitions[,x] = sample(1:fields_number, nrow(partitions))
        }
        
      }
      
    }
    print(partitions)
    ## filter by kposition and get training data
    fieldNames_traininig = fieldNames[fieldNames$id %in% subset[partitions[,run]],]
    
    
    dataTraining_per_Iteration = split(fieldNames_traininig, fieldNames_traininig$id)     
    
    train_data = do.call(rbind,sample_phendata(dataTraining_per_Iteration))
    row.names(train_data) = train_data$id
    
    ## filter by kposition and get validation data
    fieldNames_val = fieldNames[!fieldNames$id %in% subset[partitions[,run]],]
    
    
    dataTraining_per_val = split(fieldNames_val, fieldNames_val$id)     
    
    val_data = do.call(rbind,sample_phendata(dataTraining_per_val,validation = T))
    row.names(val_data) = val_data$id
    
    cat("\t training fields:",length(unique(as.character(train_data$field_name))))
    cat("training size:",nrow(train_data))
    cat("\tvalidation fields:",length(unique(as.character(val_data$field_name))),"\n")
    cat("validation size:",nrow(val_data),"\n")
    
    return(list(train_data , val_data ))
  })
  
  Iterations_training = do.call(rbind,lapply(CrossValInfo, function(sp_info)return(sp_info[[1]])))
  Iterations_training$dataset ="training"
  
  Iterations_validation = do.call(rbind,lapply(CrossValInfo, function(sp_info)return(sp_info[[2]])))
  Iterations_validation$dataset ="validation"
  
  
  
  datarice = rbind(Iterations_training,Iterations_validation)
  row.names(datarice) = paste0(row.names(datarice),"_k",run)
  datarice$id =NULL
  datarice$field_name =NULL
  datarice$date =NULL
  datarice$kfold =run
  
  return(datarice)
  
})

lapply(dataRice, function(x) table(x$class,x$dataset))


#### Organice information of Soil and Other Vegetation 



cover_types = c("Water" ,  "Urban_Zones", "other")


other_data = training_data[training_data$class%in%
                             cover_types,]

posother = which(other_data$class%in% "other")
posnoother = which(!other_data$class%in% "other")
if(length(posother)>4500){
  set.seed(1992)
  other_data= rbind(other_data[sample(posother, 3100),],other_data[posnoother,])

}
set.seed(123)
partitions = createDataPartition(y= rep("other",nrow(other_data)), 
                                 p =(training_per/100) , times = kfolds, list=F)


data_otherTypes = lapply(1:kfolds , function(run){
  
  ## create subsets for cross validations
  
  points_training = other_data[partitions[,run],]
  points_training$dataset = "training"
  
  points_validation = other_data[-partitions[,run],]
  points_validation$dataset = "validation"
  
  otherjoined = rbind(points_training, points_validation)
  otherjoined$field_name =NULL
  otherjoined$date =NULL
  otherjoined$kfold =run
  otherjoined$class = "other"
  return(otherjoined)
  
})
lapply(data_otherTypes, function(x) table(x$class,x$dataset))

### Export Data

## stack data
dataRice = do.call(rbind,dataRice)
data_otherTypes= do.call(rbind,data_otherTypes)

data_training = rbind(dataRice,data_otherTypes)

data_training$class = as.character(data_training$class)
###
table(data_training$kfold,data_training$class,data_training$dataset)

write.csv( data_training,
          paste0("model_inputs/phen_identification/optical_data_ndvi_date_deriveg_6folds_tssum.csv"))


