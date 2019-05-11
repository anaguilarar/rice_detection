##### validate data



###########
library(reshape2)
library(ggplot2)
library(stringr)
data_test = read.csv("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/test_data/valledupar/Ground_truth_data_rice_2018_valle.csv")

red_data_test = data_test[,!names(data_test)%in%c("departamento","cultivo","lote","latitud","longitud","variedad","num","id_tile","area")]

data_fields = melt(red_data_test,id.vars = "id")
data_fields$value = as.Date(data_fields$value,format = "%m/%d/%Y")
data_fields$id = paste0(data_fields$id ,"_2018")


data_growth_stages  = read.csv("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/test_data/valledupar/dates_classified_col_t2_xgboost_veg.csv",row.names = 1)

id_field = "vall_el diamante lote 25_2018"

lapply(unique(data_fields$id)[grepl(unique(data_fields$id), pattern = "diaman")], function(id_field){
  
  data_perfield = data_growth_stages[data_growth_stages$field_name %in% paste0(id_field),]
  data_perfield = data_perfield[, !names(data_perfield)%in% c("nodata" )]
  
  
  data_fieldsplot = data_fields[data_fields$id %in% id_field,] 
  # data_fieldsplot = data_fields[data_fields$id %in% "vall_el diamante lote 4a_2018",] 
  # data_fieldsplot = data_fields[data_fields$id %in% "vall_el diamante lote 5_2018",]
  
  
  data_perfield_plot= melt(data_perfield, id.vars = c("field_name","date","tile"))
  
  data_perfield_plot$date = as.Date(as.character(data_perfield_plot$date), format = "%Y%m%d")
  data_fieldsplot$percentage = max(data_perfield_plot$value)
  
  
  m = ggplot()+
    geom_col(data = data_perfield_plot,aes(date,value, fill = variable),position = "stack", width = 1.5)+
    scale_fill_manual(values = c ('harvested'='firebrick3', "ripening"="goldenrod", "other" = "gray", 
                                  "early_vegetative"="cyan4", "reproductive" = "darkgreen","vegetative" = "chartreuse3",
                                  "soil" = "saddlebrown"
    ))+geom_vline(xintercept = as.Date(as.character(data_fieldsplot$value)), color = c("blue"))
  m = m + labs(x ="date",y = "Pixels (%)", fill = "Classification")+ theme_bw() +ggtitle(str_sub(id_field,6)) +
    geom_text( data=data_fieldsplot, mapping=aes(x=as.Date(as.character(data_fieldsplot$value)), y=0, label=variable), size=3, angle=90, vjust=-0.4, hjust=0)
  print(m)
  
})



data_test = read.csv("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/test_data/valledupar/Ground_truth_data_rice_2018_valle.csv")

red_data_test = data_test[,!names(data_test)%in%c("departamento","cultivo","lote","latitud","longitud","variedad","num","id_tile","area")]

data_fields = melt(red_data_test,id.vars = "id")

data_fields$value = as.Date(data_fields$value,format = "%m/%d/%Y")
data_fields$id = paste0(data_fields$id ,"_2018")


data_growth_stages  = read.csv("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/test_data/valledupar/dates_classified_col_t1_xgboost_veg.csv",row.names = 1)


lapply(unique(data_fields$id)[!grepl(unique(data_fields$id), pattern = "diaman")], function(id_field){
  
  data_perfield = data_growth_stages[data_growth_stages$field_name %in% paste0(id_field),]
  data_perfield = data_perfield[, !names(data_perfield)%in% c("nodata" )]
  
  
  data_fieldsplot = data_fields[data_fields$id %in% id_field,] 
  # data_fieldsplot = data_fields[data_fields$id %in% "vall_el diamante lote 4a_2018",] 
  # data_fieldsplot = data_fields[data_fields$id %in% "vall_el diamante lote 5_2018",]
  
  
  data_perfield_plot= melt(data_perfield, id.vars = c("field_name","date","tile"))
  
  data_perfield_plot$date = as.Date(as.character(data_perfield_plot$date), format = "%Y%m%d")
  data_fieldsplot$percentage = max(data_perfield_plot$value)
  
  
  m = ggplot()+
    geom_col(data = data_perfield_plot,aes(date,value, fill = variable),position = "stack", width = 1.5)+
    scale_fill_manual(values = c ('harvested'='firebrick3', "ripening"="goldenrod", "other" = "gray", 
                                  "early_vegetative"="cyan4", "reproductive" = "darkgreen","vegetative" = "chartreuse3",
                                  "soil" = "saddlebrown"
    ))+geom_vline(xintercept = as.Date(as.character(data_fieldsplot$value)), color = c("blue"))
  m = m + labs(x ="date",y = "Pixels (%)", fill = "Classification")+ theme_bw() +ggtitle(str_sub(id_field,6)) +
    geom_text( data=data_fieldsplot, mapping=aes(x=as.Date(as.character(data_fieldsplot$value)), y=0, label=variable), size=3, angle=90, vjust=-0.4, hjust=0)
  print(m)
  
})



#######

ClassificationMetrics=function(conf_matrix){
  n = sum(conf_matrix) # number of instances
  nc = nrow(conf_matrix) # number of classes
  diag = diag(conf_matrix) # number of correctly classified instances per class 
  rowsums = apply(conf_matrix, 1, sum) # number of instances per class
  colsums = apply(conf_matrix, 2, sum) # number of predictions per class
  p = rowsums / n # distribution of instances over the actual classes
  q = colsums / n # distribution of instances over the predicted classes
  expAccuracy = sum(p*q)
  accuracy = sum(diag) / n
  kappa = (accuracy - expAccuracy) / (1 - expAccuracy)
  precision = diag / colsums 
  recall = diag / rowsums 
  f1 = 2 * precision * recall / (precision + recall) 
  
  return(list(kappa,accuracy,precision, recall, f1))
}

data_test = read.csv("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/test_data/valledupar/Ground_truth_data_rice_2018_valle.csv")

red_data_test = data_test[,!names(data_test)%in%c("departamento","cultivo","lote","latitud","longitud","variedad","num","id_tile","area")]

data_fields = melt(red_data_test,id.vars = "id")
data_fields$value = as.Date(data_fields$value,format = "%m/%d/%Y")

data_splitted = split(data_fields, data_fields$id)

z = 17


data_tocompare = do.call(rbind,lapply(1:length(data_splitted), function(z){
  cat(z)
  data_perfield = data_splitted[[z]]
  cat(as.character(data_perfield$id[1]))
  data_perfield = data_perfield[order(data_perfield$value),]
  
  if(!T%in% is.na(data_perfield$value)){
    data_soil = data.frame(id =data_perfield$id[1] , dates = (data_perfield$value[1]-15):
                             data_perfield$value[data_perfield$variable%in%"fecha_emergencia"],
                           stage = "soil")
    
    data_early = data.frame(id =data_perfield$id[1] , dates = (data_perfield$value[data_perfield$variable%in%"fecha_emergencia"]+1):
                              (data_perfield$value[data_perfield$variable%in%"maximo_macollamiento"]-5),
                            stage = "vegetative")
    
    
    data_late = data.frame(id =data_perfield$id[1] , 
                           dates = (data_perfield$value[data_perfield$variable%in%"maximo_macollamiento"]-4):
                             (data_perfield$value[data_perfield$variable%in%"inicio_floracion"]-10),
                           stage = "vegetative")
    
    data_repro = data.frame(id =data_perfield$id[1] , 
                            dates = (data_perfield$value[(data_perfield$variable%in%"inicio_floracion")]-9):(data_perfield$value[data_perfield$variable%in%"floracion_100"]+14),
                            stage = "reproductive")
    
    data_ripe = data.frame(id =data_perfield$id[1] , 
                           dates = (data_perfield$value[(data_perfield$variable%in%"floracion_100")]+15):(data_perfield$value[data_perfield$variable%in%"cosecha"]+3),
                           stage = "ripening")
    
    data_harvested = data.frame(id =data_perfield$id[1] , 
                                dates = (data_perfield$value[(data_perfield$variable%in%"cosecha")]+4):(data_perfield$value[data_perfield$variable%in%"cosecha"]+20),
                                stage = "harvested")
    
    
    data_tocompare = rbind(data_soil,data_early, data_late, data_repro,data_ripe,data_harvested)
    return(data_tocompare)
  }
  
  
}))




field_name = as.character(unique(data_tocompare$id))[1]
k=1
# phasestages = c("harvested","reproductive","ripening","soil","late_vegetative",
#                 "early_vegetative")
phasestages = c("vegetative", "reproductive","ripening" , "harvested" , "soil")

compare_trueaginstpredicted = function(data_tocompare, data_growth_stages, 
                                       stages =  phasestages){
  
  
  do.call(rbind,lapply(1:length(as.character(unique(data_tocompare$id))), function(k){
    
    field_name = as.character(unique(data_tocompare$id))[k]
    
    cat("\n",field_name ,k,"\n")
    if(!field_name%in%c("San Nicolás/Las 4","tranquilandia las 11")){
      field_name_year = paste0(field_name,"_2018")
      data_perfield = data_growth_stages[data_growth_stages$field_name %in% field_name_year,]
      data_perfield = data_perfield[, !names(data_perfield)%in% c("nodata" )]
      data_perfield_plot= melt(data_perfield, id.vars = c("field_name","date","tile"))
      
      grounddatafield = data_tocompare[data_tocompare$id %in% field_name,]
      stage = "ripening"
      return(do.call(rbind,lapply(stages, function(stage){
        cat(stage)
        
        
        predicted = data_perfield_plot[data_perfield_plot$variable%in%stage,]
        dates_predicted = as.Date(as.character(predicted$date[predicted$value>40]),format = "%Y%m%d")
        true_ground= grounddatafield$stage[grounddatafield$dates %in% as.numeric(dates_predicted)]  
        
        predicted = predicted$variable[predicted$value>40]
        if(length(predicted)>0 & length(true_ground)>0){
          
          if(length(predicted) != length(true_ground)){
            cat(true_ground)
            true_ground = as.character(true_ground)
            predicted = as.character(predicted)
            predicted = predicted[((length(predicted)-length(true_ground))+1):length(predicted)]
          }
          data_metric = data.frame(id =field_name,
                                   true_data = as.character(true_ground),
                                   classification = as.character(predicted))
          return(data_metric)
        }
        
        
      } )))
    }
  })) 
}

data_gs_colt1rf  = read.csv("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/test_data/valledupar/old/dates_classified_col_t1_rf_veg.csv",row.names = 1)
data_gs_colt2rf  = read.csv("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/test_data/valledupar/old/dates_classified_col_t2_rf_veg.csv",row.names = 1)

colt1_rf = compare_trueaginstpredicted(data_tocompare, data_growth_stages = data_gs_colt1rf)
colt2_rf = compare_trueaginstpredicted(data_tocompare, data_gs_colt2rf)

valleduparrf = rbind(colt1_rf, colt2_rf)

valleduparrf$true_data = factor(valleduparrf$true_data, phasestages)


cmrf = table(valleduparrf$true_data, valleduparrf$classification)


cmrf

data_gs_colt1xgboos  = read.csv("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/test_data/valledupar/old/dates_classified_col_t1_xgboost_veg.csv",row.names = 1)
data_gs_colt2xgboos  = read.csv("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/test_data/valledupar/old/dates_classified_col_t2_xgboost_veg.csv",row.names = 1)

colt1_xgboos = compare_trueaginstpredicted(data_tocompare, data_gs_colt1xgboos)
colt2_xgboos = compare_trueaginstpredicted(data_tocompare, data_gs_colt2xgboos)

valleduparxgboos = rbind(colt1_xgboos, colt2_xgboos)
valleduparxgboos$true_data = factor(valleduparxgboos$true_data,phasestages)
valleduparxgboos$classification = factor(valleduparxgboos$classification,phasestages)
cmxgboost = table(valleduparxgboos$true_data, valleduparxgboos$classification)

ClassificationMetrics(cmrf)

ClassificationMetrics(cmxgboost)



data_gs_colt1svm  = read.csv("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/test_data/valledupar/old/dates_classified_col_t1_svm_veg.csv",row.names = 1)
data_gs_colt2svm  = read.csv("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/test_data/valledupar/old/dates_classified_col_t2_svm_veg.csv",row.names = 1)

colt1_svm = compare_trueaginstpredicted(data_tocompare, data_gs_colt1svm)
colt2_svm = compare_trueaginstpredicted(data_tocompare, data_gs_colt2svm)

valleduparsvm= rbind(colt1_svm, colt2_svm)
valleduparsvm$true_data = factor(valleduparsvm$true_data,phasestages)
cmsvm = table(valleduparsvm$true_data, valleduparsvm$classification)

ClassificationMetrics(cmrf)

ClassificationMetrics(cmxgboost)


ClassificationMetrics(cmsvm)

id_field

lapply(unique(data_fields$id)[!grepl(unique(data_fields$id), pattern = "diaman")], function(id_field){
  
  data_perfield = data_growth_stages[data_growth_stages$field_name %in% paste0(id_field),]
  data_perfield = data_perfield[, !names(data_perfield)%in% c("nodata" )]
  
  
  data_fieldsplot = data_fields[data_fields$id %in% id_field,] 
  # data_fieldsplot = data_fields[data_fields$id %in% "vall_el diamante lote 4a_2018",] 
  # data_fieldsplot = data_fields[data_fields$id %in% "vall_el diamante lote 5_2018",]
  
  
  data_perfield_plot= melt(data_perfield, id.vars = c("field_name","date","tile"))
  
  data_perfield_plot$date = as.Date(as.character(data_perfield_plot$date), format = "%Y%m%d")
  data_fieldsplot$percentage = max(data_perfield_plot$value)
  
  
  m = ggplot()+
    geom_col(data = data_perfield_plot,aes(date,value, fill = variable),position = "stack", width = 1.5)+
    scale_fill_manual(values = c ('harvested'='firebrick3', "ripening"="goldenrod", "other" = "gray", 
                                  "early_vegetative"="cyan4", "reproductive" = "darkgreen","vegetative" = "chartreuse3",
                                  "soil" = "saddlebrown"
    ))+geom_vline(xintercept = as.Date(as.character(data_fieldsplot$value)), color = c("blue"))
  m = m + labs(x ="date",y = "Pixels (%)", fill = "Classification")+ theme_bw() +ggtitle(str_sub(id_field,6)) +
    geom_text( data=data_fieldsplot, mapping=aes(x=as.Date(as.character(data_fieldsplot$value)), y=0, label=variable), size=3, angle=90, vjust=-0.4, hjust=0)
  print(m)
  
})



### 
### 

read.csv("file:///D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/graphics/testmetrics.csv")

