rm(list = ls())
library(xgboost)
library(caret)
library(doParallel)
library(stringr)
library(snowfall)
ClassificationMetrics=function(conf_matrix){
  n = sum(conf_matrix) # number of instances
  nc = nrow(conf_matrix) # number of classes
  diag_cm = diag(conf_matrix) # number of correctly classified instances per class 
  rowsums = apply(conf_matrix, 1, sum) # number of instances per class
  colsums = apply(conf_matrix, 2, sum) # number of predictions per class
  p = rowsums / n # distribution of instances over the actual classes
  q = colsums / n # distribution of instances over the predicted classes
  expAccuracy = sum(p*q)
  accuracy = sum(diag_cm) / n
  kappa = (accuracy - expAccuracy) / (1 - expAccuracy)
  precision = diag_cm / colsums 
  recall = diag_cm / rowsums 
  precision_micro = sum(diag_cm) / sum(colsums)
  
  f1 = 2 * precision * recall / (precision + recall) 
  
  
  
  return(list(kappa,accuracy,precision, recall, f1, precision_micro))
}


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


######

evaluation_table= function(y_observed, y_predicted, list_params){
  
  #### organice evluation metrics
  CM=table(y_observed,y_predicted)
  list_results = ClassificationMetrics(CM)
  kappa = list_results[[1]]
  accuracy = list_results[[2]]
  preci = unlist(list_results[[3]])
  
  names(preci) = paste0("precision_",1:length(preci))
  
  recll = unlist(list_results[[4]])
  names(recll) = paste0("recall_",1:length(recll))
  
  f1sc = unlist(list_results[[5]])
  
  names(f1sc) = paste0("f1scrore_",1:length(f1sc))
  
  f1micro_ =unlist(list_results[[6]])
  

  return(data.frame(do.call(data.frame,list_params),
                           kappa, accuracy, t(preci),t(recll),t(f1sc), 
                           f1micro = f1micro_))
  
  
  
}




get_data= function(training_set, 
         iteration, 
         model = "rf"){
  
  subset = training_set[training_set$kfold %in%iteration,]
  subset$kfold = NULL
  ## 
  training_data = subset[subset$dataset %in% "training",]
  training_data$dataset = NULL
  validation_data = subset[subset$dataset %in% "validation",]
  validation_data$dataset = NULL
  
  
  #

  if(model == "svm_radial" | model == "svm_polynomial"){
    minmaxvalues = get_featuresmin_maxvals(subset[,1:(ncol(subset)-2)])
    dates_pos = grepl("Date", names(minmaxvalues))
    #
    # ###
    #
    min_dates = min(unlist(minmaxvalues[dates_pos]))
    max_dates = max(unlist(minmaxvalues[dates_pos]))
    #
    # ## calculate min max derivatives
    #
    min_der = min(unlist(minmaxvalues[!dates_pos]))
    max_der = max(unlist(minmaxvalues[!dates_pos]))
    #
    # ## feature scaling
    #

    ## scale data
    
    dates_scaled = do.call(cbind,lapply(which(dates_pos),function(x){
      feature_scaling(subset[,x], c(min_dates, max_dates))
    }))
    
    derivatives_scaled = do.call(cbind,lapply(which(!dates_pos),function(x){
      feature_scaling(subset[,x], c(min_der, max_der))
    }))
    #
    # ##
    training_iter_scaled = data.frame(dates_scaled, derivatives_scaled,
                                      class = subset$class,
                                      dataset = subset$dataset)
    
    row.names(training_iter_scaled) = row.names(subset)
    names(training_iter_scaled)[1:(ncol(subset)-2)] =
      names(subset)[1:(ncol(subset)-2)]
    
    #
    training_data = training_iter_scaled[training_iter_scaled$dataset %in% "training",]
    training_data$dataset = NULL
    validation_data = training_iter_scaled[training_iter_scaled$dataset %in% "validation",]
    validation_data$dataset = NULL
    #
    
  }
  x_training = training_data[,-ncol(training_data)]
  x_validation = validation_data[,-ncol(validation_data)]
  
  y_training = training_data[,ncol(training_data)]
  y_training = as.numeric(y_training) - 1
  y_validation = as.numeric(validation_data[,ncol(validation_data)])-1
  
  if(model == "xgboost"){
    
    x_training <- xgboost::xgb.DMatrix(as.matrix(x_training), label = y_training, missing = NA)
    x_validation <- xgboost::xgb.DMatrix(as.matrix(x_validation),  missing = NA)
  }
  
  return(list(training_set = x_training,
              validation_set = x_validation,
              y_training = y_training,
              y_validation = y_validation
              ))
  
}

read_backup=function(params,filename,default_params,path){
  if(file.exists(paste0(path,filename))){
    results = read.csv(paste0(path,filename))
    comb_done = unique(apply(results[,default_params], 1, paste, collapse="_"))
    cat("\t", comb_done , "\n these parameters have already been processed\n")
    ## leave those parameters that were not still tested
    params = params[!params%in%comb_done]
    cat("\t", params , "\n these are left to process\n")
    
  }else{
    print(paste0(filename, " does not exist"))
  }
  return(params)
}

random_forest_class = function (model_data, iterations, list_parameters_rf,
                                read_backup = F,
                     save_backup = T, seed = 1,
                     path = "/tunning_parameters/", nclusters = 2){
  
  default_params = c("mtry", "ntree")
  file_backup = "rf_tunning.csv"
  
  if(unique(is.na(list_parameters_rf))){
    list_parameters_rf = list(
      mtry = round((ncol(model_data)-1)/3),
      ntree =500
    )
    
  }
  
  if(F%in%(default_params%in%names(list_parameters_rf))){
    cat("missing paramaters!!")
    stop()
  }
  
  all_params =  apply(expand.grid(list_parameters_rf), 1, paste, collapse="_")
  
  ## filter those parameters thatwere already processed
  all_params =read_backup(all_params,file_backup, path )
  
  ## random params 
  all_params = str_split(all_params, pattern = "_")
  names(all_params) = as.character(1:length(all_params))
  all_params = all_params[sample(1:length(all_params))]
  
  Sys.time()->start
  sfInit(parallel=T,cpus=nclusters)
  sfLibrary("snowfall", character.only=TRUE)
  sfLibrary(randomForest)
  sfLibrary(stringr)
  
  sfExport("model_data");sfExport("ClassificationMetrics")
  sfExport("all_params");  sfExport("get_data")
  sfExport("save_backup");sfExport("path")
  sfExport("iterations");sfExport("evaluation_table")
  sfExport("get_featuresmin_maxvals");sfExport("feature_scaling")
  sfExport("file_backup")
  
  results = sfLapply( 1:length(all_params), function(params_i){
    results = do.call(rbind,lapply(1:iterations , function(iteration){
      ### organice data
        data_split = get_data(model_data, iteration, "rf")
        mtry_p =as.numeric( all_params[[as.character(params_i)]][1])
        ntree_p = as.numeric(all_params[[as.character(params_i)]][2])
        ## create_model
        set.seed(seed)
        modelRF = randomForest::randomForest(data_split$training_set, 
                                             as.factor(data_split$y_training), 
                                             mtry = mtry_p, 
                                             ntree = ntree_p)
        
        predictedValues = unname(predict(modelRF,data_split$validation_set))
        
        ### select params
        
        evaluation_table(y_observed = data_split$y_validation, y_predicted = predictedValues,
                         list_params = list(mtry = as.numeric(mtry_p), 
                                            ntree = as.numeric(ntree_p)))
      }))
      results$iteration = 1:iterations
      if(save_backup){
        if(file.exists(paste0(path,file_backup))){
          resultsprevious = read.csv(paste0(path,file_backup))
          results = rbind(resultsprevious,results)
        }
        write.csv(results,paste0(path,file_backup),row.names = F)
      }
        
    })
  
  print(Sys.time()-start)
  stopCluster(cl)
      
}


svm_radial_class = function (model_data, iterations, list_parameters_svm,
                                read_backup = F,
                                save_backup = T, seed = 1,
                                path = "/tunning_parameters/", nclusters = 2){
  
  default_params = c("gamma", "cost")
  file_backup = "svmradial_tunning.csv"
  
  if(unique(is.na(list_parameters_svm))){
    list_parameters_svm = list(
      gamma = round((ncol(model_data)-1)/3),
      ntree =500
    )
    
  }
  
  if(F%in%(default_params%in%names(list_parameters_svm))){
    cat("missing paramaters!!")
    stop()
  }
  
  all_params =  apply(expand.grid(list_parameters_svm), 1, paste, collapse="_")
  
  ## filter those parameters thatwere already processed
  all_params =read_backup(all_params,file_backup, path )
  
  ## random params 
  all_params = str_split(all_params, pattern = "_")
  names(all_params) = as.character(1:length(all_params))
  all_params = all_params[sample(1:length(all_params))]
  
  Sys.time()->start
  sfInit(parallel=T,cpus=nclusters)
  sfLibrary("snowfall", character.only=TRUE)
  sfLibrary(randomForest)
  sfLibrary(stringr)
  sfLibrary(xgboost)
  sfLibrary(e1071)
  
  sfExport("model_data");sfExport("ClassificationMetrics")
  sfExport("all_params");  sfExport("get_data")
  sfExport("save_backup");sfExport("path")
  sfExport("iterations");sfExport("evaluation_table")
  sfExport("get_featuresmin_maxvals");sfExport("feature_scaling")
  sfExport("file_backup")
  
  results = sfLapply( 1:length(all_params), function(params_i){
    results = do.call(rbind,lapply(1:iterations , function(iteration){
      ### organice data
        data_split = get_data(model_data, iteration, "svm_radial")
        gamma_p = as.numeric(all_params[[as.character(params_i)]][1])
        cost_p = as.numeric(all_params[[as.character(params_i)]][2])
        
        data_trainingsvm = data.frame(data_split$training_set,y = as.factor(data_split$y_training))
        ## create_model
        set.seed(seed)
        svmfit = svm(y ~ ., data = data_trainingsvm, kernel = "radial",
                     gamma = gamma_p,
                     cost = cost_p, 
                     scale = F)
        
        predictedValues = unname(predict(svmfit,data_split$validation_set))
        
        ### select params
        
        evaluation_table(y_observed = data_split$y_validation, y_predicted = predictedValues,
                         list_params = list(gamma = gamma_p, 
                                            cost = cost_p))
      }))
      results$iteration = 1:iterations
      if(save_backup){
        if(file.exists(paste0(path,file_backup))){
          resultsprevious = read.csv(paste0(path,file_backup))
          results = rbind(resultsprevious,results)
        }
        write.csv(results,paste0(path,file_backup),row.names = F)
      }
      
    })
  
  print(Sys.time()-start)
  sfStop()
  
}



svm_poly_class = function (model_data, iterations, list_parameters_svm,
                                    read_backup = F,
                                    save_backup = T, seed = 1,
                                    path = "/tunning_parameters/", nclusters = 2){
  
  default_params = c("gamma", "cost","degrees","coef0")
  file_backup = "svmpoly_tunning.csv"
  
  if(unique(is.na(list_parameters_svm))){
    list_parameters_svm = list(
      gamma = c(0.05),
      cost = c(2),
      degrees = c(2),
      coef0= c(0.05)
    )
    
  }
  
  if(F%in%(default_params%in%names(list_parameters_svm))){
    cat("missing paramaters!!")
    stop()
  }
  
  all_params =  apply(expand.grid(list_parameters_svm), 1, paste, collapse="_")
  
  ## filter those parameters thatwere already processed
  all_params =read_backup(all_params,file_backup,default_params, path )
  
  ## random params 
  all_params = str_split(all_params, pattern = "_")
  names(all_params) = as.character(1:length(all_params))
  all_params = all_params[sample(1:length(all_params))]
  
  Sys.time()->start
  sfInit(parallel=T,cpus=nclusters)
  sfLibrary("snowfall", character.only=TRUE)
  sfLibrary(randomForest)
  sfLibrary(stringr)
  sfLibrary(xgboost)
  sfLibrary(e1071)
  
  sfExport("model_data");sfExport("ClassificationMetrics")
  sfExport("all_params");  sfExport("get_data")
  sfExport("save_backup");sfExport("path")
  sfExport("iterations");sfExport("evaluation_table")
  sfExport("get_featuresmin_maxvals");sfExport("feature_scaling")
  sfExport("file_backup")
  
  results = sfLapply( 1:length(all_params), function(params_i){
    results = do.call(rbind,lapply(1:iterations , function(iteration){
      ### organice data
      data_split = get_data(model_data, iteration, "svm_polynomial")
      data_trainingsvm = data.frame(data_split$training_set,y = as.factor(data_split$y_training))
      
      ## get coefficients
      gamma_p = as.numeric(all_params[[as.character(params_i)]][1])
      cost_p = as.numeric(all_params[[as.character(params_i)]][2])
      degreep = as.numeric(all_params[[as.character(params_i)]][3])
      coef0p= as.numeric(all_params[[as.character(params_i)]][4])
      
      ## create_model
      set.seed(seed)
      svmfit = svm(y ~ ., data = data_trainingsvm, kernel = "polynomial",
                   degree = degreep,
                   gamma = gamma_p,
                   cost = cost_p,
                   coef0=coef0p, 
                   scale = F)
      
      
      predictedValues = unname(predict(svmfit,data_split$validation_set))
      
      ### select params
      table(data_split$y_validation,predictedValues)
      evaluation_table(y_observed = data_split$y_validation, y_predicted = predictedValues,
                       list_params = list(gamma = as.numeric(gamma_p), 
                                          cost = as.numeric(cost_p),
                                          degrees = as.numeric(degreep),
                                          coef0 = as.numeric(coef0p)))
    }))
    results$iteration = 1:iterations
    if(save_backup){
      if(file.exists(paste0(path,file_backup))){
        resultsprevious = read.csv(paste0(path,file_backup))
        results = rbind(resultsprevious,results)
      }
      write.csv(results,paste0(path,file_backup),row.names = F)
    }
    
  })
  
  print(Sys.time()-start)
 sfStop()
  
}



xgboost_class = function (model_data, iterations, list_parameters_xgboost,
                                  read_backup = F,
                                  save_backup = T, seed = 1,
                                  path = "/tunning_parameters/", nclusters = 2){
  
  default_params = c("eta", "max_depth","gamma","colsample","sub_sample")
  file_backup = "xgboost_tunning.csv"
  
  
  if(unique(is.na(list_parameters_xgboost))){
    list_parameters_xgboost = list(
      eta = c(0.01),
      max_depth = c(4),
      gamma = c(1),
      colsample = c(0.7),
      sub_sample = c(0.7)
    )
    
  }
  
  if(F%in%(default_params%in%names(list_parameters_xgboost))){
    cat("missing paramaters!!")
    stop()
  }
  
  all_params =  apply(expand.grid(list_parameters_xgboost), 1, paste, collapse="_")
  
  ## filter those parameters thatwere already processed
  all_params =read_backup(all_params,file_backup,default_params, path )
  
  ## random params 
  all_params = str_split(all_params, pattern = "_")
  names(all_params) = as.character(1:length(all_params))
  all_params = all_params[sample(1:length(all_params))]
  
  Sys.time()->start
  sfInit(parallel=T,cpus=nclusters)
  sfLibrary("snowfall", character.only=TRUE)
  sfLibrary(randomForest)
  sfLibrary(stringr)
  sfLibrary(xgboost)
  sfLibrary(e1071)
  
  sfExport("model_data");sfExport("ClassificationMetrics")
  sfExport("all_params");  sfExport("get_data")
  sfExport("save_backup");sfExport("path")
  sfExport("iterations");sfExport("evaluation_table")
  sfExport("get_featuresmin_maxvals");sfExport("feature_scaling")
  sfExport("file_backup")
  
  results = sfLapply( 1:length(all_params), function(params_i){
    results = do.call(rbind,lapply(1:iterations , function(iteration){
      ### organice data
      data_split = get_data(model_data, iteration, "xgboost")
      
      ## get coefficients

      eta_p = as.numeric(all_params[[as.character(params_i)]][1])
      mx_p = as.numeric(all_params[[as.character(params_i)]][2])
      gam_p = as.numeric(all_params[[as.character(params_i)]][3])
      col_p = as.numeric(all_params[[as.character(params_i)]][4])
      sub_sample = as.numeric(all_params[[as.character(params_i)]][5])

      
      ## create_model
      set.seed(seed)

      
      xgmodfit = xgboost::xgb.train(list(eta = eta_p,
                                     max_depth = mx_p,
                                     gamma = gam_p,
                                     colsample_bytree = col_p,
                                     min_child_weight = 1,
                                     subsample = sub_sample),
                                data = data_split$training_set,
                                num_class = length(unique(data_split$y_training)),
                                nrounds = 400,
                                objective = "multi:softprob")

      predictedValues = predict(xgmodfit,data_split$validation_set)
      #xgb.plot.deepness(model = out)

      out <- matrix(predictedValues, ncol = 6, byrow = TRUE)
      predictedValues <- apply(out, 1, which.max)-1
      
      ### select params
      
      evaluation_table(y_observed = data_split$y_validation, y_predicted = predictedValues,
                       list_params = list(eta = eta_p,
                                          max_depth = mx_p,
                                          gamma = gam_p,
                                          colsample = col_p,
                                          sub_sample = sub_sample))
    }))
    results$iteration = 1:iterations
    if(save_backup){
      if(file.exists(paste0(path,file_backup))){
        resultsprevious = read.csv(paste0(path,file_backup))
        results = rbind(resultsprevious,results)
      }
      write.csv(results,paste0(path,file_backup),row.names = F)
    }
    
  })
  
  print(Sys.time()-start)
  sfStop()
  
}


###########################################
###############
#

setwd("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/")


data_training=read.csv(paste0("model_inputs/phen_identification/optical_data_ndvi_date_deriveg_5folds.csv"),row.names = 1)
table(data_training$kfold,data_training$class,data_training$dataset)
dim(data_training)
table(data_training$class)




table(data_training$kfold,data_training$class,data_training$dataset)



### delete na values
levelsNA= unique(unlist(apply(data_training ,  2 , function(x) unique(which(is.na(x))))))
if(length(levelsNA) > 0){
  data_training = data_training[-levelsNA,]  
}

head(data_training)
data_training$class = as.character(data_training$class)
#data_training$class[data_training$class%in%c("early_vegetative", "late_vegetative")]="vegetative"
data_training$class = factor(data_training$class, levels = c("vegetative","reproductive",
                                                             "ripening","harvested","soil","other"))
#######################################



list_parrf = list(mtry = c(2:8),
                  ntree = c(100,200,500,800,1200,1800,2400,3000))

random_forest_class(model_data = data_training, iterations = 5,
                    list_parameters_rf = list_parrf,nclusters = 6,
                    path = "classification_models/grid_search/")

list_parsvm_radial = list(gamma = c(0.00001,0.00005,0.0001,0.0005,0.001,0.05,0.1,1,2,4),
                          cost = c(32,64,128,216,512,1048,2056,4112,8224))

svm_radial_class(model_data = data_training, iterations = 5,list_parameters_svm = list_parsvm_radial,nclusters = 6,path = "classification_models/grid_search/")


list_parsvm_polynomial = list(gamma = c(0.05,0.1,1,2,4),
                          cost = c(0.001,0.01,1,2,4,8,16,32,64,128,216,512,1048,2056,4112),
                          degrees = c(2:3),
                          coef0= c(0.05,0,2,4,8))



svm_poly_class(model_data = data_training, iterations = 5,
                        list_parameters_svm = list_parsvm_polynomial,
                        nclusters = 7,path = "classification_models/grid_search/")



list_xgboost_params = list(eta = c(0.1,0.01,0.001),
                           max_depth = c(2,4,8,16,24),
                           gamma = c(0.1,0.5,1,2,4),
                           colsample = c(0.5,0.7,0.9),
                           sub_sample = c(0.7,1))



xgboost_class(model_data = data_training, iterations = 5,
                     list_parameters_xgboost = list_xgboost_params,
                     nclusters = 7,path = "classification_models/grid_search/")




###### train final models


library(stringr)
library(dplyr)
alldata_sub = data_training

alldata_sub = alldata[alldata$kfold%in%4,]
len_id = str_sub(row.names(alldata_sub ),1,-1)

alldata = data_training
num_pos = regexpr(pattern = "_k",row.names(alldata))

id= str_sub(row.names(alldata),1, num_pos )

alldata$ID = id

alldata["dataset"] = NULL
alldata["kfold"] = NULL
dim(alldata)
alldata = data.frame(alldata %>%
                       group_by(ID) %>% filter(row_number(ID) == 1))

dim(alldata)
row.names(alldata) = alldata$ID
alldata["ID"]=NULL

y_variable = factor(as.numeric(alldata$class) - 1)
x_variables = alldata[,-ncol(alldata)]


grid_default <- expand.grid(
  nrounds = 400,
  max_depth = 16,
  eta = 0.01,
  gamma = 2,
  colsample_bytree = 0.5,
  min_child_weight = 1,
  subsample = 0.7
)

train_control <- caret::trainControl(
  method = "none",
  verboseIter = FALSE, # no training log
  allowParallel = TRUE # FALSE for reproducible results 
)

xgb_model <- caret::train(
  x = x_variables,
  y = y_variable,
  trControl = train_control,
  tuneGrid = grid_default,
  method = "xgbTree",
  verbose = TRUE
)


varxgboost <- xgboost::xgb.DMatrix(as.matrix(x_variables), 
                                   label = as.numeric(y_variable)-1, missing = NA) 

eta_p = grid_default$eta
mx_p = grid_default$max_depth
gam_p = grid_default$gamma
col_p = grid_default$colsample_bytree
nround_p = 400
mc_p = 1
sub_sample = grid_default$subsample

set.seed(1)
xgmodfit = xgboost::xgb.train(list(eta = eta_p,
                                   max_depth = mx_p,
                                   gamma = gam_p,
                                   colsample_bytree = col_p,
                                   min_child_weight = 1,
                                   subsample = sub_sample),
                              data = varxgboost,
                              num_class = length(unique(y_variable)),
                              nrounds = 400,
                              objective = "multi:softprob")

ggplot(data_testmetrics, aes(type, value, color = model))+ geom_point(size = 3 )+theme_bw() +
  labs( y = "f1 Score", x ="tipo de clasificación", colour = 'modelos')





list_parameters_rf
do.call("paste", c(list_parameters_rf, sep = "_"))


ggplot(data_training, aes(NDVI_derivative_1,NDVI_derivative_2, color = class)) + 
  geom_point()+xlim(-0.03,0.07)#+ylim(-0.03,0.02)


wrongid = row.names(data_training[(data_training$class == "reproductive") &(data_training$NDVI_derivative_3 < (-0.05)),])

stringr::str_sub(wrongid , 1,49)

############3
iteration = 1

read_backup = T

cl =  makeCluster(5)
registerDoParallel(cl)
Sys.time()->start


model_data =  data_training
quantile_imags = do.call(rbind,foreach(iteration = 1:5, 
                                       .export = c("data_training","ClassificationMetrics"),
                                       .packages=c('randomForest','stringr','xgboost','e1071')) %dopar% {
  
  subset = data_training[data_training$kfold %in%iteration,]
  subset$kfold = NULL
  ## 
  training_data = subset[subset$dataset %in% "training",]
  training_data$dataset = NULL
  validation_data = subset[subset$dataset %in% "validation",]
  validation_data$dataset = NULL
  ########### RANDOM FOREST

  #
  x_training = training_data[,-ncol(training_data)]
  x_validation = validation_data[,-ncol(validation_data)]

  y_training = training_data[,ncol(training_data)]
  y_training = as.numeric(y_training) - 1
  y_validation = as.numeric(validation_data[,ncol(validation_data)])-1

  first = T
  mtry_p = 3
  ntree_p = 200
  for(mtry_p in 3:8){
    for(ntree_p in c(200,500,800,1000,1200,1600,1800,2000,2200,2400))rr{
      cat(mtry_p, " " , ntree_p, " ")

      set.seed(1)
      modelRF = randomForest::randomForest(x_training, as.factor(y_training), 
                                           mtry = mtry_p, ntree = ntree_p)

      predictedValues = unname(predict(modelRF,x_validation))

      # testinfo = (validation_data[!y_validation ==predictedValues,])[
      # grepl(row.names(validation_data[!y_validation ==predictedValues,]), pattern = "harvested"),]
      #
      # table(stringr::str_sub(row.names(testinfo),14,49))
      # total_ = validation_data[
      # grepl(row.names(validation_data), pattern = "harvested"),]
      #    table(stringr::str_sub(row.names(total_),14,49))
      CM=table(y_validation,predictedValues)
      kappa = ClassificationMetrics(CM)[[1]]
      accuracy = ClassificationMetrics(CM)[[2]]
      preci = unlist(ClassificationMetrics(CM)[[3]])
      names(preci) = paste0("precision_",1:length(preci))

      recll = unlist(ClassificationMetrics(CM)[[4]])
      names(recll) = paste0("recall_",1:length(recll))

      f1sc = unlist(ClassificationMetrics(CM)[[5]])

      names(f1sc) = paste0("f1scrore_",1:length(f1sc))

      ClassificationMetrics(CM)
      row_results = data.frame(iteration = iteration,mtry = mtry_p ,ntree =  ntree_p,
                               kappa, accuracy, t(preci),t(recll),t(f1sc), f1micro = unlist(ClassificationMetrics(CM)[[6]]))


      if(first){
        results_table = row_results
        first = F
      }else{
        results_table = rbind(results_table,
                              row_results)
      }

      cat( " kappa_cohen:" , kappa," f1score: ", mean(f1sc),"\n")
      rm(modelRF)

    }

  }
  #
  write.csv(results_table,paste0("classification_models/grid_search/rftuning", iteration, ".csv"))
  #
  # ############# SUPPORT VECTOR MACHINE
  #
  # #
  # #
  # #
  #
  minmaxvalues = get_featuresmin_maxvals(subset[,1:9])
  dates_pos = grepl("Date", names(minmaxvalues))
  #
  # ###
  #
  min_dates = min(unlist(minmaxvalues[dates_pos]))
  max_dates = max(unlist(minmaxvalues[dates_pos]))
  #
  # ## calculate min max derivatives
  #
  min_der = min(unlist(minmaxvalues[!dates_pos]))
  max_der = max(unlist(minmaxvalues[!dates_pos]))
  #
  # ## feature scaling
  #
  first = T
  #


  ## scale data

  dates_scaled = do.call(cbind,lapply(which(dates_pos),function(x){
    feature_scaling(subset[,x], c(min_dates, max_dates))
  }))

  derivatives_scaled = do.call(cbind,lapply(which(!dates_pos),function(x){
    feature_scaling(subset[,x], c(min_der, max_der))
  }))
  #
  # ##
  training_iter_scaled = data.frame(dates_scaled, derivatives_scaled,
                                    class = subset$class,
                                    dataset = subset$dataset)
  # ## reassign names
  #
  row.names(training_iter_scaled) = row.names(subset)
  names(training_iter_scaled)[1:9] =
    names(subset)[1:9]

  #
  training_data = training_iter_scaled[training_iter_scaled$dataset %in% "training",]
  training_data$dataset = NULL
  validation_data = training_iter_scaled[training_iter_scaled$dataset %in% "validation",]
  validation_data$dataset = NULL
  #
  x_training = training_data[,-ncol(training_data)]
  x_validation = validation_data[,-ncol(validation_data)]

  y_training = training_data[,ncol(training_data)]
  y_training = as.numeric(y_training) - 1
  y_validation = as.numeric(validation_data[,ncol(validation_data)])-1

  data_trainingsvm = data.frame(x_training,y = as.factor(y_training))
  #
  #
  #
  # #### RADIAL
  first = T
  for(cost_p in c(32,64,128,216,512,1048,2056,4112,8224)){
    cat(cost_p, '\n')
    for(gamma_p in c(0.00001,0.00005,0.0001,0.0005,0.001,0.05,0.1,1,2,4)){

      set.seed(1)
      svmfit = svm(y ~ ., data = data_trainingsvm, kernel = "radial",gamma = gamma_p,cost = cost_p, scale = F)
      #svmfit = svm(y ~ ., data = data_trainingsvm, kernel = "polynomial",degree = degreep,gamma = gamma_p,cost = cost_p,coef0=coef0p, scale = F)

      predictedValues = unname(predict(svmfit,x_validation))

      CM=table(y_validation,predictedValues)
      kappa = ClassificationMetrics(CM)[[1]]
      accuracy = ClassificationMetrics(CM)[[2]]
      preci = unlist(ClassificationMetrics(CM)[[3]])
      names(preci) = paste0("precision_",1:length(preci))

      recll = unlist(ClassificationMetrics(CM)[[4]])
      names(recll) = paste0("recall_",1:length(recll))

      f1sc = unlist(ClassificationMetrics(CM)[[5]])
      names(f1sc) = paste0("f1scrore_",1:length(f1sc))

      results_scores = data.frame(iteration = iteration,gamma = gamma_p ,cost =  cost_p,kappa, accuracy, t(preci),t(recll),t(f1sc), f1micro = unlist(ClassificationMetrics(CM)[[6]]) )
      if(first){
        results_table = results_scores
        first = F
      }else{
        results_table = rbind(results_table,
                              results_scores)
      }

      cat( " kappa_cohen:" , kappa," f1score: ", mean(f1sc),"\n")
      rm(svmfit)
    }
    write.csv(results_table,paste0("classification_models/grid_search/svm_radialtuning",iteration,".csv"))
  }
  write.csv(results_table,paste0("classification_models/grid_search/svm_radialtuning",iteration,".csv"))
  #
  # ##### polynomial kernel
  #
  gammagrid  = c(0.0001,0.0005,0.001,0.05,0.1,1,2,4)
  costgrid = c(0.001,0.01,1,2,4,8,16,32,64,128,216,512,1048,2056,4112)
  degrees = c(2:3)
  coef0= c(0,0.05,0,2,4,8)
  first = T
  #
  for(degreep in degrees){
    cat(degreep, '***\n')
    for(coef0p in c(0.05,0,2,4,8,16)){
      for(cost_p in costgrid){
        cat(cost_p, '\n')
        for(gamma_p in gammagrid){

          set.seed(1)
          #svmfit = svm(y ~ ., data = data_trainingsvm, kernel = "radial basis",,gamma = gamma_p,cost = cost_p, scale = F)
          svmfit = svm(y ~ ., data = data_trainingsvm, kernel = "polynomial",degree = degreep,gamma = gamma_p,cost = cost_p,coef0=coef0p, scale = F)

          predictedValues = unname(predict(svmfit,x_validation))

          CM=table(y_validation,predictedValues)
          kappa = ClassificationMetrics(CM)[[1]]
          accuracy = ClassificationMetrics(CM)[[2]]
          preci = unlist(ClassificationMetrics(CM)[[3]])
          names(preci) = paste0("precision_",1:length(preci))

          recll = unlist(ClassificationMetrics(CM)[[4]])
          names(recll) = paste0("recall_",1:length(recll))

          f1sc = unlist(ClassificationMetrics(CM)[[5]])
          names(f1sc) = paste0("f1scrore_",1:length(f1sc))

          results_scores = data.frame(iteration = iteration,gamma = gamma_p ,cost =  cost_p, degree = degreep, coef0 = coef0p,kappa, accuracy, t(preci),t(recll),t(f1sc), f1micro = unlist(ClassificationMetrics(CM)[[6]]) )
          if(first){
            results_table = results_scores
            first = F
          }else{
            results_table = rbind(results_table,
                                  results_scores)
          }

          cat( " kappa_cohen:" , kappa," f1score: ", mean(f1sc),"\n")
          rm(svmfit)
        }
      }
    }
    write.csv(results_table,paste0("classification_models/grid_search/svm_polynomialktuning",iteration,".csv"))
  }

  write.csv(results_table,paste0("classification_models/grid_search/svm_polynomialktuning",iteration,".csv"))
  #
  # ############GBTrees
  
  # first = T
  eta_params = c(0.1,0.01,0.001)
  max_depth_params = c(2,4,8,16,24)
  gamma_params = c(0.1,0.5,1,2,4,8,16)
  colsample_params = c(0.5,0.7,0.9)
  sub_sample_params = c(0.7,1)
  nrounds_params = c(200,400,800,1200)
  # possible_comb = unlist(lapply(eta_params, function(x)
  #   paste0(x, "_",unlist(lapply(max_depth_params, function(y) paste0(y,"_", gamma_params)))))) 
  # 
  # do_process = TRUE
  # if(read_backup){
  #   file_name = paste0("classification_models/grid_search/xboostingtuning",iteration,".csv")
  #   if(file.exists(file_name)){
  #     results_table = read.csv(file_name, row.names = 1)
  #          process_done = unique(paste0(results_table$eta, "_", as.character(results_table$max_depth),"_", as.character(results_table$gamma)))
  #     
  #     if(length(possible_comb[!(possible_comb%in%process_done)])==0){
  #       do_process = FALSE
  #     }
  #     first = F    
  #   }else{
  #     first = T
  #   }
  # }
  x_training = training_data[,-ncol(training_data)]
  x_validation = validation_data[,-ncol(validation_data)]

  y_training = training_data[,ncol(training_data)]
  y_training = as.numeric(y_training) - 1
  y_validation = as.numeric(validation_data[,ncol(validation_data)])-1


  x_training <- xgboost::xgb.DMatrix(as.matrix(x_training), label = y_training, missing = NA)
  x_validation <- xgboost::xgb.DMatrix(as.matrix(x_validation),  missing = NA)
  # 
  # mc_p = 1
  # eta_p = 1; mx_p =4; col_p =0.7; gam_p = 0; nround_p = 400
  # if(do_process){
  #   Sys.time()->start
  #   for(eta_p in eta_params){
  #     
  #     for(mx_p in max_depth_params){
  #       dosubprocess = TRUE  
  #       
  #       for(gam_p in gamma_params){
  #         if(read_backup){
  #           if(paste0(as.character(eta_p),"_",as.character(mx_p),"_",as.character(gam_p))%in%process_done){
  #             dosubprocess = FALSE
  #           }
  #         }
  #         if(dosubprocess){
  #           for(col_p in colsample_params){
  #             for(sub_p in sub_sample_params){
  #               for(nround_p in nrounds_params){
  #                 cat(eta_p, " " , mx_p, " " ,gam_p, " " ,col_p, " " ,mc_p, " " ,nround_p )
  #                 
  #                 set.seed(1)
  #                 out <- xgboost::xgb.train(list(eta = eta_p,
  #                                                max_depth = mx_p,
  #                                                gamma = gam_p,
  #                                                colsample_bytree = col_p,
  #                                                min_child_weight = 1,
  #                                                subsample = 1,
  #                                                colsample_bytree = 0.5),
  #                                           data = x_training,
  #                                           num_class = length(unique(y_training)),
  #                                           nrounds = nround_p,
  #                                           objective = "multi:softprob")
  #                 
  #                 predictedValues = predict(out,x_validation)
  #                 #xgb.plot.deepness(model = out)
  #                 
  #                 out <- matrix(predictedValues, ncol = 6, byrow = TRUE)
  #                 out <- apply(out, 1, which.max)-1
  #                 
  #                 
  #                 CM=table(y_validation,out)
  #                 kappa = ClassificationMetrics(CM)[[1]]
  #                 accuracy = ClassificationMetrics(CM)[[2]]
  #                 preci = unlist(ClassificationMetrics(CM)[[3]])
  #                 names(preci) = paste0("precision_",1:length(preci))
  #                 
  #                 recll = unlist(ClassificationMetrics(CM)[[4]])
  #                 names(recll) = paste0("recall_",1:length(recll))
  #                 
  #                 f1sc = unlist(ClassificationMetrics(CM)[[5]])
  #                 names(f1sc) = paste0("f1scrore_",1:length(f1sc))
  #                 
  #                 
  #                 if(first){
  #                   results_table = data.frame(iteration = iteration,eta = eta_p ,max_depth =  mx_p, gamma = gam_p, colsample_bytree = col_p, 
  #                                              min_child_weight = mc_p, nrounds =nround_p, kappa, accuracy, t(preci),t(recll),t(f1sc), f1micro = unlist(ClassificationMetrics(CM)[[6]]) )
  #                   first = F
  #                 }else{
  #                   results_table = rbind(results_table,
  #                                         data.frame(iteration = iteration,eta = eta_p ,max_depth =  mx_p, gamma = gam_p, colsample_bytree = col_p, 
  #                                                    min_child_weight = mc_p, nrounds =nround_p, kappa, accuracy, t(preci),t(recll),t(f1sc) , f1micro = unlist(ClassificationMetrics(CM)[[6]])))
  #                 }
  #                 
  #                 cat( " kappa_cohen:" , kappa," f1score: ", mean(f1sc),"\n")
  #                 rm(out)
  #               }
  #               
  #             }
  #           }
  #           cat("\n")
  #           print(Sys.time()-start)
  #         }
  #         write.csv(results_table,paste0("classification_models/grid_search/xboostingtuning", iteration, ".csv"))
  #       }
  #     }  
  #     
  #   }
  #      
  #   
  #   write.csv(results_table,paste0("classification_models/grid_search/xboostingtuning", iteration, ".csv"))
  #   
  # }
  # 
  
})

stopCluster(cl)






################ metrics

model_toevaluate = "random_forest"
model_toevaluate = "svm"
model_toevaluate = "xgboost"

library(dplyr)

data_classes = c("vegetative","reproductive",
                 "ripening","harvested","soil","other")
length(data_classes)

##
models_toevaluate = c("random_forest","svm_radial","svm_polynomial")

model_metrics=lapply(models_toevaluate, function(model_toevaluate){
  
  iterations = 1:5
  if(model_toevaluate == "random_forest"){
    
    suffix = "rf"
    random_tunnings = organice_metrics(suffix,data_classes,files_path = "classification_models/grid_search/")
  }
  if(model_toevaluate == "svm_radial"){
    suffix = c("svm_radial")
    random_tunnings = organice_metrics(suffix,data_classes,files_path = "classification_models/grid_search/")
    
  }
  if(model_toevaluate == "svm_polynomial"){
    suffix ="svm_polynomialk"
    random_tunnings = organice_metrics(suffix,data_classes,files_path = "classification_models/grid_search/")
  }
  if(model_toevaluate == "xgboost"){
    suffix ="xboosting"
    random_tunnings = organice_metrics(suffix,data_classes,files_path = "classification_models/grid_search/")
  }
  
  random_tunnings$scores$model = model_toevaluate
  
  
  return(random_tunnings)
})


## graphic

data_results = do.call(rbind,lapply(model_metrics, function(x){
  x$scores[,-1]
}))

data_plot  =reshape2::melt(data_results)

levels(data_plot$variable)= data_classes
library(ggplot2)
ggplot(data_plot, aes(variable, value, colour = model))+ geom_boxplot()+theme_bw() +
  labs( y = "f1 Score", x ="Classes", colour = 'models')+ylim(0.75,1)


### export best models
### prepare data

library(stringr)
alldata_sub = data_training

alldata_sub = alldata[alldata$kfold%in%4,]
len_id = str_sub(row.names(alldata_sub ),1,-1)

num_pos = regexpr(pattern = "_k",row.names(alldata))

id= str_sub(row.names(alldata),1, num_pos )

alldata$ID = id

alldata["dataset"] = NULL
alldata["kfold"] = NULL

alldata = data.frame(alldata %>%
             group_by(ID) %>% filter(row_number(ID) == 1))

row.names(alldata) = alldata$ID
alldata["ID"]=NULL

y_variable = factor(as.numeric(alldata$class) - 1)
x_variables = alldata[,-ncol(alldata)]


## random forest
mtry = model_metrics[[1]]$best_tune[1,2]
ntree = model_metrics[[1]]$best_tune[1,3]
set.seed(1)
model = randomForest::randomForest(x_variables, 
                                   y_variable,
                                   mtry = mtry, ntree = ntree)

save(model, file = "D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/classification_models/rf_phen_identification_veg.RData")

## svm polynomial

minmaxvalues = get_featuresmin_maxvals(x_variables)
dates_pos = grepl("Date", names(minmaxvalues))


min_dates = min(unlist(minmaxvalues[dates_pos]))
max_dates = max(unlist(minmaxvalues[dates_pos]))
# 
# ## calculate min max derivatives
# 
min_der = min(unlist(minmaxvalues[!dates_pos]))
max_der = max(unlist(minmaxvalues[!dates_pos]))
# 
# ## feature scaling

## scale data

dates_scaled = do.call(cbind,lapply(which(dates_pos),function(x){
  feature_scaling(x_variables[,x], c(min_dates, max_dates))
}))

derivatives_scaled = do.call(cbind,lapply(which(!dates_pos),function(x){
  feature_scaling(x_variables[,x], c(min_der, max_der))
}))

## data scaled
x_variables_scaled = data.frame(cbind(dates_scaled,
                                      derivatives_scaled))
names(x_variables_scaled) = names(x_variables)

x_variables_scaled$y = y_variable

degreep = model_metrics[[3]]$best_tune[1,"degree"]
gamma_p = model_metrics[[3]]$best_tune[1,"gamma"]
cost_p = model_metrics[[3]]$best_tune[1,"cost"]
coef0p = model_metrics[[3]]$best_tune[1,"coef0"]

library(e1071)
set.seed(1)
svmfit = svm(y ~ ., data = x_variables_scaled, 
             kernel = "polynomial",
             degree = degreep,gamma = gamma_p,
             cost = cost_p,
             coef0=coef0p, 
             scale = F)

model = list(model = svmfit,
min_maxvalues = minmaxvalues)

save(model, file = "D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/classification_models/svm_polynomial_phen_identification_veg.RData")


## Gradient Boosting machine

varxgboost <- xgboost::xgb.DMatrix(as.matrix(x_variables), 
                                     label = as.numeric(y_variable)-1, missing = NA) 

eta_p = model_metrics[[4]]$best_tune[1,"eta"]
mx_p = model_metrics[[4]]$best_tune[1,"max_depth"]
gam_p = model_metrics[[4]]$best_tune[1,"gamma"]
col_p = model_metrics[[4]]$best_tune[1,"colsample_bytree"]
nround_p = model_metrics[[4]]$best_tune[1,"nrounds"]
mc_p = model_metrics[[4]]$best_tune[1,"min_child_weight"]

set.seed(1)
model <- xgboost::xgb.train(list(eta = eta_p,
                               max_depth = mx_p,
                               gamma = gam_p,
                               colsample_bytree = col_p,
                               min_child_weight = mc_p,
                               subsample = 1),
                          data = varxgboost,
                          num_class = length(unique(y_variable)),
                          nrounds = nround_p,
                          objective = "multi:softprob")

save(model, file = "D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/classification_models/xgboost_phen_identification_veg.RData")


### summarize iterations









set.seed(1)
class_model <- xgboost::xgb.train(list(eta = 0.01,
                               max_depth = 8,
                               gamma = 1,
                               colsample_bytree = 0.7,
                               min_child_weight = 1,
                               subsample = 1),
                          data = x_training,
                          num_class = length(unique(y_training)),
                          nrounds = 100,
                          objective = "multi:softprob")

predictedValues = predict(class_model,x_validation)

out <- matrix(predictedValues, ncol = 7, byrow = TRUE)
out <- apply(out, 1, which.max)-1


CM=table(y_validation,out)
ClassificationMetrics(CM)
save(class_model,file="D:/temp/xgbTree_conf_testcontrolparams.RData")

class_model

library(ggplot2)
library(dplyr)
library(reshape2)
results_table = read.csv("D:/temp/xboostingtuning.csv")

results_table = results_table[results_table$kappa>0.92,]

data_toPlot = reshape2::melt( results_table , "X")

data_f1 = data_toPlot[grepl(data_toPlot$variable, pattern = "recal"),]
data_f2 = data_toPlot[grepl(data_toPlot$variable, pattern = "precision"),]
data_f3 = data_toPlot[grepl(data_toPlot$variable, pattern = "f1scrore"),]
data_f = rbind(data_f1, data_f2, data_f3)

data_f$Grwoth_stage = ""
data_f$Grwoth_stage[grepl(data_f$variable, pattern = "_1")] = "Harvested"
data_f$Grwoth_stage[grepl(data_f$variable, pattern = "_2")] = "other"
data_f$Grwoth_stage[grepl(data_f$variable, pattern = "_3")] = "Reproductive"
data_f$Grwoth_stage[grepl(data_f$variable, pattern = "_4")] = "Ripening"
data_f$Grwoth_stage[grepl(data_f$variable, pattern = "_5")] = "Soil"
data_f$Grwoth_stage[grepl(data_f$variable, pattern = "_6")] = "Late_Vegetative"
data_f$Grwoth_stage[grepl(data_f$variable, pattern = "_7")] = "Early_Vegetative"

data_f$Metric = ""
data_f$Metric[grepl(data_f$variable, pattern = "recal")] = "recall"
data_f$Metric[grepl(data_f$variable, pattern = "precision")] = "precision"
data_f$Metric[grepl(data_f$variable, pattern = "f1scrore")] = "f1scrore"

results_table[results_table$kappa>0.9,]
ggplot(data_f, aes(Grwoth_stage, value))+geom_boxplot()+
  facet_grid(rows = vars(Metric))


##################RAndom Forest
load(file="D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/sup_class_models/phen_identification/rf_conf_8optical_ndvi_deri.RData")


lapply(ClassificationModels, function(x)x[[1]])


mtry = caret::var_seq(p = ncol(x_training), 
                      classification = is.factor(y_training), 
                      len = length(is.factor(y_training)))

round((ncol(x_training)-1)/3)

mtry_params = 3:(ncol(data_training[,-c(10:14)])-1)
first = T
iteration = 2
for(iteration in 1:3){
  
  data_training_iter = data_training[iter %in% paste0("i_",iteration),-c(11:14)]
  data_training_iter$Class = as.character(data_training_iter$Class)
  data_training_iter$Class[data_training_iter$Class%in%c("Urban_Zones",  "Roads" , "Water")] = 'other'
  data_training_iter$Class = as.factor(data_training_iter$Class)
  
  posPattern=regexpr("*_t_", row.names(data_training_iter))
  val = (str_sub(row.names(data_training_iter),posPattern+2))
  
  valpos = which(val %in% paste0("_training"))
  
  x_training = data_training_iter[valpos,-ncol(data_training_iter)]
  x_validation = data_training_iter[-valpos,-ncol(data_training_iter)]
  y_training = data_training_iter[valpos,ncol(data_training_iter)]
  y_training <- as.numeric(y_training) - 1
  y_validation = as.numeric(data_training_iter[-valpos,ncol(data_training_iter)])-1
  
  for(mtry_p in mtry_params){
    for(ntree_p in c(200,500,800,1000,1200,1600)){
      cat(mtry_p, " " , ntree_p, " ")
      
      set.seed(1)
      modelRF = randomForest::randomForest(x_training, as.factor(y_training), mtry = mtry_p, ntree = ntree_p)
      
      predictedValues = unname(predict(modelRF,x_validation))
      
      CM=table(y_validation,predictedValues)
      kappa = ClassificationMetrics(CM)[[1]]
      accuracy = ClassificationMetrics(CM)[[2]]
      preci = unlist(ClassificationMetrics(CM)[[3]])
      names(preci) = paste0("precision_",1:length(preci))
      
      recll = unlist(ClassificationMetrics(CM)[[4]])
      names(recll) = paste0("recall_",1:length(recll))
      
      f1sc = unlist(ClassificationMetrics(CM)[[5]])
      
      names(f1sc) = paste0("f1scrore_",1:length(f1sc))
      
      row_results = data.frame(iteration = iteration,mtry = mtry_p ,ntree =  ntree_p, kappa, accuracy, t(preci),t(recll),t(f1sc) )

      
      if(first){
        results_table = row_results
        first = F
      }else{
        results_table = rbind(results_table,
                              row_results)
      }
      
      cat( " kappa_cohen:" , kappa," f1score: ", mean(f1sc),"\n")
      rm(modelRF)
      
    }
  }
  
  write.csv(results_table,"D:/temp/rftingtuning.csv")
}





## save best tuning

set.seed(1)
modelRF = randomForest::randomForest(x_training, as.factor(y_training), mtry = 4, ntree = 500)

predictedValues = unname(predict(modelRF,x_validation))

CM=table(y_validation,predictedValues)

save(modelRF,file="D:/temp/rf_conf_testcontrolparams.RData")




kappa = ClassificationMetrics(CM)[[1]]
accuracy = ClassificationMetrics(CM)[[2]]
preci = unlist(ClassificationMetrics(CM)[[3]])
names(preci) = paste0("precision_",1:length(preci))

recll = unlist(ClassificationMetrics(CM)[[4]])
names(recll) = paste0("recall_",1:length(recll))

f1sc = unlist(ClassificationMetrics(CM)[[5]])
names(f1sc) = paste0("f1scrore_",1:length(f1sc))



data_metrics = read.csv("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/graphics/validationmetrics.csv")


ggplot(data_metrics, aes(type, value, color = model))+ geom_point(size = 3 )+theme_bw() +
  labs( y = "f1 Score", x ="tipo de clasificación", colour = 'modelos')


data_testmetrics = read.csv("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/graphics/testmetrics.csv")

data_testmetrics$type = factor(data_testmetrics$type, levels = c("soil", "early vegetative",
                                                                 "late vegetative","reproductive",
                                                                 "ripening","harvested"))



grid_default <- expand.grid(
  nrounds = 100,
  max_depth = 6,
  eta = 0.3,
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1
)

train_control <- caret::trainControl(
  method = "none",
  verboseIter = FALSE, # no training log
  allowParallel = TRUE # FALSE for reproducible results 
)

xgb_base <- caret::train(
  x = input_x,
  y = input_y,
  trControl = train_control,
  tuneGrid = grid_default,
  method = "xgbTree",
  verbose = TRUE
)


ggplot(data_testmetrics, aes(type, value, color = model))+ geom_point(size = 3 )+theme_bw() +
  labs( y = "f1 Score", x ="tipo de clasificación", colour = 'modelos')

######################
#####
##### Support vector machines radial kernel



# set training and validation samples
load(file="D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/sup_class_models/phen_identification/svmRadial_conf_8optical_ndvi_deri.RData")


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

feature_scaling =  function(feature, limits = c(0,1), type = "min_max"){

  ## calculate min max feature scaling
  (feature  - limits[1])/(limits[2]-limits[1])
  
}



minmaxvalues = get_featuresmin_maxvals(data_training[,1:9])
dates_pos = grepl("Date", names(minmaxvalues))

### 

min_dates = min(unlist(minmaxvalues[dates_pos]))
max_dates = max(unlist(minmaxvalues[dates_pos]))

## calculate min max derivatives

min_der = min(unlist(minmaxvalues[!dates_pos]))
max_der = max(unlist(minmaxvalues[!dates_pos]))

## feature scaling

first = T
iteration = 2
for(iteration in 1:3){
  
  data_training_iter = data_training[iter %in% paste0("i_",iteration),-c(11:14)]
  data_training_iter$Class = as.character(data_training_iter$Class)
  data_training_iter$Class[data_training_iter$Class%in%c("Urban_Zones",  "Roads" , "Water")] = 'other'
  data_training_iter$Class = as.factor(data_training_iter$Class)
  
  posPattern=regexpr("*_t_", row.names(data_training_iter))
  val = (str_sub(row.names(data_training_iter),posPattern+2))
  
  valpos = which(val %in% paste0("_training"))
  
  x_training = data_training_iter[valpos,-ncol(data_training_iter)]
  x_validation = data_training_iter[-valpos,-ncol(data_training_iter)]
  y_training = data_training_iter[valpos,ncol(data_training_iter)]
  y_training <- as.numeric(y_training) - 1
  y_validation = as.numeric(data_training_iter[-valpos,ncol(data_training_iter)])-1
  

  dates_scaled = do.call(cbind,lapply(which(dates_pos),function(x){
    feature_scaling(data_training_iter[,x], c(min_dates, max_dates))
  }))
  
  derivatives_scaled = do.call(cbind,lapply(which(!dates_pos),function(x){
    feature_scaling(data_training_iter[,x], c(min_der, max_der))
  }))

  ##
  training_iter_scaled = data.frame(dates_scaled, derivatives_scaled, Class = data_training_iter[, ncol(data_training_iter)])
  ## reassign names
  
  names(training_iter_scaled)[-ncol(training_iter_scaled)] = names(data_training_iter)[-ncol(data_training_iter)]
  
  x_training = training_iter_scaled[valpos,-ncol(training_iter_scaled)]
  
  table(training_iter_scaled$Class)
  
  x_validation = training_iter_scaled[-valpos,-ncol(training_iter_scaled)]
  y_training = training_iter_scaled[valpos,ncol(training_iter_scaled)]
  y_training <- as.numeric(y_training) - 1
  y_validation = as.numeric(training_iter_scaled[-valpos,ncol(training_iter_scaled)])-1
  
  data_trainingsvm = data.frame(x_training,y = as.factor(y_training))
  
  gammagrid  = c(0.05,0.1,1,2,4,8)
  costgrid = c(32,64,128,216,512,1048,2056,4112,8224)
  degrees = c(2:5)
  first = T
  library(e1071)
  for(degreep in degrees){
    cat(degreep, '***\n')
    for(coef0p in c(0.05,0,1,2,4,8,16)){
      for(cost_p in costgrid){
        cat(cost_p, '\n')
        for(gamma_p in gammagrid){
          
          set.seed(1)
          #svmfit = svm(y ~ ., data = data_trainingsvm, kernel = "radial basis",,gamma = gamma_p,cost = cost_p, scale = F)
          svmfit = svm(y ~ ., data = data_trainingsvm, kernel = "polynomial",degree = degreep,gamma = gamma_p,cost = cost_p,coef0=coef0p, scale = F)
          
          predictedValues = unname(predict(svmfit,x_validation))
          
          CM=table(y_validation,predictedValues)
          kappa = ClassificationMetrics(CM)[[1]]
          accuracy = ClassificationMetrics(CM)[[2]]
          preci = unlist(ClassificationMetrics(CM)[[3]])
          names(preci) = paste0("precision_",1:length(preci))
          
          recll = unlist(ClassificationMetrics(CM)[[4]])
          names(recll) = paste0("recall_",1:length(recll))
          
          f1sc = unlist(ClassificationMetrics(CM)[[5]])
          names(f1sc) = paste0("f1scrore_",1:length(f1sc))
          
          results_scores = data.frame(iteration = iteration,gamma = gamma_p ,cost =  cost_p, degree = degreep, coef0 = coef0p,kappa, accuracy, t(preci),t(recll),t(f1sc) )
          if(first){
            results_table = results_scores
            first = F
          }else{
            results_table = rbind(results_table,
                                  results_scores)
          }
          
          cat( " kappa_cohen:" , kappa," f1score: ", mean(f1sc),"\n")
          rm(svmfit)
        }
      }
    }
    
  }
}
write.csv(results_table,"D:/temp/svm_polynomialktuning.csv")

data_plot = reshape2::melt(results_table[,-c(1:4)])

newcolumns = do.call(rbind,(stringr::str_split(data_plot$variable,"_")))

data_plot[["metric"]] = newcolumns[,1]
data_plot[["class"]] = newcolumns[,2]              

### 

ggplot(data_plot, aes(class,value,color = metric)) + geom_boxplot()

###

write.csv(results_table,"D:/temp/svm_radialtuning.csv")

rf_tunning = read.csv("D:/temp/rftingtuning.csv")
