rm(list = ls())
library(xgboost)
library(caret)

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


data_training_optico=read.csv(paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/model_inputs/phen_identification/optical_data_conf8_ndvi_lswi_derndvi.csv"))
data_training = data_training_optico
row.names(data_training) = data_training$X
#row.names(data_training) = data_training[,1]
data_training = data_training[,-1]
levelsNA= unique(unlist(apply(data_training ,  2 , function(x) unique(which(is.na(x))))))
if(length(levelsNA) > 0){
  data_training = data_training[-levelsNA,]  
}
(1:ncol(data_training))[-which(names(data_training) == "Class")]
data_training = data_training[,c((1:ncol(data_training))[-which(names(data_training) == "Class")],
                                 which(names(data_training) == "Class"))]
head(data_training)
library(stringr)

posPattern=regexpr("*_i_", row.names(data_training))
iter = (str_sub(row.names(data_training),posPattern+1, posPattern+3))

data_training_iter = data_training[iter %in% paste0("i_",3),-c(8:14)]

data_training_iter$Class = as.character(data_training_iter$Class)
data_training_iter$Class[data_training_iter$Class%in%c("Urban_Zones",  "Roads" , "Water")] = 'other'
data_training_iter$Class = as.factor(data_training_iter$Class)

posPattern=regexpr("*_t_", row.names(data_training_iter))
val = (str_sub(row.names(data_training_iter),posPattern+2))

valpos = which(val %in% paste0("_training"))



head(x_training)

x_training = data_training_iter[valpos,-ncol(data_training_iter)]
x_validation = data_training_iter[-valpos,-ncol(data_training_iter)]
y_training = data_training_iter[valpos,ncol(data_training_iter)]
y_training <- as.numeric(y_training) - 1
y_validation = as.numeric(data_training_iter[-valpos,ncol(data_training_iter)])-1
x_training <- xgboost::xgb.DMatrix(as.matrix(x_training), label = y_training, missing = NA) 
x_validation <- xgboost::xgb.DMatrix(as.matrix(x_validation),  missing = NA) 

eta_params = c(1,0.1,0.01,0.001)
max_depth_params = c(2,4,8,16,32)
gamma_params = c(0,1)
colsample_params = c(0.7,1)
min_chil_params = c(1,2)
nrounds_params = c(100,400,800,1200,1600,1800)

first = T

Sys.time()->start
for(eta_p in eta_params){
  for(mx_p in max_depth_params){
    for(gam_p in gamma_params){
      for(col_p in colsample_params){
        for(mc_p in min_chil_params){
          for(nround_p in nrounds_params){
            cat(eta_p, " " , mx_p, " " ,gam_p, " " ,col_p, " " ,mc_p, " " ,nround_p )
            
            set.seed(1)
            out <- xgboost::xgb.train(list(eta = eta_p,
                                           max_depth = mx_p,
                                           gamma = gam_p,
                                           colsample_bytree = col_p,
                                           min_child_weight = mc_p,
                                           subsample = 1),
                                      data = x_training,
                                      num_class = length(unique(y_training)),
                                      nrounds = nround_p,
                                      objective = "multi:softprob")
            
            predictedValues = predict(out,x_validation)
            
            out <- matrix(predictedValues, ncol = 7, byrow = TRUE)
            out <- apply(out, 1, which.max)-1
            
            
            CM=table(y_validation,out)
            kappa = ClassificationMetrics(CM)[[1]]
            accuracy = ClassificationMetrics(CM)[[2]]
            preci = unlist(ClassificationMetrics(CM)[[3]])
            names(preci) = paste0("precision_",1:length(preci))
            
            recll = unlist(ClassificationMetrics(CM)[[4]])
            names(recll) = paste0("recall_",1:length(recll))
            
            f1sc = unlist(ClassificationMetrics(CM)[[5]])
            names(f1sc) = paste0("f1scrore_",1:length(f1sc))
            if(first){
              results_table = data.frame(eta = eta_p ,max_depth =  mx_p, gamma = gam_p, colsample_bytree = col_p, 
                         min_child_weight = mc_p, nrounds =nround_p, kappa, accuracy, t(preci),t(recll),t(f1sc) )
              first = F
            }else{
              results_table = rbind(results_table,
                                    data.frame(eta = eta_p ,max_depth =  mx_p, gamma = gam_p, colsample_bytree = col_p, 
                                                 min_child_weight = mc_p, nrounds =nround_p, kappa, accuracy, t(preci),t(recll),t(f1sc) ))
            }
            
            cat( " kappa_cohen:" , kappa," f1score: ", mean(f1sc),"\n")
            rm(out)
          }
          
        }
      }
      cat("\n")
      print(Sys.time()-start)
    }
  }
}

write.csv(results_table,"D:/temp/xboostingtuning.csv")




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

load(file="D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/sup_class_models/phen_identification/rf_conf_8optical_ndvi_deri.RData")


lapply(ClassificationModels, function(x)x[[1]])

x_training = data_training_iter[valpos,-ncol(data_training_iter)]
x_validation = data_training_iter[-valpos,-ncol(data_training_iter)]
y_training = data_training_iter[valpos,ncol(data_training_iter)]
y_training <- as.numeric(y_training) - 1
y_validation = as.numeric(data_training_iter[-valpos,ncol(data_training_iter)])-1

mtry = caret::var_seq(p = ncol(x_training), 
                      classification = is.factor(y_training), 
                      len = length(is.factor(y_training)))

round((ncol(x_training)-1)/3)

mtry_params = 3:ncol(x_training)
first = T

for(mtry_p in mtry_params){
  for(ntree_p in c(100,200,500,800,1000,1200,1600)){
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
    if(first){
      results_table = data.frame(mtry = mtry_p ,ntree =  ntree_p, kappa, accuracy, t(preci),t(recll),t(f1sc) )
      first = F
    }else{
      results_table = rbind(results_table,
                            data.frame(mtry = mtry_p ,ntree =  ntree_p, kappa, accuracy, t(preci),t(recll),t(f1sc) ))
    }
    
    cat( " kappa_cohen:" , kappa," f1score: ", mean(f1sc),"\n")
    rm(modelRF)
    
  }
}


write.csv(results_table,"D:/temp/rftingtuning.csv")

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

ggplot(data_testmetrics, aes(type, value, color = model))+ geom_point(size = 3 )+theme_bw() +
  labs( y = "f1 Score", x ="tipo de clasificación", colour = 'modelos')

######################
#####
##### Support vector machines radial kernel

library(e1071)

# set training and validation samples


x_training = data_training_iter[valpos,-ncol(data_training_iter)]
x_validation = data_training_iter[-valpos,-ncol(data_training_iter)]
y_training = data_training_iter[valpos,ncol(data_training_iter)]
y_training <- as.numeric(y_training) - 1
y_validation = as.numeric(data_training_iter[-valpos,ncol(data_training_iter)])-1

data_trainingsvm = data.frame(x_training,y = as.factor(y_training))

gammagrid  = c(0.0001,0.001,0.01,0.1,1)
costgrid = c(1,10,100,1000,10000)
first = T
for(cost_p in costgrid){
  for(gamma_p in gammagrid){
    set.seed(1)
    svmfit = svm(y ~ ., data = data_trainingsvm, kernel = "radial", gamma = gamma_p,cost = cost_p, scale = TRUE)
    
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
    if(first){
      results_table = data.frame(gamma = gamma_p ,cost =  cost_p, kappa, accuracy, t(preci),t(recll),t(f1sc) )
      first = F
    }else{
      results_table = rbind(results_table,
                            data.frame(gamma = gamma_p ,cost =  cost_p, kappa, accuracy, t(preci),t(recll),t(f1sc) ))
    }

    cat( " kappa_cohen:" , kappa," f1score: ", mean(f1sc),"\n")
    rm(svmfit)
  }
}




                  