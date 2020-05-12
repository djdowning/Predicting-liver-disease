##################################################################
#INSTALL PKGS
##################################################################

if(!require(tidyverse)) install.packages("tidyverse") 
if(!require(caret)) install.packages("caret") 
if(!require(ROCR)) install.packages("ROCR")
if(!require(PRROC)) install.packages("PRROC")
if(!require(ROSE)) install.packages("ROSE")
if(!require(xgboost)) install.packages("xgboost")
if(!require(kableExtra)) install.packages("kableExtra")
if(!require(gam)) install.packages("gam")
if(!require(splines)) install.packages("splines")
if(!require(foreach)) install.packages("foreach")
if(!require(ggthemes)) install.packages("ggthemes")
if(!require(DMwR)) install.packages("DMwR")
if(!require(caretEnsemble)) install.packages("caretEnsemble")
if(!require(reshape2)) install.packages("reshape2")

##################################################################
#LOAD LIBRARIES
##################################################################
library(tidyverse)
library(caret)
library(ROCR)
library(PRROC)
library(ROSE)
library(xgboost)
library(kableExtra)
library(gam)
library(splines)
library(foreach)
library(ggthemes)
library(DMwR)
library(caretEnsemble)
library(reshape2)

##################################################################
#LOAD DATA
##################################################################
liver <- read_csv("data/indian_liver_patient.csv")
save(liver, file = "rda/liver.rda")

##################################################################
#INSPECT AND WRANGLE DATA
#####################################################################################

#View Columns
tibble(Columns=names(liver)) %>%   
  kable() %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                position = "center",
                font_size = 10,
                full_width = FALSE)

#Inspect Outcomes
liver %>% group_by(Dataset) %>% summarize(Count=n()) %>%
  kable() %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                position = "center",
                font_size = 10,
                full_width = FALSE)

#Plot Outcomes by Gender by Age
liver %>% group_by(Gender, Age_By_Decade=cut(Age,c(0,10,20,30,40,50,60,70,80,90,100), 
                  labels = seq(10,100,10)), 
                   Dataset=ifelse(Dataset=="1","Disease","Healthy")) %>% 
  ggplot(aes(Age_By_Decade, fill=Dataset)) +
  geom_bar(position=position_dodge()) +
  facet_wrap(~ Gender) +
  ylab("Frequency")


#Check for Duplicates
liver %>% group_by(duplicated(liver)) %>% summarize(Count = n())  %>%
  kable() %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                position = "center",
                font_size = 10,
                full_width = FALSE)

liver <- unique(liver)


#Check for missing values
sapply(liver, function(x) sum(is.na(x))) %>% 
  kable(col.names = c("Missing Values")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                position = "center",
                font_size = 10,
                full_width = FALSE)

#Remove 4 missing values
liver <- na.omit(liver)

#Look at distribution of outcomes
liver %>% group_by(Dataset) %>% summarize(n=n())
sum(liver$Dataset==1) / sum(liver$Dataset==2)

#Go ahead and change the Variable we are predicting to a logical factor since 2 way classification
#and also change Gender to logical factor
liver$Dataset <- factor(ifelse(liver$Dataset==2,0,1), levels = c(1,0))
levels(liver$Dataset)
liver$Gender <- factor(ifelse(liver$Gender=="Male",1,0), levels = c(1,0))
levels(liver$Gender)

#Look at distributions of each predictor to see if any are skewed distributions that would benefit
#from a transformation
d <- melt(liver)

ggplot(d,aes(x = value)) + 
  geom_histogram(fill="blue")  +
  ylab ("Frequency") +
  facet_wrap(~variable,scales = "free_x") 


#There are 5 very skewed distributions so lets log transform them

livertrans <- liver

livertrans$Total_Bilirubin <- log(livertrans$Total_Bilirubin +1)
livertrans$Direct_Bilirubin <- log(livertrans$Direct_Bilirubin +1)
livertrans$Alkaline_Phosphotase <- log(livertrans$Alkaline_Phosphotase +1)
livertrans$Alamine_Aminotransferase <- log(livertrans$Alamine_Aminotransferase +1)
livertrans$Aspartate_Aminotransferase <- log(livertrans$Aspartate_Aminotransferase +1)

d1 <- melt(livertrans)

#View New Histograms
ggplot(d1,aes(x = value)) + 
  geom_histogram(fill="blue")  +
  ylab ("Frequency") +
  facet_wrap(~variable,scales = "free_x") 

# Correlation Matrix to see if any very highly correlated predictors
head(livertrans)
lt <- livertrans[, (3:10)]
lt <- sapply(lt, as.numeric)
lower <- function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

upper <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

reorder_cormat <- function(cormat){
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}


corr_matrix <- round(cor(lt),2)
corr_matrix <- reorder_cormat(corr_matrix)
upper_tri <- upper(corr_matrix)
melt_corr_matrix <- melt(upper_tri, na.rm = TRUE)
ggplot(melt_corr_matrix,aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 9, hjust = 1), axis.text.y = element_text(size = 9),                    axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank()) +
  coord_fixed()

#Get Values for both bilirubins since so highly ocrrelated
cor(lt[,c(1,2)]) %>% 
  kable() %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                position = "center",
                font_size = 10,
                full_width = FALSE)

#Looks at rowMeans for both bilirubins to see which is more highly correlated with other predictors
cor <- cor(lt) 
rowMeans(cor[1:2,]) %>%
  kable() %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                position = "center",
                font_size = 10,
                full_width = FALSE)


#Finalize Data set for partition
livertrans <- select (livertrans,-c(Direct_Bilirubin))



#####################################################################################
#PARTITION DATA FOR TEST AND TRAINING
#####################################################################################

set.seed(1147, sample.kind = "Rounding") #in R 3.6 or later 
test_index <- createDataPartition(y = livertrans$Dataset, times = 1,
                                  p = 0.2, list = FALSE)
train <- livertrans[-test_index,]
test <- livertrans[test_index,]

#####################################################################################
#TRAIN A VARIETY OF MODELS TO LOOK AT OVERALL ACCURACY
#####################################################################################


#####################################################################################
#JUST ASSIGN TRUE TO ALL
#####################################################################################
#Just Guess That Everyone Has the Disease Since Dataset Skewed Towards Disease=1
set.seed(1147, sample.kind = "Rounding") #in R 3.6 or later
predictions <- factor(rep(1,times=length(test$Dataset)),levels=c(1,0))
predictions
predictions1 <- prediction(as.numeric(predictions) , test$Dataset)
cm <- confusionMatrix(predictions, test$Dataset)
auc <- performance(predictions1, "auc")
AUPRC <- performance(predictions1, "prec", "rec")
aucv <- auc@y.values[[1]]
#AUPRCv <- pr.curve(scores.class0 = predictions[test$Dataset == 1],
#                   scores.class1 = predictions[test$Dataset == 0],
#                   curve=TRUE)

results <- tibble(Method="Guess Everyone",
                  Accuracy=cm$overall["Accuracy"],
                  Sensitivity= cm$byClass["Sensitivity"], 
                  Specificity= cm$byClass["Specificity"],
                  F1= cm$byClass["F1"],
                  AUC=aucv, 
                  AUPRC=0)

results %>% kable() %>%
       kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                position = "center",
                font_size = 10,
                full_width = FALSE)



#####################################################################################
#RANDOM FOREST
#####################################################################################
set.seed(1147, sample.kind = "Rounding") #in R 3.6 or later
control <- trainControl(method="repeatedcv", number = 10, repeats=10, allowParallel = TRUE)

set.seed(1147, sample.kind = "Rounding") #in R 3.6 or later
model_rf <- train(Dataset ~ ., 
        method = "rf", 
        data = train,
        trControl = control,
        preProcess = c("scale", "center"),
        verboseIter = TRUE)

predictions <- predict(model_rf, test)
predictions1 <- prediction(as.numeric(predictions) , test$Dataset)
predictions_rf <- predict(model_rf, test)

cm <- confusionMatrix(predictions, test$Dataset)
auc <- performance(predictions1, "auc")
AUPRC <- performance(predictions1, "prec", "rec")
aucv <- auc@y.values[[1]]
AUPRCv <- pr.curve(scores.class0 = predictions[test$Dataset == 1],
                   scores.class1 = predictions[test$Dataset == 0],
                   curve=TRUE)
auc_plot <- performance(predictions1, 'sens', 'spec')




results <- bind_rows(results,
                     tibble(Method="Random Forest",
                            Accuracy=cm$overall["Accuracy"],
                            Sensitivity= cm$byClass["Sensitivity"], 
                            Specificity= cm$byClass["Specificity"], 
                            F1= cm$byClass["F1"], 
                            AUC=aucv, 
                            AUPRC=AUPRCv[[3]]))
                     
results %>% 
  kable() %>%
   kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                 position = "center",
                 font_size = 10,
                 full_width = FALSE)


#####################################################################################
#RANDOM FOREST - UP
#####################################################################################
control <- trainControl(method="repeatedcv", 
                        number = 10, 
                        repeats=10, 
                        allowParallel = TRUE,
                        sampling = "up")
set.seed(1147, sample.kind = "Rounding") #in R 3.6 or later 
model_rf_up <- train(Dataset ~ ., 
                  method = "rf", 
                  data = train,
                  trControl = control,
                  preProcess = c("scale", "center"),
                  verboseIter = TRUE
                  )

predictions <- predict(model_rf_up, test)
predictions1 <- prediction(as.numeric(predictions) , test$Dataset)
predictions_up <- predict(model_rf_up, test)

cm <- confusionMatrix(predictions, test$Dataset)
auc <- performance(predictions1, "auc")
AUPRC <- performance(predictions1, "prec", "rec")
aucv <- auc@y.values[[1]]
AUPRCv <- pr.curve(scores.class0 = predictions[test$Dataset == 1],
                   scores.class1 = predictions[test$Dataset == 0],
                   curve=TRUE)


results <- bind_rows(results,
                     tibble(Method="Random Forest - Up",
                            Accuracy=cm$overall["Accuracy"],
                            Sensitivity= cm$byClass["Sensitivity"], 
                            Specificity= cm$byClass["Specificity"],
                            F1= cm$byClass["F1"], 
                            
                            AUC=aucv, 
                            AUPRC=AUPRCv[[3]]))

results %>% 
  kable() %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                position = "center",
                font_size = 10,
                full_width = FALSE)

#####################################################################################
#RANDOM FOREST - DOWN
#####################################################################################
control <- trainControl(method="repeatedcv", 
                        number = 10, 
                        repeats=10, 
                        allowParallel = TRUE,
                        sampling = "down")
set.seed(1147, sample.kind = "Rounding") #in R 3.6 or later  
model_rf_down <- train(Dataset ~ ., 
                     method = "rf", 
                     data = train,
                     trControl = control,
                     preProcess = c("scale", "center"),
                     verboseIter = TRUE
)

predictions <- predict(model_rf_down, test)
predictions1 <- prediction(as.numeric(predictions) , test$Dataset)
predictions_down <- predict(model_rf_down, test)

cm <- confusionMatrix(predictions, test$Dataset)
auc <- performance(predictions1, "auc")
AUPRC <- performance(predictions1, "prec", "rec")
aucv <- auc@y.values[[1]]
AUPRCv <- pr.curve(scores.class0 = predictions[test$Dataset == 1],
                   scores.class1 = predictions[test$Dataset == 0],
                   curve=TRUE)


results <- bind_rows(results,
                     tibble(Method="Random Forest - Down",
                            Accuracy=cm$overall["Accuracy"],
                            Sensitivity= cm$byClass["Sensitivity"], 
                            Specificity= cm$byClass["Specificity"], 
                            F1= cm$byClass["F1"], 
                            AUC=aucv, 
                            AUPRC=AUPRCv[[3]]))

results %>% 
  kable() %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                position = "center",
                font_size = 10,
                full_width = FALSE)
#####################################################################################
#RANDOM FOREST - ROSE
#####################################################################################
control <- trainControl(method="repeatedcv", 
                        number = 10, 
                        repeats=10, 
                        allowParallel = TRUE,
                        sampling = "rose")
set.seed(1147, sample.kind = "Rounding") #in R 3.6 or later  
model_rf_rose <- train(Dataset ~ ., 
                       method = "rf", 
                       data = train,
                       trControl = control,
                       preProcess = c("scale", "center"),
                       verboseIter = TRUE
                       )

predictions <- predict(model_rf_rose, test)
predictions1 <- prediction(as.numeric(predictions) , test$Dataset)
predictions_rose <- predict(model_rf_rose, test)

cm <- confusionMatrix(predictions, test$Dataset)
auc <- performance(predictions1, "auc")
AUPRC <- performance(predictions1, "prec", "rec")
aucv <- auc@y.values[[1]]
AUPRCv <- pr.curve(scores.class0 = predictions[test$Dataset == 1],
                   scores.class1 = predictions[test$Dataset == 0],
                   curve=TRUE)


results <- bind_rows(results,
                     tibble(Method="Random Forest - ROSE",
                            Accuracy=cm$overall["Accuracy"],
                            Sensitivity= cm$byClass["Sensitivity"], 
                            Specificity= cm$byClass["Specificity"], 
                            F1= cm$byClass["F1"], 
                            AUC=aucv, 
                            AUPRC=AUPRCv[[3]]))

results %>% 
  kable() %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                position = "center",
                font_size = 10,
                full_width = FALSE)

#####################################################################################
#RANDOM FOREST - SMOTE
#####################################################################################
control <- trainControl(method="repeatedcv", 
                        number = 10, 
                        repeats=10, 
                        allowParallel = TRUE,
                        sampling = "smote")
set.seed(1147, sample.kind = "Rounding") #in R 3.6 or later  
model_rf_smote <- train(Dataset ~ ., 
                       method = "rf", 
                       data = train,
                       trControl = control,
                       preProcess = c("scale", "center"),
                       verboseIter = TRUE
                       )

predictions <- predict(model_rf_smote, test)
predictions1 <- prediction(as.numeric(predictions) , test$Dataset)
predictions_smote <- predict(model_rf_smote, test)

cm <- confusionMatrix(predictions, test$Dataset)
auc <- performance(predictions1, "auc")
AUPRC <- performance(predictions1, "prec", "rec")
aucv <- auc@y.values[[1]]
AUPRCv <- pr.curve(scores.class0 = predictions[test$Dataset == 1],
                   scores.class1 = predictions[test$Dataset == 0],
                   curve=TRUE)


results <- bind_rows(results,
                     tibble(Method="Random Forest - SMOTE",
                            Accuracy=cm$overall["Accuracy"],
                            Sensitivity= cm$byClass["Sensitivity"], 
                            Specificity= cm$byClass["Specificity"], 
                            F1= cm$byClass["F1"], 
                            AUC=aucv, 
                            AUPRC=AUPRCv[[3]]))

results %>% 
  kable() %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                position = "center",
                font_size = 10,
                full_width = FALSE)



#####################################################################################
#GLM
#####################################################################################
set.seed(1147, sample.kind = "Rounding") #in R 3.6 or later
model_glm <- train(Dataset ~ ., method = "glm", data = train)

predictions <- predict(model_glm, test)
predictions1 <- prediction(as.numeric(predictions) , test$Dataset)
#Save Preds for Ensemble
predictions_glm <- predict(model_glm, test)

cm <- confusionMatrix(predictions, test$Dataset)
auc <- performance(predictions1, "auc")
AUPRC <- performance(predictions1, "prec", "rec")
aucv <- auc@y.values[[1]]
AUPRCv <- pr.curve(scores.class0 = predictions[test$Dataset == 1],
                   scores.class1 = predictions[test$Dataset == 0],
                   curve=TRUE)


results <- bind_rows(results,
                     tibble(Method="GLM",
                            Accuracy=cm$overall["Accuracy"],
                            Sensitivity= cm$byClass["Sensitivity"], 
                            Specificity= cm$byClass["Specificity"], 
                            F1= cm$byClass["F1"], 
                            AUC=aucv, 
                            AUPRC=AUPRCv[[3]]))

results %>% 
  kable() %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                position = "center",
                font_size = 10,
                full_width = FALSE)



#####################################################################################
#KNN
#####################################################################################

control <- trainControl(method ="cv", 
                        number = 10
                        )
set.seed(1147, sample.kind = "Rounding") #in R 3.6 or later 
model_knn <- train(Dataset ~ ., 
                   method = "knn", 
                   data = train,
                   trControl = control,
                   preProcess = c("scale", "center"),
                   tuneGrid = data.frame(k = seq(1,100,4))
                   )


predictions <- predict(model_knn, test)
predictions1 <- prediction(as.numeric(predictions) , test$Dataset)
#Save Preds for Ensemble
predictions_knn <- predict(model_knn, test)

cm <- confusionMatrix(predictions, test$Dataset)
auc <- performance(predictions1, "auc")
AUPRC <- performance(predictions1, "prec", "rec")
aucv <- auc@y.values[[1]]
AUPRCv <- pr.curve(scores.class0 = predictions[test$Dataset == 1],
                   scores.class1 = predictions[test$Dataset == 0],
                   curve=TRUE)


results <- bind_rows(results,
                     tibble(Method="KNN",
                            Accuracy=cm$overall["Accuracy"],
                            Sensitivity= cm$byClass["Sensitivity"], 
                            Specificity= cm$byClass["Specificity"], 
                            F1= cm$byClass["F1"], 
                            AUC=aucv, 
                            AUPRC=AUPRCv[[3]]))

results %>% 
  kable() %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                position = "center",
                font_size = 10,
                full_width = FALSE)


#####################################################################################
#LOOP THROUGH A BUNCH OF OTHER MODELS TO LOOK FOR BEST FOR ENSEMBLE
#####################################################################################
models <- c("glm", "rpart", "lda", "naive_bayes", "svmLinear", "gamLoess", "multinom", "qda", "rf", "adaboost", "xgbTree")
set.seed(1147, sample.kind = "Rounding") #in R 3.6 or later
fits <- lapply(models, function(model){ 
  print(model)
  train(Dataset ~ ., method = model, data = train)
})

names(fits) <- models

set.seed(1147, sample.kind = "Rounding") #in R 3.6 or later
resamples <- resamples(fits)

##Plot CI Accuracy 
dotplot(resamples, metric = "Accuracy", main="Ensemble Resamples by Method")


#####################################################################################
#PREDICT ALL MODELS AND LOOK AT METRICS
#####################################################################################
pred <- sapply(fits, function(object) 
  predict(object, newdata = test))

acc <- colMeans(pred == test$Dataset)

#Add Results to Table
#Get individual models predictions and metrics
predictions1 <- sapply(models, function(d) prediction(as.numeric(pred[,d]) , test$Dataset))
cm <- sapply(models, function(p) confusionMatrix(factor(pred[,p],levels=c(1,0)), test$Dataset))

#Multinom
auc <- performance(predictions1$multinom, "auc")
AUPRC <- performance(predictions1$multinom, "prec", "rec")
aucv <- auc@y.values[[1]]
AUPRCv <- pr.curve(scores.class0 = as.numeric(pred[,"multinom"])[test$Dataset == 1],
                   scores.class1 = as.numeric(pred[,"multinom"])[test$Dataset == 0],
                   curve=TRUE)



results <- bind_rows(results,
                     tibble(Method="Multinom",
                            Accuracy=cm["overall",]$multinom["Accuracy"],
                            Sensitivity= cm["byClass",]$multinom["Sensitivity"], 
                            Specificity= cm["byClass",]$multinom["Specificity"], 
                            F1= cm["byClass",]$multinom["F1"], 
                            AUC=aucv, 
                            AUPRC=AUPRCv[[3]]))


#LDA
auc <- performance(predictions1$lda, "auc")
AUPRC <- performance(predictions1$lda, "prec", "rec")
aucv <- auc@y.values[[1]]
AUPRCv <- pr.curve(scores.class0 = as.numeric(pred[,"lda"])[test$Dataset == 1],
                   scores.class1 = as.numeric(pred[,"lda"])[test$Dataset == 0],
                   curve=TRUE)


results <- bind_rows(results,
                     tibble(Method="LDA",
                            Accuracy=cm["overall",]$lda["Accuracy"],
                            Sensitivity= cm["byClass",]$lda["Sensitivity"], 
                            Specificity= cm["byClass",]$lda["Specificity"], 
                            F1= cm["byClass",]$lda["F1"], 
                            AUC=aucv, 
                            AUPRC=AUPRCv[[3]]))

#Adaboost
auc <- performance(predictions1$adaboost, "auc")
AUPRC <- performance(predictions1$adaboost, "prec", "rec")
aucv <- auc@y.values[[1]]
AUPRCv <- pr.curve(scores.class0 = as.numeric(pred[,"adaboost"])[test$Dataset == 1],
                   scores.class1 = as.numeric(pred[,"adaboost"])[test$Dataset == 0],
                   curve=TRUE)

results <- bind_rows(results,
                     tibble(Method="Adaboost",
                            Accuracy=cm["overall",]$adaboost["Accuracy"],
                            Sensitivity= cm["byClass",]$adaboost["Sensitivity"], 
                            Specificity= cm["byClass",]$adaboost["Specificity"], 
                            F1= cm["byClass",]$adaboost["F1"], 
                            AUC=aucv, 
                            AUPRC=AUPRCv[[3]]))


#SVM
auc <- performance(predictions1$svmLinear, "auc")
AUPRC <- performance(predictions1$svmLinear, "prec", "rec")
AUPRC
aucv <- auc@y.values[[1]]
#AUPRCv <- pr.curve(scores.class0 = as.numeric(pred[,"svmLinear"])[test$Dataset == 1],
#                   scores.class1 = as.numeric(pred[,"svmLinear"])[test$Dataset == 0],
#                   curve=TRUE)


results <- bind_rows(results,
                     tibble(Method="SvmLinear",
                            Accuracy=cm["overall",]$svmLinear["Accuracy"],
                            Sensitivity= cm["byClass",]$svmLinear["Sensitivity"], 
                            Specificity= cm["byClass",]$svmLinear["Specificity"], 
                            F1= cm["byClass",]$svmLinear["F1"], 
                            AUC=aucv, 
                            AUPRC=0))



#rpart
auc <- performance(predictions1$rpart, "auc")
AUPRC <- performance(predictions1$rpart, "prec", "rec")
AUPRC
aucv <- auc@y.values[[1]]
AUPRCv <- pr.curve(scores.class0 = as.numeric(pred[,"rpart"])[test$Dataset == 1],
                   scores.class1 = as.numeric(pred[,"rpart"])[test$Dataset == 0],
                   curve=TRUE)


results <- bind_rows(results,
                     tibble(Method="Rpart",
                            Accuracy=cm["overall",]$rpart["Accuracy"],
                            Sensitivity= cm["byClass",]$rpart["Sensitivity"], 
                            Specificity= cm["byClass",]$rpart["Specificity"], 
                            F1= cm["byClass",]$rpart["F1"], 
                            AUC=aucv, 
                            AUPRC=AUPRCv[[3]]))


#gamLoess
auc <- performance(predictions1$gamLoess, "auc")
AUPRC <- performance(predictions1$gamLoess, "prec", "rec")
AUPRC
aucv <- auc@y.values[[1]]
AUPRCv <- pr.curve(scores.class0 = as.numeric(pred[,"gamLoess"])[test$Dataset == 1],
                   scores.class1 = as.numeric(pred[,"gamLoess"])[test$Dataset == 0],
                   curve=TRUE)


results <- bind_rows(results,
                     tibble(Method="GamLoess",
                            Accuracy=cm["overall",]$gamLoess["Accuracy"],
                            Sensitivity= cm["byClass",]$gamLoess["Sensitivity"], 
                            Specificity= cm["byClass",]$gamLoess["Specificity"], 
                            F1= cm["byClass",]$gamLoess["F1"], 
                            AUC=aucv, 
                            AUPRC=AUPRCv[[3]]))


#qda
auc <- performance(predictions1$qda, "auc")
AUPRC <- performance(predictions1$qda, "prec", "rec")
AUPRC
aucv <- auc@y.values[[1]]
AUPRCv <- pr.curve(scores.class0 = as.numeric(pred[,"qda"])[test$Dataset == 1],
                   scores.class1 = as.numeric(pred[,"qda"])[test$Dataset == 0],
                   curve=TRUE)


results <- bind_rows(results,
                     tibble(Method="QDA",
                            Accuracy=cm["overall",]$qda["Accuracy"],
                            Sensitivity= cm["byClass",]$qda["Sensitivity"], 
                            Specificity= cm["byClass",]$qda["Specificity"], 
                            F1= cm["byClass",]$qda["F1"], 
                            AUC=aucv, 
                            AUPRC=AUPRCv[[3]]))

#xgbTree
auc <- performance(predictions1$xgbTree, "auc")
AUPRC <- performance(predictions1$xgbTree, "prec", "rec")
AUPRC
aucv <- auc@y.values[[1]]
AUPRCv <- pr.curve(scores.class0 = as.numeric(pred[,"xgbTree"])[test$Dataset == 1],
                   scores.class1 = as.numeric(pred[,"xgbTree"])[test$Dataset == 0],
                   curve=TRUE)


results <- bind_rows(results,
                     tibble(Method="XGBTree",
                            Accuracy=cm["overall",]$xgbTree["Accuracy"],
                            Sensitivity= cm["byClass",]$xgbTree["Sensitivity"], 
                            Specificity= cm["byClass",]$xgbTree["Specificity"], 
                            F1= cm["byClass",]$xgbTree["F1"], 
                            AUC=aucv, 
                            AUPRC=AUPRCv[[3]]))


#naive_bayes
auc <- performance(predictions1$naive_bayes, "auc")
AUPRC <- performance(predictions1$naive_bayes, "prec", "rec")
AUPRC
aucv <- auc@y.values[[1]]
AUPRCv <- pr.curve(scores.class0 = as.numeric(pred[,"naive_bayes"])[test$Dataset == 1],
                   scores.class1 = as.numeric(pred[,"naive_bayes"])[test$Dataset == 0],
                   curve=TRUE)


results <- bind_rows(results,
                     tibble(Method="Naive Bayes",
                            Accuracy=cm["overall",]$naive_bayes["Accuracy"],
                            Sensitivity= cm["byClass",]$naive_bayes["Sensitivity"], 
                            Specificity= cm["byClass",]$naive_bayes["Specificity"], 
                            F1= cm["byClass",]$naive_bayes["F1"], 
                            AUC=aucv, 
                            AUPRC=AUPRCv[[3]]))


results %>% 
  kable() %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                position = "center",
                font_size = 10,
                full_width = FALSE)




#Run Top 6 Models Only Now
models <- c("glm", "lda", "svmLinear",  "gamLoess", "multinom", "rf", "adaboost")
set.seed(1147, sample.kind = "Rounding") #in R 3.6 or later 
fits <- lapply(models, function(model){ 
  print(model)
  train(Dataset ~ ., method = model, data = train)
})

names(fits) <- models

#Resample
set.seed(1147, sample.kind = "Rounding") #in R 3.6 or later
resamples <- resamples(fits)

#Plot
dotplot(resamples, metric = "Accuracy")

#Model Correlations
#Examine Model Correlations
modelCor(resamples) %>% 
  kable() %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                position = "center",
                font_size = 10,
                full_width = FALSE)



#####################################################################################
#PREDICT ALL MODELS AND LOOK AT METRICS
#####################################################################################

#Get Ensemble Predictions
votes <- rowMeans(pred == "1")
y_hat <- ifelse(votes > 0.5, "1", "0")
cm <-confusionMatrix(factor(y_hat,levels=c(1,0)), test$Dataset)
predictions1 <- prediction(as.numeric(y_hat) , test$Dataset)

auc <- performance(predictions1, "auc")
AUPRC <- performance(predictions1, "prec", "rec")
aucv <- auc@y.values[[1]]
AUPRCv <- pr.curve(scores.class0 = as.numeric(y_hat)[test$Dataset == 1],
                   scores.class1 = as.numeric(y_hat)[test$Dataset == 0],
                   curve=TRUE)


results <- bind_rows(results,
                     tibble(Method="Ensemble",
                            Accuracy=cm$overall["Accuracy"],
                            Sensitivity= cm$byClass["Sensitivity"], 
                            Specificity= cm$byClass["Specificity"], 
                            F1= cm$byClass["F1"], 
                            AUC=aucv, 
                            AUPRC=AUPRCv[[3]]))


results %>% 
  kable() %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                position = "center",
                font_size = 10,
                full_width = FALSE)


results %>% select (Method, Accuracy, Sensitivity, Specificity, AUPRC) %>%
  gather("Type", "Value",-Method) %>%
  ggplot(aes(Method, Value, fill = Type)) +
  geom_col(position = "dodge") +
  theme_bw()+
  facet_wrap(~Method,scales = "free_x")


#####################################################################################
#CORRELATION MATRIX OF METRICS
#####################################################################################
results1 <-  sapply(results[2:7], as.numeric)
corr_matrix <- round(cor(results1),2)
corr_matrix <- reorder_cormat(corr_matrix)
upper_tri <- upper(corr_matrix)
melt_corr_matrix <- melt(upper_tri, na.rm = TRUE)
ggplot(melt_corr_matrix,aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 9, hjust = 1), axis.text.y = element_text(size = 9),                    axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank()) +
  coord_fixed()

#####################################################################################
#LOOK AT VARIABLE IMPORTANCE ACROSS MODELS
#####################################################################################
imp <- varImp(fits$multinom)
imp1 <- varImp(model_glm)
imp2 <- varImp(model_rf)

a <- tibble(term = rownames(imp$importance), 
           importance = imp$importance$Overall) %>%
  mutate(rank = rank(-importance)) %>% arrange(desc(importance)) 

a1 <- tibble(term = rownames(imp1$importance), 
                 importance = imp1$importance$Overall) %>%
  mutate(rank = rank(-importance)) %>% arrange(desc(importance)) 

a2 <- tibble(term = rownames(imp2$importance), 
                 importance = imp2$importance$Overall) %>%
  mutate(rank = rank(-importance)) %>% arrange(desc(importance)) 


a3 <- inner_join(a,a1, by="term")
a4 <- inner_join(a3, a2, by="term")
a4 <- a4 %>% rename(imp.multi = importance.x, imp.glm=importance.y, imp.rf=importance,
                    rank.multi = rank.x, rank.glm = rank.y, rank.rf=rank)

a4  %>% select (term, imp.multi, imp.glm, imp.rf) %>%
  gather("Type", "Value",-term) %>%
  ggplot(aes(term, Value, fill = Type)) +
  geom_col(position = "dodge") +
  theme_bw()+
  facet_wrap(~term,scales = "free_x")
  
  
a4 %>%  
  kable() %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                position = "center",
                font_size = 10,
                full_width = FALSE)



