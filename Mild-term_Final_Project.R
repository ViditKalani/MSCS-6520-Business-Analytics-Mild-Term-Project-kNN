install.packages("corrplot")
install.packages("gridExtra")
library(tidyverse)
library(caret)
library(mlbench)
library(corrplot)

bc_data<- read_csv("/Users/vidit/Desktop/Masters/Sem_3/MSCS 6520/Mid-Term_Project/breast-cancer-wisconsin.csv")
colnames(bc_data) <- c("sample_code_number", "clump_thickness", "uniformity_of_cell_size", "uniformity_of_cell_shape", "marginal_adhesion", "single_epithelial_cell_size", 
                       "bare_nuclei", "bland_chromatin", "normal_nucleoli", "mitosis", "classes")

bc_data<-as.tibble(bc_data)
head(bc_data)
bc_data$classes <- ifelse(bc_data$classes == "2", "benign",
                          ifelse(bc_data$classes == "4", "malignant", NA))
bc_data$classes <- as.factor(bc_data$classes)

bc_data[bc_data == "?"] <- NA
colnames(bc_data)[colSums(is.na(bc_data))>0]
map_int(bc_data,function(.x) sum(is.na(.x)))
# how many NAs are in the data
length(which(is.na(bc_data)))

# impute missing data
library(mice)

bc_data[,2:10] <- apply(bc_data[, 2:10], 2, function(x) as.numeric(as.character(x)))
dataset_impute <- mice(bc_data[, 2:10],  print = FALSE)
bc_data <- cbind(bc_data[, 11, drop = FALSE], mice::complete(dataset_impute, 1))

bc_data$classes <- as.factor(bc_data$classes)

#Find thecorrelation between variables
cor_1<-cor(bc_data[-1],method="pearson")
corrplot(cor_1,method="number",type="lower")
corrplot(cor_1,type="lower")

#Filter out the columns having correlation value greater than 0.9

highly_correlated<-findCorrelation(cor_1,cutoff=0.9)
print(highly_correlated)
exclude<-colnames(cor_1)[highly_correlated]

bc_data_new<-bc_data[,setdiff(colnames(bc_data), exclude)]

# Jitter plot of the remaining variables
var_list<-names(bc_data_new)

library("gridExtra")

for (i in 2:length(var_list)) {
  gs <- lapply(2:9, function(i) 
    ggplot(data=bc_data_new,aes_string(x=var_list[[1]], y=var_list[[i]])) +
      geom_jitter(mapping = aes(color=classes)))
}
grid.arrange(grobs=gs, ncol=3)

# Box plot of the variables to see the distribution of data
for (i in 2:length(var_list)) {
  gs <- lapply(2:9, function(i) 
    ggplot(data = bc_data_new,aes_string(x=var_list[[1]], y=var_list[[i]])) +
      geom_boxplot(mapping = aes(color=classes)))
  
}
grid.arrange(grobs=gs, ncol=3)

#dev.off()

#Partition the data into training and test sets
set.seed(77889)
trainIndex <- createDataPartition(bc_data_new$classes, p = 0.8, list = FALSE, times = 1)

bc_data_new[trainIndex, ]
bc_dataTrain <- bc_data_new[trainIndex, ]

bc_data_new[-trainIndex, ]
bc_dataTest <- bc_data_new[-trainIndex, ]

df_control <- trainControl(method="cv",
                           number = 15,
                           classProbs = TRUE,
                           summaryFunction = twoClassSummary)

#Standardize the data
preProcess(bc_dataTrain, method = c("center", "scale"))
scaler <- preProcess(bc_dataTrain, method = c("center", "scale"))

predict(scaler, bc_dataTrain)
bc_dataTrain <- predict(scaler, bc_dataTrain)
predict(scaler, bc_dataTest)
bc_dataTest <- predict(scaler, bc_dataTest)
head(bc_dataTrain)

#Round:- 1

knnmodel<-train(classes ~ clump_thickness,
                data = bc_dataTrain, 
                metric = "ROC", 
                trControl = df_control,
                method = "knn")
bc_dataTestPredictions <- predict(knnmodel, bc_dataTest)
confusionMatrix(bc_dataTestPredictions, bc_dataTest$classes,positive="malignant")

knnmodel<-train(classes ~ uniformity_of_cell_shape,
                data = bc_dataTrain, 
                metric = "ROC", 
                trControl = df_control,
                method = "knn")
bc_dataTestPredictions <- predict(knnmodel, bc_dataTest)
confusionMatrix(bc_dataTestPredictions, bc_dataTest$classes,positive="malignant")

knnmodel<-train(classes ~ marginal_adhesion,
                data = bc_dataTrain, 
                metric = "ROC", 
                trControl = df_control,
                method = "knn")
bc_dataTestPredictions <- predict(knnmodel, bc_dataTest)
confusionMatrix(bc_dataTestPredictions, bc_dataTest$classes,positive="malignant")

knnmodel<-train(classes ~ single_epithelial_cell_size,
                data = bc_dataTrain, 
                metric = "ROC", 
                trControl = df_control,
                method = "knn")
bc_dataTestPredictions <- predict(knnmodel, bc_dataTest)
confusionMatrix(bc_dataTestPredictions, bc_dataTest$classes,positive="malignant")

knnmodel<-train(classes ~ bare_nuclei,
                data = bc_dataTrain, 
                metric = "ROC", 
                trControl = df_control,
                method = "knn")


bc_dataTestPredictions <- predict(knnmodel, bc_dataTest)
confusionMatrix(bc_dataTestPredictions, bc_dataTest$classes,positive="malignant")

knnmodel<-train(classes ~ bland_chromatin,
                data = bc_dataTrain, 
                metric = "ROC", 
                trControl = df_control,
                method = "knn")
bc_dataTestPredictions <- predict(knnmodel, bc_dataTest)
confusionMatrix(bc_dataTestPredictions, bc_dataTest$classes,positive="malignant")

knnmodel<-train(classes ~ normal_nucleoli,
                data = bc_dataTrain, 
                metric = "ROC", 
                trControl = df_control,
                method = "knn")
bc_dataTestPredictions <- predict(knnmodel, bc_dataTest)
confusionMatrix(bc_dataTestPredictions, bc_dataTest$classes,positive="malignant")


#Round:- 2

knnmodel<-train(classes ~ single_epithelial_cell_size+
                  clump_thickness,
                data = bc_dataTrain, 
                metric = "ROC", 
                trControl = df_control,
                method = "knn")
bc_dataTestPredictions <- predict(knnmodel, bc_dataTest)
confusionMatrix(bc_dataTestPredictions, bc_dataTest$classes,positive="malignant")

knnmodel<-train(classes ~ single_epithelial_cell_size+
                  uniformity_of_cell_shape,
                data = bc_dataTrain, 
                metric = "ROC", 
                trControl = df_control,
                method = "knn")
bc_dataTestPredictions <- predict(knnmodel, bc_dataTest)
confusionMatrix(bc_dataTestPredictions, bc_dataTest$classes,positive="malignant")

knnmodel<-train(classes ~ single_epithelial_cell_size+
                  marginal_adhesion,
                data = bc_dataTrain, 
                metric = "ROC", 
                trControl = df_control,
                method = "knn")
bc_dataTestPredictions <- predict(knnmodel, bc_dataTest)
confusionMatrix(bc_dataTestPredictions, bc_dataTest$classes,positive="malignant")

knnmodel<-train(classes ~ bland_chromatin+
                  single_epithelial_cell_size,
                data = bc_dataTrain, 
                metric = "ROC", 
                trControl = df_control,
                method = "knn")
bc_dataTestPredictions <- predict(knnmodel, bc_dataTest)
confusionMatrix(bc_dataTestPredictions, bc_dataTest$classes,positive="malignant")

knnmodel<-train(classes ~ single_epithelial_cell_size+
                  bare_nuclei,
                data = bc_dataTrain,
                metric = "ROC", 
                trControl = df_control,
                method = "knn")
bc_dataTestPredictions <- predict(knnmodel, bc_dataTest)
confusionMatrix(bc_dataTestPredictions, bc_dataTest$classes,positive="malignant")

knnmodel<-train(classes ~ single_epithelial_cell_size+
                  normal_nucleoli,
                data = bc_dataTrain,
                metric = "ROC", 
                trControl = df_control,
                method = "knn")
bc_dataTestPredictions <- predict(knnmodel, bc_dataTest)
confusionMatrix(bc_dataTestPredictions, bc_dataTest$classes,positive="malignant")


#Round 3

knnmodel<-train(classes ~ bland_chromatin+
                  single_epithelial_cell_size+
                  clump_thickness,
                data = bc_dataTrain,
                metric = "ROC", 
                trControl = df_control,
                method = "knn")
bc_dataTestPredictions <- predict(knnmodel, bc_dataTest)
confusionMatrix(bc_dataTestPredictions, bc_dataTest$classes,positive="malignant")


knnmodel<-train(classes ~ bland_chromatin+
                  single_epithelial_cell_size+
                  normal_nucleoli,
                data = bc_dataTrain,
                metric = "ROC", 
                trControl = df_control,
                method = "knn")
bc_dataTestPredictions <- predict(knnmodel, bc_dataTest)
confusionMatrix(bc_dataTestPredictions, bc_dataTest$classes,positive="malignant")

knnmodel<-train(classes ~ bland_chromatin+
                  single_epithelial_cell_size+
                  uniformity_of_cell_shape,
                data = bc_dataTrain,
                metric = "ROC", 
                trControl = df_control,
                method = "knn")
bc_dataTestPredictions <- predict(knnmodel, bc_dataTest)
confusionMatrix(bc_dataTestPredictions, bc_dataTest$classes,positive="malignant")

knnmodel<-train(classes ~ bland_chromatin+
                  single_epithelial_cell_size+
                  marginal_adhesion,
                data = bc_dataTrain,
                metric = "ROC", 
                trControl = df_control,
                method = "knn")
bc_dataTestPredictions <- predict(knnmodel, bc_dataTest)
confusionMatrix(bc_dataTestPredictions, bc_dataTest$classes,positive="malignant")

knnmodel<-train(classes ~ bland_chromatin+
                  single_epithelial_cell_size+
                  bare_nuclei,
                data = bc_dataTrain,
                metric = "ROC", 
                trControl = df_control,
                method = "knn")
bc_dataTestPredictions <- predict(knnmodel, bc_dataTest)
confusionMatrix(bc_dataTestPredictions, bc_dataTest$classes,positive="malignant")



#Round4

knnmodel<-train(classes ~ bland_chromatin+
                  single_epithelial_cell_size+
                  clump_thickness+
                  normal_nucleoli,
                data = bc_dataTrain,
                metric = "ROC", 
                trControl = df_control,
                method = "knn")
bc_dataTestPredictions <- predict(knnmodel, bc_dataTest)
confusionMatrix(bc_dataTestPredictions, bc_dataTest$classes,positive="malignant")

knnmodel<-train(classes ~ bland_chromatin+
                  single_epithelial_cell_size+
                  clump_thickness+
                  bare_nuclei,
                data = bc_dataTrain,
                metric = "ROC", 
                trControl = df_control,
                method = "knn")
bc_dataTestPredictions <- predict(knnmodel, bc_dataTest)
confusionMatrix(bc_dataTestPredictions, bc_dataTest$classes,positive="malignant")

knnmodel<-train(classes ~ bland_chromatin+
                  single_epithelial_cell_size+
                  normal_nucleoli+
                  bare_nuclei,
                data = bc_dataTrain,
                metric = "ROC", 
                trControl = df_control,
                method = "knn")
bc_dataTestPredictions <- predict(knnmodel, bc_dataTest)
confusionMatrix(bc_dataTestPredictions, bc_dataTest$classes,positive="malignant")

