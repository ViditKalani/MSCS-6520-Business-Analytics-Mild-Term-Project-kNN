library(tidyverse)
library(caret)
library(lubridate)

read.csv("/Users/vidit/Desktop/Masters/Sem_3/MSCS 6520/Mid-Term_Project/breast-cancer-wisconsin.csv")
bc_data <- read.csv("/Users/vidit/Desktop/Masters/Sem_3/MSCS 6520/Mid-Term_Project/breast-cancer-wisconsin.csv")

colnames(bc_data) <- c("sample_code_number", "clump_thickness", "uniformity_of_cell_size", "uniformity_of_cell_shape", "marginal_adhesion", "single_epithelial_cell_size", 
                       "bare_nuclei", "bland_chromatin", "normal_nucleoli", "mitosis", "classes")
bc_data$classes <- ifelse(bc_data$classes == "2", "benign",
                          ifelse(bc_data$classes == "4", "malignant", NA))

bc_data[bc_data == "?"] <- NA

# how many NAs are in the data
length(which(is.na(bc_data)))

# impute missing data
library(mice)

bc_data[,2:10] <- apply(bc_data[, 2:10], 2, function(x) as.numeric(as.character(x)))
dataset_impute <- mice(bc_data[, 2:10],  print = FALSE)
bc_data <- cbind(bc_data[, 11, drop = FALSE], mice::complete(dataset_impute, 1))

bc_data$classes <- as.factor(bc_data$classes)

# how many benign and malignant cases are there?
summary(bc_data$classes)

select(bc_data,
       classes,
       clump_thickness,
       uniformity_of_cell_size,
       uniformity_of_cell_shape,
       marginal_adhesion,
       single_epithelial_cell_size,
       bare_nuclei,
       bland_chromatin,
       normal_nucleoli,
       mitosis)

bc_data <- select(bc_data,
                  classes,
                  clump_thickness,
                  uniformity_of_cell_size,
                  uniformity_of_cell_shape,
                  marginal_adhesion,
                  single_epithelial_cell_size,
                  bare_nuclei,
                  bland_chromatin,
                  normal_nucleoli,
                  mitosis)

ggplot(data = bc_data)+
  geom_boxplot(mapping = aes(x = classes, y = clump_thickness))

ggplot(data = bc_data)+
  geom_boxplot(mapping = aes(x = classes, y = uniformity_of_cell_size))

ggplot(data = bc_data)+
  geom_boxplot(mapping = aes(x = classes, y = uniformity_of_cell_shape))

ggplot(data = bc_data)+
  geom_boxplot(mapping = aes(x = classes, y = marginal_adhesion))

ggplot(data = bc_data) +
  geom_boxplot(mapping = aes(x = classes, y = single_epithelial_cell_size))

ggplot(data = bc_data) +
  geom_boxplot(mapping = aes(x = classes, y = bare_nuclei))

ggplot(data = bc_data) +
  geom_boxplot(mapping = aes(x = classes, y = bland_chromatin))

ggplot(data = bc_data) +
  geom_boxplot(mapping = aes(x = classes, y = normal_nucleoli))

ggplot(data = bc_data) +
  geom_boxplot(mapping = aes(x = classes, y = mitosis))


#geom_jitter

ggplot(data = bc_data)+
  geom_jitter(mapping = aes(x = classes, y = clump_thickness))

ggplot(data = bc_data)+
  geom_jitter(mapping = aes(x = classes, y = uniformity_of_cell_size))

ggplot(data = bc_data)+
  geom_jitter(mapping = aes(x = classes, y = uniformity_of_cell_shape))

ggplot(data = bc_data)+
  geom_jitter(mapping = aes(x = classes, y = marginal_adhesion))

ggplot(data = bc_data) +
  geom_jitter(mapping = aes(x = classes, y = single_epithelial_cell_size))

ggplot(data = bc_data) +
  geom_jitter(mapping = aes(x = classes, y = bare_nuclei))

ggplot(data = bc_data) +
  geom_jitter(mapping = aes(x = classes, y = bland_chromatin))

ggplot(data = bc_data) +
  geom_jitter(mapping = aes(x = classes, y = normal_nucleoli))

ggplot(data = bc_data) +
  geom_jitter(mapping = aes(x = classes, y = mitosis))


set.seed(77889)
trainIndex <- createDataPartition(bc_data$classes, p = 0.8, list = FALSE, times = 1)

bc_data[trainIndex, ]
bc_dataTrain <- bc_data[trainIndex, ]

bc_data[-trainIndex, ]
bc_dataTest <- bc_data[-trainIndex, ]

preProcess(bc_dataTrain, method = c("center", "scale"))
scaler <- preProcess(bc_dataTrain, method = c("center", "scale"))

predict(scaler, bc_dataTrain)
bc_dataTrain <- predict(scaler, bc_dataTrain)
predict(scaler, bc_dataTest)
bc_dataTest <- predict(scaler, bc_dataTest)
head(bc_dataTrain)


#Round:- 1

train(classes ~ clump_thickness, data = bc_dataTrain, method = "knn")
knn_model_clump_thickness <- train(classes ~ clump_thickness, data = bc_dataTrain, method = "knn")
predict(knn_model_clump_thickness, bc_dataTest)
bc_dataTestPredictions_clump_thickness <- predict(knn_model_clump_thickness, bc_dataTest)
confusionMatrix(bc_dataTestPredictions_clump_thickness, bc_dataTest$classes)

train(classes ~ uniformity_of_cell_size, data = bc_dataTrain, method = "knn")
knn_model_uniformity_of_cell_size <- train(classes ~ uniformity_of_cell_size, data = bc_dataTrain, method = "knn")
predict(knn_model_uniformity_of_cell_size, bc_dataTest)
bc_dataTestPredictions_uniformity_of_cell_size <- predict(knn_model_uniformity_of_cell_size, bc_dataTest)
confusionMatrix(bc_dataTestPredictions_uniformity_of_cell_size, bc_dataTest$classes)

train(classes ~ uniformity_of_cell_shape, data = bc_dataTrain, method = "knn")
knn_model_uniformity_of_cell_shape <- train(classes ~ uniformity_of_cell_shape, data = bc_dataTrain, method = "knn")
predict(knn_model_uniformity_of_cell_shape, bc_dataTest)
bc_dataTestPredictions_uniformity_of_cell_shape <- predict(knn_model_uniformity_of_cell_shape, bc_dataTest)
confusionMatrix(bc_dataTestPredictions_uniformity_of_cell_shape, bc_dataTest$classes)

train(classes ~ marginal_adhesion, data = bc_dataTrain, method = "knn")
knn_model_marginal_adhesion <- train(classes ~ marginal_adhesion, data = bc_dataTrain, method = "knn")
predict(knn_model_marginal_adhesion, bc_dataTest)
bc_dataTestPredictions_marginal_adhesion <- predict(knn_model_marginal_adhesion, bc_dataTest)
confusionMatrix(bc_dataTestPredictions_marginal_adhesion, bc_dataTest$classes)

train(classes ~ single_epithelial_cell_size, data = bc_dataTrain, method = "knn")
knn_model_single_epithelial_cell_size <- train(classes ~ single_epithelial_cell_size, data = bc_dataTrain, method = "knn")
predict(knn_model_single_epithelial_cell_size, bc_dataTest)
bc_dataTestPredictions_single_epithelial_cell_size <- predict(knn_model_single_epithelial_cell_size, bc_dataTest)
confusionMatrix(bc_dataTestPredictions_single_epithelial_cell_size, bc_dataTest$classes)

train(classes ~ bare_nuclei, data = bc_dataTrain, method = "knn")
knn_model_bare_nuclei <- train(classes ~ bare_nuclei, data = bc_dataTrain, method = "knn")
predict(knn_model_bare_nuclei, bc_dataTest)
bc_dataTestPredictions_bare_nuclei <- predict(knn_model_bare_nuclei, bc_dataTest)
confusionMatrix(bc_dataTestPredictions_bare_nuclei, bc_dataTest$classes)

train(classes ~ bland_chromatin, data = bc_dataTrain, method = "knn")
knn_model_bland_chromatin <- train(classes ~ bland_chromatin, data = bc_dataTrain, method = "knn")
predict(knn_model_bland_chromatin, bc_dataTest)
bc_dataTestPredictions_bland_chromatin <- predict(knn_model_bland_chromatin, bc_dataTest)
confusionMatrix(bc_dataTestPredictions_bland_chromatin,  bc_dataTest$classes)

train(classes ~ normal_nucleoli, data = bc_dataTrain, method = "knn")
knn_model_normal_nucleoli <- train(classes ~ normal_nucleoli, data = bc_dataTrain, method = "knn")
predict(knn_model_normal_nucleoli, bc_dataTest)
bc_dataTestPredictions_normal_nucleoli <- predict(knn_model_normal_nucleoli, bc_dataTest)
confusionMatrix(bc_dataTestPredictions_normal_nucleoli,  bc_dataTest$classes)

train(classes ~ mitosis, data = bc_dataTrain, method = "knn")
knn_model_mitosis <- train(classes ~ mitosis, data = bc_dataTrain, method = "knn")
predict(knn_model_mitosis, bc_dataTest)
bc_dataTestPredictions_mitosis <- predict(knn_model_mitosis, bc_dataTest)
confusionMatrix(bc_dataTestPredictions_mitosis,  bc_dataTest$classes)


#Round:- 2

train(classes ~ uniformity_of_cell_size + 
                clump_thickness, 
                data = bc_dataTrain, 
                method = "knn")
knn_model_uniformity_of_cell_size_clump_thickness <- train(classes ~ uniformity_of_cell_size + 
                                                                     clump_thickness, 
                                                                     data = bc_dataTrain, 
                                                                     method = "knn")
predict(knn_model_uniformity_of_cell_size_clump_thickness, bc_dataTest)
bc_dataTestPredictions_uniformity_of_cell_size_clump_thickness <- predict(knn_model_uniformity_of_cell_size_clump_thickness, bc_dataTest)
confusionMatrix(bc_dataTestPredictions_uniformity_of_cell_size_clump_thickness, bc_dataTest$classes)


train(classes ~ uniformity_of_cell_size + 
                uniformity_of_cell_shape, 
                data = bc_dataTrain, 
                method = "knn")
knn_model_uniformity_of_cell_size_uniformity_of_cell_shape <- train(classes ~ uniformity_of_cell_size + 
                                                                      uniformity_of_cell_shape, 
                                                                      data = bc_dataTrain, 
                                                                      method = "knn")
predict(knn_model_uniformity_of_cell_size_uniformity_of_cell_shape, bc_dataTest)
bc_dataTestPredictions_uniformity_of_cell_size_uniformity_of_cell_shape <- predict(knn_model_uniformity_of_cell_size_uniformity_of_cell_shape, bc_dataTest)
confusionMatrix(bc_dataTestPredictions_uniformity_of_cell_size_uniformity_of_cell_shape, bc_dataTest$classes)


train(classes ~ uniformity_of_cell_size + 
                marginal_adhesion , 
                data = bc_dataTrain, 
                method = "knn")
knn_model_uniformity_of_cell_size_marginal_adhesion  <- train(classes ~ uniformity_of_cell_size + 
                                                                        marginal_adhesion , 
                                                                        data = bc_dataTrain, 
                                                                        method = "knn")
predict(knn_model_uniformity_of_cell_size_marginal_adhesion , bc_dataTest)
bc_dataTestPredictions_uniformity_of_cell_size_marginal_adhesion  <- predict(knn_model_uniformity_of_cell_size_marginal_adhesion , bc_dataTest)
confusionMatrix(bc_dataTestPredictions_uniformity_of_cell_size_marginal_adhesion , bc_dataTest$classes)


train(classes ~ uniformity_of_cell_size + 
                single_epithelial_cell_size  , 
                data = bc_dataTrain, 
                method = "knn")
knn_model_uniformity_of_cell_size_single_epithelial_cell_size <- train(classes ~ uniformity_of_cell_size + 
                                                                                single_epithelial_cell_size  , 
                                                                                data = bc_dataTrain, 
                                                                                method = "knn")
predict(knn_model_uniformity_of_cell_size_single_epithelial_cell_size  , bc_dataTest)
bc_dataTestPredictions_uniformity_of_cell_size_single_epithelial_cell_size   <- predict(knn_model_uniformity_of_cell_size_single_epithelial_cell_size  , bc_dataTest)
confusionMatrix(bc_dataTestPredictions_uniformity_of_cell_size_single_epithelial_cell_size  , bc_dataTest$classes)


train(classes ~ uniformity_of_cell_size + 
                bare_nuclei, 
                data = bc_dataTrain, 
                method = "knn")
knn_model_uniformity_of_cell_size_bare_nuclei  <- train(classes ~ uniformity_of_cell_size + 
                                                                  bare_nuclei, 
                                                                  data = bc_dataTrain, 
                                                                  method = "knn")
predict(knn_model_uniformity_of_cell_size_bare_nuclei, bc_dataTest)
bc_dataTestPredictions_uniformity_of_cell_size_bare_nuclei <- predict(knn_model_uniformity_of_cell_size_bare_nuclei, bc_dataTest)
confusionMatrix(bc_dataTestPredictions_uniformity_of_cell_size_bare_nuclei, bc_dataTest$classes)


train(classes ~ uniformity_of_cell_size + 
                bland_chromatin, 
                data = bc_dataTrain, 
                method = "knn")
knn_model_uniformity_of_cell_size_bland_chromatin <- train(classes ~ uniformity_of_cell_size + 
                                                                     bland_chromatin , 
                                                                     data = bc_dataTrain, 
                                                                     method = "knn")
predict(knn_model_uniformity_of_cell_size_bland_chromatin , bc_dataTest)
bc_dataTestPredictions_uniformity_of_cell_size_bland_chromatin <- predict(knn_model_uniformity_of_cell_size_bland_chromatin , bc_dataTest)
confusionMatrix(bc_dataTestPredictions_uniformity_of_cell_size_bland_chromatin, bc_dataTest$classes)


train(classes ~ uniformity_of_cell_size + 
                normal_nucleoli , 
                data = bc_dataTrain, 
                method = "knn")
knn_model_uniformity_of_cell_size_normal_nucleoli  <- train(classes ~ uniformity_of_cell_size + 
                                                                      normal_nucleoli, 
                                                                      data = bc_dataTrain, 
                                                                      method = "knn")
predict(knn_model_uniformity_of_cell_size_normal_nucleoli, bc_dataTest)
bc_dataTestPredictions_uniformity_of_cell_size_normal_nucleoli <- predict(knn_model_uniformity_of_cell_size_normal_nucleoli, bc_dataTest)
confusionMatrix(bc_dataTestPredictions_uniformity_of_cell_size_normal_nucleoli, bc_dataTest$classes)


train(classes ~ uniformity_of_cell_size + 
                mitosis, 
                data = bc_dataTrain, 
                method = "knn")
knn_model_uniformity_of_cell_size_mitosis <- train(classes ~ uniformity_of_cell_size + 
                                                             mitosis, 
                                                             data = bc_dataTrain, 
                                                             method = "knn")
predict(knn_model_uniformity_of_cell_size_mitosis, bc_dataTest)
bc_dataTestPredictions_uniformity_of_cell_size_mitosis <- predict(knn_model_uniformity_of_cell_size_mitosis, bc_dataTest)
confusionMatrix(bc_dataTestPredictions_uniformity_of_cell_size_mitosis, bc_dataTest$classes)


#Round:- 3

train(classes ~ uniformity_of_cell_size + 
                clump_thickness +
                uniformity_of_cell_shape,
                data = bc_dataTrain, 
                method = "knn")
knn_model_uniformity_of_cell_size_clump_thickness_uniformity_of_cell_shape <- train(classes ~ uniformity_of_cell_size + 
                                                                                              clump_thickness +
                                                                                              uniformity_of_cell_shape, 
                                                                                              data = bc_dataTrain, 
                                                                                              method = "knn")
predict(knn_model_uniformity_of_cell_size_clump_thickness_uniformity_of_cell_shape, bc_dataTest)
bc_dataTestPredictions_uniformity_of_cell_size_clump_thickness_uniformity_of_cell_shape <- predict(knn_model_uniformity_of_cell_size_clump_thickness_uniformity_of_cell_shape, bc_dataTest)
confusionMatrix(bc_dataTestPredictions_uniformity_of_cell_size_clump_thickness_uniformity_of_cell_shape, bc_dataTest$classes)

train(classes ~ uniformity_of_cell_size + 
                clump_thickness +
                marginal_adhesion,
                data = bc_dataTrain, 
                method = "knn")
knn_model_uniformity_of_cell_size_clump_thickness_marginal_adhesion <- train(classes ~ uniformity_of_cell_size + 
                                                                                              clump_thickness +
                                                                                              marginal_adhesion, 
                                                                                              data = bc_dataTrain, 
                                                                                              method = "knn")
predict(knn_model_uniformity_of_cell_size_clump_thickness_marginal_adhesion, bc_dataTest)
bc_dataTestPredictions_uniformity_of_cell_size_clump_thickness_marginal_adhesion <- predict(knn_model_uniformity_of_cell_size_clump_thickness_marginal_adhesion, bc_dataTest)
confusionMatrix(bc_dataTestPredictions_uniformity_of_cell_size_clump_thickness_marginal_adhesion, bc_dataTest$classes)

train(classes ~ uniformity_of_cell_size + 
                clump_thickness +
                single_epithelial_cell_size,
                data = bc_dataTrain, 
                method = "knn")
knn_model_uniformity_of_cell_size_clump_thickness_single_epithelial_cell_size <- train(classes ~ uniformity_of_cell_size + 
                                                                                                 clump_thickness +
                                                                                                 single_epithelial_cell_size, 
                                                                                                 data = bc_dataTrain, 
                                                                                                method = "knn")
predict(knn_model_uniformity_of_cell_size_clump_thickness_single_epithelial_cell_size, bc_dataTest)
bc_dataTestPredictions_uniformity_of_cell_size_clump_thickness_single_epithelial_cell_size <- predict(knn_model_uniformity_of_cell_size_clump_thickness_single_epithelial_cell_size, bc_dataTest)
confusionMatrix(bc_dataTestPredictions_uniformity_of_cell_size_clump_thickness_single_epithelial_cell_size, bc_dataTest$classes)

train(classes ~ uniformity_of_cell_size + 
                clump_thickness +
                bare_nuclei,
                data = bc_dataTrain, 
                method = "knn")
knn_model_uniformity_of_cell_size_clump_thickness_bare_nuclei <- train(classes ~ uniformity_of_cell_size + 
                                                                                                 clump_thickness +
                                                                                                 bare_nuclei, 
                                                                                                 data = bc_dataTrain, 
                                                                                                 method = "knn")
predict(knn_model_uniformity_of_cell_size_clump_thickness_bare_nuclei, bc_dataTest)
bc_dataTestPredictions_uniformity_of_cell_size_clump_thickness_bare_nuclei <- predict(knn_model_uniformity_of_cell_size_clump_thickness_bare_nuclei, bc_dataTest)
confusionMatrix(bc_dataTestPredictions_uniformity_of_cell_size_clump_thickness_bare_nuclei, bc_dataTest$classes)

train(classes ~ uniformity_of_cell_size + 
                clump_thickness +
                bland_chromatin,
                data = bc_dataTrain, 
                method = "knn")
knn_model_uniformity_of_cell_size_clump_thickness_bland_chromatin <- train(classes ~ uniformity_of_cell_size + 
                                                                                     clump_thickness +
                                                                                     bland_chromatin, 
                                                                                     data = bc_dataTrain, 
                                                                                     method = "knn")
predict(knn_model_uniformity_of_cell_size_clump_thickness_bland_chromatin, bc_dataTest)
bc_dataTestPredictions_uniformity_of_cell_size_clump_thickness_bland_chromatin <- predict(knn_model_uniformity_of_cell_size_clump_thickness_bland_chromatin, bc_dataTest)
confusionMatrix(bc_dataTestPredictions_uniformity_of_cell_size_clump_thickness_bland_chromatin, bc_dataTest$classes)

train(classes ~ uniformity_of_cell_size + 
                clump_thickness +
                normal_nucleoli,
                data = bc_dataTrain, 
                method = "knn")
knn_model_uniformity_of_cell_size_clump_thickness_normal_nucleoli <- train(classes ~ uniformity_of_cell_size + 
                                                                                     clump_thickness +
                                                                                     normal_nucleoli, 
                                                                                     data = bc_dataTrain, 
                                                                                     method = "knn")
predict(knn_model_uniformity_of_cell_size_clump_thickness_normal_nucleoli, bc_dataTest)
bc_dataTestPredictions_uniformity_of_cell_size_clump_thickness_normal_nucleoli <- predict(knn_model_uniformity_of_cell_size_clump_thickness_normal_nucleoli, bc_dataTest)
confusionMatrix(bc_dataTestPredictions_uniformity_of_cell_size_clump_thickness_normal_nucleoli, bc_dataTest$classes)
