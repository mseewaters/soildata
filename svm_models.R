library(e1071)
library(performanceEstimation)

train <- read.csv("training.csv")
test <- read.csv("sorted_test.csv")
soil_properties <- c("Ca", "P", "pH", "SOC", "Sand")

train$Depth <- as.numeric(train$Depth)
test$Depth <- as.numeric(test$Depth)


# take the first derivatives to smoothe out the measurement noise
# training data
MIR_measurements <- train[, 2671:3579]
MIR_DER <- MIR_measurements- cbind(NA, MIR_measurements)[, -(dim(MIR_measurements)[2]+1)]
X_train <- cbind(train[, 3580:3595], MIR_DER[, -1])

# testing data
MIR_measurements <- test[, 2671:3579]
MIR_DER <- MIR_measurements- cbind(NA, MIR_measurements)[, -(dim(MIR_measurements)[2]+1)]
X_test <- cbind(test[, 3580:3595], MIR_DER[, -1])



predictions <- rep(NA, dim(test)[1])
for(soil_property in soil_properties){
  
  print(soil_property)
  sp.col <- which(names(train)==soil_property)
  
  train.d1 <- cbind(X_train,train[ ,sp.col])
  
  model.svm <- svm(train.d1[,-ncol(train.d1)], train.d1[,ncol(train.d1)], cost=100, gamma=0.0001)
  pred.svm <- predict(model.svm, X_test)
  predictions <- cbind(predictions, pred.svm)  
  
}

predictions <- predictions[,-1]
colnames(predictions) <-  soil_properties
write.csv(cbind(PIDN= as.character(test[,1]), predictions), "predictions.csv", row.names=FALSE)
