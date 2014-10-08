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
  wgts <- switch(soil_property,
                 "Ca" = c(0.5,0.2,0.3),
                 "P" = c(0.2,0.6,0.2),
                 "pH" = c(0.1,0.4,0.5),
                 "SOC" = c(0.5,0.2,0.3),
                 "Sand" = c(0.2,0.5,0.3))
  
  train.d1 <- cbind(X_train,train[ ,sp.col])
  
  model.rf <- randomForest(train.d1[,-ncol(train.d1)], train.d1[,ncol(train.d1)], ntree=500, mtry=10)
  pred.rf <- predict(model.rf, X_test)
  pred.w.rf <- pred.rf * wgts[1]
  
  model.svm1 <- svm(train.d1[,-ncol(train.d1)], train.d1[,ncol(train.d1)], cost=100, gamma=0.001)
  pred.svm1 <- predict(model.svm1, X_test)
  pred.w.svm1 <- pred.svm1 * wgts[2]
  
  model.svm2 <- svm(train.d1[,-ncol(train.d1)], train.d1[,ncol(train.d1)], cost=1000, gamma=0.0001)
  pred.svm2 <- predict(model.svm2, X_test)
  pred.w.svm2 <- pred.svm2 * wgts[3]
  
  pred.all <- cbind(pred.w.rf,pred.w.svm1,pred.w.svm2)
  pred.avg <- rowSums(pred.all)
  
  predictions <- cbind(predictions, pred.avg)  
  
}

predictions <- predictions[,-1]
colnames(predictions) <-  soil_properties
write.csv(cbind(PIDN= as.character(test[,1]), predictions), "predictions.csv", row.names=FALSE)
