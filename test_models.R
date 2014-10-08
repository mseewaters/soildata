library(e1071)
library(randomForest)
library(rpart)
library(performanceEstimation)
library(DMwR)
library(psych)
require(medley)
library(earth)

train <- read.csv("training.csv")
#test <- read.csv("sorted_test.csv")
soil_properties <- c("Ca", "P", "pH", "SOC", "Sand")

train$Depth <- as.numeric(train$Depth)
test$Depth <- as.numeric(test$Depth)

# CO2_bands <- 2656:2670
names(train)[2656:2670]
which(names(train)==soil_properties)

# take the first derivatives to smoothe out the measurement noise
# training data
MIR_measurements <- train[, 2:2655]
MIR_DER <- MIR_measurements- cbind(NA, MIR_measurements)[, -(dim(MIR_measurements)[2]+1)]
X_train <- cbind(train[, 3580:3595], MIR_DER[,-1])
MIR_measurements <- train[, 2671:3579]
MIR_DER <- MIR_measurements- cbind(NA, MIR_measurements)[, -(dim(MIR_measurements)[2]+1)]
X_train <- cbind(X_train, MIR_DER[, -1])


# training data
MIR_measurements <- train[, 2671:3579]
MIR_DER <- MIR_measurements- cbind(NA, MIR_measurements)[, -(dim(MIR_measurements)[2]+1)]
X_train1 <- cbind(train[, 3580:3595], MIR_DER[, -1])


# training data
MIR_measurements <- X_train1[, 17:924]
MIR_DER <- MIR_measurements- cbind(NA, MIR_measurements)[, -(dim(MIR_measurements)[2]+1)]
X_train2 <- cbind(train[, 3580:3595], MIR_DER[, -1])


# testing data
MIR_measurements <- test[, 2:2655]
MIR_DER <- MIR_measurements- cbind(NA, MIR_measurements)[, -(dim(MIR_measurements)[2]+1)]
X_test <- cbind(test[, 3580:3595], MIR_DER[,-1])
MIR_measurements <- test[, 2671:3579]
MIR_DER <- MIR_measurements- cbind(NA, MIR_measurements)[, -(dim(MIR_measurements)[2]+1)]
X_test <- cbind(X_test, MIR_DER[, -1])

train.base <- data.frame(matrix(ncol = 3579, nrow=1157))
#subtract baseline
for (i in 1:nrow(train))
{
  base.min[i] <- min(train[i,3500:3579])
  delta <- (base.min[i] - train[i,2671])/908
  
  for (j in 2671:3579)
  {
    train.base[i,j]=train[i,j]-delta*(j-2671)
    
  }
  print(i)
  
}
write.csv(train.base, file = "baseline.csv")
train.base <- read.csv("baseline.csv")

names(X_train[,1:16])
names(train.norm[,3550:3580])
names(train)[3600]

soil_property = "Sand"


# screening ---------------------------------------------------------------
trPerc = .7
soil_property = "Ca"
model <- earth(tr.der1[,-ncol(tr.der1)], tr.der1[,ncol(tr.der1)])
pred <- predict(model, ts.der1[,-ncol(ts.der1)])
re <- regr.eval(ts.der1[,n2], pred)
re


# medley ------------------------------------------------------------------


trPerc = .7

runlist = as.data.frame(rep(NA, 5))
for (i in 1:5)
{
  error <- NULL
  for(soil_property in soil_properties){
    idx <- sample(1:nrow(train),as.integer(trPerc*nrow(train)))
    sp.col <- which(names(train)==soil_property)
    
    
    train.der1 <- cbind(X_train1,train[ ,sp.col])
    n2 <- length(train.der1)
    colnames(train.der1)[n2] <- soil_property
    tr.der1 <- train.der1[idx,]
    ts.der1 <- train.der1[-idx,]
    
    med <- create.medley(tr.der1[,-ncol(tr.der1)], tr.der1[,ncol(tr.der1)], errfunc=rmse)
    med <- add.medley(med, svm, list(cost = 10, gamma=0.001))
    med <- add.medley(med, svm, list(cost = 100, gamma=0.001))
    med <- add.medley(med, svm, list(cost = 100, gamma=0.0001))
    med <- add.medley(med, svm, list(cost = 1000, gamma=0.0001))
    med <- add.medley(med, randomForest, list(ntree = 500, mtry=c(10,50))

    pred <- predict(med, ts.der1[,-ncol(ts.der1)])
    re <- regr.eval(ts.der1[,ncol(ts.der1)], pred)
    print(re)
    error <- rbind(error,re[2])

  }
    
  runlist <- cbind(runlist,error)
}

# svm test on datasets ----------------------------------------------------

trPerc = .7
runlist = as.data.frame(rep(NA, 40))

for (i in 1:10)
{
  error <- NULL
  for(soil_property in soil_properties){
    
    print(soil_property)
    idx <- sample(1:nrow(train),as.integer(trPerc*nrow(train)))
    sp.col <- which(names(train)==soil_property)
    
    train.norm <- cbind(train[, 3580:3595],train[, 2:2655],train[, 2671:3579],train[ ,sp.col])
    n1 <- length(train.norm)
    colnames(train.norm)[n1] <- soil_property
    tr.norm <- train.norm[idx,]
    ts.norm <- train.norm[-idx,]
    
    train.der <- cbind(X_train,train[ ,sp.col])
    n2 <- length(train.der)
    colnames(train.der)[n2] <- soil_property
    tr.der <- train.der[idx,]
    ts.der <- train.der[-idx,]
    
    train.der1 <- cbind(X_train1,train[ ,sp.col])
    n21 <- length(train.der1)
    colnames(train.der1)[n21] <- soil_property
    tr.der1 <- train.der1[idx,]
    ts.der1 <- train.der1[-idx,]
    
    train.bsln <- cbind(train[, 3580:3595], train.base[, 2671:3579],train[ ,sp.col])
    n5 <- length(train.bsln)
    colnames(train.bsln)[n5] <- soil_property
    tr.bsln <- train.bsln[idx,]
    ts.bsln <- train.bsln[-idx,]
    
    thresh <- 0.0001
    train.m <- as.matrix(X_train[,17:3577])
    train.m[abs(train.m)<thresh] = 0
    train.noise <- as.data.frame(cbind(X_train[,1:16],train.m,train[ ,sp.col]))
    n3 <- length(train.noise)
    colnames(train.noise)[n3] <- soil_property
    tr.noise <- train.noise[idx,]
    ts.noise <- train.noise[-idx,]
    
    thresh <- 0.0001
    train.m <- as.matrix(X_train1[,17:924])
    train.m[abs(train.m)<thresh] = 0
    train.noise1 <- as.data.frame(cbind(X_train1[,1:16],train.m,train[ ,sp.col]))
    n31 <- length(train.noise1)
    colnames(train.noise1)[n31] <- soil_property
    tr.noise1 <- train.noise1[idx,]
    ts.noise1 <- train.noise1[-idx,]
    
    thresh <- 0.001
    train.m <- as.matrix(X_train1[,17:924])
    train.m[abs(train.m)<thresh] = 0
    train.noise2 <- as.data.frame(cbind(X_train1[,1:16],train.m,train[ ,sp.col]))
    n32 <- length(train.noise2)
    colnames(train.noise2)[n32] <- soil_property
    tr.noise2 <- train.noise2[idx,]
    ts.noise2 <- train.noise2[-idx,]
    
    train.norm1 <- cbind(train[, 3580:3595],train[, 2671:3579],train[ ,sp.col])
    n4 <- length(train.norm1)
    colnames(train.norm1)[n4] <- soil_property
    tr.norm1 <- train.norm1[idx,]
    ts.norm1 <- train.norm1[-idx,]
    
    model.svm <- svm(tr.norm[,-ncol(tr.norm)], tr.norm[,ncol(tr.norm)], cost=100, gamma=0.0001)
    pred.svm <- predict(model.svm, ts.norm[,-ncol(ts.norm)])
    re <- regr.eval(ts.norm[,n1], pred.svm)
    print("norm")
    print(re)
    error <- rbind(error,re[2])
    
    model.svm <- svm(tr.bsln[,-ncol(tr.bsln)], tr.bsln[,ncol(tr.bsln)], cost=100, gamma=0.0001)
    pred.svm <- predict(model.svm, ts.bsln[,-ncol(ts.bsln)])
    re <- regr.eval(ts.bsln[,n5], pred.svm)
    print("bsln")
    print(re)
    error <- rbind(error,re[2])
    
    model.svm <- svm(tr.der[,-ncol(tr.der)], tr.der[,ncol(tr.der)], cost=100, gamma=0.0001)
    pred.svm <- predict(model.svm, ts.der[,-ncol(ts.der)])
    re <- regr.eval(ts.der[,n2], pred.svm)
    print("der")
    print(re)
    error <- rbind(error,re[2])
    
    model.svm <- svm(tr.der1[,-ncol(tr.der1)], tr.der1[,ncol(tr.der1)], cost=100, gamma=0.0001)
    pred.svm <- predict(model.svm, ts.der1[,-ncol(ts.der1)])
    re <- regr.eval(ts.der1[,n21], pred.svm)
    print("der1")
    print(re)
    error <- rbind(error,re[2])
    
    model.svm <- svm(tr.noise[,-ncol(tr.noise)], tr.noise[,ncol(tr.noise)], cost=100, gamma=0.0001)
    pred.svm <- predict(model.svm, ts.noise[,-ncol(ts.noise)])
    re <- regr.eval(ts.noise[,n3], pred.svm)
    print("noise")
    print(re)
    error <- rbind(error,re[2])
    
    model.svm <- svm(tr.noise1[,-ncol(tr.noise1)], tr.noise1[,ncol(tr.noise1)], cost=100, gamma=0.0001)
    pred.svm <- predict(model.svm, ts.noise1[,-ncol(ts.noise1)])
    re <- regr.eval(ts.noise1[,n31], pred.svm)
    print("noise1")
    print(re)
    error <- rbind(error,re[2])
    
    model.svm <- svm(tr.noise2[,-ncol(tr.noise2)], tr.noise2[,ncol(tr.noise2)], cost=100, gamma=0.0001)
    pred.svm <- predict(model.svm, ts.noise2[,-ncol(ts.noise2)])
    re <- regr.eval(ts.noise2[,n32], pred.svm)
    print("noise2")
    print(re)
    error <- rbind(error,re[2])
    
    model.svm <- svm(tr.norm1[,-ncol(tr.norm1)], tr.norm1[,ncol(tr.norm1)], cost=100, gamma=0.0001)
    pred.svm <- predict(model.svm, ts.norm1[,-ncol(ts.norm1)])
    re <- regr.eval(ts.norm1[,n4], pred.svm)
    print("norm1")
    print(re)
    error <- rbind(error,re[2])
  
  }
  
  runlist <- cbind(runlist,error)
  
}


trPerc = .7

runlist = as.data.frame(rep(NA, 10))
for (i in 1:5)
{
  error <- NULL
  for(soil_property in soil_properties){
    
    print(soil_property)
    idx <- sample(1:nrow(train),as.integer(trPerc*nrow(train)))
    sp.col <- which(names(train)==soil_property)
    
  
    train.der1 <- cbind(X_train1,train[ ,sp.col])
    n2 <- length(train.der1)
    colnames(train.der1)[n2] <- soil_property
    tr.der1 <- train.der1[idx,]
    ts.der1 <- train.der1[-idx,]
    
 

    soil_property = "Ca"    
    model <- earth(Ca~.,tr.der1)
    pred <- predict(model, ts.der1[,-ncol(ts.der1)])
    re <- regr.eval(ts.der1[,n2], pred)
    print("der1")
    print(re)
    error <- rbind(error,re[2])
    
    model.svm <- svm(tr.der2[,-ncol(tr.der2)], tr.der1[,ncol(tr.der2)], cost=100, gamma=0.0001)
    pred.svm <- predict(model.svm, ts.der2[,-ncol(ts.der2)])
    re <- regr.eval(ts.der2[,n3], pred.svm)
    print("der1")
    print(re)
    error <- rbind(error,re[2])
    
  }
  
  runlist <- cbind(runlist,error)
  
}

write.csv(runlist, file="errorlist2.csv")

run1 <- error
run2 <- error



# randomforest ------------------------------------------------------------


runlist = as.data.frame(rep(NA, 40))

for (i in 1:10)
{
  error <- NULL
  for(soil_property in soil_properties){
    
    print(soil_property)
    idx <- sample(1:nrow(train),as.integer(trPerc*nrow(train)))
    sp.col <- which(names(train)==soil_property)
    
    train.norm <- cbind(train[, 3580:3595],train[, 2:2655],train[, 2671:3579],train[ ,sp.col])
    n1 <- length(train.norm)
    colnames(train.norm)[n1] <- soil_property
    tr.norm <- train.norm[idx,]
    ts.norm <- train.norm[-idx,]
    
    train.der <- cbind(X_train,train[ ,sp.col])
    n2 <- length(train.der)
    colnames(train.der)[n2] <- soil_property
    tr.der <- train.der[idx,]
    ts.der <- train.der[-idx,]
    
    train.der1 <- cbind(X_train1,train[ ,sp.col])
    n21 <- length(train.der1)
    colnames(train.der1)[n21] <- soil_property
    tr.der1 <- train.der1[idx,]
    ts.der1 <- train.der1[-idx,]
    
    train.bsln <- cbind(train[, 3580:3595], train.base[, 2671:3579],train[ ,sp.col])
    n5 <- length(train.bsln)
    colnames(train.bsln)[n5] <- soil_property
    tr.bsln <- train.bsln[idx,]
    ts.bsln <- train.bsln[-idx,]
    
    thresh <- 0.0001
    train.m <- as.matrix(X_train[,17:3577])
    train.m[abs(train.m)<thresh] = 0
    train.noise <- as.data.frame(cbind(X_train[,1:16],train.m,train[ ,sp.col]))
    n3 <- length(train.noise)
    colnames(train.noise)[n3] <- soil_property
    tr.noise <- train.noise[idx,]
    ts.noise <- train.noise[-idx,]
    
    thresh <- 0.0001
    train.m <- as.matrix(X_train1[,17:924])
    train.m[abs(train.m)<thresh] = 0
    train.noise1 <- as.data.frame(cbind(X_train1[,1:16],train.m,train[ ,sp.col]))
    n31 <- length(train.noise1)
    colnames(train.noise1)[n31] <- soil_property
    tr.noise1 <- train.noise1[idx,]
    ts.noise1 <- train.noise1[-idx,]
    
    thresh <- 0.001
    train.m <- as.matrix(X_train1[,17:924])
    train.m[abs(train.m)<thresh] = 0
    train.noise2 <- as.data.frame(cbind(X_train1[,1:16],train.m,train[ ,sp.col]))
    n32 <- length(train.noise2)
    colnames(train.noise2)[n32] <- soil_property
    tr.noise2 <- train.noise2[idx,]
    ts.noise2 <- train.noise2[-idx,]
    
    train.norm1 <- cbind(train[, 3580:3595],train[, 2671:3579],train[ ,sp.col])
    n4 <- length(train.norm1)
    colnames(train.norm1)[n4] <- soil_property
    tr.norm1 <- train.norm1[idx,]
    ts.norm1 <- train.norm1[-idx,]
    
    model.rf <- randomForest(tr.norm[,-ncol(tr.norm)], tr.norm[,ncol(tr.norm)])
    pred.rf <- predict(model.rf, ts.norm[,-ncol(ts.norm)])
    re <- regr.eval(ts.norm[,n1], pred.rf)
    print("norm")
    print(re)
    error <- rbind(error,re[2])
    
    model.rf <- randomForest(tr.bsln[,-ncol(tr.bsln)], tr.bsln[,ncol(tr.bsln)])
    pred.rf <- predict(model.rf, ts.bsln[,-ncol(ts.bsln)])
    re <- regr.eval(ts.bsln[,n5], pred.rf)
    print("bsln")
    print(re)
    error <- rbind(error,re[2])
    
    model.rf <- randomForest(tr.der[,-ncol(tr.der)], tr.der[,ncol(tr.der)])
    pred.rf <- predict(model.rf, ts.der[,-ncol(ts.der)])
    re <- regr.eval(ts.der[,n2], pred.rf)
    print("der")
    print(re)
    error <- rbind(error,re[2])
    
    model.rf <- randomForest(tr.der1[,-ncol(tr.der1)], tr.der1[,ncol(tr.der1)])
    pred.rf <- predict(model.rf, ts.der1[,-ncol(ts.der1)])
    re <- regr.eval(ts.der1[,n21], pred.rf)
    print("der1")
    print(re)
    error <- rbind(error,re[2])
    
    model.rf <- randomForest(tr.noise[,-ncol(tr.noise)], tr.noise[,ncol(tr.noise)])
    pred.rf <- predict(model.rf, ts.noise[,-ncol(ts.noise)])
    re <- regr.eval(ts.noise[,n3], pred.rf)
    print("noise")
    print(re)
    error <- rbind(error,re[2])
    
    model.rf <- randomForest(tr.noise1[,-ncol(tr.noise1)], tr.noise1[,ncol(tr.noise1)])
    pred.rf <- predict(model.rf, ts.noise1[,-ncol(ts.noise1)])
    re <- regr.eval(ts.noise1[,n31], pred.rf)
    print("noise1")
    print(re)
    error <- rbind(error,re[2])
    
    model.rf <- randomForest(tr.noise2[,-ncol(tr.noise2)], tr.noise2[,ncol(tr.noise2)])
    pred.rf <- predict(model.rf, ts.noise2[,-ncol(ts.noise2)])
    re <- regr.eval(ts.noise2[,n32], pred.rf)
    print("noise2")
    print(re)
    error <- rbind(error,re[2])
    
    model.rf <- randomForest(tr.norm1[,-ncol(tr.norm1)], tr.norm1[,ncol(tr.norm1)])
    pred.rf <- predict(model.rf, ts.norm1[,-ncol(ts.norm1)])
    re <- regr.eval(ts.norm1[,n4], pred.rf)
    print("norm1")
    print(re)
    error <- rbind(error,re[2])
    
  }
  
  runlist <- cbind(runlist,error)
  
}

runlist = as.data.frame(rep(NA, 20))

for (i in 1:5)
{
  error <- NULL
  for(soil_property in soil_properties){
    
    print(soil_property)
    idx <- sample(1:nrow(train),as.integer(trPerc*nrow(train)))
    sp.col <- which(names(train)==soil_property)
    
    train.der1 <- cbind(X_train1,train[ ,sp.col])
    n21 <- length(train.der1)
    colnames(train.der1)[n21] <- soil_property
    tr.der1 <- train.der1[idx,]
    ts.der1 <- train.der1[-idx,]
    
    model.rf <- randomForest(tr.der1[,-ncol(tr.der1)], tr.der1[,ncol(tr.der1)], ntree=500, mtry=10)
    pred.rf <- predict(model.rf, ts.der1[,-ncol(ts.der1)])
    re <- regr.eval(ts.der1[,n21], pred.rf)
    print("der1")
    print(re)
    error <- rbind(error,re[2])
    
    model.rf <- randomForest(tr.der1[,-ncol(tr.der1)], tr.der1[,ncol(tr.der1)], ntree=500, mtry=25)
    pred.rf <- predict(model.rf, ts.der1[,-ncol(ts.der1)])
    re <- regr.eval(ts.der1[,n21], pred.rf)
    print("der1")
    print(re)
    error <- rbind(error,re[2])
    
    model.rf <- randomForest(tr.der1[,-ncol(tr.der1)], tr.der1[,ncol(tr.der1)],ntree=200, mtry=100)
    pred.rf <- predict(model.rf, ts.der1[,-ncol(ts.der1)])
    re <- regr.eval(ts.der1[,n21], pred.rf)
    print("der1")
    print(re)
    error <- rbind(error,re[2])
    
    model.rf <- randomForest(tr.der1[,-ncol(tr.der1)], tr.der1[,ncol(tr.der1)], ntree=500, mtry=50)
    pred.rf <- predict(model.rf, ts.der1[,-ncol(ts.der1)])
    re <- regr.eval(ts.der1[,n21], pred.rf)
    print("der1")
    print(re)
    error <- rbind(error,re[2])
    
  }
  
  runlist <- cbind(runlist,error)
  
}

write.csv(runlist, file="errorlist4.csv")
