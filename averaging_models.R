library(e1071)
library(randomForest)
library(rpart)
library(performanceEstimation)
library(DMwR)
library(psych)
require(medley)
library(earth)
library(quantreg)
library(ridge)
library(lars)
library(glmnet)


train <- read.csv("training.csv")
#test <- read.csv("sorted_test.csv")
soil_properties <- c("Ca", "P", "pH", "SOC", "Sand")

train$Depth <- as.numeric(train$Depth)

# training data
MIR_measurements <- train[, 2671:3579]
MIR_DER <- MIR_measurements- cbind(NA, MIR_measurements)[, -(dim(MIR_measurements)[2]+1)]
X_train1 <- cbind(train[, 3580:3595], MIR_DER[, -1])


trPerc = .7
soil_property = "Ca"
idx <- sample(1:nrow(train),as.integer(trPerc*nrow(train)))
sp.col <- which(names(train)==soil_property)

train.der1 <- cbind(X_train1,train[ ,sp.col])
n2 <- length(train.der1)
colnames(train.der1)[n2] <- soil_property
tr.der1 <- train.der1[idx,]
ts.der1 <- train.der1[-idx,]
mx <- as.matrix(tr.der1[,-ncol(tr.der1)])
tx <- as.matrix(ts.der1[,-ncol(ts.der1)])

# screening ---------------------------------------------------------------

model <- glmnet(mx, tr.der1[,ncol(tr.der1)])
pred <- predict(model, tx)
re <- regr.eval(ts.der1[,n2], pred)
re

model <- linearRidge(Ca~.,tr.der1, lambda=3)
pred <- predict(model, ts.der1[,-ncol(ts.der1)])
re <- regr.eval(ts.der1[,n2], pred)
re


runlist = as.data.frame(rep(NA, 50))

for (i in 1:10)
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
    pred.rf1 <- predict(model.rf, ts.der1[,-ncol(ts.der1)])
    re <- regr.eval(ts.der1[,ncol(ts.der1)], pred.rf1)
    print(re)
    error <- rbind(error,re[2])
 
    model.rf <- randomForest(tr.der1[,-ncol(tr.der1)], tr.der1[,ncol(tr.der1)], ntree=500, mtry=50)
    pred.rf2 <- predict(model.rf, ts.der1[,-ncol(ts.der1)])
    re <- regr.eval(ts.der1[,ncol(ts.der1)], pred.rf2)
    print(re)
    error <- rbind(error,re[2])
    
    model.svm <- svm(tr.der1[,-ncol(tr.der1)], tr.der1[,ncol(tr.der1)], cost=10, gamma=0.001)
    pred.svm1 <- predict(model.svm, ts.der1[,-ncol(ts.der1)])
    re <- regr.eval(ts.der1[,ncol(ts.der1)], pred.svm1)
    print(re)
    error <- rbind(error,re[2])
    
    model.svm <- svm(tr.der1[,-ncol(tr.der1)], tr.der1[,ncol(tr.der1)], cost=100, gamma=0.001)
    pred.svm2 <- predict(model.svm, ts.der1[,-ncol(ts.der1)])
    re <- regr.eval(ts.der1[,ncol(ts.der1)], pred.svm2)
    print(re)
    error <- rbind(error,re[2])
    
    model.svm <- svm(tr.der1[,-ncol(tr.der1)], tr.der1[,ncol(tr.der1)], cost=100, gamma=0.0001)
    pred.svm3 <- predict(model.svm, ts.der1[,-ncol(ts.der1)])
    re <- regr.eval(ts.der1[,ncol(ts.der1)], pred.svm3)
    print(re)
    error <- rbind(error,re[2])
    
    model.svm <- svm(tr.der1[,-ncol(tr.der1)], tr.der1[,ncol(tr.der1)], cost=1000, gamma=0.0001)
    pred.svm4 <- predict(model.svm, ts.der1[,-ncol(ts.der1)])
    re <- regr.eval(ts.der1[,ncol(ts.der1)], pred.svm4)
    print(re)
    error <- rbind(error,re[2])
    
    formula <- paste0(soil_property," ~ .")
    model.rr <- linearRidge(formula, tr.der1, lambda = 2) 
    pred.rr1 <- predict(model.rr, ts.der1[,-ncol(ts.der1)])
    re <- regr.eval(ts.der1[,ncol(ts.der1)], pred.rr1)
    print(re)
    error <- rbind(error,re[2])
    
    model.rr <- linearRidge(formula, tr.der1, lambda = 3)
    pred.rr2 <- predict(model.rr, ts.der1[,-ncol(ts.der1)])
    re <- regr.eval(ts.der1[,ncol(ts.der1)], pred.rr2)
    print(re)
    error <- rbind(error,re[2])
    
    model.rr <- linearRidge(formula, tr.der1, lambda = 4)
    pred.rr3 <- predict(model.rr, ts.der1[,-ncol(ts.der1)])
    re <- regr.eval(ts.der1[,ncol(ts.der1)], pred.rr2)
    print(re)
    error <- rbind(error,re[2])
    
    allpred <- cbind(pred.rf1, pred.rf2, pred.svm1, pred.svm2, pred.svm3, pred.svm4, pred.rr1, pred.rr2, pred.rr3)
    avgpred <- rowMeans(allpred) 
    re <- regr.eval(ts.der1[,ncol(ts.der1)], avgpred)
    print(re)
    error <- rbind(error,re[2])
  }
  
  runlist <- cbind(runlist,error)
  
}

write.csv(runlist, file="error1.csv")

cor(allpred)
