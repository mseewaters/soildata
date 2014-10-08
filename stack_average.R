library(e1071)
library(randomForest)
require(medley)
library(ridge)
library(performanceEstimation)
library(DMwR)
library(earth)


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

# stacking ---------------------------------------------------------------


runlist = as.data.frame(rep(NA, 25))

for (i in 1:10)
{
  error <- NULL
  for(soil_property in soil_properties){
    
    print(soil_property)
    idx <- sample(1:nrow(train),as.integer(trPerc*nrow(train)))
    sp.col <- which(names(train)==soil_property)
    
    wgts <- switch(soil_property,
                   "Ca" = c(0.5,0.2,0.3),
                   "P" = c(0.2,0.6,0.2),
                   "pH" = c(0.1,0.4,0.5),
                   "SOC" = c(0.5,0.2,0.3),
                   "Sand" = c(0.2,0.5,0.3))
    
    train.der1 <- cbind(X_train1,train[ ,sp.col])
    n21 <- length(train.der1)
    colnames(train.der1)[n21] <- soil_property
    tr.der1 <- train.der1[idx,]
    ts.der1 <- train.der1[-idx,]
    
    model.rf <- randomForest(tr.der1[,-ncol(tr.der1)], tr.der1[,ncol(tr.der1)], ntree=500, mtry=10)
    pred.rf1 <- predict(model.rf, ts.der1[,-ncol(ts.der1)])
    pred.w.rf1 <- pred.rf1 * wgts[1]
    pred.tr.rf1 <- predict(model.rf, tr.der1[,-ncol(tr.der1)])
    re <- regr.eval(ts.der1[,ncol(ts.der1)], pred.rf1)
    print(re)
    error <- rbind(error,re[2])

    

    model.svm <- svm(tr.der1[,-ncol(tr.der1)], tr.der1[,ncol(tr.der1)], cost=100, gamma=0.001)
    pred.svm2 <- predict(model.svm, ts.der1[,-ncol(ts.der1)])
    pred.w.svm2 <- pred.svm2 * wgts[2]
    pred.tr.svm2 <- predict(model.svm, tr.der1[,-ncol(tr.der1)])
    re <- regr.eval(ts.der1[,ncol(ts.der1)], pred.svm2)
    print(re)
    error <- rbind(error,re[2])
    
    model.svm <- svm(tr.der1[,-ncol(tr.der1)], tr.der1[,ncol(tr.der1)], cost=1000, gamma=0.0001)
    pred.svm4 <- predict(model.svm, ts.der1[,-ncol(ts.der1)])
    pred.w.svm4 <- pred.svm4 * wgts[3]
    pred.tr.svm4 <- predict(model.svm, tr.der1[,-ncol(tr.der1)])
    re <- regr.eval(ts.der1[,ncol(ts.der1)], pred.svm4)
    print(re)
    error <- rbind(error,re[2])
      
    allpred.tr <- cbind(pred.tr.rf1, pred.tr.svm2, pred.tr.svm4)
    allpred <- cbind(pred.w.rf1, pred.w.svm2, pred.w.svm4)
    
    model.stk.svm <- earth(allpred.tr, tr.der1[,ncol(tr.der1)])
    pred.stk.svm <- predict(model.stk.svm, allpred)
    
    
    avgpred <- rowSums(allpred) 
    re <- regr.eval(ts.der1[,ncol(ts.der1)], avgpred)
    re2 <- regr.eval(ts.der1[,ncol(ts.der1)], pred.stk.svm)
    print("average result")
    print(re)
    print("stack result")
    print(re2)
    error <- rbind(error,re[2])
    error <- rbind(error,re2[2])
  }
  
  runlist <- cbind(runlist,error)
  
}

write.csv(runlist, file="error4.csv")
