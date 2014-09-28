library(e1071)
library(randomForest)
library(rpart)
library(performanceEstimation)
library(DMwR)
library(psych)

train <- read.csv("training.csv")
#test <- read.csv("sorted_test.csv")
soil_properties <- c("Ca", "P", "pH", "SOC", "Sand")

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

# testing data
MIR_measurements <- test[, 2:2655]
MIR_DER <- MIR_measurements- cbind(NA, MIR_measurements)[, -(dim(MIR_measurements)[2]+1)]
X_test <- cbind(test[, 3580:3595], MIR_DER[,-1])
MIR_measurements <- test[, 2671:3579]
MIR_DER <- MIR_measurements- cbind(NA, MIR_measurements)[, -(dim(MIR_measurements)[2]+1)]
X_test <- cbind(X_test, MIR_DER[, -1])

names(X_train[,1:16])
names(train.norm[,3550:3580])

trPerc = .7
idx <- sample(1:nrow(train),as.integer(trPerc*nrow(train)))

train.norm <- cbind(train[, 3580:3595],train[, 2:2655],train[, 2671:3579],train$Ca)
colnames(train.norm)[3580] <- "Ca"
tr.norm <- train.norm[idx,]
ts.norm <- train.norm[-idx,]

train.der <- cbind(X_train,train$Ca)
colnames(train.der)[3578] <- "Ca"
tr.der <- train.der[idx,]
ts.der <- train.der[-idx,]

thresh <- 0.001
train.m <- as.matrix(X_train[,17:3577])
train.m[abs(train.m)<thresh] = 0
train.noise <- as.data.frame(cbind(X_train[,1:16],train.m,train$Ca))
colnames(train.noise)[3578] <- "Ca"
tr.noise <- train.noise[idx,]
ts.noise <- train.noise[-idx,]


x.factor <- cbind(train[, 2:2655],train[, 2671:3579])
x.factor <- train[,2:600]
FA <- factanal(x.factor, 3, rotation="varimax", scores="regression")
print(FA)

x.factor <- train[,500:1100]
FA2 <- factanal(x.factor, 3, rotation="varimax", scores="regression")
x.factor <- train[,1000:1500]
FA3 <- factanal(x.factor, 3, rotation="varimax", scores="regression")
x.factor <- train[,1500:1800]
FA4 <- factanal(x.factor, 3, rotation="varimax", scores="regression")
x.factor <- train[,1800:1900]
FA42 <- factanal(x.factor, 3, rotation="varimax", scores="regression")
x.factor <- train[,1900:1950]
FA43 <- factanal(x.factor, 3, rotation="varimax", scores="regression")

x.factor <- train[,1950:2300]
FA5 <- factanal(x.factor, 3, rotation="varimax", scores="regression")
x.factor <- train[,2300:2655]
FA6 <- factanal(x.factor, 3, rotation="varimax", scores="regression")
x.factor <- train[,2671:2900]
FA7 <- factanal(x.factor, 3, rotation="varimax", scores="regression")
x.factor <- train[,2900:3000]
FA8 <- factanal(x.factor, 3, rotation="varimax", scores="regression")
x.factor <- train[,3200:3400]
FA9 <- factanal(x.factor, 3, rotation="varimax", scores="regression") 
x.factor <- train[,3300:3579]
FA9 <- factanal(x.factor, 3, rotation="varimax", scores="regression") 

# Determine Number of Factors to Extract
library(nFactors)
ev <- eigen(cor(x.factor)) # get eigenvalues
ap <- parallel(subject=nrow(x.factor),var=ncol(x.factor),
               rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS)


model.svm <- svm(Ca ~ .,tr.noise)
pred.svm <- predict(model.svm, ts.noise)
res.svm <- cbind(ts.noise$Ca, pred.svm)
regr.eval(ts.noise$Ca, pred.svm)

model.svm <- svm(Ca ~ .,tr.norm)
pred.svm <- predict(model.svm, ts.norm)
res.svm2 <- cbind(ts.norm$Ca, pred.svm)
regr.eval(ts.norm$Ca, pred.svm)

model.svm <- svm(Ca ~ .,tr.der)
pred.svm <- predict(model.svm, ts.der)
res.svm3 <- cbind(ts.der$Ca, pred.svm)
regr.eval(ts.der$Ca, pred.svm)


model.tree <- rpartXse(Ca ~ .,tr.noise)
pred.tree <- predict(model.tree, ts.noise)
res.tree <- cbind(ts.noise$Ca, pred.tree)
regr.eval(ts.noise$Ca, pred.tree)

model.tree <- rpartXse(Ca ~ .,tr.norm)
pred.tree <- predict(model.tree, ts.norm)
res.tree2 <- cbind(ts.norm$Ca, pred.tree)
regr.eval(ts.norm$Ca, pred.tree)

model.tree <- rpartXse(Ca ~ .,tr.der)
pred.tree <- predict(model.tree, ts.der)
res.tree3 <- cbind(ts.der$Ca, pred.tree)
regr.eval(ts.der$Ca, pred.tree)


# Evaluation
res <- performanceEstimation(
  c(PredTask(Ca ~ ., train.norm),PredTask(Ca ~ ., train.der),
    PredTask(Ca ~ ., train.noise)),
  c(workflowVariants("standardWF", learner = "svm",
                     learner.pars=list(cost=c(1,10,100), gamma=c(0.1,0.01)))),
  HldSettings(nReps=5,hldSz=0.2))
  

CvSettings(nReps =1, nFolds = 10))
BootSettings(type=".632", nReps=1)

plot(res)

