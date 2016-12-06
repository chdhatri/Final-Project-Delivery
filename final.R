#############################CMPE 239 PROJECT#####################
#TEAM VOLTAGE
#Sucharu Gupta
#Dhatri CHennavajula
#Andrew Wong
#Sprush Ujjwal


#UNCOMMENT THESE WHEN YOU ARE RUNNING FIRST TIME TO INSTALL PACKAGES
#install.packages("rpart")
#install.packages("rpart.plot")
#install.packages("rpart")
#install.packages("rpart.plot")
#install.packages("doSNOW")
#install.packages("foreign")
#install.packages("nnet")
#install.packages("ggplot2")
#install.packages("reshape2")
#install.packages("caret")
#install.packages("SDMTools")
#install.packages("vcd")
#install.packages("MASS")
#install.packages("ROCR")
#install.packages("e1071")
#install.packages("randomForest")
#install.packages("rpart")
#install.packages("rpart.plot")


library(rpart)
library(rpart.plot)
library(rpart)
library(rpart.plot)
library(doSNOW)
library(foreign)
library(nnet)
library(ggplot2)
library(reshape2)
library(caret)
library(SDMTools)
library(vcd)
library(MASS)
library(ROCR)
library(e1071)
library(randomForest)
library(rpart)
library(rpart.plot)

gene <- read.delim("GENETOX_Bacterial mutagenicity_NTP.txt")

# Check the datapatterns
plot(gene$STUDY_CONCLUSION, col=c("blue","green","yellow","pink"))
gene$STUDY_CONCLUSION[which(gene$STUDY_CONCLUSION == "Weakly Positive")] = "Equivocal"
plot(droplevels(gene$STUDY_CONCLUSION), col=c("blue","green","yellow"))

# Count has NA's , for all the rows that have NA are replace with COUNT_MEAN value
gene$newCount <- ifelse(is.na(gene$COUNT), gene$COUNT_MEAN, gene$COUNT)
hist(gene$DOSE, breaks=100, xlab="DOSE", main="Histogram of DOSE")
hist(gene$newCount, breaks=50, xlab="COUNT", main="Histogram of COUNT")

# check to how dose and count are related


dose <- gene$DOSE
count <-  gene$newCount
# Check the corelation between dose and count
cor(dose,count, use="complete")
#[1] -0.1003457

par(mfrow=c(1,1))
plot(dose,count, xlab="DOSE",ylab="COUNT",main="DOSE vs COUNT")

logDose <- 1/log(dose *count)
logDosePerCount <- 1/log(dose / count)

gene$logDose <- logDose  
gene$logDosePerCount <-logDosePerCount  

cor(dose,logDose, use="complete")
plot(dose,logDose, xlab="DOSE",ylab="COUNT",main="DOSEPERCOUNT vs COUNT")
# There are 74 entries with no COUNT, COUNT_SEM, COUNT_MEAN
# Not sure if we can remove these rows as there count is not there
# for now removing these rows assuming they provide no information without COUNT, but need to look into it again
# remove rest of the rows with no count
gene <- gene[-which(is.na(count)),]
gene <- gene[-which(gene$TREATMENT_GROUP_TYPE == ""),]
gene <- gene[-which(gene$TRIAL_RESULT == "Not a Valid Test"),]
gene <- gene[-which(gene$TRIAL_RESULT == "Failed Experiment"),]
gene <- gene[-which(gene$TRIAL_RESULT == "No Call"),]

#[1] 5186   26

# Add new column if study conclusion is Poistive , its Mutage
# if its Equivocal or Negative its Non mutagen
gene$RESULT <- ifelse(gene$STUDY_CONCLUSION == "Positive", 1,0)

# convert RESULT to factor variable
gene$RESULT <- as.factor(gene$RESULT)


hist(gene$logDose, breaks=50, xlab="LOG DOSE", main="Histogram of LogDose")

par(mfrow=c(1,1))
plot(gene$STRAIN, col=c(gene$STRAIN),
     main = "STAIN")
plot(gene$TRIAL_RESULT, col=gene$TRIAL_RESULT,
     main = "TRIAL_RESULT")
plot(gene$MICROSOMAL_ACTIVATION_USED, col=gene$MICROSOMAL_ACTIVATION_USED,
     main = "MICROSOMAL_ACTIVATION_USED")
# data cleaning and preprocessing
# Convert chemical name, study title, subject name,accession number to charachers
# start date to date
gene$CHEMICAL_NAME <- as.character(gene$CHEMICAL_NAME)
gene$STUDY_TITLE <- as.character(gene$STUDY_TITLE)
gene$START_DATE <- as.Date(gene$START_DATE)
gene$SUBJECT_NAME <- as.character(gene$SUBJECT_NAME)
gene$ACCESSION_NUMBER <- as.character(gene$ACCESSION_NUMBER)
gene$DEPOSITOR_STUDY_NUMBER <- as.character(gene$DEPOSITOR_STUDY_NUMBER)
gene$ORGANIZATION_NAME <- as.character(gene$ORGANIZATION_NAME)

str(gene)
# create two new variables to classify as (non)mutagen 
mutagen <- 1
nonmutagen <- 0


# take subset of data only for year 2015
gene.subset <- droplevels(subset(gene,  format(as.Date(gene$START_DATE),"%y")== "15"  ))
dim(gene.subset)
#[1] 5260   25




# Split the data to train and test
set.seed(1234)
gene.holdout.variable <- rbinom(n = dim(gene.subset)[1], size = 1, prob = .2)

# training data
training.data <- gene.subset[gene.holdout.variable==0, ]
dim(training.data)
# test data
test.data <- gene.subset[gene.holdout.variable==1, ]
dim(test.data)

#further split train data into train and validation data
train.holdout.variable <- rbinom(n = dim(training.data)[1], size = 1, prob = .2)
train.holdout.variable

# final train data
gene.train <- training.data[train.holdout.variable==0, ]
unique(gene.train$TRIAL_RESULT)
unique(gene.train$MICROSOMAL_ACTIVATION_USED)

# write it to a file
write.csv(gene.train,"gene.train.csv")
gene.valid <- training.data[train.holdout.variable==1, ]
write.csv(gene.valid,"gene.valid.csv")
write.csv(test.data, "test.data.csv")
str(gene.subset)



#####################################################################
model8 <- glm(RESULT ~ DOSE +TRIAL_RESULT +STRAIN+ MICROSOMAL_ACTIVATION_USED,family=binomial(link='logit'),data=gene.train)
model9 <- glm(RESULT ~ DOSE + newCount +STRAIN+ TRIAL_RESULT ,family=binomial(link='logit'),data=gene.train)
model10 <- glm(RESULT ~ logDose  +DOSE+newCount+STRAIN+TRIAL_RESULT ,family=binomial(link='logit'),data=gene.train)

# Check the summary of each model #####

summary(model8) #AIC:  1371.8,  1353.8 Residual
summary(model9) #AIC:  1354,  1334.9 Residua
summary(model10) #IC: 1363.8

sd_ols <- function(object) sqrt(sum(residuals(object)^2)/df.residual(object))
sd_mle <- function(object) sqrt(mean(residuals(object)^2))


sd_ols(model8) #0.6368555
sd_ols(model9) #0.6419228
sd_ols(model10) #0.6408877

sd_mle(model8) #0.635994
sd_mle(model9) #0.6410545
sd_mle(model10) #0.6399244
# ANOVA

anova(model8, test="Chisq") 
anova(model9, test="Chisq") 
anova(model10, test="Chisq")

anova(model8, model9, test="LRT")
anova(model8, model10, test="LRT")


predict(model8,newdata=gene.valid)
fitted.results8 <- predict(model8,newdata=gene.valid,type='response')
final.fitted.results8 <- ifelse(fitted.results8 > 0.5,1,0)
misClasificError8 <- mean(final.fitted.results8 != gene.valid$RESULT)
print(paste('Accuracy',1-misClasificError8)) #0.945409429280397

predict(model9,newdata=gene.valid)
fitted.results9 <- predict(model9,newdata=gene.valid,type='response')
final.fitted.results9 <- ifelse(fitted.results9 > 0.5,1,0)
misClasificError9 <- mean(final.fitted.results9 != gene.valid$RESULT)
print(paste('Accuracy',1-misClasificError9)) #0.946650124069479


predict(model10,newdata=gene.valid)
fitted.results10 <- predict(model10,newdata=gene.valid,type='response')
final.fitted.results10 <- ifelse(fitted.results10 > 0.5,1,0)
misClasificError10 <- mean(final.fitted.results10 != gene.valid$RESULT)
print(paste('Accuracy',1-misClasificError10)) #0.946650124069479


test.results9 <- predict(model9,newdata=test.data,type='response')
test.fitted.results9 <- ifelse(test.results9 > 0.5,1,0)

# for model 9 find the accuracy , model is 94% accurate
test_misClasificError9 <- mean(test.fitted.results9 != test.data$RESULT)
print(paste('Accuracy',1-test_misClasificError9)) #0.944656488549618


test.results8 <- predict(model8,newdata=test.data,type='response')
test.fitted.results8 <- ifelse(test.results8 > 0.5,1,0)

# for model 9 find the accuracy , model is 94% accurate
test_misClasificError8 <- mean(test.fitted.results8 != test.data$RESULT)
print(paste('Accuracy',1-test_misClasificError8))

glm.p1 <-  predict(model8,newdata=test.data,type='response')
glm.p2 <-  predict(model9,newdata=test.data,type='response')
glm.p3 <-  predict(model10,newdata=test.data,type='response')

glm.pred1 <- prediction(glm.p1, test.data$RESULT)
glm.pred2 <- prediction(glm.p2, test.data$RESULT)
glm.pred3 <- prediction(glm.p3, test.data$RESULT)

glm.prf1 <- performance(glm.pred1, measure = "tpr", x.measure = "fpr")
glm.prf2 <- performance(glm.pred2, measure = "tpr", x.measure = "fpr")
glm.prf3 <- performance(glm.pred3, measure = "tpr", x.measure = "fpr")

par(mfrow=c(1,1))
plot(glm.prf1, col="blue")
plot(glm.prf2, col="red", add=T)
plot(glm.prf3, col="black", add=T)

abline(a=0, b= 1)

#### MODEL1 is best in GLM##############
####PR8 is good model####

rf1 <- randomForest(RESULT ~ DOSE +TRIAL_RESULT +STRAIN+ MICROSOMAL_ACTIVATION_USED, 
                        data = gene.train, importance=T, ntree=500)
rf1
varImpPlot(rf1)

rf2 <- randomForest(RESULT ~ DOSE + newCount +STRAIN+ TRIAL_RESULT, 
                    data = gene.train, importance=T, ntree=500)
rf2
varImpPlot(rf2)


rf3 <- randomForest(RESULT ~ logDose  +DOSE+newCount+STRAIN+TRIAL_RESULT, 
                    data = gene.train, importance=T, ntree=500)
rf3
varImpPlot(rf3)

###### RF with test data and ROC curve#####

p1 <-  predict(rf1,newdata=test.data)
table(pred=p1,true=test.data$RESULT)

p2 <-  predict(rf2,newdata=test.data)
table(pred=p2,true=test.data$RESULT)

p3 <-  predict(rf3,newdata=test.data)
table(pred=p3,true=test.data$RESULT)

p1 <-  predict(rf1,newdata=test.data,type='prob')
p2 <-  predict(rf2,newdata=test.data,type='prob')
p3 <-  predict(rf3,newdata=test.data,type='prob')

pred1 <- prediction(p1[,2], test.data$RESULT)
pred2 <- prediction(p2[,2], test.data$RESULT)
pred3 <- prediction(p3[,2], test.data$RESULT)

prf1 <- performance(pred1, measure = "tpr", x.measure = "fpr")
prf2 <- performance(pred2, measure = "tpr", x.measure = "fpr")
prf3 <- performance(pred3, measure = "tpr", x.measure = "fpr")


plot(prf1, col="blue")
plot(prf2, col="red", add=T)
plot(prf3, col="black", add=T)

legend(0.8,0.4,c("model1","model2","model3"), lty=c(1,1), col=c("blue","red","black"))
mtext("Random Forest", side=3)

# model 2 is good in RV
############### ROC CURVE on GLM vs RV #############


plot(glm.prf1, col="blue")
plot(prf2, col="red", add=T)

legend(0.8,0.4,c("glm","rf"), lty=c(1,1), col=c("blue","red"))
mtext("Random Forest vs GLM", side=3)
abline(0,1)

glm.prf1 <- performance(glm.pred1, measure = "auc")
glm.auc <- glm.prf1@y.values[[1]]
glm.auc
#[1] 0.8232026

prf2 <- performance(pred2, measure = "auc")
rf.auc <- prf2@y.values[[1]]
rf.auc
#[1] 0.8458383


opt.cut = function(perf, pred){
  cut.ind = mapply(FUN=function(x, y, p){
    d = (x - 0)^2 + (y-1)^2
    ind = which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}

#### RPART TO PLOT RF###
# Create utility function

set.seed(37596)
cv.3.folds <- createMultiFolds(rf.label, k = 3, times = 10)

ctrl.3 <- trainControl(method = "repeatedcv", number = 3, repeats = 10,
                       index = cv.3.folds)
rpart.cv <- function(seed, training, labels, ctrl) {
  cl <- makeCluster(6, type = "SOCK")
  registerDoSNOW(cl)
  
  set.seed(seed)
  # Leverage formula interface for training
  rpart.cv <- train(x = training, y = labels, method = "rpart", tuneLength = 30, 
                    trControl = ctrl)
  
  #Shutdown cluster
  stopCluster(cl)
  
  return (rpart.cv)
}
rf.label <- gene.train$RESULT
# Grab features
features <- c("DOSE", "newCount", "STRAIN","TRIAL_RESULT")
rpart.train.1 <- gene.train[, features]

# Run CV and check out results
rpart.1.cv.1 <- rpart.cv(94622, rpart.train.1, rf.label, ctrl.3)
rpart.1.cv.1

# Plot
prp(rpart.1.cv.1$finalModel, type = 0, extra = 1, under = TRUE)

############################## SECOND LEVEL ##########################

##### MULTINOMIAL REGRESSON#######
mmodel1 <- multinom(STUDY_CONCLUSION ~ DOSE + newCount +STRAIN+ TRIAL_RESULT+RESULT,data=gene.train)
summary(mmodel1)
z1 <- summary(mmodel1)$coefficients/summary(mmodel1)$standard.errors
p1 <- (1 - pnorm(abs(z1), 0, 1))*2
p1 

mmodel2 <- multinom(STUDY_CONCLUSION ~ DOSE  +STRAIN+ TRIAL_RESULT+RESULT,data=gene.train)
summary(mmodel2)
z2 <- summary(mmodel2)$coefficients/summary(mmodel2)$standard.errors
p2 <- (1 - pnorm(abs(z2), 0, 1))*2
p2 

mmodel3 <- multinom(STUDY_CONCLUSION ~ logDose  +DOSE+newCount+STRAIN+TRIAL_RESULT+RESULT,data=gene.train)
summary(mmodel3)
z3 <- summary(mmodel3)$coefficients/summary(mmodel3)$standard.errors
p3 <- (1 - pnorm(abs(z3), 0, 1))*2
p3

mmodel4 <- multinom(STUDY_CONCLUSION ~ logDose  +DOSE+newCount+STRAIN+RESULT,data=gene.train)
summary(mmodel4)
z4 <- summary(mmodel4)$coefficients/summary(mmodel4)$standard.errors
p4 <- (1 - pnorm(abs(z4), 0, 1))*2
p4

###### multinomial with validation data#####

model1.predict <- predict(mmodel1,newdata=gene.valid)
model2.predict <- predict(mmodel2,newdata=gene.valid)
model3.predict <- predict(mmodel3,newdata=gene.valid)
model4.predict <- predict(mmodel4,newdata=gene.valid)

#### confusion matrix####
cm1 <- table(pred=model1.predict, true=gene.valid$STUDY_CONCLUSION)
print(cm1)

cm2 <- table(pred=model2.predict, true=gene.valid$STUDY_CONCLUSION)
print(cm2)

cm3 <- table(pred=model3.predict, true=gene.valid$STUDY_CONCLUSION)
print(cm3)

cm4 <- table(pred=model4.predict, true=gene.valid$STUDY_CONCLUSION)
print(cm4)


###### multinomial with test data#####
tmodel1.predict <- predict(mmodel1,newdata=test.data)
tmodel2.predict <- predict(mmodel2,newdata=test.data)
tmodel3.predict <- predict(mmodel3,newdata=test.data)
tmodel4.predict <- predict(mmodel4,newdata=test.data)

#### confusion matrix####
tcm1 <- table(pred=tmodel1.predict, true=test.data$STUDY_CONCLUSION)
print(tcm1)

tcm2 <- table(pred=tmodel2.predict, true=test.data$STUDY_CONCLUSION)
print(tcm2)

tcm3 <- table(pred=tmodel3.predict, true=test.data$STUDY_CONCLUSION)
print(tcm3)


# all three models are good

########################### SVM #########################
plot(gene.train$DOSE, gene.train$newCount, col=gene.train$STUDY_CONCLUSION,
     xlab="DOSE",ylab="COUNT")
# from the plot we can differentiate between positive and negatives easily but 
# not negative and Equivocal
par(mfrow=c(1,2))
plot(gene.train$DOSE, col=gene.train$STUDY_CONCLUSION, ylab="DOSE",main="Clustering on DOSE")
plot( gene.train$newCount, col=gene.train$STUDY_CONCLUSION,ylab="COUNT",main="Clustering on Count")
par(mfrow=c(1,1))


##   svm with validation data####
set.seed(1234)
svm.model1 <- svm(STUDY_CONCLUSION ~ DOSE +newCount+TRIAL_RESULT +STRAIN+ MICROSOMAL_ACTIVATION_USED+RESULT, data = gene.train, cost = 100, gamma = .1)
summary(svm.model1)
svm.pred1 <- predict(svm.model1, gene.valid)
table(pred = svm.pred1, true = gene.valid$STUDY_CONCLUSION)

svm.model2 <- svm(STUDY_CONCLUSION ~ DOSE + newCount +TRIAL_RESULT+RESULT, data = gene.train, cost = 100, gamma = .1)
summary(svm.model2)
svm.pred2 <- predict(svm.model2, gene.valid)
table(pred = svm.pred2, true = gene.valid$STUDY_CONCLUSION)


##  svm with test data#####

test.svm.pred1 <- predict(svm.model1, test.data)
table(pred = test.svm.pred1, true = test.data$STUDY_CONCLUSION)

test.svm.pred2 <- predict(svm.model2, test.data)
table(pred = test.svm.pred2, true = test.data$STUDY_CONCLUSION)


###########Random Forest#######
rf1.two <- randomForest(STUDY_CONCLUSION ~ DOSE +newCount+TRIAL_RESULT +STRAIN+ MICROSOMAL_ACTIVATION_USED+RESULT, data = gene.train, importance=T, ntree=500)
rf1.two
varImpPlot(rf1.two)
rf1.pred1 <- predict(rf1.two, gene.valid)
table(pred = rf1.pred1, true = gene.valid$STUDY_CONCLUSION)


rf2.two <- randomForest(STUDY_CONCLUSION ~ DOSE +newCount+STRAIN+TRIAL_RESULT +RESULT, data = gene.train, importance=T, ntree=500)
rf2.two
varImpPlot(rf2.two)

rf2.pred2 <- predict(rf2.two, gene.valid)
table(pred = rf2.pred2, true = gene.valid$STUDY_CONCLUSION)


rf3.two <- randomForest(STUDY_CONCLUSION ~ logDose+newCount+DOSE+TRIAL_RESULT +RESULT, data = gene.train, importance=T, ntree=500)
rf3.two
varImpPlot(rf3.two)

rf1.pred3 <- predict(rf3.two, gene.valid)
table(pred = rf1.pred3, true = gene.valid$STUDY_CONCLUSION)

############# Test data #########
final.rf1.pred1 <- predict(rf1.two, gene.valid)
table(pred = final.rf1.pred1, true = gene.valid$STUDY_CONCLUSION)

final.rf1.pred2 <- predict(rf2.two, gene.valid)
table(pred = final.rf1.pred2, true = gene.valid$STUDY_CONCLUSION)

final.rf1.pred3 <- predict(rf3.two, gene.valid)
table(pred = final.rf1.pred3, true = gene.valid$STUDY_CONCLUSION)


####################### rpart ####################

# prepare training scheme

set.seed(1234)
rf.label <- gene.train$STUDY_CONCLUSION
# create 10 fold cross validation
cv.10.folds  <- createMultiFolds(rf.label, k =10, times=10)

#setup caret's trainControl object
control <- trainControl(method="repeatedcv", number=10, repeats=3)

#setup doSNOW package for multi-core training , This is helpful as 
# we are going to train a lot of trees

features = c("DOSE","newCount", "STRAIN","MICROSOMAL_ACTIVATION_USED","TRIAL_RESULT")
rpart.train.1 <- gene.train[,features]

set.seed(1234)
cl <- makeCluster(6, type = "SOCK")
registerDoSNOW(cl)
fir.rf <- train(STUDY_CONCLUSION ~ DOSE +newCount+STRAIN+TRIAL_RESULT +RESULT, 
                data = gene.train, method="rf", tuneLength = 3,
                trControl=control)
stopCluster(cl)

set.seed(1234)
cl <- makeCluster(6, type = "SOCK")
registerDoSNOW(cl)
fit.rpart <- train(STUDY_CONCLUSION ~ DOSE +newCount+STRAIN+TRIAL_RESULT +RESULT, data = gene.train, 
                   method="rpart", tuneLength = 3,
                   trControl=control)

stopCluster(cl)

set.seed(1234)
cl <- makeCluster(6, type = "SOCK")
registerDoSNOW(cl)
fit.svm <- train(STUDY_CONCLUSION ~ DOSE +newCount+STRAIN+TRIAL_RESULT +RESULT, data = gene.train, 
                 method="svmRadial", tuneLength = 3,
                 trControl=control)
stopCluster(cl)
results <- resamples(list(CART=fit.rpart, SVM=fit.svm, RF=fir.rf))

#results <- resamples(list(CART=fit.rpart, CART=fit.rpart1))

# compare the models
# summarize differences between modes
summary(results)


#Box and Whisker Plots
#This is a useful way to look at the spread of the estimated accuracies for different methods and how they relate.
# look at the mean and the max columns.
# box and whisker plots to compare models
scales <- list(x=list(relation="free"), y=list(relation="free"))
bwplot(results, scales=scales)

#Density Plots
#show the distribution of model accuracy as density plots. 
#This is a useful way to evaluate the overlap in the estimated behavior of algorithms.

# density plots of accuracy
scales <- list(x=list(relation="free"), y=list(relation="free"))
densityplot(results, scales=scales, pch = "|")


#Dot Plot
#These are useful plots as the show both the mean estimated accuracy as well as the 95% confidence interval 
#(e.g. the range in which 95% of observed scores fell).


# dot plots of accuracy
scales <- list(x=list(relation="free"), y=list(relation="free"))
dotplot(results, scales=scales)


#Scatterplot Matrix

#This create a scatterplot matrix of all fold-trial results for an algorithm compared 
#to the same fold-trial results for all other algorithms. All pairs are compared.
# pair-wise scatterplots of predictions to compare models
splom(results)

#For example, eye-balling the graphs it looks like RF and SVM look strongly correlated, 
# SVM and CART look weekly correlated.

#Pairwise xyPlots
#You can zoom in on one pair-wise comparison of the accuracy of trial-folds for two machine learning algorithms with an xyplot.
# xyplot plots to compare models
par(mfrow=c(2,1))
xyplot(results, models=c("RF", "CART"))
xyplot(results, models=c("RF", "SVM"))
#In this case we can see the seemingly correlated accuracy of the RF and SVM models

#Statistical Significance Tests
#calculate the significance of the differences between the metric distributions of different machine
#learning algorithms. We can summarize the results directly by calling the summary() function.
# difference in model predictions
diffs <- diff(results)
# summarize p-values for pair-wise comparisons
summary(diffs)

## calculating the values for ROC curve
#pred <- prediction(target_pred, target_class)
#perf <- performance(pred,"tpr","fpr"
###########################################FINAL DENDOGRAM ################
# Grab features
rf.final.label <- gene.train$STUDY_CONCLUSION
final.features <- c("DOSE", "newCount", "STRAIN","TRIAL_RESULT","RESULT")
rpart.final.1 <- gene.train[, final.features]

# Run CV and check out results
rpart.final.cv.1 <- rpart.cv(94622, rpart.final.1, rf.final.label, ctrl.3)
rpart.final.cv.1

# Plot
prp(rpart.final.cv.1$finalModel, type = 0, extra = 1, under = TRUE)

################FINAL MODELS ######################
###STEP 1
rf2 <- randomForest(RESULT ~ DOSE + newCount +STRAIN+ TRIAL_RESULT, 
                    data = gene.train, importance=T, ntree=500)
rf2


####STEP2 
rf2.two <- randomForest(STUDY_CONCLUSION ~ DOSE +newCount+STRAIN+TRIAL_RESULT +RESULT, data = gene.train, importance=T, ntree=500)

