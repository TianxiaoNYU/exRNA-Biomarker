
set.seed(1919)
library(dplyr)
library(randomForest)
library(caret)
library(mlbench)
library(ROCR)

raw_dat <- read.table("~/SCnorm_Combat_Exp_Matrix.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
raw_dat[,1:ncol(raw_dat)] <- lapply(raw_dat[,1:ncol(raw_dat)], as.numeric)
raw_dat<- scale(raw_dat)
raw_dat <- as.data.frame(t(raw_dat))
raw_dat$names <- rownames(raw_dat)
raw_dat <- arrange(raw_dat, names)
rownames(raw_dat) <- raw_dat$names
raw_dat <- raw_dat[,-ncol(raw_dat)]







sample_dat <- read.table("~/scirep_classes.txt", header = T, sep = ",", stringsAsFactors = F)
sample_class <- as.factor(arrange(sample_dat, sample_id)$label)

subsets = c(1:10, 20, 30, 40, 50)
ctrl= rfeControl(functions = rfFuncs, method = "cv",verbose = FALSE, returnResamp = "final")
Profile = rfe(raw_dat, sample_class, sizes = subsets, rfeControl = ctrl)

Profile

rfe_dat <- raw_dat[,Profile$optVariables]

n_samples <- nrow(rfe_dat)
n_train <- floor(n_samples * 0.8)
indices <- sample(1:n_samples)
indices <- indices[1:n_train]
train_sample <- rfe_dat[indices,]
train_sample_class <- sample_class[indices]
test_sample <- rfe_dat[-indices,]
test_sample_class <- sample_class[-indices]

test_sample_class





for(i in 1:8){
    rf_classifier = randomForest(x = train_sample, y = train_sample_class,  ntree = 400, mtry = i)
    print(mean(rf_classifier$err.rate))
}



rf_classifier = randomForest(x = train_sample, y = train_sample_class,  ntree = 400, mtry = 6)
rf_classifier

predicted_classes <- predict(rf_classifier, test_sample)

predicted_probs <- predict(rf_classifier, test_sample, type = 'prob')

# 定义versicolor为正类别
positive_class <- 'Healthy Control'

test_labels <- vector('integer', length(test_sample_class))
test_labels[test_sample_class != positive_class] <- 0
test_labels[test_sample_class == positive_class] <- 1
# 通过prediction函数，使用预测为正样本的概率和真实类别创建一个对象pred
pred <- prediction(predicted_probs[, positive_class], test_labels)

roc <- performance(pred, 'tpr', 'fpr') 
plot(roc, main = 'ROC Curve')
auc <- performance(pred, 'auc')
cat('auc =', auc@y.values[[1]], '\n')


