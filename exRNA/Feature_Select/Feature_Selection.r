
## BY ZHAO TIANXIAO

set.seed(1919)
library(dplyr)
library(ggplot2)
library(caret)
library(mlbench)
library(ROCR)
library(nnet)
library(randomForest)
library(Rtsne)
library(plotly)

## Data Import $ Process
raw_dat <- read.table("SCnorm_Combat_Exp_Matrix.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
raw_dat[,1:ncol(raw_dat)] <- lapply(raw_dat[,1:ncol(raw_dat)], as.numeric)
raw_dat<- scale(raw_dat)
raw_dat <- as.data.frame(t(raw_dat))
raw_dat$names <- rownames(raw_dat)
raw_dat <- arrange(raw_dat, names)
rownames(raw_dat) <- raw_dat$names
raw_dat <- raw_dat[,-ncol(raw_dat)]

## Label of Samples
sample_dat <- read.table("scirep_classes.txt", header = T, sep = ",", stringsAsFactors = F)
sample_class <- as.factor(arrange(sample_dat, sample_id)$label)

## RFE
subsets = c(1:10, 20, 30, 40, 50)
ctrl= rfeControl(functions = rfFuncs, method = "cv",verbose = FALSE, returnResamp = "final")
Profile = rfe(raw_dat, sample_class, sizes = subsets, rfeControl = ctrl)
rfe_dat <- raw_dat[,Profile$optVariables]

Profile
Profile$optVariables

## Data Output
write.table(rfe_dat, "RFE_data.txt", col.names = T, row.names = T, sep = "\t")

## t-SNE after Feature Selection(2 Dimension)
p <- Rtsne(rfe_dat, dims = 2, pca = FALSE)
tSNE_Vis <- as.data.frame(cbind(sample_class, p$Y))
names(tSNE_Vis) <- c("class", "X", "Y")
tSNE_Vis$class <- factor(tSNE_Vis$class)
s <- ggplot(data = tSNE_Vis, aes(x = X, y = Y, col = class)) + geom_point()
s

## t-SNE after Feature Selection(3 Dimension)
p <- Rtsne(rfe_dat, dims = 3)
p$Y <- p$Y / 20
tSNE_Vis_3 <- as.data.frame(cbind(sample_class, p$Y))
names(tSNE_Vis_3) <- c("class", "X", "Y", "Z")
tSNE_Vis_3$class <- factor(tSNE_Vis_3$class)
plot_ly(x = tSNE_Vis_3$X, y = tSNE_Vis_3$Y, z = tSNE_Vis_3$Z, type = "scatter3d", mode = "markers", color = tSNE_Vis_3$class)



## t-SNE before Feature Selection
p <- Rtsne(raw_dat, dims = 2, initial_dims = 100)
tSNE_Vis_before <- as.data.frame(cbind(sample_class, p$Y))
names(tSNE_Vis_before) <- c("class", "X", "Y")
tSNE_Vis_before$class <- factor(tSNE_Vis_before$class)
s <- ggplot(data = tSNE_Vis_before, aes(x = X, y = Y, col = class)) + geom_point(size = 3)
s

## Prepare for Logistic Regression
merged_dat <- cbind(sample_class, rfe_dat)
merged_dat$sample_class = factor(merged_dat$sample_class)
temp <- colnames(merged_dat)
temp[2] <- 'piR_hsa_23317'
colnames(merged_dat) <- temp

## Divide Train & Test in Random Forests
n_samples <- nrow(rfe_dat)
n_train <- floor(n_samples * 0.8)
indices <- sample(1:n_samples)
indices <- indices[1:n_train]
RF_train_sample <- rfe_dat[indices,]
RF_train_sample_class <- sample_class[indices]
RF_test_sample <- rfe_dat[-indices,]
RF_test_sample_class <- sample_class[-indices]

## Find Best mtry in RF
for(i in 1:8){
    rf_classifier = randomForest(x = RF_train_sample, y = RF_train_sample_class,  ntree = 400, mtry = i)
    print(mean(rf_classifier$err.rate))
}

## RF
rf_classifier = randomForest(x = RF_train_sample, y = RF_train_sample_class,  ntree = 400, mtry = 7)
rf_classifier

## Test in RF
predicted_classes <- predict(rf_classifier, RF_test_sample)
predicted_probs <- predict(rf_classifier, RF_test_sample, type = 'prob')
table(predicted_classes, RF_test_sample_class)

positive_class <- 'Healthy Control'

test_labels <- vector('integer', length(test_sample_class))
test_labels[test_sample_class != positive_class] <- 0
test_labels[test_sample_class == positive_class] <- 1
pred <- prediction(predicted_probs[, positive_class], test_labels)
roc <- performance(pred, 'tpr', 'fpr') 
plot(roc, main = 'ROC Curve')
auc <- performance(pred, 'auc')
cat('auc =', auc@y.values[[1]], '\n')

## Divide Train & Test in Logistic Regression
n_samples <- nrow(merged_dat)
n_train <- floor(n_samples * 0.8)
indices <- sample(1:n_samples)
indices <- indices[1:n_train]
train_merged_sample <- merged_dat[indices,]
#train_sample_class <- sample_class[indices]
test_merged_sample <- merged_dat[-indices,]
#test_sample_class <- sample_class[-indices]

## Multinom Logistic Regression
mult_train <- multinom(sample_class ~ piR_hsa_23317+ENST00000536684.2+ENST00000626826.1+ENST00000623130.1+ENST00000385273.1+
                       ENST00000516053.2+ENST00000384010.1+ENST00000365699.4+ENST00000362154.1+ENST00000384981.3+ENST00000600213.3+
                       ENST00000383861.1+ENST00000385012.1+ENST00000383858.1+ENST00000384898.1+ENST00000383869.1+ENST00000362162.1+
                       ENST00000384278.1+ENST00000615842.1+ENST00000613023.1, data=train_merged_sample)
step_train <- step(mult_train) 

# Test on train_data
train_pred <- predict(step_train)
table(train_merged_sample$sample_class, train_pred)

# Test on test_data
test_pred <- predict(mult_train, newdata = test_merged_sample)
table(test_merged_sample$sample_class, test_pred)
















