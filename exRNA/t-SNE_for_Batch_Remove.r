
## BY ZHAO TIANXIAO

set.seed(1919)
library(dplyr)
library(ggplot2)
library(Rtsne)
library(plotly)

SCnorm_Combat_Exp <- read.table("SCnorm_Combat_Exp_Matrix.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
temp <- as.data.frame(t(SCnorm_Combat_Exp))
temp$names <- rownames(temp)
temp <- arrange(temp, names)
rownames(temp) <- temp$names
SCnorm_Combat_Exp <- as.data.frame(temp[,-ncol(temp)])

SCnorm_RUVs_Exp <- read.table("SCnorm_RUVs_Exp_Matrix.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
temp <- as.data.frame(t(SCnorm_RUVs_Exp))
temp$names <- rownames(temp)
temp <- arrange(temp, names)
rownames(temp) <- temp$names
SCnorm_RUVs_Exp <- as.data.frame(temp[,-ncol(temp)])

CPM_Combat_Exp <- read.table("CPM_Combat_Exp_Matrix.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
temp <- as.data.frame(t(CPM_Combat_Exp))
temp$names <- rownames(temp)
temp <- arrange(temp, names)
rownames(temp) <- temp$names
CPM_Combat_Exp <- as.data.frame(temp[,-ncol(temp)])                               
                                 
                    
CPM_RUVs_Exp <- read.table("CPM_RUVs_Exp_Matrix.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
temp <- as.data.frame(t(CPM_RUVs_Exp))
temp$names <- rownames(temp)
temp <- arrange(temp, names)
rownames(temp) <- temp$names
CPM_RUVs_Exp <- as.data.frame(temp[,-ncol(temp)])                                 
                                 
batch_info <- read.table("scirep_batch.txt", sep = ",", header = T, row.names = 1, stringsAsFactors = T)

batch_info$RNA.Isolation.batch <- factor(batch_info$RNA.Isolation.batch)

p <- Rtsne(SCnorm_Combat_Exp, dims = 2)
tSNE_Vis_before <- as.data.frame(cbind(batch_info$RNA.Isolation.batch, p$Y))
names(tSNE_Vis_before) <- c("Batch", "X", "Y")
tSNE_Vis_before$Batch <- factor(tSNE_Vis_before$Batch)
s <- ggplot(data = tSNE_Vis_before, aes(x = X, y = Y, col = Batch)) + geom_point(size = 3)
s

p <- Rtsne(SCnorm_RUVs_Exp, dims = 2)
tSNE_Vis_before <- as.data.frame(cbind(batch_info$RNA.Isolation.batch, p$Y))
names(tSNE_Vis_before) <- c("Batch", "X", "Y")
tSNE_Vis_before$Batch <- factor(tSNE_Vis_before$Batch)
s <- ggplot(data = tSNE_Vis_before, aes(x = X, y = Y, col = Batch)) + geom_point(size = 3)
s

p <- Rtsne(CPM_Combat_Exp, dims = 2)
tSNE_Vis_before <- as.data.frame(cbind(batch_info$RNA.Isolation.batch, p$Y))
names(tSNE_Vis_before) <- c("Batch", "X", "Y")
tSNE_Vis_before$Batch <- factor(tSNE_Vis_before$Batch)
s <- ggplot(data = tSNE_Vis_before, aes(x = X, y = Y, col = Batch)) + geom_point(size = 3)
s

p <- Rtsne(CPM_RUVs_Exp, dims = 2)
tSNE_Vis_before <- as.data.frame(cbind(batch_info$RNA.Isolation.batch, p$Y))
names(tSNE_Vis_before) <- c("Batch", "X", "Y")
tSNE_Vis_before$Batch <- factor(tSNE_Vis_before$Batch)
s <- ggplot(data = tSNE_Vis_before, aes(x = X, y = Y, col = Batch)) + geom_point(size = 3)
s


