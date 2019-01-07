
## Import Data
raw_dat <- read.table("Matrix.csv", header = T, row.names = 1, sep = " ")
samples_scirep <- read.table("scirep_classes.txt", sep = ",", header = T)
batch_dat <- read.table("scirep_batch.txt", header = T, row.names = 1, sep = ",")

reference_gene <- c("ENST00000408438.1", "ENST00000385271.1", "ENST00000607334.3", "ENST00000385059.1", 
                    "ENST00000362134.1", "ENST00000385245.1", "ENST00000385045.1", "ENST00000362117.1", 
                    "ENST00000384832.1", "ENST00000579846.3")

## Remove NA and less-transcript gene
dat <- na.omit(raw_dat)
dat <- dat[which(rowSums(dat) > 100),]

## Remove duplicate samples
dat <- dat[,-ncol(dat)]
dat <- dat[,-191]
nrow(dat)
head(dat)



## Necessary Library
library(SCnorm)
library(EDASeq)
library(RUVSeq)
library(sva)
library(scRNA.seq.funcs)
library(dplyr)
library(ggplot2)
library(reshape2)

## CPM & RUVs
total_gene <- apply(dat, 2, sum)
CPM_dat <- t( t(dat) * 1000000 / total_gene )


scirepcpm <- log(CPM_dat + 0.001)
scIdx <- matrix(-1, ncol = max(table(samples_scirep$label)), nrow = 2)
tmp1 <- which(samples_scirep$label == "Colorectal Cancer")
scIdx[1, 1:length(tmp1)] <- tmp1
tmp2 <- which(samples_scirep$label == "Healthy Control")
scIdx[2, 1:length(tmp2)] <- tmp2
cIdx <- rownames(scirepcpm)
ruvs <- RUVs(as.matrix(scirepcpm), cIdx, k = 10, scIdx = scIdx, isLog = T)

CPM_RUV_dat <- exp(ruvs$normalizedCounts)
write.table(CPM_RUV_dat, "CPM_RUVs_Exp_Matrix.txt", col.names = T, row.names = T, sep = "\t")

## CPM & Combat

#total_gene <- apply(dat, 2, sum)
#CPM_dat <- t( t(dat) * 1000000 / total_gene )

mat <- as.data.frame(CPM_dat)
batch_info <-read.table("scirep_batch.txt",sep=',',row.names=1,header=T,check.names = FALSE)
batchname <-toString(names(batch_info)[1])
batch_info=batch_info[names(mat),]
mod <- model.matrix(~ 1, data = batch_info)
if (!dim(mat)[2]==dim(batch_info)[1])
    stop('sample numbers in batch info and expression matrix should be same')
combat <- ComBat(
    dat = log(mat+0.001),
    batch = factor(batch_info[,1]),
    mod = mod,
    par.prior = TRUE,
    prior.plots = FALSE
)
mat <- exp(combat)

write.table(mat, "CPM_Combat_Exp_Matrix.txt", col.names = T, row.names = T, sep = "\t")

## SCnorm & Combat

Conditions = rep(1, dim(dat)[2])
DataNorm <- SCnorm(Data = dat, Conditions = Conditions, PrintProgressPlots = TRUE, NCores = 4)
Normalized_Data <- results(DataNorm)

mat <- as.data.frame(Normalized_Data)
batch_info <-read.table("scirep_batch.txt",sep=',',row.names=1,header=T,check.names = FALSE)
batchname <-toString(names(batch_info)[1])
batch_info=batch_info[names(mat),]
mod <- model.matrix(~ 1, data = batch_info)
if (!dim(mat)[2]==dim(batch_info)[1])
    stop('sample numbers in batch info and expression matrix should be same')
combat <- ComBat(
    dat = log(mat+0.001),
    batch = factor(batch_info[,1]),
    mod = mod,
    par.prior = TRUE,
    prior.plots = FALSE
)
mat <- exp(combat)

write.table(mat, "SCnorm_Combat_Exp_Matrix.txt", col.names = T, row.names = T, sep = "\t")



## SCnorm & RUVs

#Conditions = rep(1, dim(dat)[2])
#DataNorm <- SCnorm(Data = dat, Conditions = Conditions, PrintProgressPlots = TRUE, NCores = 4)
#Normalized_Data <- results(DataNorm)

scirepcpm <- log(Normalized_Data + 0.001)
scIdx <- matrix(-1, ncol = max(table(samples_scirep$label)), nrow = 2)
tmp1 <- which(samples_scirep$label == "Colorectal Cancer")
scIdx[1, 1:length(tmp1)] <- tmp1
tmp2 <- which(samples_scirep$label == "Healthy Control")
scIdx[2, 1:length(tmp2)] <- tmp2
cIdx <- rownames(scirepcpm)
ruvs <- RUVs(as.matrix(scirepcpm), cIdx, k = 10, scIdx = scIdx, isLog = T)

SC_RUV_dat <- exp(ruvs$normalizedCounts)
write.table(SC_RUV_dat, "SCnrm_RUVs_Exp_Matrix.txt", col.names = T, row.names = T, sep = "\t")

## Process Batch Effect Information

temp <- as.data.frame(t(dat))
temp$sum <- apply(temp, 1, sum)
sample_dat <- temp[order(rownames(temp)), ]
batch_dat <- batch_dat[order(rownames(batch_dat)),]
batch_dat$names <- rownames(batch_dat)
sample_dat$names <-  rownames(sample_dat)

batch_info_dat <- merge(batch_dat, sample_dat, all=FALSE)
rownames(batch_info_dat) <- batch_info_dat$names
batch_info_dat = batch_info_dat[, -1]

## Visualize Batch Effect_1

batch_info_dat$RNA.Isolation.batch <- factor(batch_info_dat$RNA.Isolation.batch)

pdf("RNA_Isolation_batch.pdf", 7, 7)
ggplot(data=batch_info_dat, aes(x = RNA.Isolation.batch, y = sum, col = RNA.Isolation.batch)) + 
    geom_boxplot() + 
    labs(x = "RNA Isolation batch", y = "Total gene count")
dev.off()

## Visualize Batch Effect_2

batch_info_dat$library.prepration.day <- factor(batch_info_dat$library.prepration.day)

pdf("library_prepration_day.pdf", 7, 7)
ggplot(data=batch_info_dat, aes(x = library.prepration.day, y = sum, col = library.prepration.day)) + 
    geom_violin(data=batch_info_dat, aes(x = library.prepration.day, y = sum, col = library.prepration.day, fill = library.prepration.day)) + 
    labs(x = "library prepration day", y = "Total gene count")
dev.off()

## Visualize Batch Effect_3

batch_info_dat$gel.cut.size.selection <- factor(batch_info_dat$gel.cut.size.selection)

pdf("gel_cut_size_selection.pdf", 7, 7)
ggplot(data=batch_info_dat, aes(x = gel.cut.size.selection, y = sum, col = gel.cut.size.selection)) + 
    geom_boxplot() + 
    labs(x = "gel cut size selection", y = "Total gene count")
dev.off()



## Reference gene

reference_gene_matrix <- dat[match(reference_gene, rownames(dat)),]
reference_gene_matrix <- na.omit(reference_gene_matrix)
rownames(reference_gene_matrix) <- c('MIR1228', 'MIR16-1', 'MIR16-2', 'MIR21', 'MIR23A', 'MIR23B', 'MIR23C', 'MIR451A', 'MIR15A', 'MIR15B')

temp <- melt(t(reference_gene_matrix))
temp$Var2 <- factor(temp$Var2)
temp$value <- log(temp$value)

pdf("Reference_Gene_count(log).pdf", 7, 7)
ggplot(data=temp, aes(x = Var2, y = value, fill = Var2)) + 
    geom_violin() + 
    labs(x = "Reference gene", y = "log of gene count")
dev.off()



rowSums(reference_gene_matrix)

head(Normalized_Data)

scirepcpm <- log(Normalized_Data + 0.001)

head(SC_RUV_dat, n = 20)


