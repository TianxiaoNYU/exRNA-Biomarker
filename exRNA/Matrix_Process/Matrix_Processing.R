setwd("~/ZTX_Expression_Matrix")
dat = read.table("backup.txt", sep = ",", row.names = 1, header = T)

library(SCnorm)
Conditions = rep(1, dim(dat)[2])
DataNorm <- SCnorm(Data = dat, Conditions = Conditions, PrintProgressPlots = TRUE, NCores = 4)
NormalizedData <- results(DataNorm)
write.table(NormalizedData, file="Normalized_Matrix.txt", sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)