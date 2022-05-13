BiocManager::install(version = "3.14")
BiocManager::install(c("affy","affyPLM","sva","AnnotationDbi","hgu133plus2.db"))
library(sva)
library(affyPLM)
library(affy)
library(AnnotationDbi)
library(hgu133plus2.db)

#Reading in CEL files and producing unscaled data
Data <- ReadAffy(celfile.path="/usr4/bf527/lkzhang/lavalamp/project-1-lava-lamp-1/samples/")
rma <- rma(Data, normalize = TRUE)

Data_LavaLamp <- fitPLM(Data, normalize=TRUE, background=TRUE)
RLE(Data_LavaLamp, main="Relative Log Expression Plot")
NUSE(Data_LavaLamp, main="Normalized Unscaled Standard Error Plot")
RLE_Data <- RLE(Data_LavaLamp,type = "stats")
NUSE_Data <- NUSE(Data_LavaLamp,type ="stats")

library(tidyverse)
#Histogram Production 
hist(RLE_Data["median",], xlab = "Median RLE Values", main = "Histogram of Median RLE Values")
hist(NUSE_Data["median",], xlab = "Median NUSE Values", main = "Histogram of Median NUSE Values")

#Writing corrected batch csv 
express <- exprs(rma)
meta <- read.csv("/project/bf528/project_1/doc/proj_metadata.csv")
batch = meta$normalizationcombatbatch
model <- model.matrix(~normalizationcombatmod, data = meta)

corrected <- ComBat(express, batch, model)
write.csv(corrected, file = "corrected_batch.csv")

#PCA
trans_corrected <- t(scale(t(corrected))) #Transposed, scaled, and transposed again 
pca <- prcomp(trans_corrected, center = FALSE, scale = FALSE)
rotation <- as_tibble(pca$rotation)

#Calculate Percent Variance 
variance <- pca$sdev^2
percent_variance <- variance/sum(variance)*100

labelx = paste("PC1 with Percent Variance of", toString(percent_variance[1]), "%")
labely = paste("PC2 with Percent Variance of", toString(percent_variance[2]), "%")

#Generate Graph 
ggplot(rotation, aes(x = PC1, y = PC2)) +
  geom_point() +
  labs(title = "PC1 vs PC2", x = labelx, y = labely)
