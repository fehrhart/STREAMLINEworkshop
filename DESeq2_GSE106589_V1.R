# This script is to analyse GSE106589 data with DESeq2. 

#install required packages - needs to be done only once
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("ggplot2")
BiocManager::install("tidyr")
BiocManager::install("dplyr")
BiocManager::install("stringr")
BiocManager::install("ggfortify")

#clean workspace
rm(list=ls())

#load packages
library(DESeq2)
library(ggplot2)
library(ggfortify)
library(tidyr)
library(dplyr)
library(stringr)

#Set the working directory by copy-pasting your path
setwd("C:/...")

#load raw count data
data = read.csv(file = 'GSE106589_geneCounts.csv', row.names = 1, header = TRUE)

#_______Excursion on plotting and visualising data for quality control_________

#Some raw count data visualisation with boxplot, PCA and cluster dendrogram
#Boxplots
boxplot(log(data))
#PCA
PCA <- prcomp(t(data))
autoplot(PCA, label = TRUE, label.size = 3)

#Cluster dendrogram
# Prepare data for cluster analysis
cluster_data <- na.omit(data) # listwise deletion of missing data
cluster_data <- scale(cluster_data) # standardize variables
cluster_data <- t(cluster_data) #transpose
# Determine number of clusters (wait for the plot!)
wss <- (nrow(cluster_data)-1)*sum(apply(cluster_data,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(cluster_data,
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# Performing K-Means Cluster Analysis - using 5 clusters (see plot)
fit <- kmeans(cluster_data, 5) 
# get cluster means (takes some time, about 5min)
aggregate(cluster_data,by=list(fit$cluster),FUN=mean)
# append cluster assignment
cluster_data <- data.frame(cluster_data, fit$cluster)

# Ward Hierarchical Clustering -> creating euclidean distance matrix (will take a few minutes)
d <- dist(cluster_data,
          method = "euclidean") 
fit <- hclust(d, method="ward")

# display cluster dendogram
plot(fit) 
groups <- cutree(fit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters
rect.hclust(fit, k=5, border="red")

#________back to the DESeq protocol______________________

#Load the metadata file. As the count data and metadata are not in the same order, it is required to re-order them. Count matrix and the rows of the column data (information about samples) must be in the same order and this can be achieved by sorting them by npc_line.
metadata <- read.table(file = 'SraRunTable.txt', sep = ",", header = TRUE)
metadata <- metadata[order(metadata$npc_line),]

#Create the "dds" object from count data and metadata. The experimental design compares via diagnosis (COS and control). This would be a design if you have only one variable to compare experimental groups.
dds <- DESeqDataSetFromMatrix(countData = data,
                             colData = metadata,
                             design = ~ diagnosis, tidy = FALSE)

#Pre-filtering - technically not required, but will reduce memory size and increase speed if we remove samples with too low read counts.
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

#run the DESeq2 function
dds <- DESeq(dds)

#result output
resultsNames(dds)
DEG <- results(dds)
head(DEG)

#Visualisation and interpretation of the results and the process

#compare normalised counts vs. raw counts 
dds_n <- estimateSizeFactors(dds); 
dds_n <- counts(dds_n, normalized=TRUE)
boxplot(log(dds_n))
PCA <- prcomp(t(dds_n))
autoplot(PCA, label = TRUE, label.size = 3)

write.csv(dds_n, file = "normcounts.csv")

#Volcano plot
par(mfrow=c(1,1))
# Make a basic volcano plot
with(DEG, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-6,6)))
# Add colored points: blue if log2FC < -1 and padj<0.05 and red if log2FC > 1 and padj<0.05)
with(subset(DEG, padj<0.05 & log2FoldChange < -1), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(DEG, padj<0.05 & log2FoldChange > 1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

