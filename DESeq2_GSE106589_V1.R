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

#load raw count data and create some overview plots.
data = read.csv(file = 'GSE106589_geneCounts.csv', row.names = 1, header = TRUE)

#Some raw count data visualisation with boxplot and PCA
boxplot(log(data))
PCA <- prcomp(t(data))
autoplot(PCA, label = TRUE, label.size = 3)

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

#Volcano plot
par(mfrow=c(1,1))
# Make a basic volcano plot
with(DEG, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-6,6)))
# Add colored points: blue if log2FC < -1 and padj<0.05 and red if log2FC > 1 and padj<0.05)
with(subset(DEG, padj<0.05 & log2FoldChange < -1), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(DEG, padj<0.05 & log2FoldChange > 1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

