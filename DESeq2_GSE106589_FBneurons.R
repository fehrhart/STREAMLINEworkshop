# This script is to analyse GSE106589 data with DESeq2 - separating FB neuron data from NPC data and analysing FB neuron data. 

#install required packages - if they are not installed yet
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

#Split datafile by cell type. Columns are indicated with N or F so we can use the tidyr function "contains"
dataF <- select(data, contains("F"))

#Some raw count data visualisation with boxplot and PCA
boxplot(log(dataF))
PCA <- prcomp(t(dataF))
autoplot(PCA, label = TRUE, label.size = 3)

#Load the metadata file. As the count data and metadata are not in the same order, it is required to re-order them. Count matrix and the rows of the column data (information about samples) must be in the same order and this can be achieved by sorting them by npc_line.
metadata <- read.table(file = 'SraRunTable.txt', sep = ",", header = TRUE)
metadata <- metadata[order(metadata$npc_line),]
#Subset the metadata file by NPC and FB neurons
metadataF <- subset(metadata, cell_type=="6 wk FB neuron")

#Create the "dds" object from count data and metadata. The experimental design compares via diagnosis (COS vs control). This is a design if you have only one variable to compare experimental groups.
ddsF <- DESeqDataSetFromMatrix(countData = dataF,
                             colData = metadataF,
                             design = ~ diagnosis, tidy = FALSE)

#Pre-filtering - technically not required, but will reduce memory size and increase speed if we remove samples with too low read counts.
smallestGroupSize <- 3
keep <- rowSums(counts(ddsF) >= 10) >= smallestGroupSize
ddsF <- ddsF[keep,]

#run the DESeq2 function
ddsF <- DESeq(ddsF)

#result output
resultsNames(ddsF)
DEG <- results(ddsF)

#export the complete result file to CSV for later analysis
write.csv(DEG, file = "DEG_F.csv")
#all significantly changed genes, all significantly upregulated genes, all significantly downregulated genes
SigDEG <- DEG[DEG$pvalue < 0.05,]
write.csv(SigDEG, file = "SigDEG_F.csv")
SigUP <- DEG[DEG$pvalue < 0.05 & DEG$log2FoldChange > 1,]
write.csv(SigUP, file = "SigUP_F.csv")
SigDOWN <- DEG[DEG$pvalue < 0.05 & DEG$log2FoldChange < -1,]
write.csv(SigDOWN, file = "SigDOWN_F.csv")

#Visualisation and interpretation of the results and the process

#sort by log2FC
DEG <- DEG[order(DEG$log2FoldChange),]
head(DEG)

#Volcano plot
par(mfrow=c(1,1))
# Make a basic volcano plot
with(DEG, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-6,6)))
# Add colored points: blue if log2FC < -1 and padj<0.05 and red if log2FC > 1 and padj<0.05)
with(subset(DEG, padj<0.05 & log2FoldChange < -1), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(DEG, padj<0.05 & log2FoldChange > 1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#compare normalised counts vs. raw counts 
ddsF_n <- estimateSizeFactors(ddsF); 
ddsF_n <- counts(ddsF_n, normalized=TRUE)
boxplot(log(ddsF_n))
PCA <- prcomp(t(ddsF_n))
autoplot(PCA, label = TRUE, label.size = 3)
