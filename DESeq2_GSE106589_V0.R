# This script is to analyse GSE106589 data with DESeq2. 
#clean workspace
rm(list=ls())

#install required packages and load them
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("ggplot2")
BiocManager::install("tidyr")
BiocManager::install("dplyr")
BiocManager::install("stringr")
BiocManager::install("ggfortify")

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

#load the metadata file. As the count data and metadata are not in the same order, it is required to re-order them.
metadata <- read.table(file = 'SraRunTable.txt', sep = ",", header = TRUE)
metadata <- metadata[order(metadata$npc_line),]

#create the "dds" object from count data and metadata. The experimental design compares via diagnosis (COS and C) and cell_type (FB neuron and NPC). We are interested in comparing COS vs. C for FB neurons and NPCs.
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = metadata,
                             design = ~ diagnosis + cell_type, tidy = FALSE)
#run the DESeq2 function
dds <- DESeq(dds)

#result output
resultsNames(dds)

NPC <- results(dds, contrast=c("diagnosis", "COS", "control"))
FBneuron <- results(dds, contrast=c("diagnosis", "COS", "control"))

DEG <- results(dds)
head(results(dds, tidy=TRUE))

#export to CSV for later analysis
write.csv(DEG, file = "DEG.csv")
#all significantly changed genes, all significantly upregulated genes, all significantly downregulated genes
SigDEG <- DEG[]
write.csv(SigDEG, file = "SigDEG.csv")
SigUP <- DEG[]
write.csv(SigUP, file = "SigUP.csv")
SigDOWN <- DEG[]
write.csv(SigDOWN, file = "SigDOWN.csv")

#Visualisation and interpretation of the results and the process

#compare normalised counts vs. raw counts 
dds_n <- estimateSizeFactors(dds); 
counts(dds_n, normalized=TRUE)
boxplot(log(data))
PCA <- prcomp(t(data))
autoplot(PCA, label = TRUE, label.size = 3)

#sort by log2FC
DEG <- DEG[order(DEG$log2FoldChange),]
head(DEG)

#Volcano plot
par(mfrow=c(1,1))
# Make a basic volcano plot
with(DEG, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-6,6)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(DEG, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(DEG, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

