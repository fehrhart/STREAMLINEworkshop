#This script adds HGNC symbols to Ensembl gene identifiers

#clean workspace
rm(list=ls())

#install required packages and load them
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

#library(BiocManager)
#install("BiocFileCache", force = TRUE)
#library(BiocFileCache)
#install.packages("devtools")
#devtools::install_version("dbplyr", version = "2.3.4")
#library(dbplyr)

library(biomaRt)

#set working directory - the folder with all the files in
setwd("C:/Users/friederike.ehrhart/Documents/Teaching/Workshop materials/Biological databases and ID mapping")

#Select a database and a dataset to define your "ensembl" mart object
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

#For selection, with these queries you can look at the available databases and datasets. This part is not required if you know which to choose.
databases <- listEnsembl(ensembl)
datasets <- listDatasets(ensembl)
#With this query, you you can get a list of all available attributes (gene name, ensembl ID, HGNC etc.) you can retrieve from BioMart.
attributes = listAttributes(ensembl)

#The query function is getBM. This now gets your mapping file, a list of all genes with 2 attributes: HGNC gene symbols and Ensembl ids
Ensembl_hsa_genes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl"))

#Import the file you want to add the HGNC symbols to. 
DEG_F = read.csv("DEG_F.csv", header = TRUE)

#Change column names. The common column between two files must have the same name.
colnames(Ensembl_hsa_genes) <-c('X','HGNC')

#Merge with the DEG.csv dataset. Note that not all Ensembl IDs will get a HGNC ID.
DEG_F_HGNC <- merge(DEG_F, Ensembl_hsa_genes, by="X")

#Export to a csv
write.csv(DEG_F_HGNC, file = "DEG_F_HGNC.csv")
