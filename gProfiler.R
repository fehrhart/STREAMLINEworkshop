#This is a tutorial R script on how to use the R client of g:Profiler for gene set enrichment analysis. 
#Source: https://biit.cs.ut.ee/gprofiler/page/r 

#clean workspace
rm(list=ls())

#install required packages and load them
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("gprofiler2")

library(gprofiler2)

#Set the working directory by copy-pasting your path
setwd("C:/Users/friederike.ehrhart/Documents/Teaching/Workshop materials/Dataset GSE106589")


#import data
data = read.csv(file = 'GSE106589_geneCounts.csv', row.names = 1, header = TRUE)

gostres <- gost(query = c("X:1000:1000000", "rs17396340", "GO:0005005", "ENSG00000156103", "NLRP1"),
                organism = "hsapiens")

# The result is a named list where “result” is a data.frame with the enrichment analysis results
# and “meta” containing a named list with all the metadata for the query.
head(gostres$result)

p <- gostplot(gostres, capped = FALSE, interactive = FALSE)
p
