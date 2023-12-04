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
setwd("C:/Users/friederike.ehrhart/Documents/Teaching/Workshop materials/Dataset GSE106589 and DESeq2 tutorial")


#import data from file
data = read.csv(file = 'SigUP_F.csv', header = TRUE)

#run the function "gostres" using the first column (named "X") containing the Ensembl IDs from the dataset
gostres <- gost(query = data$X,
                organism = "hsapiens")

# The result is a named list where “result” is a data.frame with the enrichment analysis results
# and “meta” containing a named list with all the metadata for the query.
head(gostres$result)
publish_gosttable(gostres, filename = "gProfiler_result_UP_F.pdf")
publish_gostplot(gostplot(gostres, capped = FALSE, interactive = FALSE, filename = "gProfiler_plot_UP_F.pdf"))


