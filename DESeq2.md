# From raw counts to differentially expressed genes using the DESeq2 R package
In this session we will start with investigating a previously published dataset that contains of RNA sequencing data. The data from the study by Hoffmann et al. 2017 is publicly available from GEO databases under accession number GSE106589.

## Content
1. Introductory lectures about the ([dataset](https://github.com/fehrhart/STREAMLINEworkshop.github.io/blob/main/Dataset.pptx)) and data processing steps done by ([DESeq2](https://github.com/fehrhart/STREAMLINEworkshop.github.io/blob/main/DESeq2.pptx)). The raw count data can be downloaded from [here](https://github.com/fehrhart/STREAMLINEworkshop.github.io/blob/main/GSE106589_geneCounts.csv), the metadata for analysis is this file: [SraRunTable](https://github.com/fehrhart/STREAMLINEworkshop.github.io/blob/main/SraRunTable.txt), while for inspection there is a separate [Excel file](https://github.com/fehrhart/STREAMLINEworkshop.github.io/blob/main/SraRunTable_for_inspection.xlsx) prepared. (The reason is that the SraRunTable contains several columns that Excel misinterprets as dates and it becomes partly unreadable.)
2. Hands-on analysing the data to get from raw counts to differentially expressed genes using the whole dataset and exploring the basic steps ([R code_V1](https://github.com/fehrhart/STREAMLINEworkshop.github.io/blob/main/DESeq2_GSE106589_V1.R))
3. Hands-on analysing the data to get from raw counts to differentially expressed genes for the subset of neuronal cells ([R code_FBneurons](https://github.com/fehrhart/STREAMLINEworkshop.github.io/blob/main/DESeq2_GSE106589_FBneurons.R))
4. If there is time left - try to adapt the code in order to analyse the NPC subset!

## Further reading
1. Original publication on DESeq2 by [Love et al. 2014](https://doi.org/10.1186/s13059-014-0550-8)
2. Extensive [tutorial](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) for all eventual use cases by Love et al. 
