# Gene Set Enrichment Analysis (GSEA) and overrepresentation analysis (ORA)
In this session you will learn how to apply different methods of gene set enrichment and overrepresentation analyis. These methods investigate if there is a functionally related group of genes that are differentially expressed in your dataset and indicate affected biological processes or molecular pathways. For this analysis we will use the files you have created during the DESeq2 tutorial. Alternatively, you can download them also here:
* The list of all significantly up-regulated genes (p-value < 0.05 AND logFC > 1) [SigUP_F.csv](https://github.com/fehrhart/STREAMLINEworkshop.github.io/blob/main/SigUP_F.csv)
* The list of all significantly down-regulated genes (p-value < 0.05 AND logFC < -1) [SigDOWN_F.csv](https://github.com/fehrhart/STREAMLINEworkshop.github.io/blob/main/SigDOWN_F.csv)

## Content
1. Introduction lecture on GSEA and ORA [slides](https://github.com/fehrhart/STREAMLINEworkshop.github.io/blob/main/GSE%20and%20ORA.pptx).
2. [Hands on document](https://github.com/fehrhart/STREAMLINEworkshop.github.io/blob/main/Hands-on%20EnrichR%20and%20gProfiler.docx) for doing GSEA using the webtools EnrichR [webtool](https://maayanlab.cloud/Enrichr/) and g:Profiler [webtool g:GOSt](https://biit.cs.ut.ee/gprofiler/gost)
3. Hands on using the g:Profiler R package ([R code for gprofiler](https://github.com/fehrhart/STREAMLINEworkshop.github.io/blob/main/gProfiler.R))

## Further reading
1. Interpreting omics data with pathway enrichment analysis [Zhao et al. 2023](https://doi.org/10.1016/j.tig.2023.01.003)
