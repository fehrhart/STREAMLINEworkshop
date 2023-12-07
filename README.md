# STREAMLINE Bioinformatics workshop 2023
From Monday 11. - Friday 15.12.2023 online workshop

## Timetable and further information:
https://docs.google.com/document/d/1PCKnD2mtoXnuLVavG-FzF_KroXQiMBtqEmn133-RqTo/edit 

## Links to the course material:

Please make sure to prepare by installing the required software and download the larger datafiles a week before and contact the organsisers if you get problems! During the workshop we won't have time to deal with installation problems. 
* Part 1: Installation of R, R Studio, PathVisio and Cytoscape - [Part 1](https://github.com/fehrhart/STREAMLINEworkshop.github.io/blob/main/Installation%20Instructions%20part%201.docx)
* Part 2: Installation of a virtual machine for analysis of FASTQ files - [Part 2]()

### The dataset - raw and processed files
The dataset [GSE106589](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106589) from a previously published study: Hoffman GE, Hartley BJ, Flaherty E, Ladran I et al. Transcriptional signatures of schizophrenia in hiPSC-derived NPCs and neurons are concordant with post-mortem adult brains. Nat Commun 2017 Dec 20;8(1):2225. PMID: [29263384](https://doi.org/10.1038/s41467-017-02330-5). The raw count data can be downloaded from [here](https://github.com/fehrhart/STREAMLINEworkshop.github.io/blob/main/GSE106589_geneCounts.csv), the metadata for analysis is this file: [SraRunTable](https://github.com/fehrhart/STREAMLINEworkshop.github.io/blob/main/SraRunTable.txt), while for inspection there is a separate [Excel file](https://github.com/fehrhart/STREAMLINEworkshop.github.io/blob/main/SraRunTable_for_inspection.xlsx) prepared. The list of differentially expressed genes for the FB neurons you are going to create during the tutorial. If needed, the different result files can be downloaded from this list:
* The list of all genes with their expression data [DEG_F.csv](https://github.com/fehrhart/STREAMLINEworkshop.github.io/blob/main/DEG_F.csv)
* The list of significantly (p-value < 0.05) differentially expressed genes [SigDEG_F.csv](https://github.com/fehrhart/STREAMLINEworkshop.github.io/blob/main/SigDEG_F.csv)
* The list of all significantly up-regulated genes (p-value < 0.05 AND logFC > 1) [SigUP_F.csv](https://github.com/fehrhart/STREAMLINEworkshop.github.io/blob/main/SigUP_F.csv)
* The list of all significantly down-regulated genes (p-value < 0.05 AND logFC < -1) [SigDOWN_F.csv](https://github.com/fehrhart/STREAMLINEworkshop.github.io/blob/main/SigDOWN_F.csv)

### The course materials - tutorials and lectures - for the different blocks
1. [From FASTQ to raw counts part 1 and part 2](https://github.com/fehrhart/STREAMLINEworkshop.github.io/blob/main/FASTQ-%3Erawcounts.md)
2. [From raw counts to differentially expressed genes](https://github.com/fehrhart/STREAMLINEworkshop.github.io/blob/main/DESeq2.md)
3. [GSEA and overrepresentation analysis](https://github.com/fehrhart/STREAMLINEworkshop.github.io/blob/main/GSEA.md)
4. [WGCNA](https://github.com/fehrhart/STREAMLINEworkshop.github.io/blob/main/WGCNA.md)
5. [Biological databases](https://github.com/fehrhart/STREAMLINEworkshop.github.io/blob/main/BiologicalDatabases.md), and identifier mapping
6. [Pathway models](https://github.com/fehrhart/STREAMLINEworkshop.github.io/blob/main/PathwayModels.md)
7. [(Multi) omics data visualisation](https://github.com/fehrhart/STREAMLINEworkshop.github.io/blob/main/PathwayDataVisualisation.md)
8. [Networks analysis and network extension](https://github.com/fehrhart/STREAMLINEworkshop.github.io/blob/main/NetworkAnalysis.md)
