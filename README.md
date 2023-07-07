# Analysis-of-a-GEO-set-Database

# Purpose of the Code
This code performs various analyses on gene expression data to identify differentially expressed genes between two groups. The code normalizes the data, calculates correlation matrices, performs t-tests between groups, and generates plots including histograms and volcano plots. The code also uses the limma and fgsea packages to perform a differential gene expression analysis and gene set enrichment analysis (GSEA), respectively. Finally, the code generates a heatmap and boxplot of the results of the GSEA.
Required Packages
This code requires the following R packages to be installed: GEOquery, ggplot2, tidyverse, Biobase, corrplot, preprocessCore, limma, and fgsea.
# Workflow
The code first loads the required packages.
The code then sets the GEO dataset ID as a variable and retrieves the dataset using the getGEO function from the GEOquery package. The gene expression data is then extracted using the exprs function and stored as a data frame.
The code calculates the correlation matrix of the gene expression data and generates a correlation plot using the corrplot function from the corrplot package. The code also generates a boxplot of the gene expression data.
The code normalizes the gene expression data using the normalize.quantiles function from the preprocessCore package and recalculates the correlation matrix and boxplot.
The code performs a log2 transformation on the normalized gene expression data and recalculates the correlation matrix and boxplot.
The code splits the gene expression data into six groups and performs t-tests between groups.
The code calculates the log fold change and adjusted p-values for the differential gene expression analysis using the limma package. The code then generates a volcano plot of the log fold change and adjusted p-values.
The code performs a gene set enrichment analysis (GSEA) using the fgsea package. The code also generates a heatmap and boxplot of the results of the GSEA.
    
 The code performs a gene ontology (GO) enrichment analysis using the enrichGO and plotEnrichment functions from the clusterProfiler package.
 
# Output
The output of the code includes various plots, such as a correlation plot, boxplots, histograms, and volcano plots. The code also generates a table of differentially expressed genes and performs a GSEA and GO enrichment analysis. The code outputs a heatmap and boxplot of the results of the GSEA.

# Parameters Used in FGSEA
pathways: A list of gene sets or pathways to be tested. Each gene set is represented as a character vector of gene symbols.
stats: A numeric vector of gene-level statistics, typically log fold changes or t-statistics, that will be used to rank the genes.
nperm: The number of permutations used to estimate the p-values of the enrichment scores. Increasing the number of permutations can increase the accuracy of the p-values, but it also increases the computation time.
Parameters Used in enrichGO
de_genes: A vector of differentially expressed gene identifiers.
OrgDb: A character string specifying the organism database to be used for gene annotation. In this case, "org.Hs.eg.db" is the annotation package for human genome and it is specified using the OrgDb parameter.
keyType: A character string specifying the type of ID used for gene annotation. In this case, "ENSEMBL" is used to specify Ensembl gene IDs.
ont: A character string specifying the type of Gene Ontology (GO) classification to use. In this case, "BP" is used to specify Biological Process classification.
pAdjustMethod: A character string specifying the method used to adjust p-values for multiple comparisons. In this case, "BH" is used to specify the Benjamini-Hochberg method.
    

 
