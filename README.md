# BulkECexplorer
Welcome to the BulkECexplorer github. Here you can find the code related to the publication "BulkECexplorer: an online resource to survey endothelial bulk transcriptomes and predict functional transcription". As well as some files required by the scripts. These files are also provided as supplementary files in the publication.
## scripts folder
This folder contains three scripts:  
-BulkECexplorer_Bulk_database_generation.R: takes sample expression files and creates a database with expression data and the gene classifications carried out by the 0TPM and zTPM methods.  
-BulkECexplorer_GMM_EM_database_generation.R: takes sample expression files and creates a database with expression data and the gene expression classifications carried out by the GMM_EM method.  
-ConfusionMatrix_FalsePositiveRateAnalysis.R: carries out a confusion matrix analysis given an input gene list with actual positive and negative markers, as well as an analysis of false positive rates.  
