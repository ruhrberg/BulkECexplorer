# BulkECexplorer
Welcome to the BulkECexplorer github. Here you can find the code related to the app https://ruhrberglab.shinyapps.io/BulkECexplorer/ described in the publication **BulkECexplorer: an online resource to survey endothelial bulk transcriptomes and predict functional transcription**. Several files required by the scripts can be found here too. These files are also provided as supplementary files in the publication.  
## scripts folder
This folder contains four scripts:  
-**app.R**: the code used by the app to run.  
-**BulkECexplorer_Bulk_database_generation.R**: takes sample expression files and creates a database with expression data and the gene classifications carried out by the 0TPM and zTPM methods.  
-**BulkECexplorer_GMM_EM_database_generation.R**: takes sample expression files and creates a database with expression data and the gene expression classifications carried out by the GMM_EM method.  
-**ConfusionMatrix_FalsePositiveRateAnalysis.R**: carries out a confusion matrix analysis given an input gene list with actual positive and negative markers, as well as an analysis of false positive rates.  
## resources folder
Contains several files required by the aforementioned scripts. These are:  
-**align_master.csv**: align rate data output by HISAT2, used to filter out samples with low align rates.  
-**gene_list_all_runs**: gene lists with actual positive/negative markers used in the confusion matrix and false positive rate analyses.  
-**HUVEC_proteome.csv**: - HUVEC proteome data from https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD009687. Used to filter the gene list in the 2nd and 3rd rounds of confusion matrix and false positive rate analyses.  
-**hprot_coding**: human protein-coding gene list obtained from GRCh38-p13_gencode_annotation.csv (supplementary file in publication).  
-**mprot_coding**: mouse protein-coding gene list obtained from GRCm38-p6_gencode_annotation.csv (supplementary file in publication).  
