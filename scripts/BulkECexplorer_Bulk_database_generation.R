#######################################################################################
# load packages & files
library(tidyverse)
library(zFPKM)
library(ggthemes)


#######################################################################################
#### List protein coding transcripts in mouse and human genome
## Load Processed Genome annotation files

# Gencode files "GRCh38-p13_gencode_annotation.csv" and "GRCm38-p6_gencode_annotation.csv" have been provided as supplementary files
# Processed Gencode files "hprot_coding.csv" and "mprot_coding.csv" have been provided as supplementary files
# Processing of "GRCh38-p13_gencode_annotation.csv" and "GRCm38-p6_gencode_annotation.csv"
# into "hprot_coding.csv" and "mprot_coding.csv" can be found in script BulkECexplorer_Bulk_database_generation.R

hprot_coding <- read.csv("hprot_coding.csv")
mprot_coding <- read.csv("mprot_coding.csv")

# the files were processed as following:

# hgtf <- rtracklayer::import(choose.files())           # import GRCh38-p13_gencode_annotation
# hprot_coding <- as.data.frame(hgtf) %>%
#   select(,c(7, 11, 12)) %>%
#   filter(gene_type == 'protein_coding') %>%
#   filter(type == 'gene') %>% 
#   select(gene_name) %>%
#   unique()

# Murine
# mgtf <- rtracklayer::import(choose.files())           # import GRCm38-p6_gencode_annotation
# mprot_coding  <- as.data.frame(mgtf) %>%
#   select(,c(7, 11, 12)) %>%
#   filter(gene_type == 'protein_coding') %>%
#   filter(type == 'gene') %>% 
#   select(gene_name) %>%
#   unique()



#######################################################################################
# Make dataframe
Bulk <- data.frame(Gene.ID = character(0), 
                   Gene.Name = character(0),
                   FPKM = numeric(0),
                   zFPKM = numeric(0), 
                   TPM = numeric(0),
                   zTPM = numeric(0), 
                   Coverage = numeric(0),
                   cell.type = character(0),
                   run = character(0))


#######################################################################################
###### ------- Extraction and wrangle loop -------##################################### 
# Transcript abundance files are separated into directories according to their cell type
# Run the loop for each cell type, changing the working directory for each cell type.
# Loop pulls TPM values and calculates zTPM 

file_list_HUVEC = list.files(pattern = "*.tab") # tabs with transcript abundance data from Stringtie
file_list_HDMEC = list.files(pattern = "*.tab")
file_list_lung = list.files(pattern = "*.tab")
file_list_brain = list.files(pattern = "*.tab")
file_list_retina = list.files(pattern = "*.tab")

# Function to process the files with the gene TPM info from each cell type
# cell parameter: Change according to cell  - HUVEC, HDEMC, Brain, Lung, Retina
# files parameter: .tab files with transcript abundance data from Stringtie

generateBulk <- function(files, cell, bulk_df){
  
  # Loop through all files corresponding to a cell type
  for (i in 1:length(files)) {
    read <- read.delim(files[i])
    
    # Detect species and filter using corresponding coding protein df
    if(cell %in% c("HDMEC", "HUVEC")){
      read <- read %>%
        filter(Gene.Name %in% hprot_coding$gene_name)
    }else if(cell %in% c("Lung", "Brain", "Retina")){
      read <- read %>%
        filter(Gene.Name %in% mprot_coding$gene_name)
    }
    
    #zFPKM and zTPM
    read <- read %>%
      rownames_to_column()
    zread <- read %>% 
      select(FPKM) %>% 
      rename(zFPKM = FPKM) %>%
      zFPKM(FPKM) %>% 
      rownames_to_column()
    read <- left_join(read, zread, by = 'rowname')
    
    zread <- read %>% 
      select(TPM) %>% 
      rename(zTPM = TPM) %>%
      zFPKM(TPM) %>% 
      rownames_to_column()
    read <- left_join(read, zread, by = 'rowname') %>%
      select(, c(2:12))
    
    read$Gene.ID <- as.character(read$Gene.ID)
    read$Gene.Name <- as.character(read$Gene.Name)
    
    read <- read %>%
      add_column(cell.type = cell) %>%                  
      add_column(run = files[i]) %>%
      select(Gene.ID, Gene.Name, FPKM, zFPKM, TPM, zTPM, Coverage, cell.type, run)
    bulk_df <- bind_rows(bulk_df, read)
  }
  return(bulk_df)
}

Bulk <- generateBulk(file_list_HUVEC, "HUVEC", Bulk)
Bulk <- generateBulk(file_list_HDMEC, "HDMEC", Bulk)
Bulk <- generateBulk(file_list_lung, "Lung", Bulk)
Bulk <- generateBulk(file_list_brain, "Brain", Bulk)
Bulk <- generateBulk(file_list_retina, "Retina", Bulk)


#######################################################################################
# Remove ".tab" from run id & capitalize all genes

Bulk <- Bulk %>%
  mutate(run = str_replace(run, '\\.tab', '')) %>%
  mutate(Gene.Name = toupper(Gene.Name))

#######################################################################################
# Quality Control #1

# Remove samples with CDH5 < 1 TPM & KDR < 1 TPM
# Exclusion of Project PRJEB14163 (exclude)

exclude <- c('ERR1424945', 'ERR1424946', 'ERR1424947', 'ERR1424948', 'ERR1424949', 'ERR1424950', 'ERR1424951', 
             'ERR1424952', 'ERR1424953', 'ERR1424954', 'ERR1424956', 'ERR1424966', 'ERR1424967', 'ERR1424968', 
             'ERR1424969', 'ERR1424970', 'ERR1424971', 'ERR1424972')

QC <- Bulk %>%
  filter(Gene.Name %in% c("CDH5", "KDR")) %>%
  filter(TPM < 1.0) %>% 
  mutate(run = as.character(run)) %>% 
  pull(run) %>% 
  unique()

QC <- c(QC, exclude)

MegaBulk <- Bulk %>% 
  filter(!run %in% QC)

#########################################################################################
# Assigning whether a gene is considered expressed at different TPM thresholds and the zTPM threshold

MegaBulk <- MegaBulk %>%
  mutate('Expressed_0_TPM' = if_else(TPM == 0, 'Not expressed', 'Expressed')) %>%
  mutate('Expressed_1_TPM' = if_else (TPM < 1, 'Not expressed', 'Expressed')) %>%
  mutate('zTPM_threshold' = if_else(zTPM < -2.38, 'Not expressed', 'Expressed'))

######################################################################################
# Quality control #2

# Remove samples with low align rate.

# Add read information
# 'info' brings in file with align_rate extracted from HISAT2 output
info <- read.csv(choose.files()) # align_master.csv, provided as supplementary file
MegaBulk <- left_join(MegaBulk, info, by = c('run' = 'Run'))

# Filter samples with an align rate less than 60%; could imply poor quality sequencing
MegaBulk <- MegaBulk %>% filter(Align_rate >= 60)

