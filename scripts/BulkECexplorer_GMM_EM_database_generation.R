#######################################################################################
###### ------- GMM fitting by Expectation Maximisation -------#########################
#######################################################################################

# N.B. when this wan ran for the BulkECExplorer, the parameters of the GMM were first calculated
# for each run, and then these parameters were subsequently reapplied to each sample/run to calculate
# the posterior probabilities. As this really can be done in one execution, we have merged and tidied the script here.

# Load packages
library(tidyverse)
library(mclust)
library(mixtools)
library(ggpubr)
library(ggthemes)
library(cowplot)

#######################################################################################
#### List protein coding transcripts in mouse and human genome
## Load Processed Genome annotation files

# Gencode files "GRCh38-p13_gencode_annotation.csv" and "GRCm38-p6_gencode_annotation.csv" have been provided as supplementary files
# Processed Gencode files "hprot_coding.csv" and "mprot_coding.csv" have been provided as supplementary files
# Processing of "GRCh38-p13_gencode_annotation.csv" and "GRCm38-p6_gencode_annotation.csv"
# into "hprot_coding.csv" and "mprot_coding.csv" can be found in script BulkECexplorer_Bulk_database_generation.R

hprot_coding <- read.csv("hprot_coding.csv")
mprot_coding <- read.csv("mprot_coding.csv")


#######################################################################################
# Create dataframes to populate
# GMM is for storing the parameters of the fitted Gaussian Mixture Model for each run/sample

GMM  <- data.frame(Low_lambda = numeric(0),     
                   Low_mu = numeric(0),
                   Low_sigma = numeric(0),
                   High_lambda = numeric(0),
                   High_mu = numeric(0),
                   High_sigma = numeric(0),
                   run = character(0))


# EM is for storing the posterior probabilities
EM <- data.frame(Gene.ID = character(0), 
                 Gene.Name = character(0),
                 TPM = numeric(0),
                 cell.type = character(0),
                 run = character(0),
                 log2TPM = numeric(0))


# failEM is for recording run/samples for which a 2-component GMM, indicative of low vs high
# gene expression, could not be fit. 

failEM <- data.frame(EM_fail = character(0))

#######################################################################################
# function for plotting the GMM distributions for qualitative review

plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

#######################################################################################
#### Expectation Maximisation loop
### Transcript abundance files are separated into directories according to their cell type
### Run the loop for each cell type, changing the working directory for each cell type.
## The loop attempts to fit a two component GMM to each run/sample

files = list.files(pattern = "*.tab") # tabs with transcript abundance data from Stringtie
cell = 'HUVEC' # Change according to cell  - HUVEC, HDEMC, Brain, Lung, Retina

## Format data
for (i in 1:length(files)) {
  read <- read.delim(files[i]) %>%
    filter(Gene.Name %in% hprot_coding$gene_name) %>%              # Change hprot/mprot according to species
    mutate(log2TPM = log2(TPM)) %>% 
    filter(TPM != 0) %>%           
    add_column(run = files[i]) %>%
    add_column(cell.type = cell) %>%                                
    select('Gene.ID', 'Gene.Name', 'TPM', 'cell.type', 'run', 'log2TPM')
  read$Gene.ID <- as.character(read$Gene.ID)
  read$Gene.Name <- as.character(read$Gene.Name) %>% toupper()
  
## Run EM 
## When developing the code, we noticed a very small number of runs/samples would fit a 2-component GMM
## without clearly identifying the low/high gene expression distributions, but this could be rectified by
## adjusting the seed. In these scenarios, the mean of the low expression distribution was always clearly above 0.8 log2TPM.
##
## To accomadate the above described runs/samples programmatically, the below code first tries to fit the GMM using a starting seed of 3.0
## If mu[1] (mean of low-expression distribution) was > 0.8, the algorithm would repeat with a seed of 2.
  
  set.seed(3.0)
  mixmdl <- normalmixEM(read$log2TPM, k = 2)
  
  if (mixmdl$mu[1]  < 0.8) {
    
    # Store posterior probability data
    post.df <- as.data.frame(cbind('log2TPM' = mixmdl$x, mixmdl$posterior)) %>%
               rename('P.Low' = comp.1, 'P.High' = comp.2)
    
    read <- merge(read, post.df, by = 0) %>%
      select(-log2TPM.x) %>% 
      rename('log2TPM' = log2TPM.y)  
    
    EM <- bind_rows(EM, read)
    
    ## Stores GMM information
    
    Low_GMM <- as.data.frame(cbind(mixmdl$lambda[1], mixmdl$mu[1], mixmdl$sigma[1])) %>%  
               rename('Low_lambda' = V1, 'Low_mu' = V2, 'Low_sigma' = V3)
    High_GMM <- as.data.frame(cbind(mixmdl$lambda[2], mixmdl$mu[2], mixmdl$sigma[2])) %>%  
                rename('High_lambda' = V1, 'High_mu' = V2, 'High_sigma' = V3)
    GMM <- merge(Low_GMM, High_GMM, by = 'row.names') %>% select(,c(2:7)) %>%
                 add_column(logLik = mixmdl$loglik) %>%add_column(cell.type = cell)  %>% 
                 add_column(run = files[i]) %>% 
                 bind_rows(GMM) 
    
    rm(Low_GMM, High_GMM)  
    
    run_label <- GMM[1,9] %>% str_replace('\\.tab', '')
    
  ### Creates a plot of the GMM and gene expression distribution that can be used to review the GMM fitting  
    
    pl <- data.frame(x = mixmdl$x) %>%
      ggplot() +
      geom_density(aes(x, ..density..), binwidth = 0.5, colour = "blue", lwd= 1.0, 
                   fill = "white") +
      stat_function(geom = "line", fun = plot_mix_comps,
                    args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                    colour = "black", lwd = 1.5) +
      stat_function(geom = "line", fun = plot_mix_comps,
                    args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                    colour = "#E69F00", lwd = 1.5) +
      ylab("Density") +
      xlab("log2TPM") +
      xlim(-10,10) 
    
    GMM[1,9] %>% str_replace('\\.tab', '_EM.tiff') %>% ggsave(., plot = pl, device = 'tiff')
    rm(pl)
    
  } else {
    set.seed(2.0)    # For small minority of samples that did not fit
    mixmdl <- normalmixEM(read$log2TPM, k = 2)
    
    if (mixmdl$mu[1]  < 0.8){
      
      # Store posterior probability data
      post.df <- as.data.frame(cbind('log2TPM' = mixmdl$x, mixmdl$posterior)) %>%
                 rename('P.Low' = comp.1, 'P.High' = comp.2)
      
      read <- merge(read, post.df, by = 0) %>% 
        select(-log2TPM.x) %>% rename('log2TPM' = log2TPM.y) 
      
      EM <- bind_rows(EM, read)
      
      
      Low_GMM <- as.data.frame(cbind(mixmdl$lambda[1], mixmdl$mu[1], mixmdl$sigma[1])) %>%  
                 rename('Low_lambda' = V1, 'Low_mu' = V2, 'Low_sigma' = V3)
      
      High_GMM <- as.data.frame(cbind(mixmdl$lambda[2], mixmdl$mu[2], mixmdl$sigma[2])) %>%  
                  rename('High_lambda' = V1, 'High_mu' = V2, 'High_sigma' = V3)
      
      GMM <- merge(Low_GMM, High_GMM, by = 'row.names') %>% select(,c(2:7)) %>%
            add_column(run = files[i]) %>% add_column(cell.type = cell) %>% bind_rows(GMM)
      
      rm(Low_GMM, High_GMM)  
      
      run_label <- GMM[1,9] %>% str_replace('\\.tab', '')
      
 ### Creates a plot of the GMM and gene expression distribution that can be used to review the GMM fitting  
      
      pl <- data.frame(x = mixmdl$x) %>%
        ggplot() +
        geom_density(aes(x, ..density..), binwidth = 0.5, colour = "blue", lwd= 1.0, 
                     fill = "white") +
        stat_function(geom = "line", fun = plot_mix_comps,
                      args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                      colour = "black", lwd = 1.5) +
        stat_function(geom = "line", fun = plot_mix_comps,
                      args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                      colour = "#E69F00", lwd = 1.5) +
        ylab("Density") +
        xlab("log2TPM") +
        xlim(-10,10) 
      
      GMM[1,9] %>% str_replace('\\.tab', '_EM.tiff') %>% ggsave(., plot = pl, device = 'tiff')
      rm(pl)
      
### Records runs/smaples that failed      
    } else {   
      fail <- data.frame(EM_fail = files[i],cell.type = cell, RNAseq = 'Bulk') 
      failEM <- rbind(failEM, fail)
      rm(fail)
    }
  }
}

#### Remove non-fit samples
#### Qualitatively assess fitting of models, and remove sample/runs w/o a fit indicative of a low/high genes expression
# Add samples/runs to failEM

poor_fit <- as.character(c('SRR8309034.tab')) # Manually added
GMM <- GMM %>% filter(!run %in% poor_fit)
EM <- EM %>% filter(!run %in% poor_fit)

poor_fit <- data.frame(EM_fail = 'SRR8309034.tab', cell.type = 'Lung') 
failEM <- bind_rows(failEM, poor_fit) # list poor fit models into fail EM


########################################################################################################################
# QC datasets and clean up
# QC can be generated from the 'wrangle_loop_prot.R' script

failEM$EM_fail <- failEM$EM_fail  %>% str_replace('\\.tab', '')
EM$run <- EM$run %>% str_replace('\\.tab', '')
GMM$run <- GMM$run %>% str_replace('\\.tab', '')


### Remove datasets that failed QC from MegaWrangle script (Scripts were ran sequentially, so QC must be in R environment)
## QC are samples with a TPM < 1 for CDH5 and KDR

EM <- EM %>% filter(!run %in% QC)
GMM <- GMM %>% filter(!run %in% QC)
failEM <- failEM %>% filter(!EM_fail %in% QC)

## Removed 3 samples that did not have high CDH5/PECAM expression (defined as >1zTPM) or had low read alignment
## Explicitly declared samples, identified by alignQC and zQC in Mega wrangle
EM <- EM %>% filter(!run %in% c("SRR10776647.tab", "SRR10776649.tab", "SRR1303626.tab")) 
GMM <- GMM %>% filter(!run %in% c("SRR10776647.tab", "SRR10776649.tab", "SRR1303626.tab"))

######################################################################################################################
# Categorizing P.High as Unlikely (< 0.33), Undetermined (0.33 - 0.66) or Likey (> 0.66)

EM <- EM %>% mutate(cluster = if_else(P.High > 2/3, 'Active', 
                                      if_else(P.High < 1/3, 'Leaky', 'Undetermined')))


################################################################################################################

