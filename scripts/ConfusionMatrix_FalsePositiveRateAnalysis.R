###Code for building a confusion matrix using postulated EC-expressed (positive) and not-expressed (negative) genes


###############################
# Load library & declare functios
library(tidyverse)
library(gt)

`%nin%` = Negate(`%in%`)

###############################

# List of runs with unimodal distribution that are not included in zTPM data
unimodal <- c('SRR11234920', 'SRR10066475', 'SRR10066476', 'SRR10066477', 'SRR10066478', 'SRR10066479', 'SRR10066480', 'SRR10066481',
              'SRR10066482', 'SRR10066483', 'SRR8147640', 'SRR8147648', 'SRR8147649', 'SRR8147656', 'SRR8147657', 'SRR8661657',
              'SRR8661659', 'SRR2657960', 'SRR1303626', 'SRR9163308')

##### Load data
# Load data
CleanBulk <- readRDS('CleanBulk.RDS') %>%   ## Load cleanbulk.RDS containing the gene expression values (TPM and zTPM) found in the BulkECExplorer
  mutate(Expressed_1_TPM = if_else(TPM >= 1, 'Expressed', 'Not expressed')) # %>% 
  #mutate(zTPM_threshold = ifelse(run %in% unimodal, NA, zTPM_threshold))

CleanEM <- readRDS('CleanEM.RDS') ## load CleanEM.RDS containing classification data from GMMs included in the BulkExplorer

# First run: all 145 genes. 108 neg, 37 pos
genes_first_run <- read.csv('gene_list_all_runs.csv') %>%
  filter(first_run == "Y") %>% 
  dplyr::select(Gene.Name, Type)

# Second run: 145 genes minus those negatives (6 genes) that are present in the huvec proteome file = 139 genes. 102 neg, 37 pos
genes_second_run <- read.csv('gene_list_all_runs.csv') %>%
  filter(second_run == "Y") %>% 
  dplyr::select(Gene.Name, Type)

# How those 6 genes were identified:
# genes_proteome <- read.csv('HUVEC_proteome.csv', sep = ";") %>%  # Load list of genes labelled by type (negative = not thought to be expressed in EC, positive = thought to be expressed in ECs)
#   dplyr::filter(!is.na(Protein.evidence)) %>% 
#   dplyr::filter(!Protein.FDR.Confidence == "Low") %>% 
#   mutate(Type = "Positive") %>% 
#   dplyr::select(Gene.symbol, Type) %>% 
#   rename(Gene.Name = Gene.symbol) 
# 
# genes_second_run <- genes_first_run %>%
#   filter(!(Gene.Name %in% genes_proteome$Gene.Name & Type == "Negative")) %>%  # "MAP2"   "GLUL" "GAD1"  "GAPDHS" "PNPLA2" "TAGLN" removed
#   dplyr::select(Gene.Name, Type)

# Third run: 140 genes from second run minus neuronal contaminants = 89 genes. 52 neg, 37 pos
genes_third_run <- read.csv('gene_list_all_runs.csv') %>%
  filter(third_run == "Y") %>% 
  dplyr::select(Gene.Name, Type)
  
###############################
### Prepare data

# For 0TPM and 1TPM data, 240 samples
TPM_database <- CleanBulk %>%
  dplyr::select(Gene.Name, cell.type, TPM, run, Expressed_0_TPM, Expressed_1_TPM)

# For zTPM data, 220 samples, after removing unimodals
zTPM_database <- CleanBulk %>%
  filter(!run %in% unimodal) %>% 
  dplyr::select(Gene.Name, cell.type, run, zTPM_threshold)

# For GMM data, 198 samples
GMM_database <- CleanBulk %>% 
  left_join(CleanEM %>% dplyr::select(!TPM), by = c("Gene.Name", "cell.type", "run")) %>% 
  rename(GMM = cluster) %>% 
  mutate(GMM = case_when(TPM == 0 ~ 'Not expressed',
                         GMM == 'Active' ~ 'Expressed',
                         GMM == 'Leaky' ~ 'Not expressed',
                         GMM == 'Undetermined' ~ NA_character_)) %>%    # Do not consider undetermined group in confusion matrix
  filter(run %in% unique(CleanEM$run)) %>% 
  dplyr::select(Gene.Name, cell.type, run, GMM)

#################################

# Benchmark function
runBenchmark <- function(TPM_data, zTPM_data, GMM_data, genes, filePath){
  
  # Subset tested genes
  TPM_data <- inner_join(genes, TPM_data, join_by(Gene.Name == Gene.Name))
  zTPM_data <- inner_join(genes, zTPM_data, join_by(Gene.Name == Gene.Name))
  GMM_data <- inner_join(genes, GMM_data, join_by(Gene.Name == Gene.Name))

  # Assign TP/TN/FP/FN
  TPM_data <- TPM_data %>% 
    mutate(zero_class = case_when(Type == 'Negative' & Expressed_0_TPM == 'Not expressed' ~ 'TN',
                                  Type == 'Negative' & Expressed_0_TPM == 'Expressed' ~ 'FP',
                                  Type == 'Positive' & Expressed_0_TPM == 'Expressed' ~ 'TP',
                                  Type == 'Positive' & Expressed_0_TPM == 'Not expressed' ~ 'FN'),
           one_class = case_when(Type == 'Negative' & Expressed_1_TPM == 'Not expressed' ~ 'TN',
                                 Type == 'Negative' & Expressed_1_TPM == 'Expressed' ~ 'FP',
                                 Type == 'Positive' & Expressed_1_TPM == 'Expressed' ~ 'TP',
                                 Type == 'Positive' & Expressed_1_TPM == 'Not expressed' ~ 'FN'))

  zTPM_data <- zTPM_data %>% 
    mutate(zfpkm_class = case_when(Type == 'Negative' & zTPM_threshold == 'Not expressed' ~ 'TN',
                                   Type == 'Negative' & zTPM_threshold == 'Expressed' ~ 'FP',
                                   Type == 'Positive' & zTPM_threshold == 'Expressed' ~ 'TP',
                                   Type == 'Positive' & zTPM_threshold == 'Not expressed' ~ 'FN'))
  
  GMM_data <- GMM_data %>% 
    mutate(GMM_class = case_when(Type == 'Negative' & GMM == 'Not expressed' ~ 'TN',
                                 Type == 'Negative' & GMM == 'Expressed' ~ 'FP',
                                 Type == 'Positive' & GMM == 'Expressed' ~ 'TP',
                                 Type == 'Positive' & GMM == 'Not expressed' ~ 'FN'))
    
  # Merge data from all methods
  CM_data <- TPM_data %>% 
    left_join(zTPM_data, by = c("Gene.Name", "Type", "cell.type", "run")) %>% 
    left_join(GMM_data, by = c("Gene.Name", "Type", "cell.type", "run")) %>% 
    relocate("Gene.Name", "Type", "cell.type", "run", "Expressed_0_TPM",
             "Expressed_1_TPM", "zTPM_threshold", "GMM", "zero_class",
             "one_class", "zfpkm_class", "GMM_class")
  
  ##### Combine cell-agnostic and cell-specific metrics for each threshold and classification model
  ## Analysis
  cells <- c('HDMEC', 'HUVEC', 'Brain', 'Lung', 'Retina')
  
  ## First cell agnostic
  datalist <- list()
  
  datalist[['All']] <- (CM_data$zero_class %>% 
                          table() %>% 
                          as.data.frame() %>%
                          pivot_wider(names_from = '.', values_from = Freq) %>%
                          mutate(Classifier = '>0TPM', EC_type = 'All', Negative_class = 'All')) %>%
    bind_rows(
      (CM_data$one_class %>% 
         table() %>% as.data.frame() %>% 
         pivot_wider(names_from = '.', values_from = Freq) %>%
         mutate(Classifier = '>1TPM', EC_type = 'All', Negative_class = 'All'))) %>%
    bind_rows(    
      (CM_data$zfpkm_class %>% 
         table() %>% as.data.frame() %>% 
         pivot_wider(names_from = '.', values_from = Freq) %>%
         mutate(Classifier = 'zFPKM', EC_type = 'All', Negative_class = 'All'))) %>%
    bind_rows(
      CM_data$GMM_class %>% 
        table() %>% as.data.frame() %>% 
        pivot_wider(names_from = '.', values_from = Freq) %>%
        mutate(Classifier = 'GMM', EC_type = 'All', Negative_class = 'All')) %>%
    
    select(Classifier, EC_type, Negative_class, TP, TN, FP, FN)
  
  #then by cell type
  
  for (i in cells) {
    
    cell_df <- CM_data %>% filter(cell.type == i)
    
    datalist[[i]]  <-  (cell_df$zero_class %>% 
                          table() %>% as.data.frame() %>%
                          pivot_wider(names_from = '.', values_from = Freq) %>%
                          mutate(Classifier = '>0TPM', EC_type = i, Negative_class = 'All')) %>% 
      bind_rows(
        (cell_df$one_class %>% 
           table() %>% as.data.frame() %>% 
           pivot_wider(names_from = '.', values_from = Freq) %>%
           mutate(Classifier = '>1TPM', EC_type = i, Negative_class = 'All'))) %>%
      bind_rows(
        (cell_df$zfpkm_class %>% 
           table() %>% as.data.frame() %>% 
           pivot_wider(names_from = '.', values_from = Freq) %>%
           mutate(Classifier = 'zFPKM', EC_type = i, Negative_class = 'All'))) %>%
      bind_rows(
        cell_df$GMM_class %>% 
          table() %>% as.data.frame() %>% 
          pivot_wider(names_from = '.', values_from = Freq) %>%
          mutate(Classifier = 'GMM', EC_type = i, Negative_class = 'All')) %>%
      select(Classifier, EC_type, Negative_class, TP, TN, FP, FN)
  }
  
  ### Caclulate metrics
  Results <- do.call(rbind, datalist) %>%
    replace(is.na(.), 0) %>%
    mutate(Sensitivity = (TP/(FN + TP)),
           Specificity = (TN/(TN + FP))) %>% 
    mutate(across(.cols = c(Sensitivity, Specificity), ~round(.x, 2)))
  
  write.csv(Results, filePath, row.names = F)
  return(CM_data)
}

first_run <- runBenchmark(TPM_database, zTPM_database, GMM_database, genes_first_run, "results/FIRST_RUN_RESULTS.csv")
second_run <- runBenchmark(TPM_database, zTPM_database, GMM_database, genes_second_run, "results/SECOND_RUN_RESULTS.csv")
third_run <- runBenchmark(TPM_database, zTPM_database, GMM_database, genes_third_run, "results/THIRD_RUN_RESULTS.csv")

# FP analysis. Produces data and plots
run_FP_analysis <- function(outDir, runResults){
  
  false_distrib <- runResults %>% 
    dplyr::select(Gene.Name, Type, run, cell.type, TPM, zero_class, one_class, zfpkm_class, GMM_class) %>%
    filter(Type == "Negative") %>%
    mutate(zero_class_distrib = ifelse(Type == "Negative" & zero_class == "TN", "Negative_TN", "Negative_FP")) %>% 
    mutate(one_class_distrib = ifelse(Type == "Negative" & one_class == "TN", "Negative_TN", "Negative_FP")) %>% 
    mutate(zfpkm_class_distrib = ifelse(Type == "Negative" & zfpkm_class == "TN", "Negative_TN", "Negative_FP")) %>% 
    mutate(GMM_class_distrib = ifelse(Type == "Negative" & GMM_class == "TN", "Negative_TN", "Negative_FP")) %>% 
    dplyr::select(Gene.Name, cell.type, run, TPM, zero_class_distrib, one_class_distrib, zfpkm_class_distrib, GMM_class_distrib) %>% 
    mutate(across(c(zero_class_distrib, one_class_distrib, zfpkm_class_distrib, GMM_class_distrib), ~as.factor(.x)))
  
  # Vectors with methods and cell types
  method_columns <- c("zero_class_distrib", "one_class_distrib", "zfpkm_class_distrib", "GMM_class_distrib")
  cell.types <- as.character(unique(false_distrib$cell.type))
  
  # Create vector with the breaks for plotting
  xlabs <- c()
  for(i in seq(0, 100, 2)){
    if(i == 0 | i%%10 == 0){
      xlabs <- c(xlabs, i)
    }else{
      xlabs <- c(xlabs, "")
    }
  }
  
  for(i in 1:length(method_columns)){
  
    ### 1) Grouping by gene only:
    
    # Calculate false positive percentage per gene only
    FP_df_geneOnly <- false_distrib %>% 
      count(Gene.Name, !!sym(method_columns[i]), .drop = F) %>% 
      spread(!!sym(method_columns[i]), n) %>%
      mutate(sum_samples = Negative_TN + Negative_FP) %>% 
      mutate(FP_proportion = Negative_FP/(Negative_TN+Negative_FP)) %>% 
      mutate(FP_percentage = 100*FP_proportion) %>%
      arrange(FP_percentage) 
    
    # Calculate the average/median TPM per gene
    TPM_df <- false_distrib %>%
      group_by(Gene.Name) %>% 
      summarise(avg_TPM = mean(TPM), median_TPM = median(TPM)) %>% 
      mutate(avg_TPM = round(avg_TPM, 2)) %>% 
      mutate(median_TPM = round(median_TPM, 2))
    
    # Add the avg/median TPM info
    # log10_TPM for plotting text color based on avg_TPM. log10 because range is very large
    FP_df_geneOnly <- FP_df_geneOnly %>% 
      left_join(TPM_df, by = "Gene.Name") %>% 
      mutate(log10_TPM = log10(avg_TPM))
    
    # Needs to be a different command, returns error if joint to previous pipe
    FP_df_geneOnly <- FP_df_geneOnly %>% 
      mutate(log10_TPM = ifelse(is.infinite(log10_TPM),
                           min(FP_df_geneOnly %>% filter(is.infinite(log10_TPM) == F) %>% select(log10_TPM)),
                           log10_TPM)) %>% 
      mutate(is_zero_percent = ifelse(FP_percentage == 0, "0%", "")) # geom_text plotting of those with 0% FP
    
    # Needs to go in a separate command to get the proper level order
    FP_df_geneOnly <- FP_df_geneOnly %>% 
      mutate(Gene.Name = factor(Gene.Name, levels = FP_df_geneOnly$Gene.Name)) %>% 
      arrange(FP_percentage, median_TPM)
      
    #FP_df_geneOnly$Gene.Name
    
    # Save the df in a variable and write csv
    write.csv(FP_df_geneOnly, 
              paste0(outDir, method_columns[i], "/individual_gene_distributions/", method_columns[i], "_FPperc_geneOnly.csv"), 
              row.names = F)
    
    ## 1.1) Compact version of the FP% per gene
    png(paste0(outDir, method_columns[i], "/individual_gene_distributions/", method_columns[i], "_FPperc_geneOnly_compact.png"), res = 300, height = (297/210)*2000, width = 2000)
    print(
      ggplot(FP_df_geneOnly, aes(x = Gene.Name, y = FP_percentage)) +
        geom_col(aes(fill = FP_percentage), color = "orange") +
        scale_fill_gradient(low = "white", high = "red") + # Color for bar fill
        coord_flip() + # Horizontal bars
        scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), expand = c(0,1)) + # Need to increase Y axis limit so that text is not cropped
        ggtitle(paste0(method_columns[i], " False positive rate per gene [%]")) +
        xlab(NULL) + 
        ylab("False positive rate [%]") +
        ggprism::theme_prism(base_size = 9) +
        theme(panel.grid.major.x = element_line(color = "black",
                                                size = 0.75,
                                                linetype = 2)) +
        theme(panel.background = element_rect(fill = NA),
              panel.ontop = TRUE, # Grid lines on top of bars
              axis.text.y = element_blank(),
              axis.line.y = element_blank(),
              axis.title.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.title = element_text(size=14,face="bold"), 
              title = element_text(size=14,face="bold"), legend.position = "none")
    )
    dev.off()

    ## 1.2) Detailed version of the FP% per gene
    png(paste0(outDir, method_columns[i], "/individual_gene_distributions/", method_columns[i], "_FPperc_geneOnly.png"), res = 300, height = 4000, width = 1500)
    print(
      ggplot(FP_df_geneOnly, aes(x = Gene.Name, y = FP_percentage)) +
        geom_col(aes(fill = FP_percentage), color = "orange") +
        geom_text(
          aes(label = median_TPM, y  = -6),
          size = 3.5, fontface = "bold", inherit.aes = TRUE
        ) +
        geom_text(
          aes(label = is_zero_percent),
          hjust = -0.5, size = 3, position = position_dodge(width = 1),
          fontface = "bold",
          inherit.aes = TRUE
        ) +
        #scale_color_gradient2(low = "green", mid = "darkgrey", high = "black", midpoint = 0) + # Color for avg_TPM text
        scale_fill_gradient(low = "white", high = "red") + # Color for bar fill
        coord_flip() + # Horizontal bars
        scale_y_continuous(limits = c(-12, 100), breaks = c(0, 25, 50, 75, 100), expand = c(0, .7)) + # Need to increase Y axis limit so that text is not cropped
        ggtitle("False positive rate per gene [%]") +
        xlab(NULL) + 
        ylab("False positive rate [%]") +
        ggprism::theme_prism(base_size = 9) +
        theme(panel.grid.major.x = element_line(color = "black",
                                                size = 0.75,
                                                linetype = 2)) +
        theme(panel.background = element_rect(fill = NA),
              panel.ontop = TRUE, # Grid lines on top of bars
              axis.text.x = element_text(size=12),
              axis.title = element_text(size=14,face="bold"), 
              title = element_text(size=14,face="bold"), legend.position = "none")
    )
    dev.off()
    
    ## 1.3) Histogram: FP% distribution per gene
    png(paste0(outDir, method_columns[i], "/false_positive_rate_histograms/", method_columns[i], "_FPperc_geneOnly_hist.png"), res = 300, width = 4000/2, height = ((297/210)/4)*4000)
      print(
        ggplot(FP_df_geneOnly %>% filter(!FP_percentage == 0), aes(x=FP_percentage)) +
        geom_bar() +
        scale_x_binned(breaks = seq(0, 100, 2), labels = xlabs, limits = c(0, 100)) +
        ggtitle(paste0("FP rate [%] histogram per gene: ", method_columns[i])) +
        ylab("Count") +
        xlab("False positive rate [%]") +
        expand_limits(y = c(0, 35)) +
        scale_y_continuous(breaks = seq(0, 35, 5)) +
        ggprism::theme_prism(base_size = 8)
      )
    dev.off()
      
    ## 1.4) Get info about how many genes belong to each custom group
    df_groups_geneOnly <- FP_df_geneOnly %>%
      filter(!FP_percentage == 0) %>% 
      mutate(FP_group = case_when(FP_percentage >= 90 ~ "high_over_90perc",
                                  FP_percentage < 50 ~ "low_under_50perc")) %>% 
      count(FP_group) %>% 
      mutate(group_percentage = n*100/sum(n))
    
    write.csv(df_groups_geneOnly, 
              paste0(outDir, method_columns[i], "/custom_over90_under50_groups/", method_columns[i], "_FPperc_geneOnly_customGroups.csv"), 
              row.names = F)
    
    ### 2) Per gene-celltype combination
    
    # Calculate false positive percentage per gene-cellype combination
    FP_df_geneCelltype <- false_distrib %>% 
      count(Gene.Name, cell.type, !!sym(method_columns[i]), .drop = F) %>% 
      spread(!!sym(method_columns[i]), n) %>%
      mutate(sum_samples = Negative_TN + Negative_FP) %>% 
      mutate(FP_proportion = Negative_FP/(Negative_TN+Negative_FP)) %>% 
      mutate(FP_percentage = 100*FP_proportion) %>%
      arrange(FP_percentage) %>% 
      mutate(is_zero_percent = ifelse(FP_percentage == 0, "0%", "")) # geom_text plotting of those with 0% FP
    
    # Calculate the average/median TPM per gene
    TPM_geneCelltype <- false_distrib %>%
      group_by(Gene.Name, cell.type) %>% 
      summarise(avg_TPM = mean(TPM), median_TPM = median(TPM), .groups = "drop") %>% 
      mutate(avg_TPM = round(avg_TPM, 2)) %>% 
      mutate(median_TPM = round(median_TPM, 2))
    
    # Add the avg/median TPM info
    FP_df_geneCelltype <- FP_df_geneCelltype %>% 
      left_join(TPM_geneCelltype, by = c("Gene.Name", "cell.type")) %>% 
      arrange(FP_percentage, median_TPM)
    
    # Create 5 dfs, one per cell type
    FP_df_geneHUVEC <- FP_df_geneCelltype %>% filter(cell.type == "HUVEC")
    FP_df_geneHDMEC <- FP_df_geneCelltype %>% filter(cell.type == "HDMEC")
    FP_df_geneLung <- FP_df_geneCelltype %>% filter(cell.type == "Lung")
    FP_df_geneBrain <- FP_df_geneCelltype %>% filter(cell.type == "Brain")
    FP_df_geneRetina <- FP_df_geneCelltype %>% filter(cell.type == "Retina")
    
    # Proper level order for plotting genes by FP%
    makeLevelsFP <- function(df){
      df <- df %>% 
        mutate(Gene.Name = factor(Gene.Name, levels = df$Gene.Name)) %>% 
        mutate(no_data = ifelse(is.na(FP_percentage), "No gene data for this cell type", ""))
      df[is.na(df)] <- 0
      return(df)
    }
    
    # Use same order as the cell.types chr vector
    celltype_dfs <- list(FP_df_geneHDMEC, FP_df_geneHUVEC, FP_df_geneBrain, FP_df_geneLung, FP_df_geneRetina) %>% 
      lapply(makeLevelsFP)
    
    # Save the df in a variable and write csv
    for(j in 1:length(cell.types)){
      write.csv(celltype_dfs[[j]], 
                paste0(outDir, method_columns[i], "/individual_gene_distributions/", method_columns[i], "_",  cell.types[j], "_FPperc_geneCelltype.csv"), 
                row.names = F)
    }

    # 2.1 and 2.2) Gene-celltype FP% plots
    
    for(j in 1:length(cell.types)){
    
      ## 2.1) Compact version of the FP% per gene by celltype
      png(paste0(outDir, method_columns[i], "/individual_gene_distributions/", method_columns[i], "_", cell.types[j], "_FPperc_geneCelltype_compact.png"), res = 300, height = (297/210)*2000, width = 2000)
      print(
        ggplot(celltype_dfs[[j]], aes(x = Gene.Name, y = FP_percentage)) +
          geom_col(aes(fill = FP_percentage), color = "orange") +
          geom_text(
            aes(label = no_data),
            hjust = -0.5, size = 3, position = position_dodge(width = 1),
            fontface = "bold",
            inherit.aes = TRUE
          ) +
          scale_fill_gradient(low = "white", high = "red") + # Color for bar fill
          coord_flip() + # Horizontal bars
          scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), expand = c(0,1)) + # Need to increase Y axis limit so that text is not cropped
          ggtitle(paste0(method_columns[i], " False positive rate per gene-celltype [%]: ", cell.types[j])) +
          xlab(NULL) + 
          ylab("False positive rate [%]") +
          ggprism::theme_prism(base_size = 9) +
          theme(panel.grid.major.x = element_line(color = "black",
                                                  size = 0.75,
                                                  linetype = 2)) +
          theme(panel.background = element_rect(fill = NA),
                panel.ontop = TRUE, # Grid lines on top of bars
                axis.text.y = element_blank(),
                axis.line.y = element_blank(),
                axis.title.y=element_blank(),
                axis.ticks.y=element_blank(),
                axis.title = element_text(size=14,face="bold"), 
                title = element_text(size=14,face="bold"), legend.position = "none")
      )
      dev.off()
      
      ## 2.2) Detailed version of the FP% per gene by celltype
      png(paste0(outDir, method_columns[i], "/individual_gene_distributions/", method_columns[i], "_", cell.types[j], "_FPperc_geneCelltype.png"), res = 300, height = 4000, width = 1500)
      print(
        ggplot(celltype_dfs[[j]], aes(x = Gene.Name, y = FP_percentage)) +
          geom_col(aes(fill = FP_percentage), color = "orange") +
          geom_text(
            aes(label = median_TPM, y  = -6),
            size = 3.5, fontface = "bold", inherit.aes = TRUE
          ) +
          geom_text(
            aes(label = no_data),
            hjust = -0.5, size = 3, position = position_dodge(width = 1),
            fontface = "bold",
            inherit.aes = TRUE
          ) +
          geom_text(
            aes(label = is_zero_percent),
            hjust = -0.5, size = 3, position = position_dodge(width = 1),
            fontface = "bold",
            inherit.aes = TRUE
          ) +
          scale_fill_gradient(low = "white", high = "red") + # Color for bar fill
          coord_flip() + # Horizontal bars
          scale_y_continuous(limits = c(-12, 100), breaks = c(0, 25, 50, 75, 100), expand = c(0, .7)) + # Need to increase Y axis limit so that text is not cropped
          ggtitle(paste0(method_columns[i], " False positive rate per gene-celltype [%]: ", cell.types[j])) +
          xlab(NULL) + 
          ylab("False positive rate [%]") +
          ggprism::theme_prism(base_size = 9) +
          theme(panel.grid.major.x = element_line(color = "black",
                                                  size = 0.75,
                                                  linetype = 2)) +
          theme(panel.background = element_rect(fill = NA),
                panel.ontop = TRUE, # Grid lines on top of bars
                axis.text.x = element_text(size=12),
                axis.title = element_text(size=14,face="bold"), 
                title = element_text(size=14,face="bold"), legend.position = "none")
      )
      dev.off()
    }
    
    ## 2.3) Get info about how many gene-celltype combinations belong to each custom group
    
    getGroupsInfo <- function(df, celltype){
      df <- df %>% 
        filter(!FP_percentage == 0) %>% 
        mutate(FP_group = case_when(FP_percentage >= 90 ~ "high_over_90perc",
                                    FP_percentage < 50 ~ "low_under_50_perc")) %>% 
        count(FP_group) %>% 
        mutate(group_percentage = n*100/sum(n))
      
      write.csv(df, 
                paste0(outDir, method_columns[i], "/custom_over90_under50_groups/", method_columns[i], "_", celltype, "_FPperc_geneCelltype_customGroups.csv"), 
                row.names = F)
    }
    
    # All celltypes together
    getGroupsInfo(FP_df_geneCelltype, "allCellTypes")
    
    # Split by celltype (dfs must be in same order as celltypes)
    for(j in 1:length(cell.types)){
      getGroupsInfo(df = celltype_dfs[[j]], celltype = cell.types[j])
    }
    
    
    ## 2.4) Histograms: All celltypes together in the same plot
    
    FP_df_geneCelltype_noZeroes <- FP_df_geneCelltype %>% filter(!FP_percentage == 0)
    
    png(paste0(outDir, method_columns[i], "/false_positive_rate_histograms/", method_columns[i], "_allCelltypes_FPperc_geneCelltype_hist.png"), res = 300, width = 4000/2, height = ((297/210)/4)*4000)
    print(
      ggplot(FP_df_geneCelltype_noZeroes, aes(x=FP_percentage, fill = cell.type)) +
        geom_bar() +
        scale_x_binned(breaks = seq(0, 100, 2), labels = xlabs, limits = c(0, 100)) +
        ggtitle(paste0("FP rate [%] histogram per gene and cell type: ", method_columns[i])) +
        ylab("Count") +
        xlab("False positive rate [%]") +
        scale_y_continuous(breaks = function(z) seq(0, range(z)[2], by = 5)) +
        ggprism::theme_prism(base_size = 8)
    )
    dev.off()
    
    ## 2.5) Histograms: One celltype only per plot
    
    hist_geneCelltype_list <- list()
    max_counts <- 0
    
    for(j in 1:length(cell.types)){
      
      # Subset one celltype
      FP_df_celltype <- FP_df_geneCelltype_noZeroes %>%
        filter(cell.type == cell.types[j])
      
      # Store the max value of histogram, to later on make all have the same Y axis
      max_counts_temp <- max(hist(FP_df_celltype$FP_percentage, breaks = seq(0, 100, 2))$counts)
      if(max_counts_temp > max_counts){max_counts <- max_counts_temp}
      
      # Save the plot with default Y axis lim
      hist_geneCelltype_list[[j]] <- ggplot(FP_df_celltype, aes(FP_percentage)) +
        geom_bar() +
        scale_x_binned(breaks = seq(0, 100, 2), labels = xlabs, limits = c(0, 100)) +
        ggtitle(paste0("FP rate [%] histogram per gene in ",  cell.types[j], " ", method_columns[i])) +
        ylab("Count") +
        xlab("False positive rate [%]") +
        ggprism::theme_prism(base_size = 8)
    }
    
    # Make all plots have the same Y lim depending on the highest value of the histograms
    set_Y_axis <- function(fp_hist, y_lim){
      
      y_lim <- plyr::round_any(y_lim, 5, ceiling) # Round to closest higher number ending with 5 or 10
      
      fp_hist <- fp_hist + 
        scale_y_continuous(breaks = seq(0, y_lim, 5), limits = c(0, y_lim))
      
      return(fp_hist)
    }
    
    # Apply the max Y axis limit to all plots
    hist_geneCelltype_list_processed <- lapply(hist_geneCelltype_list, set_Y_axis, y_lim = max_counts)
    
    # Save the plots
    for(j in 1:length(cell.types)){
      
      png(paste0(outDir, method_columns[i], "/false_positive_rate_histograms/", method_columns[i], "_", cell.types[j],  "_FPperc_geneCelltype_hist.png"), res = 300, width = 4000/2, height = ((297/210)/4)*4000)
      print(
        hist_geneCelltype_list_processed[[j]]
      )
      dev.off()
    }
  }
}

run_FP_analysis("results/FP_analysis/FP_analysis_firstRun/", first_run)
run_FP_analysis("results/FP_analysis/FP_analysis_secondRun/", second_run)
run_FP_analysis("results/FP_analysis/FP_analysis_thirdRun/", third_run)
