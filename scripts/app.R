# BulkECexplorer

{
library(shiny)
library(tidyverse)
library(shinydashboard)
library(ggthemes)
library(gt)
library(shinyjs)
}

# Load data
CleanEM <- readRDS('CleanEM.RDS')
CleanBulk <- readRDS('CleanBulk.RDS')

# Assign correct cell type labels
CleanBulk <- CleanBulk %>%
  mutate(cell.type = as.character(cell.type)) %>%
  mutate(cell.type = replace(cell.type, cell.type == "Brain", "Brain EC")) %>%
  mutate(cell.type = replace(cell.type, cell.type == "Lung", "Lung EC")) %>%
  mutate(cell.type = replace(cell.type, cell.type == "Retina", "Retina EC")) %>%
  mutate(cell.type = factor(cell.type, levels = c("HUVEC", "HDMEC", "Lung EC", "Brain EC", "Retina EC")))

CleanEM <- CleanEM %>%
  mutate(cell.type = as.character(cell.type)) %>%
  mutate(cell.type = replace(cell.type, cell.type == "Brain", "Brain EC")) %>%
  mutate(cell.type = replace(cell.type, cell.type == "Lung", "Lung EC")) %>%
  mutate(cell.type = replace(cell.type, cell.type == "Retina", "Retina EC")) %>%
  mutate(cell.type = factor(cell.type, levels = c("HUVEC", "HDMEC", "Lung EC", "Brain EC", "Retina EC")))

unimodal <- c('SRR11234920', 'SRR10066475', 'SRR10066476', 'SRR10066477', 'SRR10066478', 'SRR10066479', 'SRR10066480', 'SRR10066481',
              'SRR10066482', 'SRR10066483', 'SRR8147640', 'SRR8147648', 'SRR8147649', 'SRR8147656', 'SRR8147657', 'SRR8661657',
              'SRR8661659', 'SRR2657960', 'SRR1303626', 'SRR9163308')

# Genes with multiple entries due to multiple ENSEMBL IDs sharing the same GENE SYMBOL
dup_genes <- c("CCDC39", "DUS4L-BCAP29", "AHRR", "MATR3", "SOD2", "POLR2J3", 
               "TBCE", "ARMCX5-GPRASP2", "TMSB15B", "HSPA14", "GOLGA8M", "GGT1", 
               "ABCF2", "ARHGAP11B", "SDHC", "PINX1", "GUCA1A", "PDE11A", "SMIM40", 
               "ZNF883", "FAM220A", "ZKSCAN7", "PPP2R3D", "FAM205A2", "IL11RA2", 
               "CCL27A", "ST6GALNAC2", "NHEJ1", "PCDHA11", "PTP4A1", "NDOR1", 
               "MIA3", "LARP1B", "4933434E20RIK", "PAKAP", "CCL21B", "CCL19", 
               "ZFP91", "ARHGEF4", "SNHG4", "CCL21A", "ALDOA", "JAKMIP3", "DDIT3", 
               "SEPTIN2", "ZC3H11A", "DPEP2", "ARHGAP26", "NNT", "MUC2", "GCAT")

################################################################################
##### Define UI for application#################################################

ui <- dashboardPage(
  dashboardHeader(title = 'BulkECexplorer',
                    titleWidth = 300),
  
  ### Side menu (dark blue area)
  dashboardSidebar(
    
    # Tab structure definition, with home and analysis tabs
    sidebarMenu(
      id = 'tabs',
      
      # Home tab selection from the side panel
      menuItem('Home',
               tabName = 'Home_tab',
               icon = icon('home')),
      
      # Analysis tab selection from the side panel
      menuItem('Analysis',
               tabName = 'Analysis_tab',
               icon = icon('chart-bar')),
      
      # Fields for gene and custom TPM inputs. Default TPM value is 1.
      menuItem(textInput(inputId = 'GeneSearch', label = 'Search gene', placeholder = 'Gene')),
      menuItem(textInput(inputId = 'CustomTPM', label = 'Custom TPM threshold', value = 1)),
      menuItem(submitButton('Search')),
      
      menuItem(htmltools::HTML("<b>Dataset types: <br/><br/>
               Human:<br/>
               -HUVEC<br/>
               -HDMEC<br/><br/>
               Mouse:<br/>
               -Lung EC<br/>
               -Brain EC<br/>
               -Retina EC</b>")),
      
      # These UI elements for downloads are within a conditionalPanel() and 
      # only appear after the button GeneSearch is pressed,
      # as defined by "condition = input.GeneSearch"
      menuItem(
        conditionalPanel(
          condition = "input.GeneSearch",
          # Table Download Button
          downloadButton("download.plots.png", "Download Plots (.png)",
                         style = "color: #fff; background-color: #27ae60; border-color: #023020; padding: 5px 14px 5px 14px; margin: 0px 20px 5px 20px; ")
          )
        ),
      menuItem(
        conditionalPanel(
          condition = "input.GeneSearch",
          # Table Download Button
          downloadButton("download.plots.tiff", "Download Plots (.tiff)",
                         style = "color: #fff; background-color: #27ae60; border-color: #023020;padding: 5px 14px 5px 14px; margin: 0px 20px 5px 20px; ")
        )
      ),
      menuItem(
        conditionalPanel(
          condition = "input.GeneSearch",
          # Table Download Button
          downloadButton("download.plots.pdf", "Download Plots (.pdf)",
                         style = "color: #fff; background-color: #27ae60; border-color: #023020;padding: 5px 14px 5px 14px; margin: 0px 20px 5px 20px; ")
        )
      ),
               width = 300)),
  
  ### UI body
  dashboardBody(
    tabItems(
      
      # Home tab
      tabItem(tabName = 'Home_tab',
              fluidPage(
                fluidRow(
                  img(src = 'Header3.png', width = '95%', align = 'center'),
                  br(),
                  br(),
                  box(title = 'Information on the Resource', width = 12, solidHeader = T, status = 'primary',
                      tags$div(
                        h4(
                          "This resource is described in the following publication: ",
                          tags$a(href = "https://www.nature.com/articles/s44161-024-00436-w",
                                 "The BulkECexplorer compiles endothelial bulk transcriptomes to predict functional versus leaky transcription")
                        )
                      ),
                    h4('The Bulk RNAseq Vascular Endothelial Cell (BulkEC) explorer compiles >200 bulk RNAseq datasets from the European Nucleotide Archive (ENA) 
                      for 5 endothelial cell (EC) subtypes, collected in 2020 as an exhaustive compilation of publicly available datasets meeting 
                      strict inclusion criteria. The EC subtypes included are: human dermal microvascular ECs (HDMECs), human umbilical vein ECs (HUVECs),
                      mouse lung ECs, mouse brain ECs and mouse retina ECs. A description of the inclusion criteria and a list of all data sets included 
                      in the BulkECexplorer can be obtained from the authors.'),
                    h4('To improve representation of the functional EC transcriptome, the BulkECexplorer incorporates published computational methods [1,2] 
                      to classify genes detected in the datasets as actively expressed or leaky. A leaky classification infers that a gene is unlikely
                      to be functional in ECs [1-3]. Nevertheless, experimental evidence should be considered before concluding that a gene classified 
                      as actively expressed is functional in ECs, and vice versa, whether a gene classified as a leaky gene has no function in ECs.'),
                    h4("Note: Datasets from freshly isolated mouse endothelial cells contain transcripts of host organ cells that interact with endothelium.
                    This app does not inform whether these transcripts reflect transcriptional mimicry and/or host cell contamination."),
                      br(),
                    h4('Correspondence: Christiana Ruhrberg at c.ruhrberg@ucl.ac.uk'),
                  )),
                box(title = 'Background',  
                    width = 12, solidHeader = T, status = 'primary',
                    h4('Bulk, scRNAseq and spatial RNAseq data provide information on which genes are expressed in the assayed cell type but 
                      are unlikely to adequately reflect a cell type’s functional transcriptome. First, mRNA transcripts for functional genes
                      may remain undetected in any individual RNAseq experiment for technical reasons [e.g. 2,4,5]. Second, bulk RNAseq 
                      typically detects lowly expressed mRNA transcripts thought to be the biologically irrelevant products of leaky
                      transcription [1,6]. Unlike more moderately expressed genes, leaky genes are defined as neither associated with active 
                      chromatin markers nor functional within the assayed cell type [1-3].'),
                    br(),
                    h4("We propose that compiling data from many bulk RNAseq datasets of vascular ECs into one compendium makes technical variables,
                      which might prevent a gene from being detected in an individual RNAseq experiment, become less impactful. To interrogate EC
                      gene expression, we have compiled data from >200 bulk EC RNAseq datasets into an in silico resource termed the BulkECexplorer.
                      By applying computational methods to the EC datasets contained in this compendium, this resource helps predict whether the
                      transcript levels of lowly expressed genes for each of five EC subtypes meet the predefined criteria of a functional gene.")
                      ),
                box(title = 'References',  
                    width = 12, solidHeader = T, status = 'primary',
                    p('1. Hebenstreit, D., M. Fang, M. Gu, V. Charoensawan, A. van Oudenaarden and S. A. Teichmann (2011). RNA sequencing reveals two major classes of gene expression levels in metazoan cells. Mol Syst Biol 7: 497. PMID 21654674'),
                    p('2. Hart, T., H. K. Komori, S. LaMere, K. Podshivalova and D. R. Salomon (2013). Finding the active genes in deep RNA-seq gene expression studies. BMC Genomics 14: 778. PMID 24215113'),
                    p('3. Nagaraj, N., J. R. Wisniewski, T. Geiger, J. Cox, M. Kircher, J. Kelso, S. Paabo and M. Mann (2011). Deep proteome and transcriptome mapping of a human cancer cell line. Mol Syst Biol 7: 548. PMID 22068331'),
                    p('4. Toung, J. M., M. Morley, M. Li and V. G. Cheung (2011). "RNA-sequence analysis of human B-cells." Genome Res 21(6): 991-998. PMID 21536721'),
                    p('5. Marinov, G. K., B. A. Williams, K. McCue, G. P. Schroth, J. Gertz, R. M. Myers and B. J. Wold (2014). From single-cell to cell-pool transcriptomes: stochasticity in gene expression and RNA splicing. Genome Res 24(3): 496-510. PMID 24299736'),
                    p('6. Gray, J. M., D. A. Harmin, S. A. Boswell, N. Cloonan, T. E. Mullen, J. J. Ling, N. Miller, S. Kuersten, Y. C. Ma, S. A. McCarroll, S. M. Grimmond and M. Springer (2014). SnapShot-Seq: a method for extracting genome-wide, in vivo mRNA dynamics from a single total RNA sample. PLoS One 9(2): e89673. PMID 24586954')
                    )
                )
              ),
      
      # Analysis tab
      tabItem(tabName = 'Analysis_tab',
              h2(textOutput('title.gene')),
              span(htmlOutput('dup_genes_text'), style="color:red;font-size:20px;font-weight:bold"), ###################################
              fluidPage(
                fluidRow(
                  box(title = 'Gene Detection and expression range', status = 'primary', solidHeader = T, width = 12,
                      h4('The panel on the left shows how frequently a gene of interest is detected in the compendium datasets that express the gene of interests with a value of transcripts per million (TPM) > 0). Each column represents a different EC subtype.'),
                      br(),
                      h4("The panel on the right shows the TPM value for each dataset, whereby each datapoint represents an individual sample. TPM units are not directly comparable between experiments, but help the user to understand the expression range of a gene.")

                        )),
                
                ### TPM > 0 Plots UI
                fluidRow(
                      box(title = textOutput('Det.gene'), status = 'primary', solidHeader = T,
                          plotOutput('Det.plot'),
                          ),
                       box(title = textOutput('Range.gene'), status = 'primary', solidHeader = T,
                               plotOutput('Range.plot', 
                                          hover = 'Range_hover'),
                       ),
                      ),
                br(),
                
                ### TPM > custom_threshold UI
                fluidRow(
                  box(title = textOutput('Det.gene.custom'), status = 'primary', solidHeader = T,
                      plotOutput('Det.plot.custom'),
                  ),
                  box(title = textOutput('Range.gene.custom'), status = 'primary', solidHeader = T, # WARNING: even if it says the same as the other range plot, it needs to be a different output to avoid errors!
                      plotOutput('Range.plot.custom',
                                 hover = 'Range_hover'),
                  ),
                ),

                br(),
                br(),
                    
                fluidRow(
                  box(title = 'Predicting leaky vs. active genes - Gaussian Mixture Models', status = 'success', solidHeader = T, width = 12,
                      h4("Several lines of evidence suggests that two classes of protein-coding transcripts exist within a bulk RNAseq dataset [1-4]: Moderate-to-highly expressed transcripts, which are thought to encode the functional proteome of the assayed cell type (active genes), and lowly expressed transcripts, which are non-functional by-products of neighbouring genomic activity (leaky genes)."),
                      br(),
                      h4('The mixture of moderate-to-highly expressed and lowly expressed transcripts within a bulk RNAseq sample produces a bimodal distribution of transcript abundance [1]. A 2-component Gaussian mixture model (GMM) can be fit to bulk RNAseq data to estimate the parameters of the active vs. leaky transcript distribution of transcripts, after which one can calculate the probability of a gene belonging to either distribution (example below) [1]. We have fit 2-component GMMs to each sample in our compendium to demonstrate how gene of interest are classified (samples that were unimodal or required 3-components to fit a GMM were excluded).'),
                      p())
                    ),
                
                ### GMM plots UI
                fluidRow(
                  box(title = textOutput('GMM.gene'), status = 'success', solidHeader = T,
                             plotOutput('GMM.plot')),
                  box(title = 'Illustrative GMM', status = 'success', solidHeader = T,
                      img(src = 'GMM3.png', width = '75%', height = '40%', style="display: block; margin-left: auto; margin-right: auto;"),
                      ),
                  br(),
                  br(),
                    
                  fluidRow(
                    box(title = 'Predicting leaky vs. active genes - zTPM expression standardisation', status = 'info', solidHeader = T, width = 12,
                        h4('Active versus. leaky genes in a given cell type have been identified with the zFPKM algorithm, which transforms gene expression values into z-scores (zTPM) ([2]. Compared to the TPM unit, the zTPM unit indicates more accurately of how highly expressed a gene of interest is, because the zTPM unit records the expression of a gene relative to the overall pattern of gene expression in that sample.'),
                        br(),
                        h4('Based on ENCODE ChIP-seq and RNAseq data, cell-specific zTPM thresholds have been proposed that can be used to classify genes as likely active versus likely leaky. In HUVECs, this threshold was -2.38 zTPM [2]. We have applied the zFPKM algorithm to each sample in our compendium.'),
                        br(),
                        h4('The left panel shows how frequently a gene of interest passed the -2.38 zTPM active expression threshold.'),
                        br(),
                        h4('The right panel shows the zTPM expression range for their gene, whereby each datapoint represents an individual sample (samples that were unimodal or required 3-components to fit a GMM were excluded).')
                        )
                    ),
                  
                  ### zTPM plots
                  fluidRow(
                    box(title = textOutput('zDet.gene'), status = 'info', solidHeader = T,
                               plotOutput('Z.plot'),
                        ),
                     box(title = textOutput('zRange.gene'), status = 'info', solidHeader = T,
                               plotOutput('Z.range'))
                    ),
                  br(),
                  br(),
                  
                  ### Summary table
                  fluidRow(
                    box(title = 'Summary by cell type', status = 'warning', solidHeader = T, width = 12,
                    gt_output('cell.table'))
                  ),
                  
                  # Table Download Button
                  downloadButton("download.cell.table", "Download Table (.csv)",
                                 style = "color: #fff; background-color: #27ae60; border-color: #fff;padding: 5px 14px 5px 14px; margin: 0px 5px 20px 5px; "),
                  br(),
                    
                  fluidRow(
                    box(title = 'References', status = 'primary', solidHeader = T, width = 12,
                        p('1. Hebenstreit, D., M. Fang, M. Gu, V. Charoensawan, A. van Oudenaarden and S. A. Teichmann (2011). RNA sequencing reveals two major classes of gene expression levels in metazoan cells. Mol Syst Biol 7: 497. PMID 21654674'),
                        p('2. Hart, T., H. K. Komori, S. LaMere, K. Podshivalova and D. R. Salomon (2013). Finding the active genes in deep RNA-seq gene expression studies. BMC Genomics 14: 778. PMID 24215113'),
                        p('3. George, N. I. and C. W. Chang (2014). "DAFS: a data-adaptive flag method for RNA-sequencing data to differentiate genes with low and high expression." BMC Bioinformatics 15: 92. PMID: 24685233'),
                        p('4. Thompson, A., May, M. R., Moore, B. R. & Kopp, A. (2020). A hierarchical Bayesian mixture model for inferring the expression state of genes in transcriptomes. Proc Natl Acad Sci U S A 117: 19339. PMID: 32709743')
                    )
                    )
                  )
                )
              )
      )
    )
  )

                     
################################################################################
#########SERVER#################################################################

server <- function(input, output, session) {
  
  # Fetch the custom TPM threshold
  custom_TPM_thresh <- reactive({ # If we do: 'eventReactive(input$GeneSearch,{', it won't be able to update the value unless the gene is also updated
    thresh <- as.numeric(input$CustomTPM)
  })
  
  # Filter CleanBulk according to the threshold
  CleanBulk_custom_thresh <- eventReactive(custom_TPM_thresh(), { # Every time custom_TPM_thresh() changes, it changes
    df <- CleanBulk %>% 
      mutate(custom_tpm_thr = ifelse(TPM > custom_TPM_thresh(),
                                      "Expressed",
                                      "Not expressed"))
    return(df)
  })
  
  # Handling genes with multiple ENSG IDs per GENE SYMBOL that give inaccurate numbers
  dup_genes_text_pre <- reactive({
    ifelse(toupper(str_trim(input$GeneSearch)) %in% dup_genes,
           'As >1 entry exists for this gene in the Ensembl reference genome,
           this application cannot provide correct classification prediction
           for this gene.',
           '') 
  })
  
  # Only show warning if problematic gene is input, also, declare bool for flow control
  output$dup_genes_text <- renderText(dup_genes_text_pre())
  
  allow_proceed <- reactive({
    ifelse(toupper(str_trim(input$GeneSearch)) %in% dup_genes,
           FALSE,
           TRUE) 
  })
  
  # App will only continue processing if allow_proceed==TRUE
  Gene <- eventReactive(input$GeneSearch,{
    req(allow_proceed())
    input$GeneSearch %>% str_trim() %>% toupper()
  })
  
  
  # Switch tabs from home to results when a gene is searched
  observeEvent(input$GeneSearch,{
    newtab <- switch(input$tabs,
                     'Home_tab' = 'Analysis_tab',
                     'Analysis_tab' = 'Analysis_tab')
    updateTabItems(session = session, inputId = 'tabs', selected = newtab)},
    ignoreInit = TRUE
  )  

#######################
#######################
    
    
  ### Data filtering, plots and tables
    
  # For TPM results, filter CleanBulk by the input gene
  Det <- eventReactive(input$GeneSearch,
                       {CleanBulk %>% filter(Gene.Name == Gene())})
  
  Det.custom <- reactive({ # eventReactive(CleanBulk_custom_thresh(), {
    CleanBulk_custom_thresh() %>% filter(Gene.Name == Gene())})
  
  # For zTPM results, filter CleanBulk by the input gene and keep only bimodal datasets
  zDet <- eventReactive(input$GeneSearch, 
                        {Det() %>% filter(!run %in% unimodal)})
  
  # For GMM results, filter CleanEM by the input gene
  EMDet <- eventReactive(input$GeneSearch,
                         {CleanEM %>% filter(Gene.Name == Gene())})
  
  # Summary table creation
  CellTab <- reactive( # eventReactive(input$GeneSearch,
                           {
                             ### TPM
                             
                             # Count total number of datasets per cell type
                             # that contain the queried gene
                             tpm.counts.total <- Det() %>% 
                               group_by(cell.type) %>%
                               summarise(count = n())
                             
                             # Count total number of datasets per cell type
                             # that contain the queried gene with TPM >0 
                             tpm.counts.detected <- Det() %>%
                               filter(TPM > 0) %>%
                               group_by(cell.type) %>%
                               summarise(filt_count = n(),
                                         Median_TPM = (median(TPM) %>% round(1)))
                             
                             # Join both dataframes and get % of datasets with
                             # the queried gene with TPM > 0
                             tpm.counts.final <- left_join(tpm.counts.total, tpm.counts.detected) %>%
                               mutate(Detection_frequency = (((filt_count/count)*100) %>% round(1))) %>%
                               select(cell.type, count, filt_count, Median_TPM, Detection_frequency)
                             
                             ### TPM Custom threshold
                             
                             # Count total number of datasets per cell type
                             # that contain the queried gene
                             tpmcustom.counts.total <- Det() %>% 
                               group_by(cell.type) %>%
                               summarise(custom_count = n())
                             
                             # Count total number of datasets per cell type
                             # that contain the queried gene with TPM >0 
                             tpmcustom.counts.detected <- Det() %>%
                               filter(TPM > custom_TPM_thresh()) %>%
                               group_by(cell.type) %>%
                               summarise(custom_filt_count = n(),
                                         custom_Median_TPM = (median(TPM) %>% round(1)))
                             
                             # Join both dataframes and get % of datasets with
                             # the queried gene with TPM > 0
                             tpmcustom.counts.final <- left_join(tpmcustom.counts.total, tpmcustom.counts.detected) %>%
                               mutate(custom_Detection_frequency = (((custom_filt_count/custom_count)*100) %>% round(1))) %>%
                               select(cell.type, custom_count, custom_filt_count, custom_Median_TPM, custom_Detection_frequency)
                             
                             ### zTPM
                             
                             # Count total number of datasets per cell type
                             # that contain the queried gene
                             ztpm.counts.total <- zDet() %>%
                               group_by(cell.type) %>%
                               summarise(z_count = n())
                             
                             # Count total number of datasets per cell type
                             # that contain the queried gene with zTPM over threshold
                             ztpm.counts.thresh <- zDet() %>%
                               filter(zTPM_threshold == 'Expressed') %>%
                               group_by(cell.type) %>%
                               summarise(zfilt_count = n())
                             
                             # Median zTPM
                             ztpm.median <- zDet() %>%
                               filter(TPM > 0) %>% 
                               group_by(cell.type) %>%
                               summarise(Median_zTPM = (median(zTPM) %>% round(3)))
                             
                             # Join all zTPM tables, get frequency of datasets
                             # with the gene zTPM over threshold
                             ztpm.final <- left_join(ztpm.counts.total, ztpm.counts.thresh) %>%  
                               left_join(ztpm.median) %>%
                               mutate(Frequency_above_zTPM = (((zfilt_count/z_count)*100) %>% round(1))) %>%
                               select(cell.type, z_count, zfilt_count, Frequency_above_zTPM, Median_zTPM)
                             
                             ### GMM
                             
                             # Count total number of datasets per cell type
                             # that contain the queried gene
                             em.counts <- EMDet() %>% 
                               group_by(cell.type) %>%
                               summarise(e_count = n())
                             
                             # Count datasets by cluster category
                             em.active <- EMDet() %>% filter(cluster == 'Active') %>% group_by(cell.type) %>% summarise('Act' = n()) 
                             em.undet <- EMDet() %>% filter(cluster == 'Undetermined') %>% group_by(cell.type) %>% summarise('Und' = n()) 
                             em.leaky <- EMDet() %>% filter(cluster == 'Leaky') %>% group_by(cell.type) %>% summarise('Leak' = n())
                             
                             # Assemble the final GMM table
                             gmm.final <- left_join(em.counts, em.active) %>% 
                               left_join(.,em.undet) %>% 
                               left_join(., em.leaky) %>%
                               mutate(Active = ((Act/e_count)*100) %>% round(1),
                                      Undetermined = ((Und/e_count)*100) %>% round(1),
                                      Leaky = ((Leak/e_count)*100) %>% round(1)) %>%
                               select(cell.type, e_count, Act, Active, Undetermined, Leaky)  %>%
                               rename(., "Active %" = Active) %>%
                               rename(., "Undetermined %" = Undetermined) %>%
                               rename(., "Leaky %" = Leaky) %>%
                               rename(., "Active" = Act)
                             
                             # Assemble the final results table
                             results.final <- left_join(tpm.counts.final, tpmcustom.counts.final) %>% 
                               left_join(gmm.final) %>%
                               left_join(ztpm.final) %>%
                               replace(is.na(.), 0) 
                             
                             return(results.final)
                           })
  
  ### >0TPM Detection plot ###
  ############################
  
  Perc_by_cellltype_0TPM <- eventReactive(input$GeneSearch,
                                          {
                                            # Create empty dataframe with all cell types (for correct plotting)
                                            TPM_zero_counts <- data.frame(cell.type = as.factor(c("HUVEC", "HDMEC", "Lung EC", "Brain EC", "Retina EC")))
                                            
                                            # Count n TPM > 0
                                            TPM_expressed <- Det() %>%
                                              filter(Expressed_0_TPM == "Expressed") %>%
                                              group_by(cell.type) %>%
                                              count() %>%
                                              rename(Detected = n)
                                            
                                            # Count n datasets TPM == 0
                                            TPM_notexpressed <- Det() %>%
                                              filter(Expressed_0_TPM == "Not expressed") %>%
                                              group_by(cell.type) %>%
                                              count() %>%
                                              rename(Not_detected = n)
                                            
                                            # Join TPM_zero_counts and TPM_expressed
                                            TPM_expressed <- left_join(TPM_zero_counts, TPM_expressed, by = "cell.type") %>%
                                              replace_na(., list(Detected = 0))
                                            
                                            # Join TPM_zero_counts and TPM_notexpressed
                                            TMP_notexpressed <- left_join(TPM_zero_counts, TPM_notexpressed, by = "cell.type") %>%
                                              replace_na(., list(Not_detected = 0))
                                            
                                            # Join expressed and not expressed, and get dataset expression percentage
                                            TPM_all <- left_join(TPM_expressed, TMP_notexpressed, "cell.type") %>%
                                              mutate(expr_perc = round((Detected/(Detected+Not_detected)*100), 0)) %>%
                                              pivot_longer(cols = 2:3, names_to = "Expression", values_to = "Count") 
                                            
                                            return(TPM_all)
                                          })
  
      output$Det.gene <- renderText({paste(Gene(), 'detection frequency - ', 
                                          ((((Det() %>% filter(TPM > 0) %>% count())/(Det() %>% count()))*100) %>% 
                                               round(digits = 1)), '%')})
      
      renderDetPlot  <- function(){
        ggplot(Perc_by_cellltype_0TPM(), 
               aes(x = factor(cell.type, levels = c("HUVEC", "HDMEC", "Lung EC", "Brain EC", "Retina EC")),
                   y = Count,
                   fill = factor(Expression, levels = c('Not_detected', 'Detected')))) +
          geom_bar(position = position_stack(), stat = "identity") +
          geom_text(aes(label=ifelse(!is.na(expr_perc),paste0(expr_perc, "%"), "")), size = 6, y=-5) +
          scale_fill_colorblind(drop = FALSE) +
          theme_bw() +
          theme(text = element_text(size = 20)) +
          theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1)) +
          theme(legend.text = element_text(size = 20)) +
          theme(legend.position = c(0.70, 0.80)) +
          theme(axis.line = element_line(colour = "black")) +
          theme(axis.text.x = element_text(color = "black")) +
          theme(axis.text.y = element_text(color = "black")) +
          theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1)) +
          labs(fill = 'Gene > 0 TPM') +
          ylab('Number of samples') +
          xlab('') +
          ylim(-3, 160)
      }
      
      output$Det.plot <- renderPlot(renderDetPlot())
  
       ### TPM Expression range
       output$Range.gene <- renderText({paste(Gene(), 'expression range')})
       
       renderRangePlot <- function(){
         ggplot(data = Det() %>% filter(TPM > 0),
                aes(x = cell.type, y = TPM)) +
           geom_boxplot(aes(fill = factor(Expressed_0_TPM, levels = c('Not expressed', 'Expressed'))),
                        width = 0.4, lwd = 1.0, outlier.shape = NA, coef = 0, alpha = 0.1) +
           geom_point(alpha = 0.6, colour = 'darkorange2', size = 1.5, position = position_jitter(width = 0.15, height = 0.0, seed = 1)) +
           theme_bw() +
           theme(text = element_text(size = 20)) +
           theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1)) +
           scale_color_colorblind() +
           scale_x_discrete(drop=FALSE) +
           theme(axis.text.x = element_text(color = "black")) +
           theme(axis.text.y = element_text(color = "black")) +
           scale_y_log10(breaks = c(0, 0.01, 0.1, 1, 10, 100, 1000), labels = c(0, 0.01, 0.1, 1, 10, 100, 1000)) +
           annotation_logticks(sides = 'l') +
           xlab('') +
           labs(fill = '') +
           geom_hline(yintercept = 1000, colour = 'black', size = 1, alpha = 0.0) +
           geom_hline(yintercept = 0.01, colour = 'black', size = 1, alpha = 0.0) +
           theme(legend.position = "none")
       }
       
       output$Range.plot <- renderPlot(renderRangePlot())
       
       
       ### Custom TPM Detection plot ###
       #################################

       Perc_by_cellltype_customTPM <- reactive( # eventReactive(Det.custom(),
                                               {
                                                 # Create empty dataframe with all cell types (for correct plotting)
                                                 TPM_custom_counts <- data.frame(cell.type = as.factor(c("HUVEC", "HDMEC", "Lung EC", "Brain EC", "Retina EC")))
                                                 
                                                 # Count n TPM > threshold
                                                 TPM_expressed_custom <- Det.custom() %>%
                                                   filter(custom_tpm_thr == "Expressed") %>%
                                                   group_by(cell.type) %>%
                                                   count() %>%
                                                   rename(Detected = n)

                                                 # Count n datasets TPM == 0
                                                 TPM_notexpressed_custom <- Det.custom() %>%
                                                   filter(custom_tpm_thr == "Not expressed") %>%
                                                   group_by(cell.type) %>%
                                                   count() %>%
                                                   rename(Not_detected = n)

                                                 # Join TPM_custom_counts and TPM_expressed_custom
                                                 TPM_expressed_custom <- left_join(TPM_custom_counts, TPM_expressed_custom, by = "cell.type") %>%
                                                   replace_na(., list(Detected = 0))

                                                 # Join TPM_custom_counts and TPM_notexpressed_custom
                                                 TMP_notexpressed_custom <- left_join(TPM_custom_counts, TPM_notexpressed_custom, by = "cell.type") %>%
                                                   replace_na(., list(Not_detected = 0))

                                                 # Join expressed and not expressed, and get dataset expression percentage
                                                 TPM_all_custom <- left_join(TPM_expressed_custom, TMP_notexpressed_custom, "cell.type") %>%
                                                   mutate(expr_perc = round((Detected/(Detected+Not_detected)*100), 0)) %>%
                                                   pivot_longer(cols = 2:3, names_to = "Expression", values_to = "Count")

                                                 return(TPM_all_custom)
                                               })

       output$Det.gene.custom <- renderText({paste0(Gene(), " detection frequency over ", custom_TPM_thresh(), " TPM - ",
                                            ((((Det.custom() %>% filter(TPM > custom_TPM_thresh()) %>% count())/(Det.custom() %>% count()))*100) %>%
                                               round(digits = 1)), '%')})

       renderDetPlotCustom  <- function(){
         ggplot(Perc_by_cellltype_customTPM(),
                aes(x = factor(cell.type, levels = c("HUVEC", "HDMEC", "Lung EC", "Brain EC", "Retina EC")),
                    y = Count,
                    fill = factor(Expression, levels = c('Not_detected', 'Detected')))) +
           geom_bar(position = position_stack(), stat = "identity") +
           geom_text(aes(label=ifelse(!is.na(expr_perc),paste0(expr_perc, "%"), "")), size = 6, y=-5) +
           scale_fill_colorblind(drop = FALSE) +
           theme_bw() +
           theme(text = element_text(size = 20)) +
           theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1)) +
           theme(legend.text = element_text(size = 20)) +
           theme(legend.position = c(0.70, 0.80)) +
           theme(axis.line = element_line(colour = "black")) +
           theme(axis.text.x = element_text(color = "black")) +
           theme(axis.text.y = element_text(color = "black")) +
           theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1)) +
           labs(fill = paste0("Gene > ", custom_TPM_thresh(), " TPM")) +
           ylab('Number of samples') +
           xlab('') +
           ylim(-3, 160)
       }

       output$Det.plot.custom <- renderPlot({renderDetPlotCustom()})

       ### TPM Expression range (custom threshold)
       output$Range.gene.custom <- renderText({paste(Gene(), 'expression range')})
       
       renderRangePlotCustom <- reactive({
         ggplot(data = Det.custom() %>% filter(TPM > 0),
                aes(x = cell.type, y = TPM)) +
           geom_boxplot(width = 0.4, lwd = 1.0, outlier.shape = NA, coef = 0, alpha = 0.1) + # aes(fill = factor(custom_tpm_thr, levels = c('Not expressed', 'Expressed'))),
           geom_point(alpha = 0.6, colour = 'darkorange2', size = 1.5, position = position_jitter(width = 0.15, height = 0.0, seed = 1)) +
           geom_hline(yintercept = custom_TPM_thresh(), color = "dark red", linetype = "longdash") +
           theme_bw() +
           theme(text = element_text(size = 20)) +
           theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1)) +
           scale_color_colorblind() +
           scale_x_discrete(drop=FALSE) +
           theme(axis.text.x = element_text(color = "black")) +
           theme(axis.text.y = element_text(color = "black")) +
           scale_y_log10(breaks = c(0, 0.01, 0.1, 1, 10, 100, 1000), labels = c(0, 0.01, 0.1, 1, 10, 100, 1000)) +
           annotation_logticks(sides = 'l') +
           xlab('') +
           labs(fill = '') +
           geom_hline(yintercept = 1000, colour = 'black', size = 1, alpha = 0.0) +
           geom_hline(yintercept = 0.01, colour = 'black', size = 1, alpha = 0.0) +
           theme(legend.position = "none")
       })
       
         output$Range.plot.custom <- renderPlot(renderRangePlotCustom())
       
       #### zTPM ####
       ##############
       Perc_by_cellltype_zTPMThresh <- eventReactive(input$GeneSearch,
                                                     {
                                                       zTPM_zero_counts <- data.frame(cell.type = as.factor(c("HUVEC", "HDMEC", "Lung EC", "Brain EC", "Retina EC")))
                                                       
                                                       perc1 <- zDet() %>%
                                                         filter(zTPM_threshold == "Expressed") %>%
                                                         group_by(cell.type) %>%
                                                         count() %>%
                                                         rename(Detected = n)
                                                       
                                                       perc2 <- zDet() %>%
                                                         filter(zTPM_threshold == "Not expressed") %>%
                                                         group_by(cell.type) %>%
                                                         count() %>%
                                                         rename(Not_detected = n)
                                                       
                                                       # Left join (zero_counts has all cell.types),
                                                       perc1 <- left_join(zTPM_zero_counts, perc1, by = "cell.type") %>%
                                                         replace_na(., list(Detected = 0))
                                                       
                                                       perc2 <- left_join(zTPM_zero_counts, perc2, by = "cell.type") %>%
                                                         replace_na(., list(Not_detected = 0))
                                                       
                                                       perc3 <- left_join(perc1, perc2, "cell.type") %>%
                                                         mutate(expr_perc = round((Detected/(Detected+Not_detected)*100), 0)) %>%
                                                         pivot_longer(cols = 2:3, names_to = "Expression", values_to = "Count")                                      
                                                     })
     
     # If all counts are 0 or not, will determine the script to generate the plot
     # with this, we create a reactive boolean for that.
     bool_all_zero <- eventReactive(input$GeneSearch,
                                    {
                                      temp <- Perc_by_cellltype_zTPMThresh()
                                      if(all(temp$Count==0)){zero_bool=TRUE}else{zero_bool=FALSE}
                                      return(zero_bool)
                                    })
     
     # To be able to print to console
     observeEvent(input$GeneSearch, {
       print(paste0("The value of bool is: ", bool_all_zero()))
     })
     
     ### zTPM Threshold detection
      output$zDet.gene <- renderText({paste(Gene(), 'detection above zTPM threshold - ',
                                          (((zDet() %>% filter(zTPM_threshold == 'Expressed') %>% count()) /
                                              (zDet() %>% count()))*100) %>% round(digits = 1),
                                          '%')})
      
      # For this plot, we need to run one option or another depending if all values are 0 or not.
      # This is determined by the bool_all_zero() eventReactive we generated before, which returns a boolean (TRUE or FALSE).
      # This selection of plot needs to be done within an eventReactive as well, and it's stored in renderZTPMPlot
      renderZTPMPlot <- eventReactive(input$GeneSearch,
                                     {
                                       if(bool_all_zero()==FALSE){
                                          Zplot <- ggplot(Perc_by_cellltype_zTPMThresh(), aes(x = factor(cell.type, levels = c("HUVEC", "HDMEC", "Lung EC", "Brain EC", "Retina EC")), 
                                                                                              y = Count,
                                                                                              fill = factor(Expression, levels = c('Not_detected', 'Detected')))) +
                                                                        geom_bar(position = position_stack(), stat = "identity") +
                                                                        scale_x_discrete(drop=FALSE) +
                                                                        geom_text(aes(label=ifelse(!is.na(expr_perc),paste0(expr_perc, "%"), "")), size = 6, y=-5) +
                                                                        scale_fill_colorblind(drop = FALSE) +
                                                                        theme_bw() +
                                                                        theme(text = element_text(size = 20)) +
                                                                        theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1)) +
                                                                        theme(legend.text = element_text(size = 20)) +
                                                                        theme(legend.position = c(0.70, 0.80)) +
                                                                        theme(axis.line = element_line(colour = "black")) +
                                                                        theme(axis.text.x = element_text(color = "black")) +
                                                                        theme(axis.text.y = element_text(color = "black")) +
                                                                        labs(fill = 'Gene > -2.38 zTPM') +
                                                                        ylab('Number of samples') +
                                                                        xlab('') +
                                                                        ylim(-3, 160)
                                        } else {
                                          Zplot <- ggplot(Perc_by_cellltype_zTPMThresh(), aes(x = factor(cell.type, levels = c("HUVEC", "HDMEC", "Lung EC", "Brain EC", "Retina EC")),
                                                                                                         y = Count,
                                                                                                         fill = factor(Expression, levels = c('Not_detected', 'Detected')))) +
                                                                        geom_bar(position = position_stack(), stat = "identity") +
                                                                        geom_text(aes(label=ifelse(!is.na(expr_perc),paste0(expr_perc, "%"), "")), size = 6, y=-5) +
                                                                        scale_fill_colorblind(drop = FALSE) +
                                                                        theme_bw() +
                                                                        theme(text = element_text(size = 20)) +
                                                                        theme(legend.text = element_text(size = 20)) +
                                                                        theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1)) +
                                                                        theme(legend.position = c(0.70, 0.80)) +
                                                                        theme(axis.line = element_line(colour = "black")) +
                                                                        theme(axis.text.x = element_text(color = "black")) +
                                                                        theme(axis.text.y = element_text(color = "black")) +
                                                                        labs(fill = 'Gene > -2.38 zTPM') +
                                                                        ylab('Number of samples') +
                                                                        xlab('') +
                                                                        ylim(-3, 160)
                                        }
                                       return(Zplot)
                                     })
      
    output$Z.plot <- renderPlot(renderZTPMPlot())
    
    
    ### zTPM expression range
    output$zRange.gene <- renderText({paste(Gene(), 'standardised expression range')})
    
    renderZRangePlot <- function(){
      ggplot(data = zDet() %>% filter(TPM > 0), aes(x = cell.type, y = zTPM))+
        geom_boxplot(aes(fill = factor(Expressed_0_TPM, levels = c('Not expressed', 'Expressed'))),
                     width = 0.4, lwd = 1.0, outlier.shape = NA, coef = 0, alpha = 0.1) +
        geom_jitter(width = 0.15, height = 0.0, alpha = 0.6, colour = 'darkorange2', size = 1.5) +
        geom_hline(yintercept = -2.38, colour = 'dark red', linetype = 'longdash', size = 1.) +
        ylim(-5, 5) +
        theme_bw() +
        theme(text = element_text(size = 20)) +
        theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1)) +
        scale_color_colorblind() +
        scale_x_discrete(drop=FALSE) +
        xlab('') +
        labs(fill = ' ') +
        theme(axis.text.x = element_text(color = "black")) +
        theme(axis.text.y = element_text(color = "black")) +
        geom_hline(yintercept = 4.0, colour = 'black', size = 1, alpha = 0.0) +
        geom_hline(yintercept = -4.0, colour = 'black', size = 1, alpha = 0.0) +
        theme(legend.position = "none")
    }
    
    output$Z.range <- renderPlot(renderZRangePlot())
    


    #### GMM plot ####
    ##################
    #Written by Guille as an alternative way of generating the GMM plot, but adding percentages
    
    output$GMM.gene <- renderText({paste(Gene(), 'GMM classification - Active:',
                                         (((EMDet() %>% filter(cluster == 'Active') %>% count()) / 
                                            (EMDet() %>% count()))*100) %>% round(digits = 1), 
                                         '%', "| Leaky:", 
                                         (((EMDet() %>% filter(cluster == 'Leaky') %>% count()) / 
                                                             (EMDet() %>% count()))*100) %>% round(digits = 1),
                                         "%")})

    # Will contain ALL cell types, allowing left join to avoid missing cell types
    GMM_zero_counts <- data.frame(cell.type = as.factor(c("HUVEC", "HDMEC", "Lung EC", "Brain EC", "Retina EC")))
    
    GMM_data_generator <- eventReactive(input$GeneSearch,
                                                  {
                                                    GMM_class_active <- EMDet() %>%
                                                      filter(cluster == "Active") %>%
                                                      group_by(cell.type) %>%
                                                      count() %>%
                                                      rename(Active = n)
                                                    
                                                    GMM_class_leaky <- EMDet() %>%
                                                      filter(cluster == "Leaky") %>%
                                                      group_by(cell.type) %>%
                                                      count() %>%
                                                      rename(Leaky = n)
                                                    
                                                    GMM_class_undet <- EMDet() %>%
                                                      filter(cluster == "Undetermined") %>%
                                                      group_by(cell.type) %>%
                                                      count() %>%
                                                      rename(Undetermined = n)
                                                    
                                                    # Left join (zero_counts has all cell.types)
                                                    GMM_class_active <- left_join(GMM_zero_counts, GMM_class_active, by = "cell.type") %>%
                                                      replace_na(., list(Active = 0))
                                                    
                                                    GMM_class_leaky <- left_join(GMM_zero_counts, GMM_class_leaky, by = "cell.type") %>%
                                                      replace_na(., list(Leaky = 0))
                                                    
                                                    GMM_class_undet <- left_join(GMM_zero_counts, GMM_class_undet, by = "cell.type") %>%
                                                      replace_na(., list(Undetermined = 0))
                                                    
                                                    #Join the 3 classifications, replace NAs with 0s, get the % of Active and pivot longer for plotting
                                                    GMM_class_perc <- left_join(GMM_class_active, GMM_class_leaky, by = "cell.type") %>%
                                                      left_join(., GMM_class_undet, by = "cell.type") %>%
                                                      mutate(Perc_active = round((Active/(Active+Leaky+Undetermined)*100), 0)) %>%
                                                      pivot_longer(cols = 2:4, names_to = "GMM_class", values_to = "Count")
      
      return(GMM_class_perc)
      
      })
    
    observeEvent(GMM_data_generator(), {
      print(paste0("GMM data:", print(GMM_data_generator())))
    })
    
    renderGMMPlot <- function(){
      ggplot(data = GMM_data_generator(), aes(x = factor(cell.type, levels = c("HUVEC", "HDMEC", "Lung EC", "Brain EC", "Retina EC")),
                                              y = Count,
                                              fill = factor(GMM_class, levels = c('Leaky', 'Undetermined', 'Active')))) +
        geom_bar(position = position_stack(), stat = "identity") +
        geom_text(aes(label=ifelse(!is.na(Perc_active), paste0(Perc_active, "%"), "")), size = 6, y=-5) +
        scale_fill_manual(values = c("black", "light blue", "#E69F00"), drop = FALSE) +
        theme_bw() +
        theme(text = element_text(size = 20)) +
        theme(legend.text = element_text(size = 20)) +
        theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1)) +
        theme(legend.position = c(0.70, 0.80)) +
        theme(axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(color = "black")) +
        theme(axis.text.y = element_text(color = "black")) +
        labs(fill = 'Classification') +
        ylab('Sample number with gene > 0 TPM') +
        xlab('') +
        ylim(-3, 160)
    }
    
    output$GMM.plot <- renderPlot(renderGMMPlot())
    


    #### Cell specific table ####
    #############################                                              
    output$cell.table <- render_gt(gt(CellTab()) %>%
                                     tab_spanner('zTPM analysis', columns = c('z_count', 'zfilt_count', 'Frequency_above_zTPM', 'Median_zTPM')) %>%
                                     tab_spanner('GMM classification', columns = c('e_count', 'Active', 'Active %', 'Undetermined %', 'Leaky %')) %>%
                                     tab_spanner('TPM analysis', columns = c('count', 'filt_count',  'Detection_frequency', 'Median_TPM')) %>%
                                     tab_spanner('Custom TPM threshold analysis', columns = c('custom_count', 'custom_filt_count', 'custom_Detection_frequency', 'custom_Median_TPM')) %>%
                                     cols_align(align = "center", columns = everything()) %>%
                                     tab_footnote(footnote = "Median calculated  only using samples where gene was detected > 0 TPM",
                                       locations = cells_column_labels(
                                         columns = c(Median_TPM, Median_zTPM)
                                       )) %>%
                                     tab_footnote(footnote = "Median calculated  only using samples where gene TPM > custom TPM threshold",
                                                  locations = cells_column_labels(
                                                    columns = c(custom_Median_TPM)
                                                  )) %>%
                                     cols_label(cell.type = '', count = 'n', filt_count = '> 0 TPM', Detection_frequency = 'Detection %', Median_TPM = 'Median TPM',
                                                z_count = 'n', zfilt_count = '> zTPM threshold', Frequency_above_zTPM = 'Detection %', Median_zTPM = 'Median zTPM',
                                                e_count = 'n', custom_count = "n", custom_filt_count = paste0("> ", custom_TPM_thresh(), " TPM"), custom_Detection_frequency = "Detection %", custom_Median_TPM = "Median TPM"))
    
    # Downloadable csv of selected dataset ----
    output$download.cell.table <- downloadHandler(
      filename = function() {
        paste0(input$GeneSearch %>% toupper(), "_BulkECexplorer_summary_table.csv")
      },
      content = function(file) {
        write.csv(
          CellTab() %>%
            select("cell.type", "count", "filt_count", "Detection_frequency", "Median_TPM",
                   "e_count", "Active", "Active %", "Undetermined %", "Leaky %",
                   "z_count", "zfilt_count", "Frequency_above_zTPM", "Median_zTPM",
                   "custom_count", "custom_filt_count", "custom_Median_TPM", "custom_Detection_frequency") %>%
            rename_with(~ c("cell_type", "n_TPM", "over_0_TPM", "Detection_pct_TPM", "Median_TPM",
                            "n_GMM", "Active_GMM", "Active_pct_GMM", "Undetermined_pct_GMM", "Leaky_pct_GMM",
                            "n_zTPM", "over_zTPM_thresh", "Detection_pct_zTPM", "Median_zTPM",
                            "custom_TPMthreshold_count", "over_custom_TPMthreshold", "over_custom_TPMthreshold_Median", "custom_TPMthreshold_Detection_freq")),
          file, row.names = FALSE)
      }
    )
    
    ########Download Plots#####
    ###########################
    
    output$download.plots.png <- downloadHandler(
      filename = function() {
        paste0(input$GeneSearch %>% toupper(), "_BulkECexplorer_plots_PNG.zip")
      },
      content = function(file) {
        
        # Set temporary working directory
        original_wd <- setwd(tempdir()) # This stores the original wd in the variable, but changes the current dir to the new one
        on.exit(setwd(original_wd)) # On exit, return to the original wd
        
        # Save the images in the temp dir
        ggsave(paste0(input$GeneSearch %>% toupper(), "_detection_freq.png"), plot = renderDetPlot(), device = "png")
        ggsave(paste0(input$GeneSearch %>% toupper(), "_expression_range.png"), plot = renderRangePlot(), device = "png")
        ggsave(paste0(input$GeneSearch %>% toupper(), "_gmm_classification.png"), plot = renderGMMPlot(), device = "png")
        ggsave(paste0(input$GeneSearch %>% toupper(), "_detection_zTPM_theshold.png"), plot = renderZTPMPlot(), device = "png")
        ggsave(paste0(input$GeneSearch %>% toupper(), "_zTPM_expression_range.png"), plot = renderZRangePlot(), device = "png")
        ggsave(paste0(input$GeneSearch %>% toupper(), "_frequency_TPM_threshold.png"), plot = renderDetPlotCustom(), device = "png")
        ggsave(paste0(input$GeneSearch %>% toupper(), "_expression_range_TPM_threshold.png"), plot = renderRangePlotCustom(), device = "png")
        
        # zip the images
        zip(file, c(paste0(input$GeneSearch %>% toupper(), "_detection_freq.png"), 
                    paste0(input$GeneSearch %>% toupper(), "_expression_range.png"), 
                    paste0(input$GeneSearch %>% toupper(), "_gmm_classification.png"), 
                    paste0(input$GeneSearch %>% toupper(), "_detection_zTPM_theshold.png"),
                    paste0(input$GeneSearch %>% toupper(), "_zTPM_expression_range.png"),
                    paste0(input$GeneSearch %>% toupper(), "_frequency_TPM_threshold.png"),
                    paste0(input$GeneSearch %>% toupper(), "_expression_range_TPM_threshold.png"))
            )
      }
    )
    
    output$download.plots.tiff <- downloadHandler(
      filename = function() {
        paste0(input$GeneSearch %>% toupper(), "_BulkECexplorer_plots_TIFF.zip")
      },
      content = function(file) {
        
        original_wd <- setwd(tempdir())
        on.exit(setwd(original_wd))
        
        ggsave(paste0(input$GeneSearch %>% toupper(), "_detection_freq.tiff"), plot = renderDetPlot(), device = "tiff")
        ggsave(paste0(input$GeneSearch %>% toupper(), "_expression_range.tiff"), plot = renderRangePlot(), device = "tiff")
        ggsave(paste0(input$GeneSearch %>% toupper(), "_gmm_classification.tiff"), plot = renderGMMPlot(), device = "tiff")
        ggsave(paste0(input$GeneSearch %>% toupper(), "_detection_zTPM_theshold.tiff"), plot = renderZTPMPlot(), device = "tiff")
        ggsave(paste0(input$GeneSearch %>% toupper(), "_zTPM_expression_range.tiff"), plot = renderZRangePlot(), device = "tiff")
        ggsave(paste0(input$GeneSearch %>% toupper(), "_frequency_TPM_threshold.tiff"), plot = renderDetPlotCustom(), device = "tiff")
        ggsave(paste0(input$GeneSearch %>% toupper(), "_expression_range_TPM_threshold.tiff"), plot = renderRangePlotCustom(), device = "tiff")
        
        zip(file, c(paste0(input$GeneSearch %>% toupper(), "_detection_freq.tiff"), 
                    paste0(input$GeneSearch %>% toupper(), "_expression_range.tiff"), 
                    paste0(input$GeneSearch %>% toupper(), "_gmm_classification.tiff"), 
                    paste0(input$GeneSearch %>% toupper(), "_detection_zTPM_theshold.tiff"),
                    paste0(input$GeneSearch %>% toupper(), "_zTPM_expression_range.tiff"),
                    paste0(input$GeneSearch %>% toupper(), "_frequency_TPM_threshold.tiff"),
                    paste0(input$GeneSearch %>% toupper(), "_expression_range_TPM_threshold.tiff"))
        )
      }
    )
    
    output$download.plots.pdf <- downloadHandler(
      filename = function() {
        paste0(input$GeneSearch %>% toupper(), "_BulkECexplorer_plots_PDF.zip")
      },
      content = function(file) {
        
        original_wd <- setwd(tempdir())
        on.exit(setwd(original_wd))
        
        ggsave(paste0(input$GeneSearch %>% toupper(), "_detection_freq.pdf"), plot = renderDetPlot(), device = "pdf")
        ggsave(paste0(input$GeneSearch %>% toupper(), "_expression_range.pdf"), plot = renderRangePlot(), device = "pdf")
        ggsave(paste0(input$GeneSearch %>% toupper(), "_gmm_classification.pdf"), plot = renderGMMPlot(), device = "pdf")
        ggsave(paste0(input$GeneSearch %>% toupper(), "_detection_zTPM_theshold.pdf"), plot = renderZTPMPlot(), device = "pdf")
        ggsave(paste0(input$GeneSearch %>% toupper(), "_zTPM_expression_range.pdf"), plot = renderZRangePlot(), device = "pdf")
        ggsave(paste0(input$GeneSearch %>% toupper(), "_frequency_TPM_threshold.pdf"), plot = renderDetPlotCustom(), device = "pdf")
        ggsave(paste0(input$GeneSearch %>% toupper(), "_expression_range_TPM_threshold.pdf"), plot = renderRangePlotCustom(), device = "pdf")
        
        zip(file, c(paste0(input$GeneSearch %>% toupper(), "_detection_freq.pdf"), 
                    paste0(input$GeneSearch %>% toupper(), "_expression_range.pdf"), 
                    paste0(input$GeneSearch %>% toupper(), "_gmm_classification.pdf"), 
                    paste0(input$GeneSearch %>% toupper(), "_detection_zTPM_theshold.pdf"),
                    paste0(input$GeneSearch %>% toupper(), "_zTPM_expression_range.pdf"),
                    paste0(input$GeneSearch %>% toupper(), "_frequency_TPM_threshold.pdf"),
                    paste0(input$GeneSearch %>% toupper(), "_expression_range_TPM_threshold.pdf"))
        )
      }
    )

}



# Run the application 
shinyApp(ui = ui, server = server)
