setwd("/data/ruoffcj/projects/drug_treatment/")
library(tidyverse)
library(GSVA)
library(readxl)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(survminer)
library(survival)
library(presto)
source("../survival_analysis/cox_regression.R")
source("source/final_scripts/drug_treatment_functions.R")
set.seed(42)

get_deconvolved_data <- function(cell_line){
  
  cancer_type <- ""
  
  if(cell_line == "A549"){
    cancer_type <- "lung_NSCLC_adenocarcinoma"
  } else if(cell_line == "K562"){
    cancer_type <- "acute_myeloid_leukaemia"
  } else {
    cancer_type <- "breast"
  }
  
  
  return(pred_TCGA_codefacs[[cancer_type]]$pred_stageIII[,,"Cancer"])
  
}

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"
plotDirectory <- "/data/ruoffcj/projects/drug_treatment/"

################################################################################
# Survival Analysis for each cell line in cancer type matched TCGA samples
################################################################################
cell_lines <- c("A549","K562","MCF7")

all_signatures <- readRDS(paste0(dataDirectory, "genesets/global_rac_signatures.rds"))

metric_to_use <- "OS"
# metric_to_use <- "PFI"

hazard_ratios <- c()
plots <- list()

curr_cell_line <- "A549"

for(curr_cell_line in cell_lines){
  cat(curr_cell_line, "\n")
  
  # Get TCGA sample count data  
  tcga_project <- get_tcga_project(curr_cell_line)
  
  # Get read count data
  tcga_matrix_vst <- get_deconvolved_data(curr_cell_line)
  
  # Score samples
  curr_signature <- all_signatures[grepl(curr_cell_line,names(all_signatures))]
  
  ssgsea_result <- gsva(tcga_matrix_vst,curr_signature,method='ssgsea')
  
  gene_set_scores_df <- as.data.frame(t(ssgsea_result)) %>%
    tibble::rownames_to_column("submitter_id")
  
  gene_set_scores_df$submitter_id <- sapply(gene_set_scores_df$submitter_id, function(x) gsub("\\.", "-",x))
  gene_set_scores_df$submitter_id <- sapply(gene_set_scores_df$submitter_id, function(x) substring(x, 1,12))
  
  # Get clinical data
  clinical_df <- get_clinical_data(tcga_project)
  
  #Run cox regression
  if("purity" %in% colnames(clinical_df)){
    covariates <- c("purity", "sex", "age")
  } else {
    covariates <- c("sex", "age")
  }
  
  
  cox_regression_info <- cox_regression(sample_survival_df = clinical_df,
                                        survival_data_id_col = "submitter_id",
                                        duration_str = paste0(metric_to_use, ".time"),
                                        feature_data_id_col = "submitter_id",
                                        model_covariates = covariates,
                                        status_str=metric_to_use,
                                        sample_features_df = gene_set_scores_df, 
                                        km_features_to_plot=names(curr_signature),
                                        low_high_percentiles = c(.5,.49))
  
  
  if(!all(is.na(cox_regression_info$km_df[[metric_to_use]]))){
    
    hazard_ratios <- append(hazard_ratios, cox_regression_info$regression_df$hazard_ratio)
    
    modform <- as.formula(paste0("Surv(",metric_to_use,".time, ",metric_to_use,") ~ ", names(curr_signature)))
    
    fit <- eval(substitute(survfit(modform, data = cox_regression_info$km_df), list(modform = modform)))
    
    
    if(metric_to_use == "OS"){
      plot_title <- paste0("\n",curr_cell_line, " (",tcga_project,")")
    } else {
      plot_title <- paste0("\n",curr_cell_line, " (",tcga_project,")")
    }
    
    
    cox_regression_info$km_df[[names(curr_signature)]] <- factor(cox_regression_info$km_df[[names(curr_signature)]], levels = c("High","Medium","Low"))
    
    p <- ggsurvplot(fit, data = cox_regression_info$km_df,
                    pval = T,
                    pval.size = 6,
                    legend.title="Expression Level:",
                    legend.labs= c("High","Low"),
                    font.title=c(28),
                    font.x=c(28),
                    font.y=c(28),
                    font.tickslab=c(20),
                    font.legned=c(50),
                    font.caption=c(30),
                    tables.theme = clean_theme())+
      xlab("")+
      ylab("")
    
    
    p$plot <- p$plot + 
      theme(legend.text = element_text(size = 24, color = "black", face = "bold"),
            legend.title = element_text(size = 28, color = "black", face = "bold"),
            legend.key.height = unit(1.5,"cm"),
            legend.key.width = unit(1.5,"cm"))
    
    
    p <- p + ggtitle(plot_title)
    
    plots <- append(plots, list(p$plot))
    
  }
}


################################################################################
# Add MCF7 plot
################################################################################

# Read in clincal data
clinical_data <- read_tsv(paste0(dataDirectory, "metabric_data/brca_metabric_clinical_data.tsv"))
patient_sample <- read.table(paste0(dataDirectory, "metabric_data/brca_metabric/data_clinical_patient.txt"),
                             fill=T)

colnames(patient_sample) <- patient_sample[1,]
patient_sample <- patient_sample[-1,]

temp <- matrix(unlist(patient_sample),nrow=nrow(patient_sample))
colnames(temp) <- colnames(patient_sample)
rownames(temp) <- rownames(patient_sample)

patient_sample <- temp

patient_sample <- as.data.frame(patient_sample) %>% 
  filter(CLAUDIN_SUBTYPE == "LumA")

# Read in expression data
microarray_data <- read.table(paste0(dataDirectory,"metabric_data/brca_metabric/data_mrna_illumina_microarray.txt"),
                              fill=T)

colnames(microarray_data) <- microarray_data[1,]
microarray_data <- microarray_data[-1,]

rownames(microarray_data) <- make.names(microarray_data[,1],unique=T)
microarray_data <- microarray_data[,-c(1,2)]

temp <- matrix(as.numeric(unlist(microarray_data)),nrow=nrow(microarray_data))
colnames(temp) <- colnames(microarray_data)
rownames(temp) <- rownames(microarray_data)

microarray_data <- temp

# Score expression data with RAC signature
curr_signature <- all_signatures["MCF7"]

ssgsea_result <- gsva(microarray_data, curr_signature, method='ssgsea')

gene_set_scores_df <- as.data.frame(t(ssgsea_result)) %>%
  tibble::rownames_to_column("PATIENT_ID")


patient_sample$OS_MONTHS <- as.numeric(patient_sample$OS_MONTHS)
patient_sample$OS_STATUS <- as.numeric(gsub(":[A-Z]*", "", patient_sample$OS_STATUS))
patient_sample$AGE_AT_DIAGNOSIS <- as.numeric(patient_sample$AGE_AT_DIAGNOSIS)

# Perform Cox regression
cox_regression_info <- cox_regression(sample_survival_df = patient_sample,
                                      survival_data_id_col = "PATIENT_ID",
                                      duration_str = "OS_MONTHS",
                                      feature_data_id_col = "PATIENT_ID",
                                      model_covariates = c("AGE_AT_DIAGNOSIS"),
                                      status_str="OS_STATUS",
                                      sample_features_df = gene_set_scores_df, 
                                      km_features_to_plot=names(curr_signature),
                                      low_high_percentiles = c(.5,.49))

hazard_ratios <- append(hazard_ratios, cox_regression_info$regression_df$hazard_ratio)

modform <- as.formula(paste0("Surv(OS_MONTHS, OS_STATUS) ~ ", names(curr_signature)))

fit <- eval(substitute(survfit(modform, data = cox_regression_info$km_df), list(modform = modform)))


cox_regression_info$km_df[[names(curr_signature)]] <- factor(cox_regression_info$km_df[[names(curr_signature)]], levels = c("High","Medium","Low"))

# Plot KM Plot
p <- ggsurvplot(fit, data = cox_regression_info$km_df,
                pval = T,
                pval.size = 6,
                legend.title="Expression Level:",
                legend.labs= c("High","Low"),
                font.title=c(28),
                font.x=c(28),
                font.y=c(28),
                font.tickslab=c(20),
                font.legned=c(50),
                font.caption=c(30),
                tables.theme = clean_theme())+
  xlab("Time")+
  ylab("Survival Probability")


p$plot <- p$plot + 
  theme(legend.text = element_text(size = 24, color = "black", face = "bold"),
        legend.title = element_text(size = 28, color = "black", face = "bold"),
        legend.key.height = unit(1.5,"cm"),
        legend.key.width = unit(1.5,"cm"))

p <- p+ggtitle("\nMCF7 (Metabric LumA Samples)")

# Add KM plot to plot list
plots <- append(plots, list(p$plot))

################################################################################
# Plot KM plots

figure <- ggarrange(plotlist = plots, ncol=3, common.legend = T, legend=c("right"))

p <- annotate_figure(figure, left = text_grob("Survival Probability", rot = 90, vjust = 1, size=35, face="bold"),
                     bottom = text_grob("Time", size=35, face="bold"),
                     top=text_grob(paste0("RAC Signatures Overall Survival in Patient Samples"), size=40, face="bold"))



png(paste0(plotDirectory,"final_figures/figure_4a.png"),
    width=30, height=10, units= "in", res = 300)

print(p)

dev.off()


