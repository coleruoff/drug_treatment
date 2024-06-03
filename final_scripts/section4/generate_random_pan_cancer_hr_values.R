setwd("/data/ruoffcj/projects/drug_treatment/")
# library(gt)
library(GSVA)
library(readxl)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(survminer)
library(survival)
library(patchwork)
source("/data/ruoffcj/projects/survival_analysis/cox_regression.R")
source("source/final_scripts/drug_treatment_functions.R")
set.seed(42)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"
plotDirectory <- "/data/ruoffcj/projects/drug_treatment/"

##################################################################################
# Get list of TCGA projects
##################################################################################

gdcprojects <- getGDCprojects()

tcga_projects <- gdcprojects[grepl("TCGA", gdcprojects$id),]$id

tcga_projects <- tcga_projects

#################################################################################
# Generate KM plot data and Cox hazard ratios for each TCGA project
#################################################################################
supercluster_signatures <- readRDS(paste0(dataDirectory, "genesets/rac_supercluster_signatures.rds"))

all_signatures <- supercluster_signatures
metric_to_use <- "OS"
hazard_ratio_df <- list()
plot_data <- list()

for(curr_project in tcga_projects){
  
  cat(curr_project, "\n")
  
  #Get read count data
  # tcga_matrix_vst <- get_read_count_data(curr_project)
  tcga_matrix_vst <- readRDS(paste0(dataDirectory, "../../../TCGA_processed_data/",curr_project,"_data.rds"))
  
  #Filter out genes with low/zero variance
  row_vars <- rowVars(tcga_matrix_vst)
  tcga_matrix_vst <- tcga_matrix_vst[which(row_vars > 1e-10),]
  
  # Score samples
  gsvapar <- ssgseaParam(tcga_matrix_vst, all_signatures)
  ssgsea_result <- gsva(gsvapar)
  
  gene_set_scores_df <- as.data.frame(t(ssgsea_result)) %>%
    tibble::rownames_to_column("submitter_id")
  
  gene_set_scores_df$submitter_id <- sapply(gene_set_scores_df$submitter_id, function(x) gsub("\\.", "-",x))
  gene_set_scores_df$submitter_id <- sapply(gene_set_scores_df$submitter_id, function(x) substring(x, 1,12))
  
  # Get clinical data
  clinical_df <- get_clinical_data(curr_project)
  
  # Randomize OS data
  random_index <- sample(nrow(clinical_df))
  
  if(metric_to_use == "OS"){
    clinical_df$OS <- clinical_df$OS[random_index]
    clinical_df$OS.time <- clinical_df$OS.time[random_index]
  } else{
    clinical_df$PFI <- clinical_df$PFI[random_index]
    clinical_df$PFI.time <- clinical_df$PFI.time[random_index]
  }
  
  
  #Run cox regression for KM plots
  if("purity" %in% colnames(clinical_df)){
    covariates <- c("purity", "sex", "age")
  } else {
    covariates <- c("sex", "age")
  }
  
  
  clinical_df$sex <- factor(clinical_df$sex, levels = c("FEMALE","MALE"))
  
  curr_project_plots <- list()
  for(curr_signature in names(all_signatures)){
    
    # Run cox regression model again to get HRs, confidence intervals, and p values
    regression_str <- paste0("Surv(",metric_to_use,".time, ",metric_to_use,") ~ ", curr_signature, " + ", paste(covariates,collapse = " + "))
    data_df <- as.data.frame(merge(clinical_df,gene_set_scores_df, by='submitter_id'))
    
    data_df <- data_df %>%
      dplyr::select(covariates,names(all_signatures), metric_to_use, paste0(metric_to_use,".time"))
    
    data_df$sex <- factor(data_df$sex, levels = c("FEMALE","MALE"))
    
    model <- coxph(as.formula(regression_str), data = data_df)
    model_sum <- summary(model)
    
    hazard_ratio_df[["project"]] <- append(hazard_ratio_df[["project"]],curr_project)
    hazard_ratio_df[["hazard_ratio"]] <- append(hazard_ratio_df[["hazard_ratio"]], model_sum$conf.int[curr_signature,"exp(coef)"])
    hazard_ratio_df[["hazard_ratio.high"]] <- append(hazard_ratio_df[["hazard_ratio.high"]],model_sum$conf.int[curr_signature,"upper .95"])
    hazard_ratio_df[["hazard_ratio.low"]] <- append(hazard_ratio_df[["hazard_ratio.low"]],model_sum$conf.int[curr_signature,"lower .95"])
    hazard_ratio_df[["p_value"]] <- append(hazard_ratio_df[["p_value"]], model_sum$coefficients[curr_signature, "Pr(>|z|)"])
    hazard_ratio_df[["signature"]] <- append(hazard_ratio_df[["signature"]], paste0("Supercluster ",which(curr_signature==names(all_signatures))))
  }
}

saveRDS(hazard_ratio_df, paste0(dataDirectory, "survival_analysis_hazard_ratios/pan_cancer_hazard_ratios_OS_random.rds"))

