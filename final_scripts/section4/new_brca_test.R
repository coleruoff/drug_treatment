setwd("/data/ruoffcj/projects/drug_treatment/")
library(tidyverse)
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

################################################################################

supercluster_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_supercluster_signatures.rds")

##################################################################################
# Get list of projects
##################################################################################

gdcprojects <- getGDCprojects()

tcga_projects <- gdcprojects[grepl("TCGA", gdcprojects$id),]$id

tcga_projects <- tcga_projects

##################################################################################

all_signatures <- supercluster_signatures

##################################################################################

metric_to_use <- "OS"

hazard_ratio_df <- list()

all_plots <- list()

plot_data <- list()

curr_project <- tcga_projects[1]

for(curr_project in tcga_projects){
  
  cat(curr_project, "\n")
  
  if(curr_project == "TCGA-BRCA"){
    
    # Read in clincal data
    # clinical_data <- read_tsv(paste0(dataDirectory, "metabric_data/brca_metabric_clinical_data.tsv"))
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
    
    
    # total_result <- list()
    # for(i in 1:3){
    #   curr_signature <- all_signatures[i]
    #   ssgsea_result <- gsva(microarray_data, curr_signature, method = "ssgsea")
    #   
    #   total_result <- append(total_result, list(ssgsea_result))
    # }
    
    # temp <- do.call(rbind, total_result)
    
    
    
    ssgsea_result <- gsva(microarray_data, all_signatures, method = "ssgsea")
    
    names(all_signatures)
    
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
                                          km_features_to_plot=names(all_signatures),
                                          low_high_percentiles = c(.5,.5))
    
  } else{
    #Get read count data
    # tcga_matrix_vst <- get_read_count_data(curr_project)
    tcga_matrix_vst <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/TCGA_processed_data/",curr_project,"_data.rds"))
    
    # Score samples
    ssgsea_result <- gsva(tcga_matrix_vst,all_signatures,method='ssgsea')
    
    gene_set_scores_df <- as.data.frame(t(ssgsea_result)) %>%
      tibble::rownames_to_column("submitter_id")
    
    gene_set_scores_df$submitter_id <- sapply(gene_set_scores_df$submitter_id, function(x) gsub("\\.", "-",x))
    gene_set_scores_df$submitter_id <- sapply(gene_set_scores_df$submitter_id, function(x) substring(x, 1,12))
    
    # Get clinical data
    clinical_df <- get_clinical_data(curr_project)
    
    
    
    #Run cox regression
    if("purity" %in% colnames(clinical_df)){
      covariates <- c("purity", "sex", "age")
    } else {
      covariates <- c("sex", "age")
    }
    
    
    clinical_df$sex <- factor(clinical_df$sex, levels = c("FEMALE","MALE"))
    
    
    cox_regression_info <- cox_regression(sample_survival_df = clinical_df,
                                          survival_data_id_col = "submitter_id",
                                          duration_str = paste0(metric_to_use, ".time"),
                                          feature_data_id_col = "submitter_id",
                                          model_covariates = covariates,
                                          status_str=metric_to_use,
                                          sample_features_df = gene_set_scores_df, 
                                          km_features_to_plot=names(all_signatures),
                                          low_high_percentiles = c(.5,.5))
  }
  
  
  
  plots <- list()
  for(curr_signature in 1:nrow(cox_regression_info$regression_df)){
    
    signature_to_use <- names(all_signatures)[curr_signature]
    signatures_title <- gsub("_"," ",signature_to_use)
    signatures_title <- gsub("type1 ","",signatures_title)
    signatures_title <- gsub(" signature","",signatures_title)
    signatures_title <- gsub("supercluster","Supercluster ",signatures_title)
    
    
    # Cox regression
    if(curr_project == "TCGA-BRCA"){
      regression_str <- paste0("Surv(",metric_to_use,".time, ",metric_to_use,") ~ ", signature_to_use, " + ", paste(covariates,collapse = " + "))
      data_df <- as.data.frame(merge(clinical_df,gene_set_scores_df, by='submitter_id'))
      
      data_df <- data_df %>%
        dplyr::select(covariates,names(all_signatures), metric_to_use, paste0(metric_to_use,".time"))
      
      data_df$sex <- factor(data_df$sex, levels = c("FEMALE","MALE"))
      
    } else {
      regression_str <- paste0("Surv(",metric_to_use,".time, ",metric_to_use,") ~ ", signature_to_use, " + ", paste(covariates,collapse = " + "))
      data_df <- as.data.frame(merge(clinical_df,gene_set_scores_df, by='submitter_id'))
      
      data_df <- data_df %>%
        dplyr::select(covariates,names(all_signatures), metric_to_use, paste0(metric_to_use,".time"))
      
      data_df$sex <- factor(data_df$sex, levels = c("FEMALE","MALE"))
    }
    
    
    model <- coxph(as.formula(regression_str), data = data_df)
    model_sum <- summary(model)
    
    hazard_ratio_df[["project"]] <- append(hazard_ratio_df[["project"]],curr_project)
    hazard_ratio_df[["hazard_ratio"]] <- append(hazard_ratio_df[["hazard_ratio"]], model_sum$conf.int[signature_to_use,"exp(coef)"])
    hazard_ratio_df[["hazard_ratio.high"]] <- append(hazard_ratio_df[["hazard_ratio.high"]],model_sum$conf.int[signature_to_use,"upper .95"])
    hazard_ratio_df[["hazard_ratio.low"]] <- append(hazard_ratio_df[["hazard_ratio.low"]],model_sum$conf.int[signature_to_use,"lower .95"])
    hazard_ratio_df[["p_value"]] <- append(hazard_ratio_df[["p_value"]], model_sum$coefficients[signature_to_use, "Pr(>|z|)"])
    hazard_ratio_df[["signature"]] <- append(hazard_ratio_df[["signature"]], paste0("Supercluster ",curr_signature))
    
    
    # 
    # #KM Plots
    # modform <- as.formula(paste0("Surv(",metric_to_use,".time, ",metric_to_use,") ~ ", signature_to_use))
    # 
    # fit <- eval(substitute(survfit(modform, data = cox_regression_info$km_df), list(modform = modform)))
    # 
    # 
    # cox_regression_info$km_df[[signature_to_use]] <- factor(cox_regression_info$km_df[[signature_to_use]], levels = c("High","Low"))
    # 
    # 
    # plot_data[[curr_project]][["fit"]][[curr_signature]] <- fit
    # plot_data[[curr_project]][["data"]][[curr_signature]] <- cox_regression_info$km_df
    # 
    # p <- ggsurvplot(fit, data = cox_regression_info$km_df,
    #                 pval = T, 
    #                 legend.title="Expression Level:",
    #                 legend.labs= c("High", "Low"),
    #                 pval.size = 6,
    #                 font.title=c(28),
    #                 font.x=c(28),
    #                 font.y=c(28),
    #                 font.tickslab=c(20),
    #                 font.legend=c(50),
    #                 font.caption=c(30),
    #                 tables.theme = clean_theme())
    # 
    # p$plot <- p$plot + 
    #   theme(legend.text = element_text(size = 14, color = "black", face = "bold"),
    #         legend.title = element_text(size = 14, color = "black", face = "bold"))
    # 
    # p <- p + ggtitle(signatures_title)+
    #   xlab("")+
    #   ylab("")
    # theme(title = element_text(size=5))
    # 
    # plots <- append(plots, list(p$plot))
    
  }
  # 
  # if(metric_to_use == "OS"){
  #   main_title <- curr_project
  #   file_name <- paste0(plotDirectory, "final_figures/survival_plots/pan_cancer/",curr_project,"_OS.png")
  # } else {
  #   main_title <- paste0("RAC Supercluster Signatures - Progression Free Survival (",curr_project,")")
  #   file_name <- paste0(plotDirectory, "final_figures/survival_plots/pan_cancer/",curr_project,"_PFS.png")
  # }
  # 
  # 
  # project_plot <- ggarrange(plotlist = plots,ncol = 2,nrow=1, common.legend = T)
  # 
  # project_plot <- annotate_figure(project_plot, top = text_grob(main_title, face = "bold", size = 14))
  # 
  # all_plots <- append(all_plots, list(project_plot))
}

saveRDS(hazard_ratio_df, paste0(dataDirectory, "survival_analysis_hazard_ratios/pan_cancer_hazard_ratios_OS_brca.rds"))
