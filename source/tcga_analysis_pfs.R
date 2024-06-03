library(dplyr)
library(data.table)
library(tidyverse)
library(survival)
library(survminer)
library(fgsea)
library(GSVA)
source("../survival_analysis/score_gene_expression.R")
source("../survival_analysis/cox_regression.R")

##################################################################################
# RACs Signatures
##################################################################################

rac_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_signatures.rds")

##################################################################################
# Global Type 1 vs inactive Signature
##################################################################################

rac_type_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/global_rac_type_signatures.rds")

##################################################################################
# Type 1 Individual signatures
##################################################################################

rac_type1_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_type1_signatures.rds")

##################################################################################
# supercluster signature
##################################################################################

rac_supercluster_consensus_signature <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_supercluster_consensus_signature.rds")

##################################################################################
# type 1 supercluster signatures
##################################################################################

type1_supercluster_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/type1_supercluster_signatures.rds")

##################################################################################
# Create list of gene lists
##################################################################################

all_gene_lists <- list(rac_signatures,rac_type_signatures,rac_type1_signatures,rac_supercluster_consensus_signature,type1_supercluster_signatures)

names(all_gene_lists) <- c("RAC Signatures", "Global RAC Type Signatures" ,"RAC Type 1 Signatures",
                           "Supercluster Signature", "RAC Type 1 Superclusters Signatures")


################################################################################

tcga_log_fpkm_mat <- fread("../survival_analysis/TCGA_BRCA_LUAD_Combined.csv.gz") %>% 
  tibble::column_to_rownames("gene") %>% as.matrix

tcga_survival_df <- fread("../survival_analysis/TCGA_Survival.csv")


projects <- c("BRCA","LUAD")


sample_list <- list()

project <- projects[1]
for (project in projects) {

  cat(project, "\n")  
  #Create list of sample names for current cancer type (project)
  sample_list[[project]] <- dplyr::filter(tcga_survival_df,TCGA_project == project) %>%
    pull(submitter_id)
  
  #Select expression and survival data for current samples
  cancer_specific_log_fpkm_mat <- tcga_log_fpkm_mat[,sample_list[[project]]]
  cancer_specific_survival_df <- tcga_survival_df %>% dplyr::filter(TCGA_project == project)
  
  j <- 1
  for(j in 1:length(all_gene_lists)){
    
    curr_gene_list <- all_gene_lists[[j]]
    #Score samples for given genesets
    # cancer_gene_set_score_info <- compute_bulk_normalized_gene_set_scores(gene_exp_mat = cancer_specific_log_fpkm_mat,gene_sets = cluster_signatures, num_controls = 100, num_bins = 10, q_thresh=0.95, gene_universe=NULL, use_median=F)
    # 
    # #Subtract background geneset scores from foreground geneset scores
    # normalized_gene_set_score_mat <- cancer_gene_set_score_info$fg - cancer_gene_set_score_info$bg
    
    ssgsea_result <- gsva(cancer_specific_log_fpkm_mat,curr_gene_list, method='ssgsea')
    
    # ssgsea_result_norm <- t(scale(t(ssgsea_result)))
    
    gene_set_scores_df <- as.data.frame(t(ssgsea_result_norm)) %>%
      tibble::rownames_to_column("submitter_id")
    
    cox_regression_info <- cox_regression(sample_survival_df = cancer_specific_survival_df,survival_data_id_col = "submitter_id",duration_str = "PFS",
                                          feature_data_id_col = "submitter_id",
                                          model_covariates = c("Gender","purity"),
                                          status_str="PFS_Flag",
                                          sample_features_df = gene_set_scores_df, 
                                          km_features_to_plot=names(curr_gene_list),
                                          low_high_percentiles=c(0.5,0.5))
    
    
    # Fitting survival curve for each signature
    surv_plots <- list()
    for(curr_signature in colnames(gene_set_scores_df)[-1]){
      
      modform <- as.formula(paste0("Surv(PFS, PFS_Flag) ~ ", curr_signature))
      
      fit <- eval(substitute(survfit(modform, data = cox_regression_info$km_df), list(modform = modform)))
      
      p <- ggsurvplot(fit,
                      data = cox_regression_info$km_df,
                      pval = T,
                      risk.table = F)+
        ggtitle(paste0(curr_signature))
      
      surv_plots <- append(surv_plots, list(p$plot))
      
    }
    
    # Plotting all signature KM plots on same figure
    cancer_type_name <- gsub("-","_",str_to_lower(project))
    gene_signatures_folder_name <- gsub(" ", "_", tolower(names(all_gene_lists[j])))
    
    if(!file.exists(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/figures/survival_plots/",gene_signatures_folder_name))){
      dir.create(file.path(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/figures/survival_plots/",gene_signatures_folder_name)))
    } 
    
    
    if(length(curr_gene_list) == 1){
      plt <- egg::ggarrange(plots = surv_plots, ncol = 1)
      png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/figures/survival_plots/",gene_signatures_folder_name,"/",gene_signatures_folder_name,"_",cancer_type_name,"_PFS.png"), width = 1000,height = 1000)
    } else if(length(curr_gene_list) == 3){
      plt <- egg::ggarrange(plots = surv_plots, ncol = 3)
      png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/figures/survival_plots/",gene_signatures_folder_name,"/",gene_signatures_folder_name,"_",cancer_type_name,"_PFS.png"), width = 3000,height = 1000)
    } else if(length(curr_gene_list) == 5){
      plt <- egg::ggarrange(plots = surv_plots, ncol = 3)
      png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/figures/survival_plots/",gene_signatures_folder_name,"/",gene_signatures_folder_name,"_",cancer_type_name,"_PFS.png"), width = 3000,height = 1000)
    } else {
      plt <- egg::ggarrange(plots = surv_plots, ncol = 5)
      png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/figures/survival_plots/",gene_signatures_folder_name,"/",gene_signatures_folder_name,"_",cancer_type_name,"_PFS.png"), width = 3000,height = 2000)
    }
    
    main_title <- paste0(names(all_gene_lists)[j]," Survival Plots (",project,")")
    plt <- annotate_figure(plt, top = text_grob(main_title, color = "black", face = "bold", size = 35))
    plot(plt)
    
    dev.off()
  }
}
