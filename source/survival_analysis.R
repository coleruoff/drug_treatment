library(dplyr)
library(data.table)
library(tidyr)
library(tidyverse)
library(survival)
library(survminer)
library(fgsea)
source("../survival_analysis/score_gene_expression.R")
source("../survival_analysis/cox_regression.R")

tcga_log_fpkm_mat <- fread("../survival_analysis/TCGA_BRCA_LUAD_Combined.csv.gz") %>% 
  tibble::column_to_rownames("gene") %>% as.matrix

tcga_survival_df <- fread("../survival_analysis/TCGA_Survival.csv")


# raj_common_resistance_signature <- list("raj_resistance" = readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_common_resistance_signature.rds"))


curr_cell_line <- "MCF7"
signature_length <- 200

cluster_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_cluster_signatures_",signature_length,".rds"))

genesets_to_use <- cluster_signatures

projects <- c("BRCA","LUAD")


sample_list <- list()
gene_set_scores_df_list <- list()
cox_regression_list <- list()

for (project in projects) {
  
  #Get samples in current project
  sample_list[[project]] <- dplyr::filter(tcga_survival_df,TCGA_project == project) %>%
    pull(submitter_id)
  
  #Get the log fpkm matrix and survival df for current samples
  cancer_specific_log_fpkm_mat <- tcga_log_fpkm_mat[,sample_list[[project]]]
  cancer_specific_survival_df <- tcga_survival_df %>% dplyr::filter(TCGA_project == project)
  
  #Score samples for genesets
  cancer_gene_set_score_info <- compute_bulk_normalized_gene_set_scores(gene_exp_mat = cancer_specific_log_fpkm_mat,gene_sets = genesets_to_use, num_controls = 100, num_bins = 10, q_thresh=0.95, use_median=F)
  normalized_gene_set_score_mat <- cancer_gene_set_score_info$fg - cancer_gene_set_score_info$bg
  
  
  gene_set_scores_df <- as.data.frame(t(normalized_gene_set_score_mat)) %>%
    tibble::rownames_to_column("submitter_id")
  
  gene_set_scores_df_list[[project]] <- gene_set_scores_df
  
  cox_regression_info <- cox_regression(sample_survival_df = cancer_specific_survival_df,survival_data_id_col = "submitter_id",duration_str = "OS",
                                        feature_data_id_col = "submitter_id",
                                        model_covariates = c("Gender","purity"),
                                        status_str="OS_Flag",
                                        sample_features_df = gene_set_scores_df, 
                                        km_features_to_plot=c(names(genesets_to_use)))
  
  cox_regression_list[[project]] <- cox_regression_info
}



for(project in projects){
  for(curr_geneset in names(genesets_to_use)){
    
    cat(curr_geneset,"\n")

    modform <- as.formula(paste0("Surv(OS, OS_Flag) ~ ", curr_geneset))

    fit <- eval(substitute(survfit(modform, data = cox_regression_list[[project]]$km_df), list(modform = modform)))

    curr_geneset_name <- str_to_title(gsub("cluster","cluster ", gsub("_"," ",curr_geneset)))

    file_name <- paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/figures/cluster_signature_survival_plots/", curr_cell_line, "/", curr_geneset, "_survial_plot_", project, ".png")
    png(file_name, width= 18, height= 15, units="in", res = 1200)

    p <- ggsurvplot(fit, data = cox_regression_list[[project]]$km_df, pval = T,pval.size = 10) +
      ggtitle(paste0("Survival of ", curr_geneset_name, " Expression Groups in TCGA ",project, " Samples"))
    
    
    plot(ggpar(p$plot, 
               font.main = c(25, "bold"),
               font.x = c(16, "bold"),
               font.y = c(16, "bold"),
               font.caption = c(46, "bold"), 
               font.legend = c(20, "bold"), 
               font.tickslab = c(16, "bold"),
               font.text = c(30)) +
           theme(plot.title = element_text(hjust = 0.5)))
    
    dev.off()
    
  }
}






