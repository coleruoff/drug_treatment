setwd("/data/ruoffcj/projects/drug_treatment/")
library(tidyverse)
library(GSVA)
library(readxl)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(survminer)
library(survival)
source("../survival_analysis/cox_regression.R")

# Read in deconvolved data
if(!exists("pred_TCGA_codefacs")){
  load("/data/CDSL_hannenhalli/Cole/data/deconvolved_data/pred_TCGA_codefacs_12.1.RData")
  load("/data/CDSL_hannenhalli/Cole/data/deconvolved_data/prob.TCGA.surv.strata.RData")
}


get_tcga_project <- function(cell_line){
  
  if(cell_line == "A549"){
    tcga_project <- "TCGA-LUAD"
  } else if(cell_line == "K562"){
    tcga_project <- "TCGA-LAML"
  } else{
    tcga_project <- "TCGA-BRCA"
  }
  
  return(tcga_project)
}


get_clinical_data <- function(project){
  
  clinical_data <- read_xlsx("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/TCGA_clinical_data.xlsx")
  colnames(clinical_data)[2] <- "submitter_id"
  
  clinical_data$age <- clinical_data$age_at_initial_pathologic_diagnosis
  clinical_data$sex <- clinical_data$gender
  
  tcga_type <- strsplit(project, "-")[[1]][2]
  if(tcga_type %in% Tumor.purity$Cancer.type){
    
    purity_data <- Tumor.purity
    
    colnames(purity_data)[1] <- "submitter_id"
    
    purity_data$submitter_id <- sapply(purity_data$submitter_id, function(x) substring(x, 1,12))
    
    clinical_data <- merge(clinical_data,purity_data,by = "submitter_id")
    
    clinical_data$CPE <- as.numeric(sapply(clinical_data$CPE, function(x) gsub(",",".",x)))
    
    clinical_data$purity <- clinical_data$CPE
  }
  
  return(clinical_data)
}

#################################################################################

cell_lines <- c("A549","K562","MCF7")

#################################################################################

rac_type_down_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/global_rac_type_down_signatures.rds")
rac_type_down_signatures <- rac_type_down_signatures[grepl("type1",names(rac_type_down_signatures))]

type_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/global_rac_type_signatures.rds")

type_signatures <- type_signatures[grepl("type1",names(type_signatures))]

all_signatures <- type_signatures
#################################################################################

deconvolved_tcga <- c("TCGA-LAML","TCGA-BRCA","TCGA-LUAD")

names(deconvolved_tcga) <- c("acute_myeloid_leukaemia","breast","lung_NSCLC_adenocarcinoma")

metric_to_use <- "OS"
# metric_to_use <- "PFI"

hazard_ratios <- c()

curr_cell_line <- cell_lines[1]

plots <- list()
for(curr_cell_line in cell_lines){
  cat(curr_cell_line, "\n")
  
  # Get TCGA sample expr data  
  tcga_project <- get_tcga_project(curr_cell_line)
  cancer_type <- names(deconvolved_tcga)[tcga_project == deconvolved_tcga]
  expr_mat <- pred_TCGA_codefacs[[cancer_type]]$pred_stageIII[,,"Cancer"]
  
  # Score samples
  curr_signature <- all_signatures[grepl(curr_cell_line,names(all_signatures))]
  
  ssgsea_result <- gsva(expr_mat,curr_signature,method='ssgsea')
  
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
  
  
  # if(metric_to_use %in% colnames(clinical_df))
  
  cox_regression_info <- cox_regression(sample_survival_df = clinical_df,
                                        survival_data_id_col = "submitter_id",
                                        duration_str = paste0(metric_to_use, ".time"),
                                        feature_data_id_col = "submitter_id",
                                        model_covariates = covariates,
                                        status_str=metric_to_use,
                                        sample_features_df = gene_set_scores_df, 
                                        km_features_to_plot=names(curr_signature))
  
  
  if(!all(is.na(cox_regression_info$km_df[[metric_to_use]]))){
    
    hazard_ratios <- append(hazard_ratios, cox_regression_info$regression_df$hazard_ratio)
    
    modform <- as.formula(paste0("Surv(",metric_to_use,".time, ",metric_to_use,") ~ ", names(curr_signature)))
    
    fit <- eval(substitute(survfit(modform, data = cox_regression_info$km_df), list(modform = modform)))
    
    
    if(metric_to_use == "OS"){
      plot_title <- paste0("\n",curr_cell_line, " RAC Type 1 Signature Overall Survival (",tcga_project,")")
      # file_name <- paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/survival_plots/cancer_type_matched/",curr_cell_line,"_type1_OS.png")
    } else {
      plot_title <- paste0("\n",curr_cell_line, " RAC Type 1 Signature Progression Free Survival (",tcga_project,")")
      # file_name <- paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/survival_plots/cancer_type_matched/",curr_cell_line,"_type1_PFS.png")
    }
    
    # png(file_name, height = 1000, width = 1000)
    
    cox_regression_info$km_df[[names(curr_signature)]] <- factor(cox_regression_info$km_df[[names(curr_signature)]], levels = c("High","Medium","Low"))
    
    p <- ggsurvplot(fit, data = cox_regression_info$km_df,
                    pval = T,
                    pval.size = 6,
                    legend.title="Expression Level:",
                    legend.labs= c("High","Medium","Low"),
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
    
    
    # dev.off()
  }
}



png(paste0("/data/ruoffcj/projects/drug_treatment/figures/deconvolved_km_plots.png"),
    width=40, height=12, units= "in", res = 300)


figure <- ggarrange(plotlist = plots, ncol=3, common.legend = T, legend=c("right"))

p <- annotate_figure(figure, left = text_grob("Survival Probability", rot = 90, vjust = 1, size=35, face="bold"),
                     bottom = text_grob("Time", size=35, face="bold"),
                     top=text_grob(paste0("RAC Type 1 Signatures Overall Survival in Deconvolved TCGA Samples"), size=40, face="bold"))


print(p)

dev.off()







##################################################################################
# RACs Signatures
##################################################################################

rac_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_signatures.rds")

##################################################################################
# Global active vs inactive Signature
##################################################################################

resistant_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/global_resistant_signatures.rds")

##################################################################################
# DE within RACs Signatures
##################################################################################

within_rac_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/within_rac_signatures.rds")

##################################################################################
# Active subpopulation signatures
##################################################################################

all_active_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/all_active_signatures.rds")

##################################################################################
# supercluster signature
##################################################################################

supercluster_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_supercluster_signature.rds")

##################################################################################
# supercluster consensus signature
##################################################################################

supercluster_consensus_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_supercluster_consensus_signature.rds")

##################################################################################
# Active supercluster signatures
##################################################################################

active_supercluster_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/active_supercluster_signatures.rds")

##################################################################################
# Create list of gene lists
##################################################################################

all_gene_lists <- list(rac_signatures,resistant_signatures,within_rac_signatures,all_active_signatures,supercluster_signatures,supercluster_consensus_signatures,active_supercluster_signatures)

names(all_gene_lists) <- c("RAC Signatures", "Global Active vs Inactive Signatures", 
                           "Within RACs Active vs Inactive Signatures","RACs Active Subpopulation Signatures",
                           "Supercluster Signature", "Supercluster Consensus Signature", "Active Superclusters Signatures")


################################################################################
# Read in deconvolved TCGA data
##################################################################################

if(!exists("pred_TCGA_codefacs")){
  load("/data/ruoffcj/deconvolved_data/pred_TCGA_codefacs_12.1.RData")
  load("/data/ruoffcj/deconvolved_data/prob.TCGA.surv.strata.RData")
}

####################################################

cancer_types <- names(pred_TCGA_codefacs)

gdcprojects <- getGDCprojects()

tcga_projects <- gdcprojects[grepl("TCGA", gdcprojects$id),]$id

# tcga_projects <- tcga_projects[1:20]

deconvolved_tcga <- c("TCGA-LAML","TCGA-DLBC","TCGA-BLCA","TCGA-BRCA","TCGA-UCEC","TCGA-LGG","TCGA-HNSC","TCGA-KICH","TCGA-KIRC","TCGA-KIRP","TCGA-STAD","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-SKCM","TCGA-UVM","TCGA-MESO")

names(deconvolved_tcga) <- c("acute_myeloid_leukaemia","B_cell_lymphoma","bladder","breast","endometrium","Glioma","head.and.neck","kidney","kidney","kidney","large_intestine","liver","lung_NSCLC_adenocarcinoma","lung_NSCLC_squamous_cell_carcinoma","melanoma","melanoma","mesothelioma")

for(i in 1:length(deconvolved_tcga)){
  
  tcga_sample <- deconvolved_tcga[i]
  cancer_type <- names(deconvolved_tcga[i])
  
  cat(tcga_sample,": ")
  cat(cancer_type,"\n")
  
  curr_clinical <- GDCquery_clinic(tcga_sample)
  
  curr_clinical$deceased <- ifelse(curr_clinical$vital_status == "Alive", F, T)
  
  curr_clinical$OS <- ifelse(curr_clinical$vital_status == "Alive", 
                             curr_clinical$days_to_last_follow_up,
                             curr_clinical$days_to_death)
  
  
  curr_expr_mat <- pred_TCGA_codefacs[[cancer_type]]$pred_stageIII[,,1]
  
  colnames(curr_expr_mat) <- sapply(colnames(curr_expr_mat),FUN = function(x) return(substring(x,1,12)))
  
  curr_clinical <- curr_clinical %>% 
    filter(submitter_id %in% colnames(curr_expr_mat))
  
  if(nrow(curr_clinical) > 0){
    ##################################################################################
    # j <- 1
    for(j in 1:length(all_gene_lists)){
      
      curr_gene_list <- all_gene_lists[[j]]
      
      ssgsea_result <- ssgsea(curr_expr_mat,curr_gene_list)
      
      ssgsea_result_norm <- t(scale(t(ssgsea_result)))
      
      scores_mat <- matrix(NA, nrow=nrow(ssgsea_result_norm),ncol=ncol(ssgsea_result_norm))
      colnames(scores_mat) <- colnames(ssgsea_result_norm)
      rownames(scores_mat) <- rownames(ssgsea_result_norm)
      
      for(p in 1:nrow(ssgsea_result_norm)){
        # scores_mat[i,] <- ifelse(ssgsea_result_norm[i,] < quantile(ssgsea_result_norm[i,],prob=.25), 'low', ifelse(ssgsea_result_norm[i,] > quantile(ssgsea_result_norm[i,],prob=.75), 'high', "med"))
        
        scores_mat[p,] <- ifelse(ssgsea_result_norm[p,] > quantile(ssgsea_result_norm[p,],prob=.5), 'high', 'low')
      }
      
      scores_mat <- t(scores_mat) %>% 
        as.data.frame() %>% 
        rownames_to_column("submitter_id")
      
      scores_mat$submitter_id <- gsub("\\.", "-", scores_mat$submitter_id)
      
      ################################################################################
      
      final_df <- merge(scores_mat,curr_clinical,by.x='submitter_id')
      
      # Fitting survival curve
      surv_plots <- list()
      for(curr_signature in colnames(scores_mat)[-1]){
        cat(curr_signature,"\n")
        fit <- survfit(as.formula(paste0("Surv(OS, deceased) ~ ", curr_signature)), data = final_df)
        
        
        p <- ggsurvplot(fit,
                        data = final_df,
                        pval = T,
                        risk.table = F)+
          ggtitle(paste0(curr_signature))
        
        surv_plots <- append(surv_plots, list(p$plot))
      }
      
      cancer_type_name <- gsub("-","_",str_to_lower(tcga_sample))
      gene_signatures_folder_name <- gsub(" ", "_", tolower(names(all_gene_lists[j])))
      
      
      if(!file.exists(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/figures/deconvolved_survival_plots/",gene_signatures_folder_name))){
        dir.create(file.path(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/figures/deconvolved_survival_plots/",gene_signatures_folder_name)))
      } 
      
      if(length(curr_gene_list) == 1){
        plt <- egg::ggarrange(plots = surv_plots, ncol = 1)
        png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/figures/deconvolved_survival_plots/",gene_signatures_folder_name,"/",gene_signatures_folder_name,"_",cancer_type_name,".png"), width = 1400,height = 1200)
      } else if(length(curr_gene_list) == 3){
        plt <- egg::ggarrange(plots = surv_plots, ncol = 3)
        png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/figures/deconvolved_survival_plots/",gene_signatures_folder_name,"/",gene_signatures_folder_name,"_",cancer_type_name,".png"), width = 3000,height = 1000)
      } else if(length(curr_gene_list) == 5){
        plt <- egg::ggarrange(plots = surv_plots, ncol = 3)
        png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/figures/deconvolved_survival_plots/",gene_signatures_folder_name,"/",gene_signatures_folder_name,"_",cancer_type_name,".png"), width = 3000,height = 1000)
      } else {
        plt <- egg::ggarrange(plots = surv_plots, ncol = 5)
        png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/figures/deconvolved_survival_plots/",gene_signatures_folder_name,"/",gene_signatures_folder_name,"_",cancer_type_name,".png"), width = 3000,height = 2000)
      }
      
      main_title <- paste0(names(all_gene_lists)[j]," Survival Plots (",tcga_sample," Deconvolved)")
      plt <- annotate_figure(plt, top = text_grob(main_title, color = "black", face = "bold", size = 35))
      plot(plt)
      
      dev.off()
      
    }
  }
}

