source("source/final_scripts/drug_treatment_functions.R")
source("/data/ruoffcj/projects/survival_analysis/cox_regression.R")
library(GSVA)
library(tidyverse)
library(matrixStats)
library(readxl)
library(survminer)
library(survival)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"
# plotDirectory <- "/data/ruoffcj/projects/drug_treatment/"
args = commandArgs(trailingOnly=TRUE)
dataDirectory <- args[1]
plotDirectory <- args[2]

cluster_signatures <- readRDS(paste0(dataDirectory, "genesets/cluster_signatures.rds"))


cell_lines <- c("A549","K562","MCF7")
supercluster_components <- list()
supercluster_components[["sc1"]]  <- c(9,5,8)
names(supercluster_components[["sc1"]]) <- cell_lines
supercluster_components[["sc2"]]  <- c(14,9,13)
names(supercluster_components[["sc2"]]) <- cell_lines


curr_cell_line <- cell_lines[1]
curr_sc <- 1

metric_to_use <- "OS"
hazard_ratio_df <- list()
plot_data <- list()

gdcprojects <- getGDCprojects()

for(curr_cell_line in cell_lines[1:2]){
  
  cat(curr_cell_line, "\n")
  
  # Get TCGA project
  curr_project <- get_tcga_project(curr_cell_line)
  
  #Get read count data
  # tcga_matrix_vst <- get_read_count_data(curr_project)
  tcga_matrix_vst <- readRDS(paste0(dataDirectory, "patient_data/TCGA_processed_data/",tcga_project,"_data.rds"))
  
  #Filter out genes with low/zero variance
  row_vars <- rowVars(tcga_matrix_vst)
  tcga_matrix_vst <- tcga_matrix_vst[which(row_vars > 1e-10),]
  
  
  all_signatures <- list()
  for(curr_sc in 1:2){
    
    curr_component <- supercluster_components[[curr_sc]][[curr_cell_line]]
    
    curr_signature <- cluster_signatures[[curr_cell_line]][[curr_component]]
    
    all_signatures[[paste0(curr_cell_line,"_",curr_component)]] <- curr_signature
  }
  
  # Score samples
  gsvapar <- ssgseaParam(tcga_matrix_vst, all_signatures)
  ssgsea_result <- gsva(gsvapar)
  
  gene_set_scores_df <- as.data.frame(t(ssgsea_result)) %>%
    tibble::rownames_to_column("submitter_id")
  
  gene_set_scores_df$submitter_id <- sapply(gene_set_scores_df$submitter_id, function(x) gsub("\\.", "-",x))
  gene_set_scores_df$submitter_id <- sapply(gene_set_scores_df$submitter_id, function(x) substring(x, 1,12))
  
  # Get clinical data
  clinical_df <- get_clinical_data(curr_project)
  
  #Run cox regression for KM plots
  if("purity" %in% colnames(clinical_df)){
    covariates <- c("purity", "sex", "age")
  } else {
    covariates <- c("sex", "age")
  }
  
  
  clinical_df$sex <- factor(clinical_df$sex, levels = c("FEMALE","MALE"))
  
  
  # cox_regression_info <- cox_regression(sample_survival_df = clinical_df,
  #                                       survival_data_id_col = "submitter_id",
  #                                       duration_str = paste0(metric_to_use, ".time"),
  #                                       feature_data_id_col = "submitter_id",
  #                                       model_covariates = covariates,
  #                                       status_str=metric_to_use,
  #                                       sample_features_df = gene_set_scores_df,
  #                                       km_features_to_plot=names(all_signatures),
  #                                       low_high_percentiles = c(.5,.49))
  
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
    
    #Create KM Plots Formula and Data
    
    # Classify samples as High or Low
    median_score <- quantile(data_df[[curr_signature]],probs=.5)
    km_data <- data_df %>% 
      mutate("score" = ifelse(data_df[[curr_signature]] > median_score, "High","Low"))
    
    modform <- as.formula(paste0("Surv(",metric_to_use,".time, ",metric_to_use,") ~ score"))
    
    fit <- eval(substitute(survfit(modform, data = km_data), list(modform = modform)))
    
    # cox_regression_info$km_df[[curr_signature]] <- factor(cox_regression_info$km_df[[curr_signature]], levels = c("High","Low"))
    
    # Save data for KM plot for current project
    plot_data[[curr_cell_line]][["fit"]][[curr_signature]] <- fit
    plot_data[[curr_cell_line]][["data"]][[curr_signature]] <- km_data
  }
}

plot(survfit(model))
################################################################################
# Add MCF7 plots
################################################################################
curr_cell_line <- "MCF7"
curr_project <- "METABRIC"

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

patient_sample$OS_MONTHS <- as.numeric(patient_sample$OS_MONTHS)
patient_sample$OS_STATUS <- as.numeric(gsub(":[A-Z]*", "", patient_sample$OS_STATUS))
patient_sample$AGE_AT_DIAGNOSIS <- as.numeric(patient_sample$AGE_AT_DIAGNOSIS)

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

# Set up signatures to score data with
all_signatures <- list()
for(curr_sc in 1:2){
  
  curr_component <- supercluster_components[[curr_sc]][[curr_cell_line]]
  
  curr_signature <- cluster_signatures[[curr_cell_line]][[curr_component]]
  
  all_signatures[[paste0(curr_cell_line,"_",curr_component)]] <- curr_signature
}


#Filter out genes with low/zero variance
row_vars <- rowVars(microarray_data)
microarray_data <- microarray_data[which(row_vars > 1e-10),]

# Score samples
gsvapar <- ssgseaParam(microarray_data, all_signatures)
ssgsea_result <- gsva(gsvapar)

gene_set_scores_df <- as.data.frame(t(ssgsea_result)) %>%
  tibble::rownames_to_column("PATIENT_ID")

# Create Cox regression model and save data from KM plots

covariates <- c("SEX", "AGE_AT_DIAGNOSIS")
curr_project_plots <- list()
for(curr_signature in names(all_signatures)){
  
  # Run cox regression model again to get HRs, confidence intervals, and p values
  regression_str <- paste0("Surv(OS_MONTHS, OS_STATUS) ~ ", curr_signature, " + ", paste(covariates,collapse = " + "))
  data_df <- as.data.frame(merge(patient_sample,gene_set_scores_df, by='PATIENT_ID'))
  
  data_df <- data_df %>%
    dplyr::select(covariates,names(all_signatures), "OS_MONTHS", "OS_STATUS") 
    
  
  data_df$SEX <- factor(data_df$SEX, levels = c("Female","Male"))
  
  model <- coxph(as.formula(regression_str), data = data_df)
  model_sum <- summary(model)
  
  hazard_ratio_df[["project"]] <- append(hazard_ratio_df[["project"]],curr_project)
  hazard_ratio_df[["hazard_ratio"]] <- append(hazard_ratio_df[["hazard_ratio"]], model_sum$conf.int[curr_signature,"exp(coef)"])
  hazard_ratio_df[["hazard_ratio.high"]] <- append(hazard_ratio_df[["hazard_ratio.high"]],model_sum$conf.int[curr_signature,"upper .95"])
  hazard_ratio_df[["hazard_ratio.low"]] <- append(hazard_ratio_df[["hazard_ratio.low"]],model_sum$conf.int[curr_signature,"lower .95"])
  hazard_ratio_df[["p_value"]] <- append(hazard_ratio_df[["p_value"]], model_sum$coefficients[curr_signature, "Pr(>|z|)"])
  hazard_ratio_df[["signature"]] <- append(hazard_ratio_df[["signature"]], paste0("Supercluster ",which(curr_signature==names(all_signatures))))
  
  #Create KM Plots Formula and Data
  
  # Classify samples as High or Low
  median_score <- quantile(data_df[[curr_signature]],probs=.5)
  km_data <- data_df %>% 
    mutate("score" = ifelse(data_df[[curr_signature]] > median_score, "High","Low"))
  
  modform <- as.formula(paste0("Surv(OS_MONTHS, OS_STATUS) ~ score"))
  
  fit <- eval(substitute(survfit(modform, data = km_data), list(modform = modform)))
  
  # cox_regression_info$km_df[[curr_signature]] <- factor(cox_regression_info$km_df[[curr_signature]], levels = c("High","Low"))
  
  # Save data for KM plot for current project
  plot_data[[curr_cell_line]][["fit"]][[curr_signature]] <- fit
  plot_data[[curr_cell_line]][["data"]][[curr_signature]] <- km_data
}

#####################################################
curr_cell_line <- cell_lines[1]

all_plots <- list()
for(curr_cell_line in cell_lines){
  
  cat(curr_cell_line, "\n")
  
  curr_project <- get_tcga_project(curr_cell_line)
  
  curr_project_plots <- list()
  
  # Create KM plot for each signature
  for(curr_signature in 1:2){
    
    # Clean up signature name for plot title
    signatures_title <- gsub("_"," ",curr_signature)
    signatures_title <- gsub(" signature","",signatures_title)
    signatures_title <- gsub("supercluster","Supercluster ",signatures_title)
    
    signatures_title <- names(plot_data[[curr_cell_line]]$fit)[curr_signature]
    
    # KM plot
    p <- ggsurvplot(plot_data[[curr_cell_line]][["fit"]][[signatures_title]], data = plot_data[[curr_cell_line]][["data"]][[signatures_title]],
                    pval = T,
                    legend.title="Expression Level:",
                    legend.labs= c("High", "Low"),
                    pval.size = 6,
                    font.title=c(24),
                    font.x=c(12),
                    font.y=c(28),
                    font.tickslab.y=c(20),
                    font.tickslab.x=c(12),
                    font.legend=c(50),
                    font.caption=c(30),
                    tables.theme = clean_theme())
    
    
    
    p$plot <- p$plot +
      theme(legend.text = element_text(size = 20, color = "black", face = "bold"),
            legend.title = element_text(size = 20, color = "black", face = "bold"))+
      NoLegend()
    
    p <- p + ggtitle(signatures_title)+
      xlab("")+
      ylab("")
    
    curr_project_plots <- append(curr_project_plots, list(p$plot))
  }
  
  # Combine all KM plots for current TCGA project
  project_plot <- ggarrange(plotlist = curr_project_plots, ncol = 2,nrow=1, common.legend = T)
  
  project_plot <- annotate_figure(project_plot, top = text_grob(curr_project, face = "bold", size = 34))
  
  all_plots <- append(all_plots, list(project_plot))
}

names(all_plots) <- cell_lines

# saveRDS(all_plots, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/figure_data/figure_4b_plotlist.rds")
# all_plots <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/figure_data/figure_4b_plotlist.rds")

# Combine all TCGA projects plots
final_plot <- ggarrange(plotlist = all_plots, ncol = 1,nrow=3, common.legend = T)

final_title <- "Supercluster Component Signatures - Overall Survival\n"

final_plot <- annotate_figure(final_plot, top = text_grob(final_title, face = "bold", size = 60),
                              left = text_grob("Survival Probability", rot = 90, vjust = 1, size=45, face="bold"),
                              bottom = text_grob("Time", size=45, face="bold"))

# png(paste0("/data/ruoffcj/projects/drug_treatment/final_figures/figure_4b.png"),
#     width=45, height=60, units= "in", res=300)

print(final_plot)

# dev.off()


data.frame(hazard_ratio_df)





