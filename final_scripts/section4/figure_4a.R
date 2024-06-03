args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])
library(tidyverse)
library(GSVA)
library(readxl)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(survminer)
library(survival)
library(presto)
# source("../survival_analysis/cox_regression.R")
source("final_scripts/drug_treatment_functions.R")
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"
# plotDirectory <- "/data/ruoffcj/projects/drug_treatment/"

################################################################################
# Survival Analysis for each cell line in cancer type matched TCGA samples
################################################################################
cell_lines <- c("A549","K562","MCF7")

gdcprojects <- getGDCprojects()

all_signatures <- readRDS(paste0(dataDirectory, "genesets/global_rac_signatures.rds"))

metric_to_use <- "OS"
# metric_to_use <- "PFI"

hazard_ratio_df <- list()
plot_data <- list()

curr_cell_line <- "A549"
for(curr_cell_line in cell_lines[1:2]){
  cat(curr_cell_line, "\n")
  
  curr_signature <- all_signatures[grepl(curr_cell_line,names(all_signatures))]
  curr_signature_name <- names(curr_signature)
  
  # Get TCGA project
  tcga_project <- get_tcga_project(curr_cell_line)
  
  # Get read count data
  tcga_matrix_vst <- readRDS(paste0(dataDirectory, "patient_data/TCGA_processed_data/",tcga_project,"_data.rds"))
  
  #Filter out genes with low/zero variance
  row_vars <- rowVars(tcga_matrix_vst)
  tcga_matrix_vst <- tcga_matrix_vst[which(row_vars > 1e-10),]
  
  # Score samples
  gsvapar <- ssgseaParam(tcga_matrix_vst, curr_signature)
  ssgsea_result <- gsva(gsvapar)
  
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

  # Run cox regression model again to get HRs, confidence intervals, and p values
  regression_str <- paste0("Surv(",metric_to_use,".time, ",metric_to_use,") ~ ", curr_signature_name, " + ", paste(covariates,collapse = " + "))
  data_df <- as.data.frame(merge(clinical_df,gene_set_scores_df, by='submitter_id'))
  
  data_df <- data_df %>%
    dplyr::select(covariates, curr_signature_name, metric_to_use, paste0(metric_to_use,".time"))
  
  data_df$sex <- factor(data_df$sex, levels = c("FEMALE","MALE"))
  
  model <- coxph(as.formula(regression_str), data = data_df)
  model_sum <- summary(model)
  
  hazard_ratio_df[["project"]] <- append(hazard_ratio_df[["project"]],tcga_project)
  hazard_ratio_df[["hazard_ratio"]] <- append(hazard_ratio_df[["hazard_ratio"]], model_sum$conf.int[curr_signature_name,"exp(coef)"])
  hazard_ratio_df[["hazard_ratio.high"]] <- append(hazard_ratio_df[["hazard_ratio.high"]],model_sum$conf.int[curr_signature_name,"upper .95"])
  hazard_ratio_df[["hazard_ratio.low"]] <- append(hazard_ratio_df[["hazard_ratio.low"]],model_sum$conf.int[curr_signature_name,"lower .95"])
  hazard_ratio_df[["p_value"]] <- append(hazard_ratio_df[["p_value"]], model_sum$coefficients[curr_signature_name, "Pr(>|z|)"])
  hazard_ratio_df[["signature"]] <- append(hazard_ratio_df[["signature"]], curr_signature_name)
  
  #Create KM Plots Formula and Data
  
  # Classify samples as High or Low
  median_score <- quantile(data_df[[curr_signature_name]],probs=.5)
  km_data <- data_df %>% 
    mutate("score" = ifelse(data_df[[curr_signature_name]] > median_score, "High","Low"))
  
  modform <- as.formula(paste0("Surv(",metric_to_use,".time, ",metric_to_use,") ~ score"))
  
  fit <- eval(substitute(survfit(modform, data = km_data), list(modform = modform)))
  
  # Save data for KM plot for current project
  plot_data[[curr_signature_name]][["fit"]][[curr_signature_name]] <- fit
  plot_data[[curr_signature_name]][["data"]][[curr_signature_name]] <- km_data
}

################################################################################
# Add MCF7 plot
################################################################################
curr_project <- "METABRIC"

# Read in clincal data
clinical_data <- read_tsv(paste0(dataDirectory, "patient_data/metabric_data/brca_metabric_clinical_data.tsv"))
patient_sample <- read.table(paste0(dataDirectory, "patient_data/metabric_data/brca_metabric/data_clinical_patient.txt"),
                             fill=T)

colnames(patient_sample) <- patient_sample[1,]
patient_sample <- patient_sample[-1,]

temp <- matrix(unlist(patient_sample),nrow=nrow(patient_sample))
colnames(temp) <- colnames(patient_sample)
rownames(temp) <- rownames(patient_sample)

patient_sample <- temp

patient_sample <- as.data.frame(patient_sample) %>% 
  filter(CLAUDIN_SUBTYPE == "LumA")

patient_sample$OS_MONTHS <- as.numeric(patient_sample$OS_MONTHS)*30
patient_sample$OS_STATUS <- as.numeric(gsub(":[A-Z]*", "", patient_sample$OS_STATUS))
patient_sample$AGE_AT_DIAGNOSIS <- as.numeric(patient_sample$AGE_AT_DIAGNOSIS)

# Read in expression data
microarray_data <- read.table(paste0(dataDirectory, "patient_data/metabric_data/brca_metabric/data_mrna_illumina_microarray.txt"),
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
curr_signature_name <- names(curr_signature)

#Filter out genes with low/zero variance
row_vars <- rowVars(microarray_data)
microarray_data <- microarray_data[which(row_vars > 1e-10),]

# Score samples
gsvapar <- ssgseaParam(microarray_data, curr_signature)
ssgsea_result <- gsva(gsvapar)

gene_set_scores_df <- as.data.frame(t(ssgsea_result)) %>%
  tibble::rownames_to_column("PATIENT_ID")

# Create Cox regression model and save data from KM plots
covariates <- c("SEX", "AGE_AT_DIAGNOSIS")

# Run cox regression model again to get HRs, confidence intervals, and p values
regression_str <- paste0("Surv(OS_MONTHS, OS_STATUS) ~ ", curr_signature_name, " + ", paste(covariates,collapse = " + "))
data_df <- as.data.frame(merge(patient_sample,gene_set_scores_df, by='PATIENT_ID'))

data_df <- data_df %>%
  dplyr::select(covariates,curr_signature_name, "OS_MONTHS", "OS_STATUS") 

data_df$SEX <- factor(data_df$SEX, levels = c("Female","Male"))

model <- coxph(as.formula(regression_str), data = data_df)
model_sum <- summary(model)

hazard_ratio_df[["project"]] <- append(hazard_ratio_df[["project"]],curr_project)
hazard_ratio_df[["hazard_ratio"]] <- append(hazard_ratio_df[["hazard_ratio"]], model_sum$conf.int[curr_signature_name,"exp(coef)"])
hazard_ratio_df[["hazard_ratio.high"]] <- append(hazard_ratio_df[["hazard_ratio.high"]],model_sum$conf.int[curr_signature_name,"upper .95"])
hazard_ratio_df[["hazard_ratio.low"]] <- append(hazard_ratio_df[["hazard_ratio.low"]],model_sum$conf.int[curr_signature_name,"lower .95"])
hazard_ratio_df[["p_value"]] <- append(hazard_ratio_df[["p_value"]], model_sum$coefficients[curr_signature_name, "Pr(>|z|)"])
hazard_ratio_df[["signature"]] <- append(hazard_ratio_df[["signature"]], curr_signature_name)

#Create KM Plots Formula and Data

# Classify samples as High or Low
median_score <- quantile(data_df[[curr_signature_name]],probs=.5)
km_data <- data_df %>% 
  mutate("score" = ifelse(data_df[[curr_signature_name]] > median_score, "High","Low"))

modform <- as.formula(paste0("Surv(OS_MONTHS, OS_STATUS) ~ score"))

fit <- eval(substitute(survfit(modform, data = km_data), list(modform = modform)))

# Save data for KM plot for current project
plot_data[[curr_signature_name]][["fit"]][[curr_signature_name]] <- fit
plot_data[[curr_signature_name]][["data"]][[curr_signature_name]] <- km_data


if(!file.exists(paste0(dataDirectory,"survival_analysis_hazard_ratios/"))){
  dir.create(paste0(dataDirectory,"survival_analysis_hazard_ratios/")) 
}

saveRDS(hazard_ratio_df, paste0(dataDirectory, "survival_analysis_hazard_ratios/figure_4a_hazard_ratios_OS.rds"))

if(!file.exists(paste0(dataDirectory,"figure_data/"))){
  dir.create(paste0(dataDirectory,"figure_data/")) 
}

saveRDS(plot_data, paste0(dataDirectory,"figure_data/figure_4a_plot_data.rds"))

#################################################################################
# Read in plot data and plot KM plots for each TCGA project
# (plotting is done separately so plots can be customized without having to score samples all over again)
#################################################################################

hazard_ratio_df <- data.frame(hazard_ratio_df)

plot_data <- readRDS(paste0(dataDirectory,"figure_data/figure_4a_plot_data.rds"))

all_plots <- list()
for(curr_cell_line in cell_lines){
  
  cat(curr_cell_line, "\n")
  
  if(curr_cell_line == "MCF7"){
    curr_plot_title <- paste0("\n", curr_cell_line, " (METABRIC LumA Samples)")
  } else {
    curr_plot_title <- paste0("\n", curr_cell_line, " (",get_tcga_project(curr_cell_line),")")
  }
  
  # KM plot
  p <- ggsurvplot(plot_data[[curr_cell_line]][["fit"]][[curr_cell_line]], data = plot_data[[curr_cell_line]][["data"]][[curr_cell_line]],
                  pval = T,
                  legend.title="Expression Level:",
                  legend.labs= c("High", "Low"),
                  pval.size = 10,
                  font.tickslab = c(20),
                  font.legend=c(50),
                  font.caption=c(30),
                  tables.theme = clean_theme())
  
  p <- p$plot
  
  p <- p + ggtitle(curr_plot_title)+
    xlab("")+
    ylab("")+
    theme(legend.text = element_text(size = 20, color = "black", face = "bold"),
          legend.title = element_text(size = 20, color = "black", face = "bold"),
          plot.title = element_text(size=40),
          axis.text = element_text(size=20))+
    NoLegend()
  
  #Add cox regression HR and p-value
  cox_hr <- hazard_ratio_df %>% 
    filter(signature == curr_cell_line) %>% 
    pull(hazard_ratio) %>% 
    sprintf("%.1f",.)
  
  cox_pval <- hazard_ratio_df %>% 
    filter(signature == curr_cell_line) %>% 
    pull(p_value)
  
  if(cox_pval < .0001){
    cox_pval <- "p < .0001"
  } else {
    cox_pval <- paste0("p = ", sprintf("%.3f",cox_pval))
  }
    
  p <- p+annotate(
    "text",
    x = Inf, y = Inf,
    vjust = 1.5, hjust = 1,
    label = paste0("Cox Regression:\nHR = ", cox_hr, "\n ",cox_pval),
    size = 10)
  
  all_plots <- append(all_plots, list(p))
}

names(all_plots) <- cell_lines

final_plot <- ggarrange(plotlist = all_plots, ncol = 3,nrow=1, common.legend = F)

final_title <- "Global RAC Signatures - Overall Survival"

final_plot <- annotate_figure(final_plot, top = text_grob(final_title, face = "bold", size = 60),
                              left = text_grob("Survival Probability", rot = 90, vjust = 1, size=45, face="bold"),
                              bottom = text_grob("Time", size=45, face="bold"))

png(paste0(plotDirectory,"figure_4a.png"),
    width=30, height=10, units= "in", res = 300)

print(final_plot)

dev.off()
