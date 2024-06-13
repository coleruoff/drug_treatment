args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])
# library(gt)
library(GSVA)
library(readxl)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(survminer)
library(survival)
library(patchwork)
# source("/data/ruoffcj/projects/survival_analysis/cox_regression.R")
source("final_scripts/drug_treatment_functions.R")
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"
# plotDirectory <- "/data/ruoffcj/projects/drug_treatment/final_figures/"

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
  tcga_matrix_vst <- readRDS(paste0(dataDirectory, "patient_data/TCGA_processed_data/",curr_project,"_data.rds"))
  
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
    
    #Create KM Plots Formula and Data
    
    # Classify samples as High or Low
    median_score <- quantile(data_df[[curr_signature]],probs=.5)
    km_data <- data_df %>% 
      mutate("score" = ifelse(data_df[[curr_signature]] > median_score, "High","Low"))
    
    modform <- as.formula(paste0("Surv(",metric_to_use,".time, ",metric_to_use,") ~ score"))
    
    fit <- eval(substitute(survfit(modform, data = km_data), list(modform = modform)))
    
    # Save data for KM plot for current project
    plot_data[[curr_project]][["fit"]][[curr_signature]] <- fit
    plot_data[[curr_project]][["data"]][[curr_signature]] <- km_data
  }
}

saveRDS(hazard_ratio_df, paste0(dataDirectory, "survival_analysis_hazard_ratios/pan_cancer_hazard_ratios_OS.rds"))
saveRDS(plot_data, paste0(dataDirectory, "figure_data/figure_4b_plot_data.rds"))

#################################################################################
# Read in plot data and plot KM plots for each TCGA project
# (plotting is done separately so plots can be customized without having to score samples all over again)
#################################################################################

plot_data <- readRDS(paste0(dataDirectory, "figure_data/figure_4b_plot_data.rds"))

all_plots <- list()
for(curr_project in tcga_projects){

  cat(curr_project, "\n")

  curr_project_plots <- list()

  # Create KM plot for each signature
  for(curr_signature in names(all_signatures)){
    
    # Clean up signature name for plot title
    signatures_title <- gsub("_"," ",curr_signature)
    signatures_title <- gsub(" signature","",signatures_title)
    signatures_title <- gsub("supercluster","Supercluster ",signatures_title)

    # KM plot
    p <- ggsurvplot(plot_data[[curr_project]][["fit"]][[curr_signature]], data = plot_data[[curr_project]][["data"]][[curr_signature]],
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
  project_plot <- ggarrange(plotlist = curr_project_plots,ncol = 2,nrow=1, common.legend = T)

  project_plot <- annotate_figure(project_plot, top = text_grob(curr_project, face = "bold", size = 34))

  all_plots <- append(all_plots, list(project_plot))
}

names(all_plots) <- tcga_projects

# Combine all TCGA projects plots
final_plot <- ggarrange(plotlist = all_plots, ncol = 3,nrow=11, common.legend = T)

final_plot <- annotate_figure(final_plot, 
                              left = text_grob("Survival Probability", rot = 90, vjust = 1, size=45, face="bold"),
                              bottom = text_grob("Time", size=45, face="bold"))

png(paste0(plotDirectory, "figure_S4a.png"),
    width=45, height=60, units= "in", res=300)

print(final_plot)

dev.off()

