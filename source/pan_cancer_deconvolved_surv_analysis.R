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

################################################################################
# Read in deconvolved TCGA data
##################################################################################
if(!exists("pred_TCGA_codefacs")){
  load("/data/CDSL_hannenhalli/Cole/data/deconvolved_data/pred_TCGA_codefacs_12.1.RData")
  load("/data/CDSL_hannenhalli/Cole/data/deconvolved_data/prob.TCGA.surv.strata.RData")
}

################################################################################
# Get list of projects
##################################################################################
deconvolved_tcga <- c("TCGA-LAML","TCGA-DLBC","TCGA-BLCA","TCGA-BRCA","TCGA-UCEC","TCGA-LGG","TCGA-HNSC","TCGA-KICH","TCGA-KIRC","TCGA-KIRP","TCGA-STAD","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-SKCM","TCGA-UVM","TCGA-MESO")
names(deconvolved_tcga) <- c("acute_myeloid_leukaemia","B_cell_lymphoma","bladder","breast","endometrium","Glioma","head.and.neck","kidney","kidney","kidney","large_intestine","liver","lung_NSCLC_adenocarcinoma","lung_NSCLC_squamous_cell_carcinoma","melanoma","melanoma","mesothelioma")

gdcprojects <- getGDCprojects()

tcga_projects <- gdcprojects[grepl("TCGA", gdcprojects$id),]$id

tcga_projects <- tcga_projects[tcga_projects%in%deconvolved_tcga]

##################################################################################

type1_supercluster_down_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/type1_supercluster_down_signatures.rds")

type1_supercluster_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/type1_supercluster_signatures.rds")

all_signatures <- type1_supercluster_signatures

##################################################################################

metric_to_use <- "OS"

hazard_ratio_df <- list()

curr_project <- tcga_projects[1]

curr_project <- "TCGA-UCEC"

all_plots <- list()

plot_data <- list()

for(curr_project in tcga_projects){
  
  cat(curr_project, "\n")
  
  #Get read count data
  cancer_type <- names(deconvolved_tcga)[curr_project == deconvolved_tcga]
  expr_mat <- pred_TCGA_codefacs[[cancer_type]]$pred_stageIII[,,"Cancer"]  
  
  # Score samples
  ssgsea_result <- gsva(expr_mat,all_signatures,method='ssgsea')
  
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
                                        km_features_to_plot=names(all_signatures))
  
  # supercluster1_hazard_ratios <- append(supercluster1_hazard_ratios, cox_regression_info$regression_df$hazard_ratio[1])
  # supercluster2_hazard_ratios <- append(supercluster2_hazard_ratios, cox_regression_info$regression_df$hazard_ratio[2])
  
  
  plots <- list()
  for(curr_signature in 1:nrow(cox_regression_info$regression_df)){
    
    signature_to_use <- names(all_signatures)[curr_signature]
    signatures_title <- gsub("_"," ",signature_to_use)
    signatures_title <- gsub("type1 ","",signatures_title)
    signatures_title <- gsub(" signature","",signatures_title)
    signatures_title <- gsub("supercluster","Supercluster ",signatures_title)
    
    
    # Cox regression
    regression_str <- paste0("Surv(",metric_to_use,".time, ",metric_to_use,") ~ ", signature_to_use, " + ", paste(covariates,collapse = " + "))
    data_df <- as.data.frame(merge(clinical_df,gene_set_scores_df, by='submitter_id'))
    
    data_df <- data_df %>%
      dplyr::select(covariates,names(all_signatures), metric_to_use, paste0(metric_to_use,".time"))
    
    data_df$sex <- factor(data_df$sex, levels = c("FEMALE","MALE"))
    
    model <- coxph(as.formula(regression_str), data = data_df)
    model_sum <- summary(model)
    
    hazard_ratio_df[["project"]] <- append(hazard_ratio_df[["project"]],curr_project)
    hazard_ratio_df[["hazard_ratio"]] <- append(hazard_ratio_df[["hazard_ratio"]], model_sum$conf.int[signature_to_use,"exp(coef)"])
    hazard_ratio_df[["hazard_ratio.high"]] <- append(hazard_ratio_df[["hazard_ratio.high"]],model_sum$conf.int[signature_to_use,"upper .95"])
    hazard_ratio_df[["hazard_ratio.low"]] <- append(hazard_ratio_df[["hazard_ratio.low"]],model_sum$conf.int[signature_to_use,"lower .95"])
    hazard_ratio_df[["p_value"]] <- append(hazard_ratio_df[["p_value"]], model_sum$coefficients[signature_to_use, "Pr(>|z|)"])
    hazard_ratio_df[["signature"]] <- append(hazard_ratio_df[["signature"]], paste0("Supercluster ",curr_signature))
    
    
    
    #KM Plots
    modform <- as.formula(paste0("Surv(",metric_to_use,".time, ",metric_to_use,") ~ ", signature_to_use))
    
    fit <- eval(substitute(survfit(modform, data = cox_regression_info$km_df), list(modform = modform)))
    
    
    cox_regression_info$km_df[[signature_to_use]] <- factor(cox_regression_info$km_df[[signature_to_use]], levels = c("High","Medium","Low"))
    
    
    plot_data[[curr_project]][["fit"]][[curr_signature]] <- fit
    plot_data[[curr_project]][["data"]][[curr_signature]] <- cox_regression_info$km_df
    
    p <- ggsurvplot(fit, data = cox_regression_info$km_df,
                    pval = T, 
                    legend.title="Expression Level:",
                    legend.labs= c("High","Medium","Low"),
                    pval.size = 6,
                    font.title=c(28),
                    font.x=c(28),
                    font.y=c(28),
                    font.tickslab=c(20),
                    font.legend=c(50),
                    font.caption=c(30),
                    tables.theme = clean_theme())
    
    p$plot <- p$plot + 
      theme(legend.text = element_text(size = 14, color = "black", face = "bold"),
            legend.title = element_text(size = 14, color = "black", face = "bold"))
    
    p <- p + ggtitle(signatures_title)+
      xlab("")+
      ylab("")
    theme(title = element_text(size=5))
    
    plots <- append(plots, list(p$plot))
    
  }
  
  if(metric_to_use == "OS"){
    main_title <- curr_project
    file_name <- paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/survival_plots/pan_cancer/",curr_project,"_OS.png")
  } else {
    main_title <- paste0("RAC Type 1 Supercluster Signatures - Progression Free Survival (",curr_project,")")
    file_name <- paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/survival_plots/pan_cancer/",curr_project,"_PFS.png")
  }
  
  
  project_plot <- ggarrange(plotlist = plots,ncol = 2,nrow=1, common.legend = T)
  
  project_plot <- annotate_figure(project_plot, top = text_grob(main_title, face = "bold", size = 14))
  
  all_plots <- append(all_plots, list(project_plot))
  
  
  
}

hazard_ratio_df <- data.frame(hazard_ratio_df)

saveRDS(hazard_ratio_df,"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/survival_analysis_hazard_ratios/pan_cancer_deconvolved_hazard_ratios_OS.rds")

# saveRDS(plot_data, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/figure_data/figure_4b_plot_data.rds")
# saveRDS(all_plots, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/figure_data/figure_4b_plotlist.rds")
# all_plots <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/figure_data/figure_4b_plotlist.rds")

# plot_data <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/figure_data/figure_4b_plot_data.rds")


signatures_titles <- c("Supercluster 1","Supercluster 2")


all_plots <- list()
for(curr_project in tcga_projects){
  
  cat(curr_project, "\n")
  
  plots <- list()
  for(curr_signature in 1:2){
    
    p <- ggsurvplot(plot_data[[curr_project]][["fit"]][[curr_signature]], data = plot_data[[curr_project]][["data"]][[curr_signature]],
                    pval = T, 
                    legend.title="Expression Level:",
                    legend.labs= c("High","Medium","Low"),
                    pval.size = 6,
                    font.title=c(28),
                    font.x=c(28),
                    font.y=c(28),
                    font.tickslab=c(20),
                    font.legend=c(50),
                    font.caption=c(30),
                    tables.theme = clean_theme())
    
    p$plot <- p$plot + 
      theme(legend.text = element_text(size = 14, color = "black", face = "bold"),
            legend.title = element_text(size = 14, color = "black", face = "bold"))
    
    p <- p + ggtitle(signatures_titles[curr_signature])+
      xlab("")+
      ylab("")
    theme(plot.title = element_text(size=18))
    
    plots <- append(plots, list(p$plot))
    
    
    
  }
  
  project_plot <- ggarrange(plotlist = plots,ncol = 2,nrow=1, common.legend = T)
  
  project_plot <- annotate_figure(project_plot, top = text_grob(curr_project, face = "bold", size = 24))
  
  all_plots <- append(all_plots, list(project_plot))
  
}



final_plot <- ggarrange(plotlist = all_plots,ncol = 3,nrow=6, common.legend = T)

final_title <- "RAC Type 1 Supercluster Signatures in Deconvolved TCGA Data - Overall Survival\n"

final_plot <- annotate_figure(final_plot, top = text_grob(final_title, face = "bold", size = 60),
                              left = text_grob("Survival Probability", rot = 90, vjust = 1, size=45, face="bold"),
                              bottom = text_grob("Time", size=45, face="bold"),)

png(paste0("/data/ruoffcj/projects/drug_treatment/figures/pan_cancer_deconvolved_km_plots.png"),
    width=40, height=60, units= "in", res=300)

print(final_plot)

dev.off()

