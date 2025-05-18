args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

source("final_scripts/drug_treatment_functions.R")
source("/data/CDSL_hannenhalli/Cole/projects/survival_analysis/cox_regression.R")
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"
# plotDirectory <- "/data/ruoffcj/projects/drug_treatment/"

################################################################################

supercluster_signatures <- readRDS(paste0(dataDirectory, "genesets/supercluster_up_signatures.rds"))

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

# curr_project <- "TCGA-GBM"

# tcga_projects <- tcga_projects[-which(tcga_projects == "TCGA-LAML")]

for(curr_project in tcga_projects){
  
  cat(curr_project, "\n")
  
  #Get read count data
  # tcga_matrix_vst <- get_read_count_data(curr_project)

  tcga_matrix_vst <- readRDS(paste0(dataDirectory, "patient_data/TCGA_processed_data/",curr_project,"_data.rds"))
  
  rownames(tcga_matrix_vst) <- make.names(rownames(tcga_matrix_vst),unique = T)
  
  # Score samples
  gsvaPar <- ssgseaParam(tcga_matrix_vst, all_signatures)
  gsva.es <- gsva(gsvaPar, verbose=T)
  
  
  gene_set_scores_df <- as.data.frame(t(gsva.es)) %>%
    tibble::rownames_to_column("submitter_id")
  
  # gene_set_scores_df <- as.data.frame(t(ssgsea_result)) %>%
  #   tibble::rownames_to_column("submitter_id")
  
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
  
  
  
  
  plots <- list()
  for(curr_signature in 1:nrow(cox_regression_info$regression_df)){
    signature_to_use <- names(all_signatures)[curr_signature]
    signatures_title <- gsub("_"," ",signature_to_use)
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


    cox_regression_info$km_df[[signature_to_use]] <- factor(cox_regression_info$km_df[[signature_to_use]], levels = c("High","Low"))


    plot_data[[curr_project]][["fit"]][[curr_signature]] <- fit
    plot_data[[curr_project]][["data"]][[curr_signature]] <- cox_regression_info$km_df

    p <- ggsurvplot(fit, data = cox_regression_info$km_df,
                    pval = T,
                    legend.title="Expression Level:",
                    legend.labs= c("High", "Low"),
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
    # theme(title = element_text(size=5))

    plots <- append(plots, list(p$plot))
    
  }
  
  if(metric_to_use == "OS"){
    main_title <- curr_project
    # file_name <- paste0(plotDirectory, "final_figures/survival_plots/pan_cancer/",curr_project,"_OS.png")
  } else {
    main_title <- paste0("RAC Supercluster Signatures - Progression Free Survival (",curr_project,")")
    # file_name <- paste0(plotDirectory, "final_figures/survival_plots/pan_cancer/",curr_project,"_PFS.png")
  }
  
  
  project_plot <- ggarrange(plotlist = plots,ncol = 3,nrow=1, common.legend = T)
  
  project_plot <- annotate_figure(p = project_plot, top = text_grob(main_title, face = "bold", size = 14))
  
  all_plots <- append(all_plots, list(project_plot))
}

names(all_plots) <- tcga_projects

saveRDS(all_plots, paste0(dataDirectory, "survival_analysis_hazard_ratios/pan_cancer_supercluster_km_plots_",metric_to_use, ".rds"))

hazard_ratio_df <- as.data.frame(hazard_ratio_df)

saveRDS(hazard_ratio_df, paste0(dataDirectory, "survival_analysis_hazard_ratios/pan_cancer_hazard_ratios_",metric_to_use, ".rds"))
