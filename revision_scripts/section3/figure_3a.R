args = commandArgs(trailingOnly=TRUE)
# dataDirectory <- paste0(args[1],"final_data/")
# plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

# source("../survival_analysis/cox_regression.R")

source("revision_scripts/drug_treatment_functions.R")
source("/data/CDSL_hannenhalli/Cole/projects/survival_analysis/score_gene_expression.R")
dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"
plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_figures/"

################################################################################
# Survival Analysis for each cell line in cancer type matched TCGA samples
################################################################################

################################################################################
# Global RAC Signatures

cell_lines <- c("A549","K562","MCF7")
# cell_lines <- c("A549","MCF7")

gdcprojects <- getGDCprojects()

all_signatures <- readRDS(paste0(dataDirectory, "genesets/global_rac_signatures.rds"))

metric_to_use <- "OS"
# metric_to_use <- "DFI"

hazard_ratio_df <- list()
plot_data <- list()

curr_cell_line <- "MCF7"
for(curr_cell_line in cell_lines){
  cat(curr_cell_line, "\n")
  
  
  curr_signature <- all_signatures[grepl(curr_cell_line,names(all_signatures))]
  curr_signature_name <- names(curr_signature)
  
  # Get TCGA project
  tcga_project <- get_tcga_project(curr_cell_line)
  
  # Get read count data
  tcga_matrix_vst <- readRDS(paste0(dataDirectory, "patient_data/TCGA_processed_data/",tcga_project,"_data.rds"))
  
  rownames(tcga_matrix_vst) <- make.names(rownames(tcga_matrix_vst), unique = T)
  
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
  
  temp <- gene_set_scores_df
  
  # Run cox regression model again to get HRs, confidence intervals, and p values
  regression_str <- paste0("Surv(",metric_to_use,".time, ",metric_to_use,") ~ ", curr_signature_name, " + ", paste(covariates,collapse = " + "))
  data_df <- as.data.frame(merge(clinical_df,gene_set_scores_df, by='submitter_id'))
  
  data_df <- data_df %>%
    dplyr::select(covariates, curr_signature_name, metric_to_use, paste0(metric_to_use,".time")) 
  
  if("purity" %in% covariates){
    data_df <- data_df %>% 
      filter(!is.nan(purity))
  }
  
  data_df$sex <- factor(data_df$sex, levels = c("FEMALE","MALE"))
  
  # data_df %>%
  #   filter(!is.na(data_df$DFI))
  
  
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
  # top_score <- quantile(data_df[[curr_signature_name]],probs=.75)
  # bottom_score <- quantile(data_df[[curr_signature_name]],probs=.25)
   
  
  
  km_data <- data_df %>%
    mutate("score" = ifelse(data_df[[curr_signature_name]] > median_score, "High","Low"))
  
  # km_data <- data_df %>%
  #   mutate("score" = ifelse(data_df[[curr_signature_name]] > top_score, "High", ifelse(data_df[[curr_signature_name]] < bottom_score, "Low","Med"))) %>%
  #   filter(!score == "Med")
  
  
  
  modform <- as.formula(paste0("Surv(",metric_to_use,".time, ",metric_to_use,") ~ score"))
  
  fit <- eval(substitute(survfit(modform, data = km_data), list(modform = modform)))
  
  # Save data for KM plot for current project
  plot_data[[curr_signature_name]][["fit"]][[curr_signature_name]] <- fit
  plot_data[[curr_signature_name]][["data"]][[curr_signature_name]] <- km_data
}



# saveRDS(plot_data, paste0(dataDirectory,"figure_data/figure_3a_plot_data.rds"))

#################################################################################
# Read in plot data and plot KM plots for each TCGA project
# (plotting is done separately so plots can be customized without having to score samples all over again)
#################################################################################

hazard_ratio_df <- data.frame(hazard_ratio_df)

# plot_data <- readRDS(paste0(dataDirectory,"figure_data/figure_3a_plot_data.rds"))

all_plots <- list()
curr_cell_line <- "MCF7"
for(curr_cell_line in cell_lines){
  
  cat(curr_cell_line, "\n")
  
  curr_plot_title <- paste0("\n", curr_cell_line, " (",get_tcga_project(curr_cell_line),")")
  
  # KM plot
  p <- ggsurvplot(plot_data[[curr_cell_line]][["fit"]][[curr_cell_line]], data = plot_data[[curr_cell_line]][["data"]][[curr_cell_line]],
                  pval = F,
                  legend.title="Expression Level:",
                  legend.labs= c("High", "Low"),
                  font.tickslab = c(10),
                  font.legend=c(10),
                  font.caption=c(10),
                  tables.theme = clean_theme())
  
  p <- p$plot
  
  p <- p + ggtitle(curr_plot_title)+
    xlab("")+
    ylab("")+
    theme(legend.text = element_text(size = 6, color = "black"),
          legend.title = element_text(size = 6, color = "black"),
          plot.title = element_text(size=8),
          axis.text = element_text(size=6))+
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
    cox_pval <- paste0("p = ", sprintf("%.4f",cox_pval))
  }
  
  p <- p+annotate(
    "text",
    x = Inf, y = Inf,
    vjust = 1.5, hjust = 1,
    label = paste0("Cox Regression:\nHR = ", cox_hr, "\n ",cox_pval),
    size = 3)
  
  all_plots <- append(all_plots, list(p))
}

names(all_plots) <- cell_lines

final_plot <- ggarrange(plotlist = all_plots, ncol = 3, nrow=1, common.legend = F)


final_plot <- annotate_figure(final_plot,
                              left = text_grob("Survival Probability", rot = 90, vjust = 1, size=8),
                              bottom = text_grob("Time", size=8))



final_plot

jpeg(paste0(plotDirectory,"figure_3a.jpg"), width=200, height = 120, units = "mm", res = 1000)
print(final_plot)
dev.off()


