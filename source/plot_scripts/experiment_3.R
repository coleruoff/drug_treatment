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
source("source/cole_functions.R")

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

get_read_count_data <- function(i){
  ################################################################################
  # Get read count data
  ################################################################################
  
  query_tcga <- GDCquery(project = i,
                         data.category = 'Transcriptome Profiling',
                         data.type = 'Gene Expression Quantification',
                         experimental.strategy = "RNA-Seq",
                         workflow.type = 'STAR - Counts',
                         access = 'open')
  
  output_query_tcga <- getResults(query_tcga)
  
  
  all_samples <- output_query_tcga[grepl("Primary",output_query_tcga$sample_type), "cases"]
  
  query_tcga <- GDCquery(project = i,
                         data.category = 'Transcriptome Profiling',
                         data.type = 'Gene Expression Quantification',
                         experimental.strategy = "RNA-Seq",
                         workflow.type = 'STAR - Counts',
                         access = 'open',
                         barcode = all_samples)
  
  
  output_query_tcga <- getResults(query_tcga)
  
  # Download data - GDCdownload
  GDCdownload(query_tcga,
              directory = "/data/CDSL_hannenhalli/Cole/GDCdata/")
  
  # Prepare data
  tcga_data <- GDCprepare(query_tcga, summarizedExperiment = T,
                          directory = "/data/CDSL_hannenhalli/Cole/GDCdata/")
  
  tcga_matrix <- assay(tcga_data, 'unstranded')
  
  
  # extract gene and sample metadata from summarizedExperiment object
  gene_metadata <- as.data.frame(rowData(tcga_data))
  coldata <- as.data.frame(colData(tcga_data))
  
  # Normalize data with DESeq
  dds <- DESeqDataSetFromMatrix(countData = tcga_matrix,
                                colData = coldata,
                                design = ~ 1)
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  # variance stabalizing transformation
  vsd <- vst(dds, blind = F)
  tcga_matrix_vst <- assay(vsd)
  
  
  temp <- gene_metadata %>% 
    filter(gene_id %in% rownames(tcga_matrix_vst)) %>% 
    dplyr::select(gene_id,gene_name)
  
  
  if(all(temp$gene_id == rownames(tcga_matrix_vst))){
    rownames(tcga_matrix_vst) <- temp$gene_name
  }
  
  
  
  return(tcga_matrix_vst)
  
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

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

################################################################################
# Create RAC supercluster signature
################################################################################

# Components of RAC supercluster
supercluster_components <- list()
supercluster_components[["A549"]] <- c(4,9,13)
supercluster_components[["K562"]] <- c(5)
supercluster_components[["MCF7"]] <-  c(5,8,17)

cell_lines <- c("A549","K562","MCF7")

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")

#################################################################################
# RAC Signatures

component_signatures <- list()

for(curr_cell_line in cell_lines){
  
  cat(curr_cell_line,"\n")
  #Read in cell line data
  # data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))
  # 
  # #read in DR signature scores and set active cells
  # scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
  # threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
  # active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
  # clusters_of_interest <- RACs[[curr_cell_line]]
  # 
  # #Add metadata for RAC and Cell Group
  # data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
  # data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, "1", ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), 2, 0)), col.name = "cell_group")
  # data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, paste0(data$Cluster, "_1"), ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), paste0(data$Cluster, "_2"), paste0(data$Cluster, "_0"))), col.name = "cell_cluster_group")
  # 
  # Idents(data) <- data$Cluster
  
  # de_results <- FindAllMarkers(data)
  # saveRDS(de_results, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/", curr_cell_line, "_cluster_all_markers_de.rds"))
  
  de_results <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/", curr_cell_line, "_cluster_all_markers_de.rds"))
  
  curr_rac_signatures <- list()
  for(i in levels(data)){
    
    curr_signature <- de_results %>%
      filter(p_val_adj < 0.05 & avg_log2FC > 0 & cluster == i) %>%
      arrange(desc(avg_log2FC)) %>%
      pull(gene)
    
    if(length(curr_signature) > 200){
      curr_signature <- curr_signature[1:200]
    }
    
    curr_rac_signatures <- append(curr_rac_signatures, list(curr_signature))
    
  }
  
  names(curr_rac_signatures) <- paste0(curr_cell_line, "_", levels(data))
  curr_rac_signatures <- curr_rac_signatures[supercluster_components[[curr_cell_line]]]
  
  component_signatures <- append(component_signatures,curr_rac_signatures)
}


supercluster_genes <- find_consensus_geneset(component_signatures,2)

supercluster_genes <- supercluster_genes[!grepl("^MT",supercluster_genes)]

supercluster_signature <- list("supercluster_signature" = supercluster_genes)

# saveRDS(supercluster_signature, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_supercluster_signature.rds")

##################################################################################
# Get list of projects
##################################################################################

gdcprojects <- getGDCprojects()

tcga_projects <- gdcprojects[grepl("TCGA", gdcprojects$id),]$id

tcga_projects <- tcga_projects1[1:3]

##################################################################################

all_signatures <- supercluster_signature

##################################################################################

metric_to_use <- "OS"

hazard_ratio_df <- list()

all_plots <- list()

plot_data <- list()

for(curr_project in tcga_projects){
  
  cat(curr_project, "\n")
  
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
                                        km_features_to_plot=names(all_signatures))
  
  
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


# saveRDS(plot_data, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/figure_data/figure_4b_plot_data.rds")
# saveRDS(all_plots, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/figure_data/figure_4b_plotlist.rds")
# all_plots <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/figure_data/figure_4b_plotlist.rds")

# plot_data <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/figure_data/figure_4b_plot_data.rds")


gdcprojects <- getGDCprojects()

tcga_projects <- gdcprojects[grepl("TCGA", gdcprojects$id),]$id

tcga_projects <- tcga_projects
signatures_titles <- c("Supercluster 1","Supercluster 2")


all_plots <- list()
for(curr_project in tcga_projects){
  
  cat(curr_project, "\n")
  
  plots <- list()
  
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
  
  project_plot <- ggarrange(plotlist = plots,ncol = 2,nrow=1, common.legend = T)
  
  project_plot <- annotate_figure(project_plot, top = text_grob(curr_project, face = "bold", size = 24))
  
  all_plots <- append(all_plots, list(project_plot))
}



final_plot <- ggarrange(plotlist = all_plots,ncol = 3,nrow=11, common.legend = T)

final_title <- "RAC Type 1 Supercluster Signatures - Overall Survival\n"

final_plot <- annotate_figure(final_plot, top = text_grob(final_title, face = "bold", size = 60),
                              left = text_grob("Survival Probability", rot = 90, vjust = 1, size=45, face="bold"),
                              bottom = text_grob("Time", size=45, face="bold"),)

png(paste0("/data/ruoffcj/projects/drug_treatment/final_figures/figure_4b.png"),
    width=35, height=60, units= "in", res=300)

print(final_plot)

dev.off()




















