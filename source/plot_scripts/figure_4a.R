setwd("/data/ruoffcj/projects/drug_treatment/")
library(tidyverse)
library(GSVA)
library(readxl)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(survminer)
library(survival)
library(presto)
source("../survival_analysis/cox_regression.R")
source("source/cole_functions.R")

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

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
  GDCdownload(query_tcga)
  
  # Prepare data
  tcga_data <- GDCprepare(query_tcga, summarizedExperiment = T)
  
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
    select(gene_id,gene_name)
  
  
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

################################################################################
cell_lines <- c("A549","K562","MCF7")

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- cell_lines

# Create Global RAC Signatures for each cell line
rac_signatures <- list()
for(curr_cell_line in cell_lines){
  cat(curr_cell_line,"\n")
  
  #Read in cell line data
  data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))
  
  #read in DR signature scores and set active cells
  scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
  threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
  active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
  clusters_of_interest <- RACs[[curr_cell_line]]
  
  #Add metadata for RAC and Cell Group
  data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
  data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, "1", ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), 2, 0)), col.name = "cell_group")
  data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, paste0(data$Cluster, "_1"), ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), paste0(data$Cluster, "_2"), paste0(data$Cluster, "_0"))), col.name = "cell_cluster_group")
  
  Idents(data) <- data$rac
  
  #Find global RAC markers for current cell line
  de_res <- FindAllMarkers(data)
  
  curr_rac_up_genes <- de_res %>% 
    filter(cluster == "rac" & p_val_adj < 0.05 & avg_log2FC > 0) %>% 
    arrange(desc(avg_log2FC)) %>% 
    pull(gene)
  
  rac_signatures[[curr_cell_line]] <- curr_rac_up_genes[1:200]
}

saveRDS(rac_signatures, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/global_rac_signatures.rds")

################################################################################
# Survival Analysis for each cell line in cancer type matched TCGA samples

rac_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/global_rac_signatures.rds")

all_signatures <- post_rac_signatures

metric_to_use <- "OS"
# metric_to_use <- "PFI"

hazard_ratios <- c()

curr_cell_line <- cell_lines[1]

plots <- list()
for(curr_cell_line in cell_lines){
  cat(curr_cell_line, "\n")
  
  # Get TCGA sample count data  
  tcga_project <- get_tcga_project(curr_cell_line)
  
  # Get read count data
  # tcga_matrix_vst <- get_read_count_data(tcga_project)
  tcga_matrix_vst <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/TCGA_processed_data/",tcga_project,"_data.rds"))
  
  # Score samples
  curr_signature <- all_signatures[grepl(curr_cell_line,names(all_signatures))]
  
  ssgsea_result <- gsva(tcga_matrix_vst,curr_signature,method='ssgsea')
  
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
  
  
  cox_regression_info <- cox_regression(sample_survival_df = clinical_df,
                                        survival_data_id_col = "submitter_id",
                                        duration_str = paste0(metric_to_use, ".time"),
                                        feature_data_id_col = "submitter_id",
                                        model_covariates = covariates,
                                        status_str=metric_to_use,
                                        sample_features_df = gene_set_scores_df, 
                                        km_features_to_plot=names(curr_signature),
                                        low_high_percentiles = c(.5,.49))
  
  
  if(!all(is.na(cox_regression_info$km_df[[metric_to_use]]))){
    
    hazard_ratios <- append(hazard_ratios, cox_regression_info$regression_df$hazard_ratio)
    
    modform <- as.formula(paste0("Surv(",metric_to_use,".time, ",metric_to_use,") ~ ", names(curr_signature)))
    
    fit <- eval(substitute(survfit(modform, data = cox_regression_info$km_df), list(modform = modform)))
    
    
    if(metric_to_use == "OS"){
      plot_title <- paste0("\n",curr_cell_line, " (",tcga_project,")")
    } else {
      plot_title <- paste0("\n",curr_cell_line, " (",tcga_project,")")
    }
    
    
    cox_regression_info$km_df[[names(curr_signature)]] <- factor(cox_regression_info$km_df[[names(curr_signature)]], levels = c("High","Medium","Low"))
    
    p <- ggsurvplot(fit, data = cox_regression_info$km_df,
                    pval = T,
                    pval.size = 6,
                    legend.title="Expression Level:",
                    legend.labs= c("High","Low"),
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
    
  }
}

################################################################################
# Plot KM plots
figure <- ggarrange(plotlist = plots, ncol=3, common.legend = T, legend=c("right"))

p <- annotate_figure(figure, left = text_grob("Survival Probability", rot = 90, vjust = 1, size=35, face="bold"),
                     bottom = text_grob("Time", size=35, face="bold"),
                     top=text_grob(paste0("RAC Signatures Overall Survival in TCGA Samples"), size=40, face="bold"))



png(paste0("/data/ruoffcj/projects/drug_treatment/final_figures/figure_4a.png"),
    width=30, height=10, units= "in", res = 300)

print(p)

dev.off()

hazard_ratios

