library(TCGAbiolinks)
library(tidyverse)
library(SummarizedExperiment)
library(DESeq2)
library(survminer)
library(survival)
library(gtools)

ssgsea = function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
  row_names = rownames(X)
  num_genes = nrow(X)
  gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})
  
  # Ranks for genes
  R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')
  
  # Calculate enrichment score (es) for each sample (column)
  es = apply(R, 2, function(R_col) {
    gene_ranks = order(R_col, decreasing = TRUE)
    
    # Calc es for each gene set
    es_sample = sapply(gene_sets, function(gene_set_idx) {
      # pos: match (within the gene set)
      # neg: non-match (outside the gene set)
      indicator_pos = gene_ranks %in% gene_set_idx
      indicator_neg = !indicator_pos
      
      rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha
      
      step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
      step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)
      
      step_cdf_diff = step_cdf_pos - step_cdf_neg
      
      # Normalize by gene number
      if (scale) step_cdf_diff = step_cdf_diff / num_genes
      
      # Use ssGSEA or not
      if (single) {
        sum(step_cdf_diff)
      } else {
        step_cdf_diff[which.max(abs(step_cdf_diff))]
      }
    })
    unlist(es_sample)
  })
  
  if (length(gene_sets) == 1) es = matrix(es, nrow = 1)
  
  # Normalize by absolute diff between max and min
  if (norm) es = es / diff(range(es))
  
  # Prepare output
  rownames(es) = names(gene_sets)
  colnames(es) = colnames(X)
  return(es)
}

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
# Get list of projects
##################################################################################

gdcprojects <- getGDCprojects()

tcga_projects <- gdcprojects[grepl("TCGA", gdcprojects$id),]$id

# tcga_projects <- tcga_projects[1:20]
tcga_projects <- tcga_projects[1:20]

which(tcga_projects == "TCGA-LAML")
################################################################################
# Get clinical data from GDC
##################################################################################
i <- tcga_projects[15]

for(i in tcga_projects){
  cat(i,"\n")
  
  #Get clincal datqa then create new columns (deceased and OS)
  curr_clinical <- GDCquery_clinic(i)
  
  curr_clinical$deceased <- ifelse(curr_clinical$vital_status == "Alive", F, T)
  
  curr_clinical$OS <- ifelse(curr_clinical$vital_status == "Alive", 
                             curr_clinical$days_to_last_follow_up,
                             curr_clinical$days_to_death)
  
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
  
  # all_samples <- output_query_tcga[output_query_tcga$sample_type == "Primary Tumor", "cases"][1:50]
  
  
  
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
  
  ##################################################################################
  for(j in 1:length(all_gene_lists)){
    
    curr_gene_list <- all_gene_lists[[j]]
    
    ssgsea_result <- ssgsea(tcga_matrix_vst,curr_gene_list)
    
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
    
    # scores_mat$submitter_id <- gsub("-01.*", "", scores_mat$submitter_id)
    
    
    scores_mat$submitter_id <- sapply(scores_mat$submitter_id, FUN=function(x) return(substring(x,1,12)))
    
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
    
    
    cancer_type_name <- gsub("-","_",str_to_lower(i))
    gene_signatures_folder_name <- gsub(" ", "_", tolower(names(all_gene_lists[j])))
    
    
    if(!file.exists(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/figures/survival_plots/",gene_signatures_folder_name))){
      dir.create(file.path(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/figures/survival_plots/",gene_signatures_folder_name)))
    } 
    
    
    if(length(curr_gene_list) == 1){
      plt <- egg::ggarrange(plots = surv_plots, ncol = 1)
      png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/figures/survival_plots/",gene_signatures_folder_name,"/",gene_signatures_folder_name,"_",cancer_type_name,".png"), width = 1000,height = 1000)
    } else if(length(curr_gene_list) == 3){
      plt <- egg::ggarrange(plots = surv_plots, ncol = 3)
      png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/figures/survival_plots/",gene_signatures_folder_name,"/",gene_signatures_folder_name,"_",cancer_type_name,".png"), width = 3000,height = 1000)
    } else if(length(curr_gene_list) == 5){
      plt <- egg::ggarrange(plots = surv_plots, ncol = 3)
      png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/figures/survival_plots/",gene_signatures_folder_name,"/",gene_signatures_folder_name,"_",cancer_type_name,".png"), width = 3000,height = 1000)
    } else {
      plt <- egg::ggarrange(plots = surv_plots, ncol = 5)
      png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/figures/survival_plots/",gene_signatures_folder_name,"/",gene_signatures_folder_name,"_",cancer_type_name,".png"), width = 3000,height = 2000)
    }
    
    main_title <- paste0(names(all_gene_lists)[j]," Survival Plots (",i,")")
    plt <- annotate_figure(plt, top = text_grob(main_title, color = "black", face = "bold", size = 35))
    plot(plt)
    
    dev.off()
    
  }
}


