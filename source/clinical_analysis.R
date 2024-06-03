library(tidyverse)
library(ggpubr)
library(DESeq2)

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

drug_response_function <- function(all_signatures){
  
  drug_classes <- list()
  drug_classes <- append(drug_classes, list(c("MGH_Alpelisib","Sorafenib","Sorafenib_2","MK2206","Tipifarnib_1","Tipifarnib_2","Selinexor","Vismodegib")))
  drug_classes <- append(drug_classes, list(c("Anti-PD1","Anti-PD1 +- Anti-CTLA4", "Anti-PD1_2", "Anti-PD1_3","Anti-PD1_4")))
  drug_classes <- append(drug_classes, list(c("Bevacizumab","Bevacizumab_2","Bevacizumab_3","Bevacizumab_4","Trastuzumab","Trastuzumab_2","Trastuzumab_3","Trastuzumab_4","Trastuzumab_5","Cetuximab","Rituximab")))
  names(drug_classes) <- c("small_molecules","immunotherapy","mAb")
  
  enlight_response_data <- read.csv("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/enlight_data/drug_response_classifications_with_type.csv")
  
  rna_seq_datasets <- c("Trastuzumab_5","Bevacizumab_3", "Selinexor", "Vismodegib", "MGH_Ribociclib", "MGH_Alpelisib", "BRAFi_1", "Anti-PD1", "Anti-PD1_2", "Anti-PD1_5")
  
  enlight_response_data <- enlight_response_data %>%
    filter(Dataset %in% rna_seq_datasets)#&Dataset %in% drug_classes$small_molecules)
  
  # enlight_response_data <- enlight_response_data %>%
  #   filter(Dataset %in% drug_classes$small_molecules)
  
  all_drugs <- unique(enlight_response_data$Dataset)
  
  # all_drugs <- all_drugs[-which(all_drugs == "Anti-PD1 +- Anti-CTLA4")]
  
  all_sample_scores <- matrix(NA,nrow=0,ncol=(length(all_signatures)+1))
  
  # curr_drug <- all_drugs[1]
  
  for(curr_drug in all_drugs){
    cat(curr_drug, "\n")
    
    curr_data <- read.csv(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/enlight_data/", curr_drug,".csv"), row.names = 1)
    
    #Remove rows that are NAs
    curr_data <- curr_data[!is.na(curr_data[,1]),]
    
    #Normalize data if data are integers
    if(curr_data[which(curr_data != 0)[1],1]%%1 == 0){
      
      meta <- enlight_response_data %>% 
        filter(Dataset == curr_drug) %>% 
        select(Sample.ID,Response) %>% 
        column_to_rownames("Sample.ID")
      
      #reorder columns in data to match metadata
      curr_data <- curr_data[, rownames(meta)]
      
      all(colnames(curr_data) %in% rownames(meta))
      all(colnames(curr_data) == rownames(meta))
      
      curr_data <- round(curr_data)
      
      dds <- DESeqDataSetFromMatrix(countData = curr_data, colData = meta, design = ~ Response)
      
      dds <- estimateSizeFactors(dds)
      
      curr_data <- counts(dds, normalized=TRUE) 
    }
    
    ssgsea_res <- ssgsea(as.matrix(curr_data), all_signatures, scale = TRUE)
    ssgsea_res_scale <- t(scale(t(ssgsea_res)))
    
    colnames(ssgsea_res) <- colnames(curr_data)
    temp <- data.frame(t(ssgsea_res)) %>% 
      rownames_to_column()
    # cancer_gene_set_score_info <- compute_bulk_normalized_gene_set_scores(gene_exp_mat = curr_data, gene_sets = all_signatures, num_controls = 100, num_bins = 10, q_thresh=0.95, gene_universe=NULL, use_median=F)
    #Subtract background geneset scores from foreground geneset scores
    # normalized_gene_set_score_mat <- cancer_gene_set_score_info$fg - cancer_gene_set_score_info$bg
    # temp <- data.frame(t(normalized_gene_set_score_mat)) %>% 
    #   rownames_to_column()
    
    dim(all_sample_scores)
    
    dim(temp)
    if(ncol(temp) == ncol(all_sample_scores)){
      all_sample_scores <- rbind(all_sample_scores,temp)  
      cat("  added\n")
    }
    
  }
  
  colnames(all_sample_scores)[1] <- "Sample.ID"
  
  enlight_with_scores <- merge(enlight_response_data, all_sample_scores, by="Sample.ID")
  
  
  df <- enlight_with_scores %>% 
    pivot_longer(!c(Sample.ID,Dataset,Response,cancer_type), names_to = "geneset",values_to = "score")
  
  return(df)
}

##################################################################################
# RACs Signatures
cell_lines <- c("A549","K562","MCF7")

all_signatures <- list()

for(curr_cell_line in cell_lines){
  
  curr_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_cluster_signatures_200.rds"))
  
  RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
  names(RACs) <- c("A549","K562","MCF7")
  clusters_of_interest <- RACs[[curr_cell_line]]
  
  curr_signatures <- curr_signatures[clusters_of_interest]
  
  for(i in 1:length(curr_signatures)){
    if(length(curr_signatures[[i]]) > 10){
      curr_signatures[[i]] <- curr_signatures[[i]][1:10]  
    }
    
  }
  
  
  all_signatures <- append(all_signatures, curr_signatures)
}

rac_signatures <- all_signatures

##################################################################################
# Global active vs inactive Signature
cell_lines <- c("A549","K562","MCF7")

all_signatures <- list()

for(curr_cell_line in cell_lines){
  
  de_result <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/",curr_cell_line,"_all_resistant_cells_de.rds"))
  
  
  curr_resistant_genes <- de_result %>% 
    filter(cluster == "resistant" & p_val_adj < 0.05 & avg_log2FC > 0) %>% 
    arrange(desc(avg_log2FC)) %>% 
    pull(gene)
  
  
  if(length(curr_resistant_genes) > 100){
    curr_resistant_genes <- curr_resistant_genes[1:100]  
  }
    
  
  
  all_signatures <- append(all_signatures, list(curr_resistant_genes))
}

names(all_signatures) <- cell_lines

resistant_signatures <- all_signatures

##################################################################################
# DE within RACs Signatures

cell_lines <- c("A549","K562","MCF7")

all_signatures <- list()

for(curr_cell_line in cell_lines){
  
  curr_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_rac_within_cluster_de_signatures.rds"))
  
  
  for(i in 1:length(curr_signatures)){
    if(length(curr_signatures[[i]]) > 100){
      curr_signatures[[i]] <- curr_signatures[[i]][1:100]  
    }
  }

  names(curr_signatures) <- paste0(curr_cell_line, "_", names(curr_signatures)) 
  
  all_signatures <- append(all_signatures, curr_signatures)
}


within_rac_signatures <- all_signatures


#################################################################################

df <- drug_response_function(rac_signatures)

p <- ggboxplot(df, x = "Response", y = "score",
               color = "Response", palette = "jco",
               facet.by = "geneset", short.panel.labs = FALSE,
               add="")

# Use only p.format as label. Remove method name.
p + stat_compare_means(label = "p.format", method = "wilcox")+
  ggtitle(paste0("RAC Signatures ssGSEA Scores"))+
  xlab("")+
  ylab("Score")+
  theme(legend.position="right",
        legend.title = element_blank(),
        strip.text = element_text(size=20),
        title = element_text(size=30))

################################################################################
df <- drug_response_function(resistant_signatures)

p <- ggboxplot(df, x = "Response", y = "score",
               color = "Response", palette = "jco",
               facet.by = "geneset", short.panel.labs = FALSE,
               add="")

# Use only p.format as label. Remove method name.
p + stat_compare_means(label = "p.format", method = "wilcox")+
  ggtitle(paste0("Global Active vs Inactive Signatures ssGSEA Scores"))+
  xlab("")+
  ylab("Score")+
  theme(legend.position="right",
        legend.title = element_blank(),
        strip.text = element_text(size=20),
        title = element_text(size=30))

################################################################################
df <- drug_response_function(within_rac_signatures)

p <- ggboxplot(df, x = "Response", y = "score",
               color = "Response", palette = "jco",
               facet.by = "geneset", short.panel.labs = FALSE,
               add="")

# Use only p.format as label. Remove method name.
p + stat_compare_means(label = "p.format", method = "wilcox")+
  ggtitle(paste0("RAC Active vs Inactive Signatures ssGSEA Scores"))+
  xlab("")+
  ylab("Score")+
  theme(legend.position="right",
        legend.title = element_blank(),
        strip.text = element_text(size=20),
        title = element_text(size=30))

################################################################################
supercluster_signature <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/supercluster_signature.rds")

df <- drug_response_function(supercluster_signature)

p <- ggboxplot(df, x = "Response", y = "score",
               color = "Response", palette = "jco",
               facet.by = "geneset", short.panel.labs = FALSE,
               add="")

# Use only p.format as label. Remove method name.
p + stat_compare_means(label = "p.format", method = "wilcox")+
  ggtitle(paste0("Supercluster Signature ssGSEA Scores"))+
  xlab("")+
  ylab("Score")+
  theme(legend.position="right",
        legend.title = element_blank(),
        strip.text = element_text(size=20),
        title = element_text(size=30))





