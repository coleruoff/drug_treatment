library(tidyverse)
source("/data/ruoffcj/projects/survival_analysis/score_gene_expression.R")


#################################################################################

precompute_scores <- function(data_list, genesets_to_use){
  
  all_scores <- data.frame(matrix(NA, nrow=length(genesets_to_use), ncol=0))
  
  for(i in 1:length(data_list)){
    cat(names(data_list)[i], "\n")
    
    #Select current experiemnt data, drug, drug class, and cancer type
    curr_data <- data_list[[i]]
    curr_drug <- names(data_list)[i]
    
    
    #Score all samples for each geneset
    cancer_gene_set_score_info <- compute_bulk_normalized_gene_set_scores(gene_exp_mat = curr_data, gene_sets = genesets_to_use, num_controls = 100, num_bins = 10, q_thresh=0.95, gene_universe=NULL, use_median=F)
    
    #Subtract background geneset scores from foreground geneset scores
    normalized_gene_set_score_mat <- cancer_gene_set_score_info$fg - cancer_gene_set_score_info$bg
    
  
    all_scores <- cbind(all_scores, normalized_gene_set_score_mat)
  
  }
  
 
  return(all_scores) 
}
  
#################################################################################
enlight_response_data <- read.csv("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/enlight_data/drug_response_classifications_with_type.csv")

all_drugs <- unique(enlight_response_data$Dataset)

data_list <- list()
for(curr_drug in all_drugs){
  data_olist[[curr_drug]] <- read.csv(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/enlight_data/", curr_drug,".csv"), row.names = 1)
}

#################################################################################

A549_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/A549_cluster_signatures_200.rds")

K562_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/K562_cluster_signatures_200.rds")

MCF7_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/MCF7_cluster_signatures_200.rds")

supercluster_signature <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/supercluster_signature.rds"))

raj_resistance_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_resistance_signatures.rds")

all_signatures <- c(A549_signatures,K562_signatures,MCF7_signatures,supercluster_signature,raj_resistance_signatures)

#################################################################################
  
  
scores <- precompute_scores(data_list,all_signatures)


saveRDS(scores, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/enlight_data/precomputed_signature_scores.rds")


