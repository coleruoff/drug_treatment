library(tidyverse)
source("../survival_analysis/score_gene_expression.R")

#################################################################################

enlight_response_data <- read.csv("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/enlight_data/drug_response_classifications_with_type.csv")

all_drugs <- unique(enlight_response_data$Dataset)

data_list <- list()
for(curr_drug in all_drugs){
  data_list[[curr_drug]] <- read.csv(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/enlight_data/", curr_drug,".csv"), row.names = 1)
}

#################################################################################

signature_length <- 200

curr_cell_line <- "A549"
A549_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_cluster_signatures_",signature_length,".rds"))

curr_cell_line <- "K562"
K562_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_cluster_signatures_",signature_length,".rds"))

curr_cell_line <- "MCF7"
MCF7_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_cluster_signatures_",signature_length,".rds"))

supercluster_signature <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/supercluster_signature.rds"))

raj_resistance_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_resistance_signatures.rds")


all_signatures <- c(A549_signatures, K562_signatures, MCF7_signatures,supercluster_signature,raj_resistance_signatures)
#################################################################################

precompute_geneset_scores <- function(data_list, genesets_to_use){
  
  all_scores <- matrix(NA, ncol=0,nrow=length(genesets_to_use))
  
  for(i in 1:length(data_list)){
    cat(names(data_list)[i], "\n")
    
    #Select current experiemnt data, drug, drug class, and cancer type
    curr_data <- data_list[[i]]
    curr_drug <- names(data_list)[i]
    curr_cancer_type <- enlight_response_data %>% 
      filter(Dataset == curr_drug) %>% 
      pull(cancer_type) %>% 
      unique()
    
    if(curr_drug %in% drug_classes[[1]]){
      curr_class <- names(drug_classes)[1]
    } else if(curr_drug %in% drug_classes[[2]]){
      curr_class <- names(drug_classes)[2]
    } else if(curr_drug %in% drug_classes[[3]]){
      curr_class <- names(drug_classes)[3]      
    } else{
      next
    }
    
    #Remove genesets that have no overlap with current data row names
    # temp <- lapply(genesets_to_use, FUN = function(x){ if(length(intersect(x, rownames(curr_data))) == 0){
    #   return(x)
    # }})
    
    #Score all samples for each geneset
    cancer_gene_set_score_info <- compute_bulk_normalized_gene_set_scores(gene_exp_mat = curr_data, gene_sets = genesets_to_use, num_controls = 5, num_bins = 10, q_thresh=0.95, gene_universe=NULL, use_median=F)
    
    #Subtract background geneset scores from foreground geneset scores
    normalized_gene_set_score_mat <- cancer_gene_set_score_info$fg - cancer_gene_set_score_info$bg
    
    all_scores <- cbind(all_scores, normalized_gene_set_score_mat)
    
  }
  
  return(all_scores)
}

#################################################################################

scores <- precompute_geneset_scores(data_list, all_signatures)

saveRDS(scores, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/enlight_data/precompute_signature_scores.rds")



