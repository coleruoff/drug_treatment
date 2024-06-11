library(tidyverse)
library(xlsx)
source("source/final_scripts/drug_treatment_functions.R")

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

all_treatments <- c("RKO 1.9 uM Oxaliplatin 10 days","SW620 0.38 uM Oxaliplatin 21 days","SW620 0.5 uM Oxaliplatin 21 days","OVCAR8 IC25 Prexasertib 20 days","OVCAR8 IC50 Prexasertib 20 days","OVCAR8 IC75 Prexasertib 20 days")

final_df <- readRDS(paste0(dataDirectory, "supercluster_up_down_crispr_ranks_data.rds"))

top_genes <- list()
list_names <- c()
for(experiment_number in 1:6){
  if(experiment_number %in% c(2,3)){
    sc_to_use <- "Supercluster 1"
  } else {
    sc_to_use <- "Supercluster 3"
  }
  
  curr_ranks <- final_df %>% 
    filter(geneset == sc_to_use & up == "Up" & treatment == all_treatments[experiment_number]) %>% 
    pull(rank)
  
  crispr_ko_data <- read_xlsx(paste0(dataDirectory, "gottesman_crispr_data/CRISPR_KO_summary_negative.xlsx"), sheet = experiment_number)
  
  crispr_ko_ranks <- crispr_ko_data %>% 
    arrange(`neg|fdr`) %>% 
    mutate(ranks = 1:nrow(.)) %>% 
    mutate(adj_p_value = `neg|fdr`) %>% 
    select(id,ranks,adj_p_value)
  
  curr_genes <- crispr_ko_ranks %>% 
    filter(ranks %in% curr_ranks & adj_p_value < 0.05) %>%
    arrange(ranks) %>%
    pull(id) 
  # %>% sort()
  
  top_genes <- append(top_genes, list(curr_genes))
  
  list_names <- append(list_names,paste0("expr",experiment_number,"_",gsub(" ","_",sc_to_use)))
}

top_genes <- append(top_genes, list(c("")))
top_genes <- append(top_genes, list(sort(find_consensus_geneset(top_genes,4))))
top_genes <- append(top_genes, list(sort(find_consensus_geneset(top_genes,3))))
top_genes <- append(top_genes, list(sort(find_consensus_geneset(top_genes,2))))

list_names <- c(list_names,"___","Four Screens","Three Screens","Two Screens")

names(top_genes) <- list_names

max_len <- max(lengths(top_genes))

top_genes <- lapply(top_genes, FUN = function(x) append(x,rep("",max_len-length(x))))

top_genes <- data.frame(top_genes)

write.xlsx(top_genes, file="/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/supercluster_genes_in_crispr.xlsx", sheetName="supercluster_up_genes", row.names=FALSE)

################################################################################

final_df <- readRDS(paste0(dataDirectory, "supercluster_tf_top_bottom_crispr_ranks_data.rds"))

top_genes <- list()
list_names <- c()
for(experiment_number in 5){
  if(experiment_number %in% c(2,3)){
    sc_to_use <- "Supercluster 1"
  } else {
    sc_to_use <- "Supercluster 3"
  }
  
  curr_ranks <- final_df %>% 
    filter(geneset == sc_to_use & top == "Top" & treatment == all_treatments[experiment_number]) %>% 
    pull(rank)
  
  crispr_ko_data <- read_xlsx(paste0(dataDirectory, "gottesman_crispr_data/CRISPR_KO_summary_negative.xlsx"), sheet = experiment_number)
  
  crispr_ko_ranks <- crispr_ko_data %>% 
    arrange(`neg|fdr`) %>% 
    mutate(ranks = 1:nrow(.)) %>% 
    mutate(adj_p_value = `neg|fdr`) %>% 
    select(id,ranks,adj_p_value)
  
  curr_genes <- crispr_ko_ranks %>% 
    filter(ranks %in% curr_ranks & adj_p_value < 0.05) %>%
    arrange(ranks) %>%
    pull(id) 
  # %>% sort()
  
  top_genes <- append(top_genes, list(curr_genes))
  
  list_names <- append(list_names,paste0("expr",experiment_number,"_",gsub(" ","_",sc_to_use)))
}

top_genes <- append(top_genes, list(c("")))
top_genes <- append(top_genes, list(sort(find_consensus_geneset(top_genes,4))))
top_genes <- append(top_genes, list(sort(find_consensus_geneset(top_genes,3))))
top_genes <- append(top_genes, list(sort(find_consensus_geneset(top_genes,2))))

list_names <- c(list_names,"___","Four Screens","Three Screens","Two Screens")

names(top_genes) <- list_names

max_len <- max(lengths(top_genes))

top_genes <- lapply(top_genes, FUN = function(x) append(x,rep("",max_len-length(x))))

top_genes <- data.frame(top_genes)

write.xlsx(top_genes, file="/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/supercluster_genes_in_crispr.xlsx", sheetName="supercluster_top_tfs", row.names=FALSE, append = T)

