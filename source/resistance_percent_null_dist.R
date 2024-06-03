source("/data/ruoffcj/projects/aucell_scoring/aucell_thresholding.R")
library(Seurat)
library(AUCell)
library(fgsea)
library(GSEABase)
library(data.table)
library(tidyverse)

curr_cell_line <- "A549"

data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

################################################################################

genesets_filename <- "raj_watermelon_resistance_signature"

genesets <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", genesets_filename, ".rds"))

if(class(genesets) == "character"){
  genesets <- list(genesets)
  names(genesets) <- genesets_title
}

genesets_names <- names(genesets)

################################################################################

active_percent_runs <- c()
samples <- c()

for(i in 1:100){
  set.seed(i)
  shuffled_data <- data[sample(nrow(data)),]
  auc_obj <- compute_AUCell_scores(shuffled_data, genesets, compute_thresholds=F, assay_to_use = "RNA")
  
  #Compute separate thresholds across time-points
  computed_thresholds_df <- compute_shuffled_gene_set_AUCell_scores(shuffled_data, gene_sets=genesets, do_sample_wise=F, q_thresh=0.95, num_controls=100, assay_to_use = "RNA")
  
  curr_percent <- sum(auc_obj$auc_mat[,1] > computed_thresholds_df$threshold[1])/ncol(shuffled_data)
  
  active_percent_runs <- append(active_percent_runs,curr_percent)
  samples <- append(samples, list(sample(nrow(data))))
}


# plot(density(active_percent_runs))

saveRDS(active_percent_runs, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/resistance_percent_null_dist.rds")


