args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
setwd(args[1])

source("final_scripts/aucell_thresholding.R")
library(Seurat)
library(AUCell)
library(fgsea)
library(GSEABase)
library(data.table)
library(tidyverse)
library(foreach)
library(doMC)

# 1) -c cellset_filename (seurat object)
# 2) -a assay_to_use
# 3) -g genesets_filename
# 4) CPUS_PER_TASK
# 5) -o output_file

################################################################################
#Read in cells to score
################################################################################

cellset_filename <- args[1]
# Read in cell set data
data_to_use <- readRDS(paste0(dataDirectory, "processed_data/", cellset_filename, ".rds"))

cellset_title <- strsplit(cellset_filename, "/")[[1]][-1]

assay_to_use <- args[2]

################################################################################
#Read in genesets to use
#Named list of genesets
################################################################################

genesets_filename <- args[4]

genesets_title <- genesets_filename

genesets <- readRDS(paste0(dataDirectory, "genesets/", genesets_filename, ".rds"))

if(class(genesets) == "character"){
  genesets <- list(genesets)
  names(genesets) <- genesets_title
}

genesets_names <- names(genesets)

# Remove genesets that have no genes in the gene universe
gene_universe <- rownames(data_to_use)
genesets_to_keep <- c()
for(idx in 1:length(genesets)){
  
  if(sum(genesets[[idx]] %in% gene_universe) > 0){
    genesets_to_keep <- append(genesets_to_keep, idx)
  }
}

genesets <- genesets[genesets_to_keep]
################################################################################
#Run parallel AUCell Scoring
################################################################################

n_cpus <- as.numeric(args[5])

#Make sure that NormalizeData has been run on the Seurat Object
set.seed(42) #Make sure that a random seed is set to a fixed number. This ensures reproducibility across runs.
auc_obj <- compute_AUCell_scores(data_to_use, genesets, compute_thresholds=F, nCores = n_cpus, assay_to_use = assay_to_use)

# met_auc_obj$auc_mat contains raw AUCell score matrix

obj_name <- paste0(cellset_title, "_", genesets_title, "_AUCell_scores")
assign(obj_name, auc_obj$auc_mat)

file_name <- paste0(dataDirectory, "aucell_score_objects/", cellset_title,"_", genesets_title, "_aucell_scores.rds")
saveRDS(get(obj_name), file=file_name)


