args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"revision_data/")
setwd(args[1])

source("revision_scripts/aucell_thresholding.R")
source("revision_scripts/drug_treatment_functions.R")

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"
################################################################################

if(!file.exists(paste0(dataDirectory,"aucell_score_objects/"))){
  dir.create(paste0(dataDirectory,"aucell_score_objects/")) 
}

# 1) -c cellset_filename (seurat object)
# 2) -a assay_to_use
# 3) -g genesets_filename
# 4) CPUS_PER_TASK

################################################################################
#Read in cells to score
################################################################################

cellset_filename <- args[2]

cat(paste0(dataDirectory, "processed_data/", cellset_filename, ".rds\n\n"))

# Read in cell set data
data_to_use <- readRDS(paste0(dataDirectory, "processed_data/", cellset_filename, ".rds"))

cellset_title <- strsplit(cellset_filename, "/")[[1]][-1]

assay_to_use <- args[3]

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

obj_name <- paste0(cellset_title, "_", genesets_title, "_aucell_scores")
assign(obj_name, auc_obj$auc_mat)

file_name <- paste0(dataDirectory, "aucell_score_objects/", cellset_title,"_", genesets_title, "_aucell_scores.rds")
saveRDS(get(obj_name), file=file_name)


#Compute separate thresholds across time-points
# debugonce(compute_shuffled_gene_set_AUCell_scores)
computed_thresholds_df <- compute_shuffled_gene_set_AUCell_scores(data_to_use, gene_sets=genesets, nCores=n_cpus, do_sample_wise=F, q_thresh=0.95, num_controls=100, assay_to_use = assay_to_use)

# print(head(computed_thresholds_df)) #Each row is a threshold for a gene set.

obj_name <- paste0(cellset_title, "_", genesets_title, "_aucell_thresholds")
assign(obj_name, computed_thresholds_df)


file_name <- paste0(dataDirectory, "aucell_score_objects/", cellset_title,"_", genesets_title, "_aucell_thresholds.rds")
saveRDS(get(obj_name), file=file_name)




