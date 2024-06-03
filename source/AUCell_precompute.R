library(Seurat)
library(AUCell)
library(fgsea)
library(GSEABase)
library(data.table)
library(tidyverse)
library(foreach)
library(doMC)

args <- commandArgs(trailingOnly=TRUE)

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
data_to_use <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/", cellset_filename, ".rds"))

cellset_title <- strsplit(cellset_filename, "/")[[1]][-1]

assay_to_use <- args[2]

################################################################################
#Read in genesets to use
#Named list of genesets
################################################################################

genesets_filename <- args[3]

genesets_title <- genesets_filename

genesets <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", genesets_filename, ".rds"))

if(class(genesets) == "character"){
  genesets <- list(genesets)
  names(genesets) <- genesets_title
}

genesets_names <- names(genesets)

################################################################################
#Run parallel AUCell Scoring
################################################################################

n <- as.numeric(args[4])

genesets <- split(genesets, ceiling(seq_along(genesets) / (length(genesets)/n)))

registerDoMC(cores=future::availableCores())

x <- foreach (t = 1:n, .errorhandling="pass", .verbose=T, .combine=rbind) %dopar% {

  curr_geneset <- genesets[[t]]
  
  cat(names(curr_geneset), "\n")
  
  cells_AUC <- AUCell_run(as.matrix(data_to_use[[assay_to_use]]@data), curr_geneset)
  
  cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 
  
  raw_auc <- cells_AUC@assays@data$AUC
  
  return(list("scores" = raw_auc, "thresholds" = cells_assignment))
}




obj_name <- paste0(cellset_title, "_", genesets_title, "_aucell_result")
assign(obj_name, x)

file_name <- paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", obj_name, ".rds")
saveRDS(get(obj_name), file=file_name)




