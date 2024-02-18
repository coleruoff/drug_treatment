source("/data/ruoffcj/projects/drug_treatment/source/cole_functions.R")
source("/data/ruoffcj/projects/aucell_scoring/aucell_thresholding.R")

cell_lines <- c("A549","K562","MCF7")

################################################################################

all_data <- list()
all_data[["A549"]] <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/A549_processed_filtered.rds")
all_data[["K562"]] <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/K562_processed_filtered.rds")
all_data[["MCF7"]] <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/MCF7_processed_filtered.rds")

################################################################################

# Using signatures for type 1 cells (even though its called refined, its not)
refined_type1_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_type1_signatures.rds")

names1 <- gsub(pattern = "_rac_subcluster","_", names(refined_type1_signatures))
names(refined_type1_signatures) <- gsub(pattern = "_1_signature","", names1)

supercluster1_signatures <- refined_type1_signatures[c(paste0("A549_",c(4,9)),paste0("K562_",c(4,9,5)),paste0("MCF7_",c(8,12)))]
supercluster2_signatures <- refined_type1_signatures[c(paste0("A549_",c(19)),paste0("K562_",c(11)),paste0("MCF7_",c(5)))]

all_control_signatures <- list()
for(i in 1:500){
  cat(i, "\n")
  ##############
  # Find supercluster 1 control signature
  
  control_supercluster1_signatures <- supercluster1_signatures
  for(curr_cell_line in cell_lines){
    
    cat(curr_cell_line, "\n")
    
    curr_signatures <- supercluster1_signatures[grepl(curr_cell_line, names(supercluster1_signatures))]
    curr_signatures_names <- names(supercluster1_signatures)[grepl(curr_cell_line, names(supercluster1_signatures))]
    
    
    control_supercluster1_signatures[curr_signatures_names] <- find_control_gene_sets(all_data[[curr_cell_line]][["RNA"]]$counts, curr_signatures, num_bins=10)
  }
  
  
  control_supercluster1_signature <- list("type1_supercluster1_signature" = find_consensus_geneset(control_supercluster1_signatures,2))
  
  
  ##############
  # Find supercluster 2 control signature
  
  control_supercluster2_signatures <- supercluster2_signatures
  for(curr_cell_line in cell_lines){
    
    cat(curr_cell_line, "\n")
    
    curr_signatures <- supercluster2_signatures[grepl(curr_cell_line, names(supercluster2_signatures))]
    curr_signatures_names <- names(supercluster2_signatures)[grepl(curr_cell_line, names(supercluster2_signatures))]
    
    
    control_supercluster2_signatures[curr_signatures_names] <- find_control_gene_sets(all_data[[curr_cell_line]][["RNA"]]$counts, curr_signatures, num_bins=10)
  }
  
  
  control_supercluster2_signature <- list("type1_supercluster2_signature" = find_consensus_geneset(control_supercluster2_signatures,2))
  
  control_type1_supercluster_signatures <- c(control_supercluster1_signature,control_supercluster2_signature)
  
  all_control_signatures <- append(all_control_signatures, list(control_type1_supercluster_signatures))
}


saveRDS(all_control_signatures, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/control_supercluster_signatures.rds")



