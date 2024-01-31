source("source/cole_functions.R")

curr_cell_line <- "A549"

cluster_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/",curr_cell_line,"_cluster_signatures_200.rds"))

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")
clusters_of_interest <- paste0(curr_cell_line,"_cluster",RACs[[curr_cell_line]],"_signature")

RAC_signatures <- cluster_signatures[clusters_of_interest]

#################

ecoli_orthologs <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ecoli_AMR_genesets_orthologs.rds")

RAC_overlap_with_ecoli_signatures <- list()
for(i in 1:length(RAC_signatures)){
  
  temp <- c()
  for(j in 1:length(ecoli_orthologs)){
    
    temp <- append(temp, intersect(RAC_signatures[[i]],ecoli_orthologs[[j]]))
  }
  RAC_overlap_with_ecoli_signatures <- append(RAC_overlap_with_ecoli_signatures, list(unique(temp)))
  
}

names(RAC_overlap_with_ecoli_signatures) <- paste0("cluster_",RACs[[curr_cell_line]],"_overlap_with_ecoli")


#################
yeast_orthologs <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/yeast_upregulated_orthologs.rds")

RAC_overlap_with_yeast_signature <- list()
for(i in 1:length(RAC_signatures)){
  
  
  temp <- append(temp, intersect(RAC_signatures[[i]],yeast_orthologs))
  
  RAC_overlap_with_yeast_signature <- append(RAC_overlap_with_yeast_signature, list(unique(temp)))
  
}

names(RAC_overlap_with_yeast_signature) <- paste0("cluster_",RACs[[curr_cell_line]],"_overlap_with_yeast")




Heatmap(calc_jaccard_matrix(RAC_overlap_with_yeast_signature,RAC_overlap_with_yeast_signature))
