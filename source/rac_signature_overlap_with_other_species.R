source("source/cole_functions.R")


RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")
clusters_of_interest <- paste0(curr_cell_line,"_cluster",RACs[[curr_cell_line]],"_signature")



#################

upregulated_ecoli_human_orthologs_table <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/orthologs/upregulated_ecoli_human_orthologs_table.rds")


human_orthologs <- unique(upregulated_ecoli_human_orthologs_table$HUMAN_SYMBOL[!is.na(upregulated_ecoli_human_orthologs_table$HUMAN_SYMBOL)])




global_rac_type_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/global_rac_type_signatures.rds")
rac_type1_signatures <- global_rac_type_signatures[grepl("type1", names(global_rac_type_signatures))]

shared_genes <- list()
for(curr_rac_signature in rac_type1_signatures){
  
  
  shared_genes <- append(shared_genes, list(intersect(curr_rac_signature, human_orthologs)))
  
  
}

unique(unlist(shared_genes))





sort(human_orthologs)
supercluster_sigantures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/type1_supercluster_signatures.rds")



unique(intersect(supercluster_sigantures[[1]],human_orthologs),intersect(supercluster_sigantures[[2]],human_orthologs))












curr_cell_line <- "A549"

cluster_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/",curr_cell_line,"_cluster_signatures_200.rds"))
RAC_signatures <- cluster_signatures[clusters_of_interest]

RAC_overlap_with_ecoli_signatures <- list()
for(i in 1:length(RAC_signatures)){
  
  temp <- intersect(RAC_signatures[[i]],human_orthologs)
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
