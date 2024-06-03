source("source/cole_functions.R")

cell_lines <-  c("A549","K562","MCF7")

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- cell_lines

rac_type1_upregulated_genesets <- list()
rac_type2_upregulated_genesets <- list()

for(curr_cell_line in cell_lines){

  de_results <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/", curr_cell_line, "_cell_cluster_group_de.rds"))
  
  clusters <- as.character(unique(de_results$cluster))
  
  clusters <- clusters[!grepl("_0",clusters)]
  
  for(curr_cluster in clusters){
    curr_upreg_genes <- de_results %>% 
      filter(cluster == curr_cluster & p_val_adj < 0.05 & avg_log2FC > 0) %>% 
      pull(gene)
    
    
    if(grepl("_1",curr_cluster)){
      rac_type1_upregulated_genesets[[paste0(curr_cell_line,"_",curr_cluster)]] <- curr_upreg_genes  
    } else {
      rac_type2_upregulated_genesets[[paste0(curr_cell_line,"_",curr_cluster)]] <- curr_upreg_genes  
    }
  }
}




names(rac_type1_upregulated_genesets) <- gsub("_1$", "",names(rac_type1_upregulated_genesets))

supercluster1_upregulated_genesets <- rac_type1_upregulated_genesets[c(paste0("A549_",c(4,9)),paste0("K562_",c(4,9,5)),paste0("MCF7_",c(8,12)))]

supercluster1_upregulated_genes <- list("type1_supercluster1_upregulated" = find_consensus_geneset(supercluster1_upregulated_genesets,2))


##############

supercluster2_upregulated_genesets <- rac_type1_upregulated_genesets[c(paste0("A549_",c(19)),paste0("K562_",c(11)),paste0("MCF7_",c(5)))]

supercluster2_upregulated_genes <- list("type1_supercluster2_signature" = find_consensus_geneset(supercluster2_upregulated_genesets,2))

type1_supercluster_upregulated_genesets <- c(supercluster1_upregulated_genes,supercluster2_upregulated_genes)


saveRDS(type1_supercluster_upregulated_genesets,"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/type1_superclusters_upregulated_genesets.rds")
