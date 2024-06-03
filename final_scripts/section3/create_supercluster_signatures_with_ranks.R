source("source/final_scripts/drug_treatment_functions.R")

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

create_supercluster_signature <- function(supercluster_ranks, num_genes){
  supercluster_genes <- list(names(supercluster_ranks[[1]]),names(supercluster_ranks[[2]]),names(supercluster_ranks[[3]]))
  
  # Get common genes
  supercluster_genes <- find_consensus_geneset(supercluster_genes,3)
  
  # Retain common genes and sort alphabetically 
  for(i in 1:3){
    supercluster_ranks[[i]] <- supercluster_ranks[[i]][names(supercluster_ranks[[i]]) %in% supercluster_genes]
    
    supercluster_ranks[[i]] <- supercluster_ranks[[i]][sort(names(supercluster_ranks[[i]]))]
  }
  
  #calculate avg FC values for each gene
  avg_values <- c()
  for(i in 1:length(supercluster_ranks[[1]])){
    
    avg_values <- append(avg_values,mean(c(supercluster_ranks[[1]][i],supercluster_ranks[[2]][i],supercluster_ranks[[3]][i])))
    
  }
  
  names(avg_values) <- names(supercluster_ranks[[1]])
  
  # Remove MT genes
  avg_values <- avg_values[!grepl("MT-",names(avg_values))]
  
  #Select top 200 genes
  names(sort(avg_values, decreasing = T)[1:num_genes])
}


all_ranks <- readRDS(paste0(dataDirectory, "genesets/cluster_ranks.rds"))

# sc1
supercluster1_ranks <- c(list(all_ranks[["A549"]][[9]]),list(all_ranks[["K562"]][[5]]),list(all_ranks[["MCF7"]][[8]]))

supercluster1_signature <- create_supercluster_signature(supercluster1_ranks, 200)

# sc2
supercluster2_ranks <- c(list(all_ranks[["A549"]][[19]]),list(all_ranks[["K562"]][[11]]),list(all_ranks[["MCF7"]][[5]]))

supercluster2_signature <- create_supercluster_signature(supercluster2_ranks, 200)

# sc3
supercluster3_ranks <- c(list(all_ranks[["A549"]][[14]]),list(all_ranks[["K562"]][[9]]),list(all_ranks[["MCF7"]][[13]]))

supercluster3_signature <- create_supercluster_signature(supercluster3_ranks, 200)

###################

supercluster_signatures <- list(supercluster1_signature,supercluster2_signature,supercluster3_signature)

names(supercluster_signatures) <- c("supercluster1_signature","supercluster2_signature","supercluster3_signature")

saveRDS(supercluster_signatures, paste0(dataDirectory, "genesets/rac_supercluster_signatures.rds"))

