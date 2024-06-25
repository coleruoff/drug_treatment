args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])
source("final_scripts/drug_treatment_functions.R")

source("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_scripts/drug_treatment_functions.R")
# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

#################################################################################

top_tfs <- readRDS(paste0(dataDirectory, "genesets/all_cluster_top_tfs.rds"))

supercluster_components <- readRDS(paste0(dataDirectory, "processed_data/supercluster_components.rds"))

for(i in names(top_tfs)){
  for(j in names(top_tfs[[i]])){
    top_tfs[[i]][[j]] <- names(top_tfs[[i]][[j]])
  }
}

supercluster1_tf_list <- c(list(top_tfs[["A549"]][[as.character(supercluster_components[[1]][["A549"]])]]),
                           list(top_tfs[["K562"]][[as.character(supercluster_components[[1]][["K562"]])]]),
                           list(top_tfs[["MCF7"]][[as.character(supercluster_components[[1]][["MCF7"]])]]))


supercluster1_signature <- list("supercluster1_signature" = find_consensus_geneset(supercluster1_tf_list,3))

##############

supercluster2_tf_list <- c(list(top_tfs[["A549"]][[as.character(supercluster_components[[2]][["A549"]])]]),
                           list(top_tfs[["K562"]][[as.character(supercluster_components[[2]][["K562"]])]]),
                           list(top_tfs[["MCF7"]][[as.character(supercluster_components[[2]][["MCF7"]])]]))


supercluster2_signature <- list("supercluster2_signature" = find_consensus_geneset(supercluster2_tf_list,3))

supercluster_tf_list <- c(supercluster1_signature,supercluster2_signature)

saveRDS(supercluster_tf_list, paste0(dataDirectory, "genesets/rac_supercluster_top_tf_list.rds"))

#################################################################################

bottom_tfs <- readRDS(paste0(dataDirectory, "genesets/all_cluster_bottom_tfs.rds"))

supercluster_components <- readRDS(paste0(dataDirectory, "processed_data/supercluster_components.rds"))

for(i in names(bottom_tfs)){
  for(j in names(bottom_tfs[[i]])){
    bottom_tfs[[i]][[j]] <- names(bottom_tfs[[i]][[j]])
  }
}

supercluster1_tf_list <- c(list(bottom_tfs[["A549"]][[as.character(supercluster_components[[1]][["A549"]])]]),
                           list(bottom_tfs[["K562"]][[as.character(supercluster_components[[1]][["K562"]])]]),
                           list(bottom_tfs[["MCF7"]][[as.character(supercluster_components[[1]][["MCF7"]])]]))


supercluster1_signature <- list("supercluster1_signature" = find_consensus_geneset(supercluster1_tf_list,3))

##############

supercluster2_tf_list <- c(list(bottom_tfs[["A549"]][[as.character(supercluster_components[[2]][["A549"]])]]),
                           list(bottom_tfs[["K562"]][[as.character(supercluster_components[[2]][["K562"]])]]),
                           list(bottom_tfs[["MCF7"]][[as.character(supercluster_components[[2]][["MCF7"]])]]))


supercluster2_signature <- list("supercluster2_signature" = find_consensus_geneset(supercluster2_tf_list,3))

supercluster_tf_list <- c(supercluster1_signature,supercluster2_signature)

saveRDS(supercluster_tf_list, paste0(dataDirectory, "genesets/rac_supercluster_bottom_tf_list.rds"))