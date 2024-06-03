dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

top_tfs <- readRDS(paste0(dataDirectory, "genesets/all_cluster_top_tfs.rds"))

# top_tfs[["A549"]] <- sapply(top_tfs[["A549"]], FUN = function(x) names(x))
# top_tfs[["K562"]] <- sapply(top_tfs[["K562"]], FUN = function(x) names(x))
# top_tfs[["MCF7"]] <- sapply(top_tfs[["MCF7"]], FUN = function(x) names(x))

for(i in 1:length(top_tfs)){
  for(j in 1:length(top_tfs[[i]])){
    top_tfs[[i]][[j]] <- names(top_tfs[[i]][[j]])
  }
}

supercluster1_tf_list <- c(list(top_tfs[["A549"]][[9]]),list(top_tfs[["K562"]][[5]]),list(top_tfs[["MCF7"]][[8]]))

supercluster1_signature <- list("supercluster1_signature" = find_consensus_geneset(supercluster1_tf_list,3))

##############

supercluster2_tf_list <- c(list(top_tfs[["A549"]][[14]]),list(top_tfs[["K562"]][[9]]),list(top_tfs[["MCF7"]][[13]]))

supercluster2_signature <- list("supercluster2_signature" = find_consensus_geneset(supercluster2_tf_list,3))

supercluster_tf_list <- c(supercluster1_signature,supercluster2_signature)

saveRDS(supercluster_tf_list, paste0(dataDirectory, "genesets/rac_supercluster_top_tf_list.rds"))

#################################################################################

bottom_tfs <- readRDS(paste0(dataDirectory, "genesets/all_cluster_bottom_tfs.rds"))

# bottom_tfs[["A549"]] <- sapply(bottom_tfs[["A549"]], FUN = function(x) names(x))
# bottom_tfs[["K562"]] <- sapply(bottom_tfs[["K562"]], FUN = function(x) names(x))
# bottom_tfs[["MCF7"]] <- sapply(bottom_tfs[["MCF7"]], FUN = function(x) names(x))

for(i in 1:length(bottom_tfs)){
  for(j in 1:length(bottom_tfs[[i]])){
    bottom_tfs[[i]][[j]] <- names(bottom_tfs[[i]][[j]])
  }
}

supercluster1_tf_list <- c(list(bottom_tfs[["A549"]][[9]]),list(bottom_tfs[["K562"]][[5]]),list(bottom_tfs[["MCF7"]][[8]]))

supercluster1_signature <- list("supercluster1_signature" = find_consensus_geneset(supercluster1_tf_list,3))

##############

supercluster2_tf_list <- c(list(bottom_tfs[["A549"]][[14]]),list(bottom_tfs[["K562"]][[9]]),list(bottom_tfs[["MCF7"]][[13]]))

supercluster2_signature <- list("supercluster2_signature" = find_consensus_geneset(supercluster2_tf_list,3))

supercluster_tf_list <- c(supercluster1_signature,supercluster2_signature)

saveRDS(supercluster_tf_list, paste0(dataDirectory, "genesets/rac_supercluster_bottom_tf_list.rds"))


