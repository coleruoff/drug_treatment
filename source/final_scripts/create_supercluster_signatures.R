source("source/cole_functions.R")

cell_lines <- c("A549","K562","MCF7")

################################################################################
# Create RAC supercluster signature
################################################################################

# Components of RAC supercluster
A549_active <- c(4,9,13)
K562_active <- c(5)
MCF7_active <-  c(5,8,17)

# supercluster_components <- list(A549_active, K562_active, MCF7_active)
# names(supercluster_components) <- cell_lines

signature_length <- 200

# Read in cluster signatures
A549_cluster_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/A549_cluster_signatures_",signature_length,".rds"))
K562_cluster_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/K562_cluster_signatures_",signature_length,".rds"))
MCF7_cluster_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/MCF7_cluster_signatures_",signature_length,".rds"))

# Unlist all genes from RAC supercluster component RACs
A549_unlist <- unlist(A549_cluster_signatures[A549_active])
K562_unlist <- unlist(K562_cluster_signatures[K562_active])
MCF7_unlist <- unlist(MCF7_cluster_signatures[MCF7_active])

# Get intersection
supercluster_genes <- Reduce(intersect, list(A549_unlist,K562_unlist,MCF7_unlist))

supercluster_genes <- supercluster_genes[!grepl("^MT",supercluster_genes)]

supercluster_signature <- list("supercluster_signature" = supercluster_genes)

saveRDS(supercluster_signature, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_supercluster_signature.rds")

################################################################################
# Create RAC supercluster consensus signature
################################################################################

# Number of signatures a gene must appear in to be added to consensus
k <- 2
supercluster_consensus_signature <- find_consensus_geneset(c(A549_cluster_signatures[A549_active],K562_cluster_signatures[K562_active],MCF7_cluster_signatures[MCF7_active]),k)

# supercluster_consensus_signature <- supercluster_consensus_signature[!grepl("^MT",supercluster_consensus_signature)]

supercluster_consensus_signature <- list("supercluster_consensus_signature" = supercluster_consensus_signature)

saveRDS(supercluster_consensus_signature, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_supercluster_consensus_signature.rds")

################################################################################
# Create RAC type 1 supercluster consensus signatures
################################################################################

refined_type1_signatures <- list()

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")

rac_type1_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_type1_signatures.rds")
rac_type2_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_type2_signatures.rds")

for(curr_cell_line in cell_lines){
  
  cat(curr_cell_line, "\n")
  
  curr_rac_type1_signatures <- rac_type1_signatures[grepl(curr_cell_line, names(rac_type1_signatures))]
  curr_rac_type2_signatures <- rac_type2_signatures[grepl(curr_cell_line, names(rac_type2_signatures))]
  
  # Intracluster DE Upregulated genes
  curr_intra_rac_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_rac_within_cluster_de_signatures.rds"))
  
  clusters_of_interest <- RACs[[curr_cell_line]]
  
  curr_type1_signatures <- list()
  
  for(curr_cluster in clusters_of_interest){
    
    cat(curr_cluster, "\n")
    
    type1_signature <- curr_rac_type1_signatures[[paste0(curr_cell_line,"_rac_subcluster",curr_cluster,"_1_signature")]]
    
    type2_signature <- curr_rac_type2_signatures[[paste0(curr_cell_line,"_rac_subcluster",curr_cluster,"_2_signature")]]
    
    #Find shared genes between active and inactive components
    shared_genes <- intersect(type1_signature, type2_signature)
    
    #Remove shared genes
    type1_signature <- type1_signature[!type1_signature %in% shared_genes]
    
    #Add genes from DE between active and inactive within this cluster
    
    type1_signature <- union(type1_signature, curr_intra_rac_signatures[[paste0(curr_cluster)]])
    
    curr_type1_signatures <- append(curr_type1_signatures, list(type1_signature))
    
  }
  
  names(curr_type1_signatures) <- paste0(curr_cell_line, "_", clusters_of_interest)
  
  refined_type1_signatures <- append(refined_type1_signatures, curr_type1_signatures)
}

saveRDS(refined_type1_signatures, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/refined_type1_signatures.rds")

##############

# Using signatures for type 1 cells (even though its called refined, its not)
rac_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_signatures.rds")

names1 <- gsub(pattern = "_signature","", names(rac_signatures))
names(rac_signatures) <- gsub(pattern = "cluster","", names1)

supercluster1_signatures <- rac_signatures[c(paste0("A549_",c(9)),paste0("K562_",c(5)),paste0("MCF7_",c(8)))]

supercluster1_signature <- list("supercluster1_signature" = find_consensus_geneset(supercluster1_signatures,1))

# saveRDS(supercluster1_signature, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/active_supercluster1_signature.rds")

##############

supercluster2_signatures <- rac_signatures[c(paste0("A549_",c(14)),paste0("K562_",c(9)),paste0("MCF7_",c(13)))]

supercluster2_signature <- list("supercluster2_signature" = find_consensus_geneset(supercluster2_signatures,1))

# saveRDS(supercluster2_signature, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/active_supercluster2_signature.rds")

##############

supercluster3_signatures <- rac_signatures[c(paste0("A549_",c(19)),paste0("K562_",c(19)),paste0("MCF7_",c(5)))]

supercluster3_signature <- list("supercluster3_signature" = find_consensus_geneset(supercluster3_signatures,2))
# 
# # saveRDS(supercluster3_signature, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/active_supercluster3_signature.rds")
# 
# ##############
# 
# supercluster4_signatures <- refined_type1_signatures[c(paste0("A549_",c(13,17)),paste0("MCF7_",c(17)))]
# 
# supercluster4_signature <- list("type1_supercluster4_signature" = find_consensus_geneset(supercluster4_signatures,2))
# 
# # saveRDS(supercluster4_signature, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/active_supercluster4_signature.rds")
# 
# ##############
# 
# supercluster5_signatures <- refined_type1_signatures[c(paste0("A549_",c(12,16)))]
# 
# supercluster5_signature <- list("type1_supercluster5_signature" = find_consensus_geneset(supercluster5_signatures,2))

# saveRDS(supercluster5_signature, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/active_supercluster5_signature.rds")

##############

# type1_supercluster_signatures <- c(supercluster1_signature,supercluster2_signature,supercluster3_signature,supercluster4_signature,supercluster5_signature)
supercluster_signatures <- c(supercluster1_signature,supercluster2_signature)

saveRDS(supercluster_signatures, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_supercluster_signatures.rds")


################################################################################
# Consensus for downregulated genesets
rac_type1_down_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_type1_down_signatures.rds")

names1 <- gsub(pattern = "_rac_subcluster","_", names(rac_type1_down_signatures))
names(rac_type1_down_signatures) <- gsub(pattern = "_1_down_signature","", names1)

supercluster1_signatures <- rac_type1_down_signatures[c(paste0("A549_",c(4,9)),paste0("K562_",c(4,9,5)),paste0("MCF7_",c(8,12)))]

supercluster1_signature <- list("type1_supercluster1_signature" = find_consensus_geneset(supercluster1_signatures,2))

##############

supercluster2_signatures <- rac_type1_down_signatures[c(paste0("A549_",c(19)),paste0("K562_",c(11)),paste0("MCF7_",c(5)))]

supercluster2_signature <- list("type1_supercluster2_signature" = find_consensus_geneset(supercluster2_signatures,2))

type1_supercluster_down_signatures <- c(supercluster1_signature,supercluster2_signature)

saveRDS(type1_supercluster_down_signatures, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/type1_supercluster_down_signatures.rds")


