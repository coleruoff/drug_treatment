source("source/human_to_organism_conversion_function.R")
source("source/cole_functions.R")

#################################################################################

ecoli_stress_genesets <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/ecoli_data/ecoli_various_antibiotic_stress.rds")

ecoli_stress_orthologs <- list()
for(i in ecoli_stress_genesets){
  ecoli_stress_orthologs <- append(ecoli_stress_orthologs, list(convert_ecoli_to_human_genes(i)))
}

names(ecoli_stress_orthologs) <- paste0(names(ecoli_stress_genesets),"_ortholog")

saveRDS(ecoli_stress_orthologs, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ecoli_stress_orthologs.rds")

#################################################################################

ecoli_AMR_genesets <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/ecoli_data/ecoli_AMR_genesets.rds")

ecoli_AMR_genesets_orthologs <- list()
for(i in ecoli_AMR_genesets){
  ecoli_AMR_genesets_orthologs <- append(ecoli_AMR_genesets_orthologs, list(convert_ecoli_to_human_genes(i)))
  
}

names(ecoli_AMR_genesets_orthologs) <- paste0(names(ecoli_AMR_genesets),"_orthologs")

saveRDS(ecoli_AMR_genesets_orthologs, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ecoli_AMR_genesets_orthologs.rds")



#################################################################################
# Create AMR consensus genesets
ecoli_AMR_up <- ecoli_AMR_genesets[grepl("up",names(ecoli_AMR_genesets))]
ecoli_AMR_up_consensus <- find_consensus_geneset(ecoli_AMR_up, 4)

ecoli_AMR_down <- ecoli_AMR_genesets[grepl("down",names(ecoli_AMR_genesets))]
ecoli_AMR_down_consensus <- find_consensus_geneset(ecoli_AMR_down, 4)

#Convert to human orthologs
ecoli_AMR_up_consensus_orthologs <- convert_ecoli_to_human_genes(ecoli_AMR_up_consensus)
ecoli_AMR_down_consensus_orthologs <- convert_ecoli_to_human_genes(ecoli_AMR_down_consensus)

#Remove shared genes from up and down genesets
shared_genes <- intersect(ecoli_AMR_up_consensus_orthologs,ecoli_AMR_down_consensus_orthologs)
ecoli_AMR_up_consensus_orthologs <- ecoli_AMR_up_consensus_orthologs[!ecoli_AMR_up_consensus_orthologs %in% shared_genes]
ecoli_AMR_down_consensus_orthologs <- ecoli_AMR_down_consensus_orthologs[!ecoli_AMR_down_consensus_orthologs %in% shared_genes]


ecoli_AMR_consensus_orthologs <- list("ecoli_AMR_up_consensus_orthologs" = ecoli_AMR_up_consensus_orthologs, "ecoli_AMR_down_consensus_orthologs" = ecoli_AMR_down_consensus_orthologs)

saveRDS(ecoli_AMR_consensus_orthologs, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ecoli_AMR_consensus_orthologs.rds")
#################################################################################
# Create stress consensus genesets
ecoli_stress_up <- ecoli_stress_genesets[grepl("up",names(ecoli_stress_genesets))]
ecoli_stress_up_consensus <- find_consensus_geneset(ecoli_stress_up, 4)

ecoli_stress_down <- ecoli_stress_genesets[grepl("down",names(ecoli_stress_genesets))]
ecoli_stress_down_consensus <- find_consensus_geneset(ecoli_stress_down, 4)

#Convert to human orthologs
ecoli_stress_up_consensus_orthologs <- convert_ecoli_to_human_genes(ecoli_stress_up_consensus)
ecoli_stress_down_consensus_orthologs <- convert_ecoli_to_human_genes(ecoli_stress_down_consensus)

#Remove shared genes from up and down genesets
shared_genes <- intersect(ecoli_stress_up_consensus_orthologs,ecoli_stress_down_consensus_orthologs)
ecoli_stress_up_consensus_orthologs <- ecoli_stress_up_consensus_orthologs[!ecoli_stress_up_consensus_orthologs %in% shared_genes]
ecoli_stress_down_consensus_orthologs <- ecoli_stress_down_consensus_orthologs[!ecoli_stress_down_consensus_orthologs %in% shared_genes]

ecoli_stress_consensus_orthologs <- list("ecoli_stress_up_consensus_orthologs" = ecoli_stress_up_consensus_orthologs, "ecoli_stress_down_consensus_orthologs" = ecoli_stress_down_consensus_orthologs)

saveRDS(ecoli_stress_consensus_orthologs, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ecoli_stress_consensus_orthologs.rds")






