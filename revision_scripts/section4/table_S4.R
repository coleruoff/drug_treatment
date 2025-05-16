args = commandArgs(trailingOnly=TRUE)
# dataDirectory <- paste0(args[1],"final_data/")
setwd(args[1])

source("revision_scripts/drug_treatment_functions.R")

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"
plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_figures/"
################################################################################

df <- list()

yeast_human_orthologs_up <- readRDS(paste0(dataDirectory, "genesets/yeast_human_orthologs_up.rds"))

ecoli_human_orthologs_up <- readRDS(paste0(dataDirectory, "genesets/ecoli_human_orthologs_up.rds"))

df[["yeast_upregulated_orthologs"]] <- yeast_human_orthologs_up$yeast_human_orthologs_up
df[["ecoli_upregulated_orthologs"]] <- ecoli_human_orthologs_up$ecoli_human_orthologs_up

#################################################################################
# Yeast
yeast_overlap <- readRDS(paste0(dataDirectory, "genesets/yeast_sc_overlap.rds"))
df[["yeast_supercluster2_overlap"]] <- yeast_overlap[[1]]
df[["yeast_supercluster3_overlap"]] <- yeast_overlap[[2]]

#############################
# E. coli

ecoli_overlap <- readRDS(paste0(dataDirectory, "genesets/ecoli_sc_overlap.rds"))
df[["ecoli_supercluster2_overlap"]] <- ecoli_overlap[[1]]


write.xlsx(genelist_to_table(df), file=paste0(dataDirectory, "supplementary_tables/table_S4.xlsx"), sheetName="supercluster_ortholog_genes", append=F, row.names=FALSE)
