args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
setwd(args[1])
library(org.EcK12.eg.db)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(tidyverse)
library(readxl)
source("final_scripts/drug_treatment_functions.R")
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"

################################################################################

AMR_de <- read_xlsx(paste0(dataDirectory, "raw_data/ecoli_data/AMR_de_genes.xlsx"))

colnames(AMR_de) <- AMR_de[1,]
AMR_de <- AMR_de[-1,]

experiments <- c("TETvs","MMCvs","IPMvs","CAZvs","KANvs","CIPvs","PMEvs","ERYvs","CHLvs")

ecoli_AMR_genesets_up <- list()
ecoli_AMR_genesets_down <- list()
for(curr_experiment in experiments){
  curr_experiment_data <- AMR_de[,grepl(curr_experiment, colnames(AMR_de)) | grepl("Gene",colnames(AMR_de))]
  
  curr_experiment_data <- curr_experiment_data[,-c(2,3)]
  
  colnames(curr_experiment_data) <- c("gene_id","log2FC","pval","padj","significant","gene_name")
  
  curr_experiment_genes_up <- curr_experiment_data %>% 
    filter(as.numeric(padj) < 0.05 & as.numeric(log2FC) > 0) %>% 
    pull(gene_id)
  
  ecoli_AMR_genesets_up <- append(ecoli_AMR_genesets_up, list(curr_experiment_genes_up))
  
  curr_experiment_genes_down <- curr_experiment_data %>% 
    filter(as.numeric(padj) < 0.05 & as.numeric(log2FC) < 0) %>% 
    pull(gene_id)
  
  ecoli_AMR_genesets_down <- append(ecoli_AMR_genesets_down, list(curr_experiment_genes_down))
}

names(ecoli_AMR_genesets_up) <- gsub("vs","",experiments)

names(ecoli_AMR_genesets_down) <- gsub("vs","",experiments)


saveRDS(ecoli_AMR_genesets_up, paste0(dataDirectory, "genesets/ecoli_AMR_genesets_up.rds"))
saveRDS(ecoli_AMR_genesets_down, paste0(dataDirectory, "genesets/ecoli_AMR_genesets_down.rds"))

ecoli_AMR_genesets_up <- readRDS(paste0(dataDirectory, "genesets/ecoli_AMR_genesets_up.rds"))
ecoli_AMR_genesets_down <- readRDS(paste0(dataDirectory, "genesets/ecoli_AMR_genesets_down.rds"))

ecoli_AMR_consensus_geneset_up <- find_consensus_geneset(ecoli_AMR_genesets_up,2)
ecoli_AMR_consensus_geneset_down <- find_consensus_geneset(ecoli_AMR_genesets_down,2)


#Read in conversion table from OMA browser
human_ecoli_conversion <- read.table(paste0(dataDirectory, "raw_data/ecoli_data/human_ecoli_source_id.txt"), fill = T, sep='\t')[,1:3]
colnames(human_ecoli_conversion) <- c("HUMAN", "ECOLI", "MAPPING")
human_ecoli_conversion$HUMAN <- gsub("\\..*", "", human_ecoli_conversion$HUMAN)


upregulated_ecoli_human_orthologs_table <- human_ecoli_conversion %>%
  filter(ECOLI %in% ecoli_AMR_consensus_geneset_up) 

downregulated_ecoli_human_orthologs_table <- human_ecoli_conversion %>%
  filter(ECOLI %in% ecoli_AMR_consensus_geneset_down) 



#Add ECOLI id to symbol conversion
ecoli_genename_geneid_conversion_table <- read.table(paste0(dataDirectory, "raw_data/ecoli_data/ecoli_genename_geneid_conversion_table.txt"), fill = T)
colnames(ecoli_genename_geneid_conversion_table) <- ecoli_genename_geneid_conversion_table[1,]
ecoli_genename_geneid_conversion_table <- ecoli_genename_geneid_conversion_table[-1,]

upregulated_ecoli_human_orthologs_table <- merge(upregulated_ecoli_human_orthologs_table, ecoli_genename_geneid_conversion_table, by.x="ECOLI",by.y="gene_id")
downregulated_ecoli_human_orthologs_table <- merge(downregulated_ecoli_human_orthologs_table, ecoli_genename_geneid_conversion_table, by.x="ECOLI",by.y="gene_id")


# Add HUMAN ensembl to symbol conversion
human_ensembl_symbol_conversion <- AnnotationDbi::select(org.Hs.eg.db, 
                                                        keys = unique(c(upregulated_ecoli_human_orthologs_table$HUMAN,downregulated_ecoli_human_orthologs_table$HUMAN)),
                                                        columns = c("ENSEMBL", "SYMBOL"),
                                                        keytype = "ENSEMBL")



upregulated_ecoli_human_orthologs_table <- distinct(merge(upregulated_ecoli_human_orthologs_table, human_ensembl_symbol_conversion, by.x="HUMAN",by.y="ENSEMBL"))
downregulated_ecoli_human_orthologs_table <- distinct(merge(downregulated_ecoli_human_orthologs_table, human_ensembl_symbol_conversion, by.x="HUMAN",by.y="ENSEMBL"))



upregulated_ecoli_human_orthologs_table <- as.data.frame(upregulated_ecoli_human_orthologs_table) %>% 
  dplyr::select(HUMAN,SYMBOL,MAPPING,ECOLI,gene_name)

downregulated_ecoli_human_orthologs_table <- as.data.frame(downregulated_ecoli_human_orthologs_table) %>% 
  dplyr::select(HUMAN,SYMBOL,MAPPING,ECOLI,gene_name)


colnames(upregulated_ecoli_human_orthologs_table) <- c("HUMAN_ENSEMBL","HUMAN_SYMBOL","MAPPING","ECOLI_GENEID","ECOLI_SYMBOL")
colnames(downregulated_ecoli_human_orthologs_table) <- c("HUMAN_ENSEMBL","HUMAN_SYMBOL","MAPPING","ECOLI_GENEID","ECOLI_SYMBOL")

saveRDS(upregulated_ecoli_human_orthologs_table, paste0(dataDirectory, "genesets/upregulated_ecoli_human_orthologs_table.rds"))
saveRDS(downregulated_ecoli_human_orthologs_table, paste0(dataDirectory, "genesets/downregulated_ecoli_human_orthologs_table.rds"))


#read in ecoli-human orthologs and save as RDS
ecoli_human_orthologs_up <- unique(upregulated_ecoli_human_orthologs_table$HUMAN_SYMBOL[!is.na(upregulated_ecoli_human_orthologs_table$HUMAN_SYMBOL)])
ecoli_human_orthologs_up <- list("ecoli_human_orthologs_up"=ecoli_human_orthologs_up)


ecoli_human_orthologs_down <- unique(downregulated_ecoli_human_orthologs_table$HUMAN_SYMBOL[!is.na(downregulated_ecoli_human_orthologs_table$HUMAN_SYMBOL)])
ecoli_human_orthologs_down <- list("ecoli_human_orthologs_down"=ecoli_human_orthologs_down)


#Find and remove intersecting genes
ecoli_human_orthologs_intersect <- intersect(ecoli_human_orthologs_up$ecoli_human_orthologs_up, ecoli_human_orthologs_down$ecoli_human_orthologs_down)

ecoli_human_orthologs_up$ecoli_human_orthologs_up <- ecoli_human_orthologs_up$ecoli_human_orthologs_up[!ecoli_human_orthologs_up$ecoli_human_orthologs_up %in% ecoli_human_orthologs_intersect]
ecoli_human_orthologs_down$ecoli_human_orthologs_down <- ecoli_human_orthologs_down$ecoli_human_orthologs_down[!ecoli_human_orthologs_down$ecoli_human_orthologs_down %in% ecoli_human_orthologs_intersect]

ecoli_human_orthologs_intersect <- list("ecoli_human_orthologs_intersect"=ecoli_human_orthologs_intersect)

saveRDS(ecoli_human_orthologs_intersect, paste0(dataDirectory, "genesets/ecoli_human_orthologs_intersect.rds"))
saveRDS(ecoli_human_orthologs_up, paste0(dataDirectory, "genesets/ecoli_human_orthologs_up.rds"))
saveRDS(ecoli_human_orthologs_down, paste0(dataDirectory, "genesets/ecoli_human_orthologs_down.rds"))

################################################################################
# Create table of all ecoli-ortholog genes, regardless if up/down regulated
################################################################################
human_ecoli_conversion <- read.table(paste0(dataDirectory, "raw_data/ecoli_data/human_ecoli_source_id.txt"), fill = T, sep='\t')[,1:3]
colnames(human_ecoli_conversion) <- c("HUMAN", "ECOLI", "MAPPING")
human_ecoli_conversion$HUMAN <- gsub("\\..*", "", human_ecoli_conversion$HUMAN)


ecoli_genename_geneid_conversion_table <- read.table(paste0(dataDirectory, "raw_data/ecoli_data/ecoli_genename_geneid_conversion_table.txt"), fill = T)
colnames(ecoli_genename_geneid_conversion_table) <- ecoli_genename_geneid_conversion_table[1,]
ecoli_genename_geneid_conversion_table <- ecoli_genename_geneid_conversion_table[-1,]

ecoli_human_orthologs_table <- merge(human_ecoli_conversion, ecoli_genename_geneid_conversion_table, by.x="ECOLI",by.y="gene_id")

human_ensembl_symbol_conversion <- AnnotationDbi::select(org.Hs.eg.db, 
                                                         keys = ecoli_human_orthologs_table$HUMAN,
                                                         columns = c("ENSEMBL", "SYMBOL"),
                                                         keytype = "ENSEMBL")



ecoli_human_orthologs_table <- distinct(merge(ecoli_human_orthologs_table, human_ensembl_symbol_conversion, by.x="HUMAN",by.y="ENSEMBL"))


ecoli_human_orthologs_table <- as.data.frame(ecoli_human_orthologs_table) %>% 
  dplyr::select(HUMAN,SYMBOL,MAPPING,ECOLI,gene_name)


colnames(ecoli_human_orthologs_table) <- c("HUMAN_ENSEMBL","HUMAN_SYMBOL","MAPPING","ECOLI_GENEID","ECOLI_SYMBOL")

saveRDS(ecoli_human_orthologs_table, paste0(dataDirectory, "genesets/ecoli_human_orthologs_table.rds"))



