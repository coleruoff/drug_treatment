source("source/cole_functions.R")
library(org.EcK12.eg.db)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(tidyverse)

ecoli_AMR_genesets_up <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/ecoli_data/ecoli_AMR_genesets_up.rds")
ecoli_AMR_genesets_down <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/ecoli_data/ecoli_AMR_genesets_down.rds")

ecoli_AMR_consensus_geneset_up <- find_consensus_geneset(ecoli_AMR_genesets_up,2)
ecoli_AMR_consensus_geneset_down <- find_consensus_geneset(ecoli_AMR_genesets_down,2)


#Read in conversion table from OMA browser
human_ecoli_conversion <- read.table("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/ortholog_mapping/human_ecoli_source_id.txt", fill = T, sep='\t')[,1:3]
colnames(human_ecoli_conversion) <- c("HUMAN", "ECOLI", "MAPPING")
human_ecoli_conversion$HUMAN <- gsub("\\..*", "", human_ecoli_conversion$HUMAN)


upregulated_ecoli_human_orthologs_table <- human_ecoli_conversion %>%
  filter(ECOLI %in% ecoli_AMR_consensus_geneset_up) 

downregulated_ecoli_human_orthologs_table <- human_ecoli_conversion %>%
  filter(ECOLI %in% ecoli_AMR_consensus_geneset_down) 



#Add ECOLI id to symbol conversion
ecoli_genename_geneid_conversion_table <- read.table("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/ecoli_data/ecoli_genename_geneid_conversion_table.txt", fill = T)
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

saveRDS(upregulated_ecoli_human_orthologs_table, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/orthologs/upregulated_ecoli_human_orthologs_table.rds")
saveRDS(downregulated_ecoli_human_orthologs_table, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/orthologs/downregulated_ecoli_human_orthologs_table.rds")






# 
# # Merge all conversion tables
# colnames(human_ecoli_conversion) <- paste0(colnames(human_ecoli_conversion),"_ENTREZID")
# colnames(human_entrez_symbol_conversion) <- paste0("HUMAN_",colnames(human_entrez_symbol_conversion))
# colnames(ecoli_entrez_symbol_conversion) <- paste0("ECOLI_",colnames(ecoli_entrez_symbol_conversion))
# 
# temp <- merge(human_ecoli_conversion,human_entrez_symbol_conversion, by="HUMAN_ENTREZID")
# temp2 <- merge(temp,ecoli_entrez_symbol_conversion, by="ECOLI_ENTREZID")
# 
# ecoli_human_symbol_conversion <- temp2 %>% 
#   select("ECOLI_SYMBOL","HUMAN_SYMBOL")
# 
# # Remove rows with any NAs
# rows_to_keep <- apply(ecoli_human_symbol_conversion, 1, function(x) sum(is.na(x)) == 0)
# ecoli_human_symbol_conversion <- ecoli_human_symbol_conversion[rows_to_keep,]
# 
# 
# ecoli_human_symbol_conversion <- ecoli_human_symbol_conversion[!duplicated(ecoli_human_symbol_conversion$HUMAN_SYMBOL),]
# 
# 
# upregulated_ecoli_human_orthologs_symbol <- ecoli_human_symbol_conversion$HUMAN_SYMBOL
# names(upregulated_ecoli_human_orthologs_symbol) <- ecoli_human_symbol_conversion$ECOLI_SYMBOL
# 
# length(upregulated_ecoli_human_orthologs_entrez)
# length(upregulated_ecoli_human_orthologs_symbol)
# 
# length(unique(upregulated_ecoli_human_orthologs_symbol))
# 
# 
# 
# 
# saveRDS(upregulated_ecoli_human_orthologs_symbol, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/orthologs/upregulated_ecoli_human_orthologs_symbol.rds")
# 
# 
# 
# length(upregulated_ecoli_human_orthologs_symbol)






