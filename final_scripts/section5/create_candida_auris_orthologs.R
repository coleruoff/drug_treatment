args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
setwd(args[1])
library(readxl)
library(tidyverse)
library(org.Hs.eg.db)
library(org.Sc.sgd.db)
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

################################################################################

yeast_de <- read_xlsx(paste0(dataDirectory, "raw_data/yeast_data/candida_auris_de.XLSX"), col_names = T)
colnames(yeast_de) <- yeast_de[1,]
yeast_de <- yeast_de[-1,]

upregulated_genes <- yeast_de %>% 
  filter(as.numeric(padj) < 0.05 & up_down == "Up") %>% 
  pull(gene_id)

downregulated_genes <- yeast_de %>% 
  filter(as.numeric(padj) < 0.05 & up_down == "Down") %>% 
  pull(gene_id)


length(upregulated_genes)

#Read in conversion table from OMA browser
human_yeast_conversion <- read.table(paste0(dataDirectory, "raw_data/yeast_data/human_candida_auris_source_id.txt"), fill = T, sep='\t')[,1:3]
colnames(human_yeast_conversion) <- c("HUMAN","YEAST", "MAPPING")

# Trim ensembl gene names
human_yeast_conversion$HUMAN <- gsub("\\..*", "", human_yeast_conversion$HUMAN)

#Select genes that map to one yeast gene and those that are present in yeast upregulated list
upregulated_yeast_human_orthologs_table <- human_yeast_conversion %>% 
  filter(YEAST %in% upregulated_genes)

downregulated_yeast_human_orthologs_table <- human_yeast_conversion %>% 
  filter(YEAST %in% downregulated_genes)


# Add HUMAN ensembl to symbol conversion
human_ensembl_symbol_conversion <- AnnotationDbi::select(org.Hs.eg.db, 
                                                         keys = unique(c(upregulated_yeast_human_orthologs_table$HUMAN,downregulated_yeast_human_orthologs_table$HUMAN)),
                                                         columns = c("ENSEMBL", "SYMBOL"),
                                                         keytype = "ENSEMBL")


upregulated_yeast_human_orthologs_table <- distinct(merge(upregulated_yeast_human_orthologs_table, human_ensembl_symbol_conversion, by.x="HUMAN",by.y="ENSEMBL"))

downregulated_yeast_human_orthologs_table <- distinct(merge(downregulated_yeast_human_orthologs_table, human_ensembl_symbol_conversion, by.x="HUMAN",by.y="ENSEMBL"))


upregulated_yeast_human_orthologs_table <- upregulated_yeast_human_orthologs_table %>% 
  dplyr::select(HUMAN,SYMBOL,MAPPING,YEAST)

downregulated_yeast_human_orthologs_table <- downregulated_yeast_human_orthologs_table %>% 
  dplyr::select(HUMAN,SYMBOL,MAPPING,YEAST)

colnames(upregulated_yeast_human_orthologs_table) <- c("HUMAN_ENSEMBL","HUMAN_SYMBOL","MAPPING","YEAST_GENEID")
colnames(downregulated_yeast_human_orthologs_table) <- c("HUMAN_ENSEMBL","HUMAN_SYMBOL","MAPPING","YEAST_GENEID")


saveRDS(upregulated_yeast_human_orthologs_table, paste0(dataDirectory, "genesets/upregulated_yeast_human_orthologs_table.rds"))
saveRDS(downregulated_yeast_human_orthologs_table, paste0(dataDirectory, "genesets/downregulated_yeast_human_orthologs_table.rds"))


#read in yeast-human orthologs and save as RDS
yeast_human_orthologs_up <- unique(upregulated_yeast_human_orthologs_table$HUMAN_SYMBOL[!is.na(upregulated_yeast_human_orthologs_table$HUMAN_SYMBOL)])
yeast_human_orthologs_up <- list("yeast_human_orthologs_up"=yeast_human_orthologs_up)
saveRDS(yeast_human_orthologs_up, paste0(dataDirectory, "genesets/yeast_human_orthologs_up.rds"))

yeast_human_orthologs_down <- unique(downregulated_yeast_human_orthologs_table$HUMAN_SYMBOL[!is.na(downregulated_yeast_human_orthologs_table$HUMAN_SYMBOL)])
yeast_human_orthologs_down <- list("yeast_human_orthologs_down"=yeast_human_orthologs_down)
saveRDS(yeast_human_orthologs_down, paste0(dataDirectory, "genesets/yeast_human_orthologs_down.rds"))




















# 
# 
# ################################################################################
# # Create list of 1:m genes
# upregulated_yeast_human_orthologs_ensembl_many <- human_yeast_conversion %>% 
#   filter((MAPPING == "1:m") & YEAST %in% upregulated_genes) %>% 
#   pull(HUMAN)
# 
# 
# names(upregulated_yeast_human_orthologs_ensembl_many) <- human_yeast_conversion %>% 
#   filter((MAPPING == "1:m") & YEAST %in% upregulated_genes) %>% 
#   pull(YEAST)
# 
# # Trim ensembl gene names
# upregulated_yeast_human_orthologs_ensembl_many <- gsub("\\..*", "", upregulated_yeast_human_orthologs_ensembl_many)
# 
# many_list <- c()
# many_list_names <- c()
# for(i in unique(upregulated_yeast_human_orthologs_ensembl_many)){
#   cat(i, "\n")
#   
#   many_list <- append(many_list,i)
#   
#   many_list_names <- append(many_list_names,paste(names(upregulated_yeast_human_orthologs_ensembl_many)[upregulated_yeast_human_orthologs_ensembl_many == i],
#         collapse = "/"))
#   
# }
# 
# names(many_list) <- many_list_names
# 
# upregulated_yeast_human_orthologs_ensembl <- c(upregulated_yeast_human_orthologs_ensembl,many_list)
# ################################################################################
# 
# length(upregulated_yeast_human_orthologs_ensembl)
# length(unique(upregulated_yeast_human_orthologs_ensembl))
# 
# length(names(upregulated_yeast_human_orthologs_ensembl))
# length(unique(names(upregulated_yeast_human_orthologs_ensembl)))
# 
# 
# saveRDS(upregulated_yeast_human_orthologs_ensembl, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/orthologs/upregulated_yeast_human_orthologs_ensembl.rds")
# 
# 
# keytypes <- keytypes(org.Sc.sgd.db)
# 
# for(curr in keytypes){
#   cat(keys(org.Sc.sgd.db, keytype = curr)[1:5], "\n")
# }
# keys(org.Sc.sgd.db, keytype = "ENSEMBL")
# 
# human_ensembl_symbol_conversion <- AnnotationDbi::select(org.Hs.eg.db, 
#                                               keys = upregulated_yeast_human_orthologs_ensembl,
#                                               columns = c("ENSEMBL", "SYMBOL"),
#                                               keytype = "ENSEMBL")
# 
# 
# 
# human_yeast_conversion_ensembl <- cbind(upregulated_yeast_human_orthologs_ensembl,names(upregulated_yeast_human_orthologs_ensembl))
# colnames(human_yeast_conversion_ensembl) <- c("HUMAN_ENSEMBL","YEAST_ENSEMBL")
# 
# colnames(human_ensembl_symbol_conversion) <- paste0("HUMAN_",colnames(human_ensembl_symbol_conversion))
# 
# 
# # Merge conversion tables and make distinct
# temp <- merge(human_yeast_conversion_ensembl,human_ensembl_symbol_conversion, by="HUMAN_ENSEMBL")
# temp <- temp[!duplicated(temp$HUMAN_SYMBOL),]
# 
# 
# # Create named list for gene symbols
# upregulated_yeast_human_orthologs_symbol <- temp$HUMAN_SYMBOL
# names(upregulated_yeast_human_orthologs_symbol) <- temp$YEAST_ENSEMBL
# 
# 
# #Remove NAs
# upregulated_yeast_human_orthologs_symbol <- upregulated_yeast_human_orthologs_symbol[!is.na(upregulated_yeast_human_orthologs_symbol)]
# 
# saveRDS(upregulated_yeast_human_orthologs_symbol, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/orthologs/upregulated_yeast_human_orthologs_symbol.rds")
# 
# 



