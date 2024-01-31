if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.EcK12.eg.db")
BiocManager::install("org.Sc.sgd.db")

library(org.EcK12.eg.db)
library(org.Hs.eg.db)
library(org.Sc.sgd.db)
library(tidyverse)

convert_ecoli_to_human_genes <- function(genes){
  
  
  #Convert ecoli genes to Entrez IDs
  genes_conversion_ecoli <- AnnotationDbi::select(org.EcK12.eg.db, 
                             keys = genes,
                             columns = c("ENTREZID", "SYMBOL"),
                             keytype = "SYMBOL")
  
  ecoli_entrez <- genes_conversion_ecoli$ENTREZID
  
  ecoli_entrez <- ecoli_entrez[!is.na(ecoli_entrez)]
  
  num_genes_lost <- length(genes) - length(ecoli_entrez)
  
  #Select corresponding human entrez IDs
  human_to_ecoli_conversion_table <- read.table("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/ecoli_data/human_ecoli_entrez.txt", fill = T, sep='\t')[,1:3]
  colnames(human_to_ecoli_conversion_table) <- c("HUMAN", "ECOLI", "MAPPING")
  
  human_entrez <- human_to_ecoli_conversion_table %>% 
    filter(ECOLI %in% ecoli_entrez) %>% 
    pull(HUMAN)
  
  #Convert human entrez IDs to SYMBOLs
  genes_conversion_human <- AnnotationDbi::select(org.Hs.eg.db, 
                             keys = human_entrez,
                             columns = c("ENTREZID", "SYMBOL"),
                             keytype = "ENTREZID")
  #Return human symbols
  human_symbols <- genes_conversion_human$SYMBOL
  
  human_symbols <- (human_symbols[!is.na(human_symbols)])
  
  return(human_symbols)
}


convert_human_to_ecoli_genes <- function(genes){
  
  
  #Convert human genes to Entrez IDs
  genes_conversion_human <- AnnotationDbi::select(org.Hs.eg.db, 
                                                  keys = genes,
                                                  columns = c("ENTREZID", "SYMBOL"),
                                                  keytype = "SYMBOL")
  
  human_entrez <- genes_conversion_human$ENTREZID
  
  human_entrez <- human_entrez[!is.na(human_entrez)]
  
  num_genes_lost <- length(genes) - length(human_entrez)
  
  #Select corresponding ecoli entrez IDs
  human_to_ecoli_conversion_table <- read.table("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/ecoli_data/human_ecoli_entrez.txt", fill = T)
  colnames(human_to_ecoli_conversion_table) <- c("HUMAN", "ECOLI", "MAPPING")
  
  ecoli_entrez <- human_to_ecoli_conversion_table %>% 
    filter(HUMAN %in% human_entrez) %>% 
    pull(ECOLI)
  
  #Convert ecoli entrez IDs to SYMBOLs
  genes_conversion_ecoli <- AnnotationDbi::select(org.EcK12.eg.db, 
                                                  keys = ecoli_entrez,
                                                  columns = c("ENTREZID", "SYMBOL"),
                                                  keytype = "ENTREZID")
  #Return human symbols
  ecoli_symbols <- genes_conversion_ecoli$SYMBOL
  
  ecoli_symbols <- unique(ecoli_symbols[!is.na(ecoli_symbols)])
  
  return(ecoli_symbols)
}


convert_yeast_to_human_genes <- function(genes){
  
  #Convert yeast genes to Entrez IDs
  human_yeast_conversion <- read.table("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/yeast_data/human_candida_auris_source_id.txt", fill = T, sep='\t')[,1:3]
  colnames(human_yeast_conversion) <- c("HUMAN","YEAST", "MAPPING")
  
  
  human_genes <- human_yeast_conversion %>% 
    filter(YEAST %in% genes) %>% 
    pull(HUMAN)
  
  
  human_genes <- sapply(human_genes, FUN = function(x) strsplit(x, "\\.")[[1]][1])
  
  keytypes(org.Hs.eg.db)
  keys(org.Hs.eg.db, keytype = "ENSEMBL")
  
  human_ensembl_symbol <- AnnotationDbi::select(org.Hs.eg.db, 
                                                keys = human_genes,
                                                columns = c("ENSEMBL", "SYMBOL"),
                                                keytype = "ENSEMBL")
  
  
  return(human_ensembl_symbol$SYMBOL)
}


convert_human_to_yeast_genes <- function(genes){
  
  
  #Convert human genes to Entrez IDs
  genes_conversion_human <- AnnotationDbi::select(org.Hs.eg.db, 
                                                  keys = genes,
                                                  columns = c("ENTREZID", "SYMBOL"),
                                                  keytype = "SYMBOL")
  
  human_entrez <- genes_conversion_human$ENTREZID
  
  human_entrez <- human_entrez[!is.na(human_entrez)]
  
  num_genes_lost <- length(genes) - length(human_entrez)
  
  #Select corresponding yeast entrez IDs
  human_to_yeast_conversion_table <- read.table("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/yeast_data/human_yeast_entrez.txt", fill = T)
  colnames(human_to_yeast_conversion_table) <- c("HUMAN", "YEAST", "MAPPING")
  
  yeast_entrez <- human_to_yeast_conversion_table %>% 
    filter(HUMAN %in% human_entrez) %>% 
    pull(YEAST)
  
  #Convert yeast entrez IDs to SYMBOLs
  genes_conversion_yeast <- AnnotationDbi::select(org.Sc.sgd.db, 
                                                  keys = yeast_entrez,
                                                  columns = c("ENTREZID", "SYMBOL"),
                                                  keytype = "ENTREZID")
  #Return human symbols
  yeast_symbols <- genes_conversion_yeast$SYMBOL
  
  yeast_symbols <- unique(yeast_symbols[!is.na(yeast_symbols)])
  
  return(yeast_symbols)
}


