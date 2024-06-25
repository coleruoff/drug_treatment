args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])
library(readxl)
library(xlsx)
library(tidyverse)

supercluster_top_tf_list <- readRDS(paste0(dataDirectory, "genesets/rac_supercluster_top_tf_list.rds"))
supercluster_bottom_tf_list <- readRDS(paste0(dataDirectory, "genesets/rac_supercluster_bottom_tf_list.rds"))

df <- list()
for(i in 1:2){
  df[[paste0("supercluster",i,"_top_tfs")]] <- supercluster_top_tf_list[[i]]
  df[[paste0("supercluster",i,"_bottom_tfs")]] <- supercluster_bottom_tf_list[[i]]
  
}

max_len <- max(lengths(df))

df <- lapply(df, FUN = function(x) append(x,rep("",max_len-length(x))))

df <- data.frame(df)

write.xlsx(df, file=paste0(dataDirectory, "supplementary_tables/table_S3.xlsx"), sheetName="supercluster_tfs", row.names=FALSE)

################################################################################
crispr_ko_data <- read_xlsx("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/gottesman_crispr_data/CRISPR_KO_summary_negative.xlsx", sheet = 5)

top_genes <- list()
for(j in 1:2){
  
  curr_top_tfs <- crispr_ko_data %>% 
    filter(`neg|score` < quantile(crispr_ko_data$`neg|score`, probs = .25) & id %in% supercluster_top_tfs[[j]]) %>% 
    pull(id)
  
  top_genes <- append(top_genes, list(curr_top_tfs))
}

names(top_genes) <- paste0("supercluster", 1:2)

max_len <- max(lengths(top_genes))

top_genes <- lapply(top_genes, FUN = function(x) append(x,rep("",max_len-length(x))))

top_genes <- data.frame(top_genes)

write.xlsx(top_genes, file=paste0(dataDirectory, "supplementary_tables/table_S3.xlsx"), sheetName="supercluster_CRISPR_top_tfs", row.names=FALSE, append=T)
