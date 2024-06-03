setwd("/data/ruoffcj/projects/drug_treatment/")
library(tidyverse)
library(Seurat)
set.seed(42)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

#################################################################################

files <- list.files(paste0(dataDirectory, "processed_data/watermelon_data/"))

files <- files[1:8]

# Remove cell lines without cycling annotations
files <- files[-which(files %in% c("BT474_BREAST_processed.rds", "EFM192A_BREAST_processed.rds","MMACSF_SKIN_processed.rds"))]

all_de_list <- list()

for(curr_file in files){
  cat(curr_file, "\n")
  watermelon_data <- readRDS(paste0(dataDirectory, "processed_data/watermelon_data/", curr_file))
  
  watermelon_data$cycling <- ifelse(grepl("Non",watermelon_data$sample_pool),"noncycling",
                                    ifelse(grepl("Cyc",watermelon_data$sample_pool), "cycling","pre-treatment"))
  
  watermelon_data$time_point <- ifelse(grepl("_10_", watermelon_data$sample_pool), 1, 0)
  
  Idents(watermelon_data) <- watermelon_data$cycling
  
  de_res <- FindMarkers(watermelon_data, ident.1 = "cycling",ident.2 = "noncycling", test.use = "MAST")
  
  all_de_list <- append(all_de_list, list(de_res))
}

names(all_de_list) <- mapply(files, FUN = function(x) gsub("_processed.rds","",x))


saveRDS(all_de_list, paste0(dataDirectory, "de_results/oren_cell_lines_cycling_vs_non_de.rds"))

all_de_list <- readRDS(paste0(dataDirectory,"de_results/oren_cell_lines_cycling_vs_non_de.rds"))


cycling_de_genes <- list()
for(i in all_de_list){
  temp <- i %>% 
    filter(p_val_adj < 0.05 & avg_log2FC > 0) %>% 
    arrange(desc(avg_log2FC)) %>% 
    rownames_to_column("gene") %>% 
    pull(gene)
  
  cycling_de_genes <- append(cycling_de_genes,list(temp))
  
}

consensus_cycling_signature <- list("cycling_signature"=find_consensus_geneset(cycling_de_genes, 4))

data <- all_data[["A549"]]

cells_AUC <- AUCell_run(data@assays$RNA$data, consensus_cycling_signature)


saveRDS(consensus_cycling_signature, paste0(dataDirectory, "genesets/consensus_cycling_signature.rds"))


