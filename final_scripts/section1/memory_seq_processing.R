args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")

library(Seurat)
library(tidyverse)
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"


if(!file.exists(paste0(dataDirectory,"processed_data/memorySeq_data/"))){
  dir.create(paste0(dataDirectory,"processed_data/"))
  dir.create(paste0(dataDirectory,"processed_data/memorySeq_data/")) 
}

################################################################################

MDA_table <- read.csv(paste0(dataDirectory, "raw_data/memorySeq_data/GSE151375_full_MDA.csv"))

MDA.data <- MDA_table %>% 
  filter(controls == "Controls") %>% 
  select(gene_id, GeneSymbol, sampleID, total_counts)   %>% 
  pivot_wider(names_from = sampleID, values_from = total_counts)    %>% 
  select(-gene_id)  %>% 
  group_by(GeneSymbol)


MDA.data <- aggregate(. ~ GeneSymbol, MDA.data, sum) %>% 
  column_to_rownames("GeneSymbol")

saveRDS(MDA.data, paste0(dataDirectory,"processed_data/memorySeq_data/MDA_control_counts.rds"))


################################################################################

WM989_table <- read.csv(paste0(dataDirectory,"raw_data/memorySeq_data/GSE151375_full_WM989.csv"))

WM989.data <- WM989_table %>% 
  filter(controls == "Controls") %>% 
  select(gene_id, GeneSymbol, sampleID, total_counts)   %>% 
  pivot_wider(names_from = sampleID, values_from = total_counts)    %>% 
  select(-gene_id)  %>% 
  group_by(GeneSymbol)


WM989.data <- aggregate(. ~ GeneSymbol, WM989.data, sum) %>% 
  column_to_rownames("GeneSymbol")


saveRDS(WM989.data, paste0(dataDirectory,"processed_data/memorySeq_data/WM989_control_counts.rds"))