library(Seurat)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"

################################################################################

watermelon_data <- readRDS(paste0(dataDirectory,"processed_data/watermelon_data/watermelon_pc9_processed.rds"))

watermelon_data_old <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/watermelon_data/watermelon_pc9_processed.rds"))

Idents(watermelon_data) <- watermelon_data$time_point

watermelon_de <- FindAllMarkers(watermelon_data)

saveRDS(watermelon_de, paste0(dataDirectory,"de_results/watermelon_time_point_de.rds"))

dim(watermelon_de)

watermelon_data[["RNA3"]] <- as(object = watermelon_data[["RNA"]], Class = "Assay")

watermelon_de_og %>% 
  filter(cluster == 14 & p_val_adj < 0.05 & avg_log2FC > 0) %>% 
  pull(gene)


watermelon_de %>% 
  filter(cluster == 14 & p_val_adj < 0.05 & avg_log2FC > 0) %>% 
  pull(gene)

