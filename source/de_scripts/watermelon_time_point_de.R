
library(Seurat)

watermelon_data <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/watermelon_data/watermelon_pc9_processed.rds")


Idents(watermelon_data) <- watermelon_data$time_point


watermelon_de <- FindAllMarkers(watermelon_data)


saveRDS(watermelon_de, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/watermelon_time_point_de.rds")