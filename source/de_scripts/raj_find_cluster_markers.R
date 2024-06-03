library(Seurat)

raj_data <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/diverse_clonal_fates_data/raj_resistant_breast_processed.rds")

de_results <- FindAllMarkers(raj_data)

saveRDS(de_results, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/raj_breast_cluster_de.rds")


