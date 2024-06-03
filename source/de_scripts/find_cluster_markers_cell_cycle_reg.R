library(Seurat)

curr_cell_line <- "A549"

data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/",curr_cell_line,"_processed_filtered_cell_cycle.rds"))

Idents(data) <- data$Cluster

de.results <- FindAllMarkers(data)

saveRDS(de.results, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/",curr_cell_line,"_cluster_de_cell_cycle.rds"))


