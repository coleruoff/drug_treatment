library(Seurat)


A549.data <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/A549_processed_filtered.rds")

K562.data <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/K562_processed_filtered.rds")

MCF7.data <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/MCF7_processed_filtered.rds")

cell_lines <- c("A549","K562","MCF7")

merged_trapnell <- merge(A549.data, y = c(K562.data, MCF7.data),
                         add.cell.ids = cell_lines,
                         project = "merged_trapnell")



temp <- merge(A549.data,K562.data)

temp2 <- merge(temp,MCF7.data)




ncol(A549.data)+ncol(K562.data)+ncol(MCF7.data)
ncol(merged_trapnell)


saveRDS(merged_trapnell, "/data/CDSL_hannenhalli/Cole/drug_treatment/processed_data/sciPlex_data/merged_trapnell.rds")

Idents(merged_trapnell) <- merged_trapnell$cell_type

cell_line.de <- FindAllMarkers(merged_trapnell)

saveRDS(cell_line.de, "/data/CDSL_hannenhalli/Cole/drug_treatment/processed_data/sciPlex_data/cell_line_de.rds")