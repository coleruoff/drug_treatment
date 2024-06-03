library(Seurat)

curr_cell_line <- "MCF7"

data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

curr_markers <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/", curr_cell_line, "_cluster_all_markers_de.rds"))

curr_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

data <- ScaleData(data)


DoHeatmap(subset(data, downsample = 100), 
          features = top10$gene, 
          assay = "RNA")

