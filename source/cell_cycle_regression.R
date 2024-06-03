library(Seurat)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

args <- commandArgs(trailingOnly=TRUE)
curr_cell_line <- args[1]

data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/",curr_cell_line, "_processed_filtered.rds"))

data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

data <- ScaleData(data, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(data))

data <- FindVariableFeatures(data)

data <- RunPCA(data, features = VariableFeatures(data))

# When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
data <- RunPCA(data, features = c(s.genes, g2m.genes))
# DimPlot(marrow)


saveRDS(data, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/",curr_cell_line,"_processed_filtered_cell_cycle.rds"))