


pan_stem_genes <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/pan_stem_genes.rds")




csc_table <- read.csv("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/csc_table.csv")


csc_table <- csc_table %>% 
  select(Gene,Marker.Type) %>% 
  filter(Gene %in% rownames(data))

all_csc_markers <- csc_table$Gene



curr_cell_line <- "A549"

data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

lung_csc_markers <- csc_table %>% 
  filter(grepl("Lung", Marker.Type) & Gene %in% rownames(data)) %>% 
  pull(Gene)


curr_cell_line <- "MCF7"

data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))


breast_csc_markers <- csc_table %>% 
  filter(grepl("Breast", Marker.Type) & Gene %in% rownames(data)) %>% 
  pull(Gene)


curr_cell_line <- "K562"

data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))


blood_csc_markers <- csc_table %>% 
  filter(grepl("Blood", Marker.Type) | grepl("Leukemia", Marker.Type) & Gene %in% rownames(data)) %>% 
  pull(Gene)



saveRDS(lung_csc_markers, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/lung_csc_markers.rds")
saveRDS(breast_csc_markers, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/breast_csc_markers.rds")
saveRDS(blood_csc_markers, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/blood_csc_markers.rds")






