



cell_lines <- c("A549","K562","MCF7")

#################################################################################
# Read in Trapnell data
#################################################################################
curr_cell_line <- cell_lines[1]

data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

data_ad <- Convert(from = data, to = "anndata", filename = paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/",curr_cell_line,"_processed_filtered.h5ad"))
data_ad