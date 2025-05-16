args = commandArgs(trailingOnly=TRUE)
# dataDirectory <- paste0(args[1],"final_data/")
setwd(args[1])

source("revision_scripts/drug_treatment_functions.R")

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"
################################################################################

cell_lines <- c("A549","K562","MCF7")
# cell_lines <- c("MCF7")

all_data <- readRDS(paste0(dataDirectory, "processed_data/sciplex_data/all_cell_lines_data.rds"))

for(curr_cell_line in cell_lines){
  
  cat(curr_cell_line,"\n")
  
  #Read in cell line data
  data <- all_data[[curr_cell_line]]
  
  data <- data[,data$treatment_stage=="post"]
  
  Idents(data) <- data$Cluster
  
  de_res <- FindAllMarkers(data,
                           logfc.threshold = 0,
                           return.thresh = 1,
                           only.pos = FALSE)
  
  saveRDS(de_res, paste0(dataDirectory, "de_results/", curr_cell_line, "_cluster_de.rds"))
}


# curr_cell_line <- "A549"
# 
# de_res <- readRDS(paste0(dataDirectory, "de_results/", curr_cell_line, "_cluster_de.rds"))
# 
# 
# 
