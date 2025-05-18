args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

source("final_scripts/drug_treatment_functions.R")

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"

################################################################################

cell_lines <- c("A549","K562","MCF7")

all_data <- readRDS(paste0(dataDirectory, "processed_data/sciplex_data/all_cell_lines_data.rds"))

for(curr_cell_line in cell_lines) {
  
  cat(curr_cell_line,"\n")
  
  #Read in cell line data
  data <- all_data[[curr_cell_line]]
  
  data <- data[,data$treatment_stage=="post"]
  
  Idents(data) <- data$rac
  
  de_res <- FindAllMarkers(data)
  
  saveRDS(de_res, paste0(dataDirectory, "de_results/", curr_cell_line, "_global_rac_de.rds"))
}







