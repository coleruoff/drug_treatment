args = commandArgs(trailingOnly=TRUE)
# dataDirectory <- paste0(args[1],"final_data/")
setwd(args[1])

source("revision_scripts/drug_treatment_functions.R")

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"
################################################################################

cell_lines <- readRDS(paste0(dataDirectory, "oren_cell_lines_names.rds"))

# cell_lines <- cell_lines[-1]
# cell_lines <- cell_lines[-7]

# Run DE between day 0 and day 10 for each cell line
all_cell_line_de_res <- list()
for(curr_cell_line in cell_lines){
  file_name <- paste0(dataDirectory,"processed_data/oren_data/", curr_cell_line, "_processed.rds")
  
  data <- readRDS(file_name)

  data$timepoint <- ifelse(data$sample_type == "naive", 0, 10)

  Idents(data) <- data$timepoint
  
  dim(data)

  curr_de_res <- FindAllMarkers(data,
                                logfc.threshold = 0,
                                return.thresh = 1,
                                only.pos = FALSE)

  all_cell_line_de_res <- append(all_cell_line_de_res, list(curr_de_res))

}



# name each DE result based on cell line
names(all_cell_line_de_res) <- paste0(cell_lines, "_de_res")

saveRDS(all_cell_line_de_res, paste0(dataDirectory, "de_results/oren_all_cell_line_de_res.rds"))



