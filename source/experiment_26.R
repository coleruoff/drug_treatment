source("source/read_in_all_cell_lines.R")
library(Seurat)

cell_lines <- c("A549","K562","MCF7")

supercluster1_components <- c(9,5,8)
supercluster2_components <- c(14,9,13)

###############################################################################

curr_cell_line <- cell_lines[1]


supercluster1_unique_genes <- list()
supercluster2_unique_genes <- list()
for(curr_cell_line in cell_lines){
  
  data <- all_data[[curr_cell_line]]
  
  
  
  
}
