args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
setwd(args[1])
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"
# plotDirectory <- "/data/ruoffcj/projects/drug_treatment/final_figures"

cell_lines <- c("A549","K562","MCF7")

cell_line_universes <- list()
for(curr_cell_line in cell_lines){

  cat(curr_cell_line, "\n")

  #find genes expressed in > 1% of cells
  data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/",curr_cell_line,"_processed_filtered.rds"))
  gene_universe <- names(apply(data@assays$RNA@data, 1, FUN = function(x) sum(x>0) >.01*ncol(data)))

  cell_line_universes <- append(cell_line_universes, list(gene_universe))
}

names(cell_line_universes) <- cell_lines

gene_universe <- intersect(cell_line_universes[[1]], cell_line_universes[[2]])
gene_universe <- intersect(gene_universe, cell_line_universes[[3]])

saveRDS(gene_universe, paste0(dataDirectory, "cell_line_gene_universe_intersection.rds"))
saveRDS(cell_line_universes, paste0(dataDirectory, "cell_line_universes.rds"))