args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

source("final_scripts/drug_treatment_functions.R")
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"
# plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/"

################################################################################

supercluster_signatures <- readRDS(paste0(dataDirectory, "genesets/supercluster_up_signatures.rds"))

yeast_human_orthologs_up <- readRDS(paste0(dataDirectory, "genesets/yeast_human_orthologs_up.rds"))

all_plots <- list()
overlap_signatures <- list()
for(i in c(2,3)){
  cat(i,"\n")
  curr_geneset <- intersect(supercluster_signatures[[i]], yeast_human_orthologs_up$yeast_human_orthologs_up)
  
  overlap_signatures <- append(overlap_signatures,list(curr_geneset))
  
}


names(overlap_signatures) <- c("Supercluster 2","Supercluster 3")

final_plot <- plot_enrichment_ora(overlap_signatures)

jpeg(paste0(plotDirectory,"figure_4b.jpg"), width=250, height = 140, units = "mm", res = 600)
print(final_plot)
dev.off()






