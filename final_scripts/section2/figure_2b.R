args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

source("final_scripts/drug_treatment_functions.R")

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"
# plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/"
################################################################################
# Plotting for supercluster enrichment

supercluster_signatures <- readRDS(paste0(dataDirectory, "genesets/supercluster_up_signatures.rds"))

names(supercluster_signatures) <- c("Supercluster 1","Supercluster 2","Supercluster 3")

final_plot <- plot_enrichment_ora(supercluster_signatures)


jpeg(paste0(plotDirectory,"figure_2b.jpg"), width=300, height = 180, units = "mm", res = 600)
print(final_plot)
dev.off()


################################################################################



