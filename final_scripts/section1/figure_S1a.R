args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

source("final_scripts/drug_treatment_functions.R")

set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"
# plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/"

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
#################################################################################

resistance_signature <- readRDS(paste0(dataDirectory, "genesets/resistance_signature.rds"))

names(resistance_signature) <- "Resistance Signature"

final_plot <- plot_enrichment_ora(resistance_signature)

jpeg(paste0(plotDirectory,"figure_S1a.jpg"), width=250, height=80, units = "mm", res = 600)
print(final_plot)
dev.off()



