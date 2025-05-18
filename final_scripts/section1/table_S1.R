args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

source("final_scripts/drug_treatment_functions.R")

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"



if (file.exists(paste0(dataDirectory, "supplementary_tables/table_S1.xlsx"))) {
  #Delete file if it exists
  file.remove(paste0(dataDirectory, "supplementary_tables/table_S1.xlsx"))
}

if(!file.exists(paste0(dataDirectory,"supplementary_tables/"))){
  dir.create(paste0(dataDirectory,"supplementary_tables/")) 
}

################################################################################
# Create supplementary tables

# drug resistance signatures
drug_resistance_signature <- readRDS(paste0(dataDirectory, "genesets/resistance_signature.rds"))

# global rac signatures
all_signatures <- readRDS(paste0(dataDirectory, "genesets/global_rac_signatures.rds"))
names(all_signatures) <- paste0(names(all_signatures), "_global_rac_signature")

# supercluster signatures
supercluster_signatures <- readRDS(paste0(dataDirectory, "genesets/supercluster_up_signatures.rds"))
names(supercluster_signatures) <- paste0(names(supercluster_signatures), "_up")

supercluster_signatures_down <- readRDS(paste0(dataDirectory, "genesets/supercluster_down_signatures.rds"))
names(supercluster_signatures_down) <- paste0(names(supercluster_signatures_down), "_down")

all_signatures_final <- c(drug_resistance_signature,all_signatures,supercluster_signatures,supercluster_signatures_down)

write.xlsx(genelist_to_table(all_signatures_final), file=paste0(dataDirectory, "supplementary_tables/table_S1.xlsx"), sheetName="gene_signatures", append=F, row.names=FALSE)


