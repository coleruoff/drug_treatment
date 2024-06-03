args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
setwd(args[1])

source("final_scripts/drug_treatment_functions.R")
library(xlsx)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

if(!file.exists(paste0(dataDirectory,"supplementary_tables/"))){
  dir.create(paste0(dataDirectory,"supplementary_tables/")) 
}

if (file.exists(paste0(dataDirectory, "supplementary_tables/supplementary_table1.xlsx"))) {
  #Delete file if it exists
  file.remove(paste0(dataDirectory, "supplementary_tables/supplementary_table1.xlsx"))
}

################################################################################
# Create supplementary tables

# All resistant genesets
all_resistant_genesets <- readRDS(paste0(dataDirectory, "genesets/all_resistant_genesets.rds"))

write.xlsx(genelist_to_table(all_resistant_genesets), file=paste0(dataDirectory, "supplementary_tables/supplementary_table1.xlsx"), sheetName="all_resistant_genesets", row.names=FALSE)

# Common drug resistance signatures
drug_resistance_signatures <- readRDS(paste0(dataDirectory, "genesets/drug_resistance_signatures.rds"))

write.xlsx(genelist_to_table(drug_resistance_signatures), file=paste0(dataDirectory, "supplementary_tables/supplementary_table1.xlsx"), sheetName="drug_resistance_signatures", append=T, row.names=FALSE)

#hallmarks
hallmarks <- readRDS(paste0(dataDirectory, "genesets/hallmarks.rds"))

write.xlsx(genelist_to_table(hallmarks), file=paste0(dataDirectory, "supplementary_tables/supplementary_table1.xlsx"), sheetName="cancer_hallmarks", append=T, row.names=FALSE)

#ITH MPs
MPs <- readRDS(paste0(dataDirectory, "genesets/ITH_meta_programs.rds"))

write.xlsx(genelist_to_table(MPs), file=paste0(dataDirectory, "supplementary_tables/supplementary_table1.xlsx"), sheetName="ITH_meta_programs", append=T, row.names=FALSE)
