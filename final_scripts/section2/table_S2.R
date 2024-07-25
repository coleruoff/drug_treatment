args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
setwd(args[1])

source("final_scripts/drug_treatment_functions.R")
library(xlsx)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"
################################################################################

all_signatures <- readRDS(paste0(dataDirectory, "genesets/cluster_signatures.rds"))
RACs <- readRDS(paste0(dataDirectory, "processed_data/all_RACs.rds"))

cell_lines <- c("A549","K562","MCF7")

all_rac_signatures <- list()

for(curr_cell_line in cell_lines){
  curr_signatures <- all_signatures[[curr_cell_line]]
  
  curr_rac_signatures <- curr_signatures[RACs[[curr_cell_line]]]
  
  names(curr_rac_signatures) <- paste0(curr_cell_line, "_", names(curr_rac_signatures))
  
  all_rac_signatures <- append(all_rac_signatures,curr_rac_signatures)
}

max_len <- max(lengths(all_rac_signatures))

df <- lapply(all_rac_signatures, FUN = function(x) append(x,rep("",max_len-length(x))))

df <- data.frame(df)

write.xlsx(df, file=paste0(dataDirectory, "supplementary_tables/table_S2.xlsx"), sheetName="rac_signatures", row.names=FALSE)
################################################################################

supercluster_signatures <- readRDS(paste0(dataDirectory, "genesets/rac_supercluster_signatures.rds"))

max_len <- max(lengths(supercluster_signatures))

df <- lapply(supercluster_signatures, FUN = function(x) append(x,rep("",max_len-length(x))))

df <- data.frame(df)

write.xlsx(df, file=paste0(dataDirectory, "supplementary_tables/table_S2.xlsx"), sheetName="supercluster_signatures", row.names=FALSE, append = T)

