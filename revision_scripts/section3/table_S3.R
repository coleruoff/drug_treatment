args = commandArgs(trailingOnly=TRUE)
# dataDirectory <- paste0(args[1],"final_data/")
setwd(args[1])

source("revision_scripts/drug_treatment_functions.R")

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"
plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_figures/"
################################################################################


df <- readRDS(paste0(dataDirectory, "dinstag_data_scores.rds"))


write.xlsx(df, file=paste0(dataDirectory, "supplementary_tables/table_S3.xlsx"), sheetName="dinstag_data_scores", append=F, row.names=FALSE)
