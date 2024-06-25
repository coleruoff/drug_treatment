args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
setwd(args[1])
library(xlsx)

#################################################################################

df <- readRDS(paste0(dataDirectory, "dinstag_data_scores.rds"))

write.xlsx(df, file=paste0(dataDirectory, "supplementary_tables/table_S4.xlsx"), sheetName="dinstag_sample_scores", row.names=FALSE)



