library(DESeq2)
library(tidyverse)
library(affy)

enlight_response_data <- read.csv("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/enlight_data/drug_response_classifications_with_type.csv")

enlight_datsets <- unique(enlight_response_data$Dataset)

enlight_data_types <- c("TPM", "TPM","FPKM","Microarray",
  "Microarray","Microarray","Microarray","readcounts",
  "Microarray","Microarray","Microarray","Microarray",
  "Microarray","readcounts","Microarray","Microarray",
  "FPKM","Microarray","TPM","Microarray",
  "Microarray","Microarray","Microarray","Microarray",
  "Microarray","Microarray","Microarray","readcounts",
  "FPKM","TPM")

enlight_data_types <- cbind(enlight_datsets,enlight_data_types)


sample_responses <- enlight_response_data %>% 
  select(Sample.ID,Response) 

all_drugs <- unique(enlight_response_data$Dataset)

data_list <- list()
curr_drug <- all_drugs[4]

for(curr_drug in all_drugs){
  
  cat(curr_drug, "\n")
  curr_data <- read.csv(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/enlight_data/", curr_drug,".csv"), row.names = 1)
  
  head(curr_data)
  
  row_to_remove <- lapply(rownames(curr_data), FUN = function(x){
    if(sum(is.na(curr_data[x,])) > 0){
      return(T)
    } else{
      return(F)
    }
  })
  
  rows_to_remove <- unlist(row_to_remove)
  
  curr_data <- curr_data[!rows_to_remove,]
  
  # curr_data[curr_data < 0] <- 0
  
  #Removes samples from data matrix that arent in enlight response data table
  curr_data <- curr_data[,gsub("\\.","-", colnames(curr_data)) %in% sample_responses$Sample.ID]
  
  curr_sample_responses <- sample_responses %>% 
    filter(Sample.ID %in% gsub("\\.","-", colnames(curr_data)))
  
  
  dds <- DESeqDataSetFromMatrix(countData = round(curr_data), colData = curr_sample_responses, design = ~ Response)
  dds <- estimateSizeFactors(dds)
  
  normalized_counts <- counts(dds, normalized=TRUE)
  
  data_list[[curr_drug]] <- normalized_counts
}


saveRDS(data_list, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/enlight_data/normalized_data_list.rds")




