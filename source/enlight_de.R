library(DESeq2)


enlight_response_data <- read.csv("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/enlight_data/drug_response_classifications_with_type.csv")

# rename enlight samples for ease of use
rownames(enlight_response_data) <- make.names(enlight_response_data$Sample.ID, unique = TRUE)

all_drugs <- unique(enlight_response_data$Dataset)

all_non_responder_up_genes <- list()

for(curr_drug in all_drugs){
  cat(curr_drug, "\n")
  curr_data <- read.csv(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/enlight_data/", curr_drug,".csv"), row.names = 1)
  cat(paste(curr_data[1,]))
}

for(curr_drug in all_drugs){
  cat(curr_drug, "\n")
  curr_data <- read.csv(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/enlight_data/", curr_drug,".csv"), row.names = 1)
  
  sum(is.na(curr_data))
  
  dim(curr_data)
  
  curr_data <- curr_data[!is.na(curr_data[,1]),]
  
  colnames(curr_data) <- make.names(colnames(curr_data), unique = TRUE)
  
  # get names of responders and non-responders
  coldata <- enlight_response_data %>% 
    filter(Dataset == curr_drug) %>% 
    select(Response)
  
  coldata <- data.frame(coldata)
  
  coldata$Response <- factor(coldata$Response)
  
  # coldata
  
  dds <- DESeqDataSetFromMatrix(countData = round(curr_data),
                                colData = coldata,
                                design = ~ Response)
  
  # all(colnames(curr_data) == rownames(coldata))
  
  
  
  dds
  
  keep <- rowSums(counts(dds)) >= 10
  
  dds <- dds[keep,]
  
  dds$condition <- relevel(dds$Response, ref="Responder")
  
  dds <- DESeq(dds)
  
  res <- results(dds)
  # res
  
  res <- as.data.frame(as.matrix(res))
  # do DE between two types of samples
  
  
  # Get upregulated in non-responder genes and save as vector. Add to list of vectors
  curr_upreg_genes <- res %>% 
    filter(log2FoldChange > 1 & padj < 0.05) %>% 
    arrange(desc(log2FoldChange)) %>% 
    rownames_to_column() %>% 
    pull(rowname)
  
  all_non_responder_up_genes <- append(all_non_responder_up_genes, list(curr_upreg_genes))
  
}



# Find common genes across all drugs' non-responders


# See if this signature is more highly expressed in trapnell resistant cells


