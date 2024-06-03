source("/data/ruoffcj/projects/aucell_scoring/aucell_thresholding.R")
library(DESeq2)
library(tidyverse)
library(ggpubr)
library(GSVA)
library(scales)

# install.packages('corto')

#################################################################################
all_control_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/control_rac_supercluster_signatures.rds")

drug_classes <- list()
drug_classes <- append(drug_classes, list(c("BRAFi_3", "BRAFi_2","BRAFi_1","MGH_Ribociclib","MGH_Alpelisib","Sorafenib","Sorafenib_2","MK2206","Tipifarnib_1","Tipifarnib_2","Selinexor","Vismodegib","Lapatinib")))
drug_classes <- append(drug_classes, list(c("Anti-PD1","Anti-PD1 +- Anti-CTLA4", "Anti-PD1_2", "Anti-PD1_3","Anti-PD1_4")))
drug_classes <- append(drug_classes, list(c("Bevacizumab","Bevacizumab_2","Bevacizumab_3","Bevacizumab_4","Trastuzumab","Trastuzumab_2","Trastuzumab_3","Trastuzumab_4","Trastuzumab_5","Cetuximab","Rituximab")))
names(drug_classes) <- c("small_molecules","immunotherapy","mAb")

enlight_response_data <- read.csv("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/enlight_data/drug_response_classifications_with_type.csv")

rna_seq_datasets <- c("Trastuzumab_5","Bevacizumab_3", "Selinexor", "Vismodegib", "MGH_Ribociclib", "MGH_Alpelisib", "BRAFi_1", "Anti-PD1", "Anti-PD1_2", "Anti-PD1_5")

# enlight_response_data <- enlight_response_data %>%
#   filter(Dataset %in% rna_seq_datasets)#!Dataset %in% drug_classes$immunotherapy)

enlight_response_data <- enlight_response_data %>%
  filter(Dataset %in% drug_classes$small_molecules)

all_drugs <- unique(enlight_response_data$Dataset)

if("Anti-PD1 +- Anti-CTLA4" %in% all_drugs){
  all_drugs <- all_drugs[-which(all_drugs == "Anti-PD1 +- Anti-CTLA4")]  
}


# Pre load all data
all_data <- list()
for(curr_drug in all_drugs){
  
  curr_data <- read.csv(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/enlight_data/", curr_drug,".csv"), row.names = 1)
  
  #Remove rows that are NAs
  curr_data <- curr_data[!is.na(curr_data[,1]),]
  
  #Normalize data if data are integers
  if(curr_data[which(curr_data != 0)[1],1]%%1 == 0){
    
    meta <- enlight_response_data %>% 
      filter(Dataset == curr_drug) %>% 
      dplyr::select(Sample.ID,Response) %>% 
      column_to_rownames("Sample.ID")
    
    #reorder columns in data to match metadata
    curr_data <- curr_data[, rownames(meta)]
    
    all(colnames(curr_data) %in% rownames(meta))
    all(colnames(curr_data) == rownames(meta))
    
    curr_data <- round(curr_data)
    
    dds <- DESeqDataSetFromMatrix(countData = curr_data, colData = meta, design = ~ Response)
    
    dds <- estimateSizeFactors(dds)
    
    curr_data <- counts(dds, normalized=TRUE) 
    
  }
  all_data[[curr_drug]] <- curr_data
}



all_counts_vector <- c()
for(curr_control_signatures in all_control_signatures){
  
  all_signatures <- list(curr_control_signatures)
  
  all_sample_scores <- matrix(NA,nrow=0,ncol=(length(all_signatures)+1))
  
  for(curr_drug in all_drugs){
    cat(curr_drug, "\n")
    
    curr_data <- all_data[[curr_drug]]
    
    ssgsea_res <- gsva(as.matrix(curr_data), all_signatures, method="ssgsea")
    
    ssgsea_res <- t(apply(ssgsea_res, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))
    
    colnames(ssgsea_res) <- colnames(curr_data)
    temp <- data.frame(t(ssgsea_res)) %>% 
      rownames_to_column()
    
    
    if(ncol(temp) == ncol(all_sample_scores)){
      all_sample_scores <- rbind(all_sample_scores,temp)  
      cat("  added\n")
    }
    
  }
  
  colnames(all_sample_scores)[1] <- "Sample.ID"
  
  enlight_with_scores <- merge(enlight_response_data, all_sample_scores, by="Sample.ID")
  
  df <- enlight_with_scores %>% 
    pivot_longer(!c(Sample.ID,Dataset,Response,cancer_type), names_to = "geneset",values_to = "score")
  
  names1 <- sapply(as.character(df$geneset), FUN=function(x) gsub("_", " ",x))
  names2 <- gsub("type1 ","",names1)
  names3 <- gsub("supercluster","supercluster ",names2)
  df$geneset <- str_to_title(names3)
  
  df$drug <- sapply(df$Dataset, function(x) strsplit(x, "_")[[1]][1])
  
  
  R_df <- df %>% 
    dplyr::select(cancer_type, drug,geneset,score, Response) %>% 
    mutate(cohort = paste0(cancer_type, " + ", drug), cancer_type=NULL, drug=NULL) %>% 
    group_by(cohort,Response,geneset) %>% 
    summarise(mean_R = mean(score), .groups="drop") %>% 
    filter(Response == "Responder") %>% 
    select(cohort,geneset,mean_R)
  
  NR_df <- df %>% 
    dplyr::select(cancer_type, drug,geneset,score, Response) %>% 
    mutate(cohort = paste0(cancer_type, " + ", drug), cancer_type=NULL, drug=NULL) %>% 
    group_by(cohort,Response,geneset) %>% 
    summarise(mean_NR = mean(score), .groups="drop") %>% 
    filter(Response == "Non-responder") %>% 
    select(cohort,geneset,mean_NR)
  
  df <- df %>% 
    dplyr::select(cancer_type, drug,geneset,score, Response) %>% 
    mutate(cohort = paste0(cancer_type, " + ", drug), cancer_type=NULL, drug=NULL) %>% 
    group_by(cohort,Response,geneset) %>% 
    summarise(mean = mean(score),.groups= "drop") %>% 
    select(cohort,geneset) %>% 
    distinct()
  
  
  df <- merge(df,R_df, by=c("cohort","geneset"))
  df <- merge(df,NR_df, by=c("cohort","geneset"))
  
  temp <- df %>% 
    group_by(cohort) %>% 
    mutate("residual"=mean_NR - mean_R) %>% 
    ungroup() %>% 
    select(cohort,residual) %>% 
    distinct()
  
  
  
  count <- 0
  group1 <- c()
  group2 <- c()
  for(i in unique(temp$cohort)){
    
    curr_residuals <- temp %>% 
      filter(cohort == i) %>% 
      pull(residual)
    
    if(sum(curr_residuals > 0) > 0){
      count <- count + 1
      group1 <- append(group1,i)
    } else{
      group2 <- append(group2,i)
    }
  }
  
  
  all_counts_vector <- append(all_counts_vector, count)
}

saveRDS(all_counts_vector, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/ENLIGHT_scatter_null_distribution_rac_supercluster.rds")



all_counts_vector <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/ENLIGHT_scatter_null_distribution_rac_supercluster.rds")


plot(density(all_counts_vector))


sum(all_counts_vector>=6)/length(all_control_signatures)

