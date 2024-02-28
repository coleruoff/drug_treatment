library(DESeq2)
library(tidyverse)
library(ggpubr)
library(GSVA)
library(scales)

# install.packages('corto')

#################################################################################
rac_type_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/global_rac_type_signatures.rds")
rac_type_signatures <- rac_type_signatures[grepl("type1",names(rac_type_signatures))]

rac_supercluster_consensus_signature <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_supercluster_consensus_signature.rds")


supercluster_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_supercluster_signatures.rds")

all_signatures <- supercluster_signatures
geneset_title <- "RAC Type 1 Supercluster Signatures"
#################################################################################
drug_classes <- list()
drug_classes <- append(drug_classes, list(c("BRAFi_3", "BRAFi_2","BRAFi_1","MGH_Ribociclib","MGH_Alpelisib","Sorafenib","Sorafenib_2","MK2206","Tipifarnib_1","Tipifarnib_2","Selinexor","Vismodegib","Lapatinib")))
drug_classes <- append(drug_classes, list(c("Anti-PD1","Anti-PD1 +- Anti-CTLA4", "Anti-PD1_2", "Anti-PD1_3","Anti-PD1_4")))
drug_classes <- append(drug_classes, list(c("Bevacizumab","Bevacizumab_2","Bevacizumab_3","Bevacizumab_4","Trastuzumab","Trastuzumab_2","Trastuzumab_3","Trastuzumab_4","Trastuzumab_5","Cetuximab","Rituximab")))
names(drug_classes) <- c("small_molecules","immunotherapy","mAb")

drug_class_df <- cbind(names(drug_classes[1]),drug_classes[[1]])
drug_class_df <- rbind(drug_class_df,cbind(names(drug_classes[2]),drug_classes[[2]]))
drug_class_df <- as.data.frame(rbind(drug_class_df,cbind(names(drug_classes[3]),drug_classes[[3]])))

colnames(drug_class_df) <- c("drug_class","Dataset")

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

all_sample_scores <- matrix(NA,nrow=0,ncol=(length(all_signatures)+1))



drug_class_df <- drug_class_df %>% 
  filter(Dataset %in% all_drugs)

for(curr_drug in all_drugs){
  cat(curr_drug, "\n")
  
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
  
  ssgsea_res <- gsva(as.matrix(curr_data), all_signatures, method="ssgsea")
  
  ssgsea_res <- t(apply(ssgsea_res, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))
  # ssgsea_res <- t(scale(t(ssgsea_res)))
  
  
  colnames(ssgsea_res) <- colnames(curr_data)
  temp <- data.frame(t(ssgsea_res)) %>% 
    rownames_to_column()
  # cancer_gene_set_score_info <- compute_bulk_normalized_gene_set_scores(gene_exp_mat = curr_data, gene_sets = all_signatures, num_controls = 100, num_bins = 10, q_thresh=0.95, gene_universe=NULL, use_median=F)
  #Subtract background geneset scores from foreground geneset scores
  # normalized_gene_set_score_mat <- cancer_gene_set_score_info$fg - cancer_gene_set_score_info$bg
  # temp <- data.frame(t(normalized_gene_set_score_mat)) %>% 
  #   rownames_to_column()
  
  dim(all_sample_scores)
  
  dim(temp)
  if(ncol(temp) == ncol(all_sample_scores)){
    all_sample_scores <- rbind(all_sample_scores,temp)  
    cat("  added\n")
  }
  
}

colnames(all_sample_scores)[1] <- "Sample.ID"

enlight_with_scores <- merge(enlight_response_data, all_sample_scores, by="Sample.ID")

df <- enlight_with_scores %>% 
  pivot_longer(!c(Sample.ID,Dataset,Response,cancer_type), names_to = "geneset",values_to = "score")

df <- merge(df,drug_class_df, by="Dataset",all.x=T)

names1 <- sapply(as.character(df$geneset), FUN=function(x) gsub("_", " ",x))
names2 <- gsub("type1 ","",names1)
names3 <- gsub("supercluster","supercluster ",names2)
df$geneset <- str_to_title(names3)

df$drug <- sapply(df$Dataset, function(x) strsplit(x, "_[0-9]")[[1]][1])


R_df <- df %>% 
  dplyr::select(cancer_type, drug,geneset,score, Response,drug_class) %>% 
  mutate(cohort = paste0(cancer_type, " + ", drug), cancer_type=NULL, drug=NULL) %>% 
  group_by(cohort,Response,geneset,drug_class) %>% 
  summarise(mean_R = mean(score), .groups="drop") %>% 
  filter(Response == "Responder") %>% 
  dplyr::select(cohort,geneset,mean_R,drug_class)

NR_df <- df %>% 
  dplyr::select(cancer_type, drug,geneset,score, Response,drug_class) %>% 
  mutate(cohort = paste0(cancer_type, " + ", drug), cancer_type=NULL, drug=NULL) %>% 
  group_by(cohort,Response,geneset,drug_class) %>% 
  summarise(mean_NR = mean(score), .groups="drop") %>% 
  filter(Response == "Non-responder") %>% 
  dplyr::select(cohort,geneset,mean_NR,drug_class)

df <- df %>% 
  dplyr::select(cancer_type, drug,geneset,score, Response,drug_class) %>% 
  mutate(cohort = paste0(cancer_type, " + ", drug), cancer_type=NULL, drug=NULL) %>% 
  group_by(cohort,Response,geneset,drug_class) %>% 
  summarise(mean = mean(score),.groups= "drop") %>% 
  dplyr::select(cohort,geneset,drug_class) %>% 
  distinct()


df <- merge(df,R_df, by=c("cohort","geneset","drug_class"))
df <- merge(df,NR_df, by=c("cohort","geneset","drug_class"))



library(RColorBrewer)
n <- length(unique(df$cohort))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))


colors_to_use <- sample(col_vector,n)


p <- ggplot(df)+
  geom_point(aes(x=mean_R, y=mean_NR,color=cohort,shape=geneset), size=4)+
  facet_wrap(~drug_class)+
  scale_color_manual(values=colors_to_use,name="Cohort")+
  scale_shape_manual(values=c(15, 16), name="Signature")+
  xlim(0, 1)+
  ylim(0, 1)+
  geom_abline(intercept = 0,slope = 1)+
  xlab("Mean Responder Score")+
  ylab("Mean Non-Responder Score")+
  ggtitle("Mean RAC Type 1 Supercluster Signature Scores Across Cancer Treatment Cohorts")
  


p
png("/data/ruoffcj/projects/drug_treatment/figures/ENLIGHT_scatter.png",
    width=10, height = 8, units="in",res=300)


print(p)

dev.off()




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

count/17
