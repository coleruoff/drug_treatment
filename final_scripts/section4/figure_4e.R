args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])
library(DESeq2)
library(tidyverse)
library(ggpubr)
library(GSVA)
library(scales)
set.seed(42)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"
# plotDirectory <- "/data/ruoffcj/projects/drug_treatment/"

#################################################################################

supercluster_signatures <- readRDS(paste0(dataDirectory,"genesets/rac_supercluster_signatures.rds"))

all_signatures <- supercluster_signatures
geneset_title <- "RAC Supercluster Signatures"

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

enlight_response_data <- read.csv(paste0(dataDirectory, "raw_data/enlight_data/drug_response_classifications_with_type.csv"))

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
  
  curr_data <- read.csv(paste0(dataDirectory,"raw_data/enlight_data/", curr_drug,".csv"), row.names = 1)
  
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
  
  #Filter out genes with low/zero variance
  data_to_use <- as.matrix(curr_data)
  row_vars <- rowVars(data_to_use)
  data_to_use <- data_to_use[which(row_vars > 1e-10),]
  
  # Score samples
  gsvapar <- ssgseaParam(data_to_use, all_signatures)
  ssgsea_res <- gsva(gsvapar)
  
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

regression_df <- df

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



# library(RColorBrewer)
# n <- length(unique(df$cohort))
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# pie(rep(1,n), col=sample(col_vector, n))
# 
# 
# colors_to_use <- sample(col_vector,n)

saveRDS(df, paste0(dataDirectory, "enlight_scatter_df.rds"))

p <- ggplot(df)+
  geom_point(aes(x=mean_R, y=mean_NR,color=cohort,shape=geneset), size=4)+
  scale_shape_manual(values=c(15, 16,17), name="Signature")+
  scale_x_continuous(expand=c(0,0),limits = c(0,1))+
  scale_y_continuous(expand=c(0,0),limits = c(0,1))+
  geom_abline(intercept = 0,slope = 1)+
  xlab("Mean Responder Score")+
  ylab("Mean Non-Responder Score")+
  ggtitle("Mean Supercluster Signature Scores Across Cancer Treatment Cohorts")+
  guides(color=guide_legend(title="Cohort"))
  
# scale_color_manual(values=colors_to_use)+

ann_text<-data.frame( 
  x = 4, y = 20, 
  label = "geeks for geeks") 

p+geom_text(data=ann_text,aes(x=x,y=y,label=label),size=10)

ann_text<-data.frame( 
    x = 4, y = 20, 
    label = "geeks for geeks"
) 


png(paste0(plotDirectory,"figure_4e.png"),
    width=10, height = 8, units="in",res=300)


print(p)

dev.off()

df <-  regression_df %>% 
  dplyr::select(cancer_type, drug,geneset,score, Response) %>% 
  mutate(cohort = paste0(cancer_type, " + ", drug), cancer_type=NULL, drug=NULL)


for(curr_signature in unique(df$geneset)){
  temp <- df %>% 
    filter(geneset == curr_signature)
  
  temp$Response <- ifelse(temp$Response == "Non-responder", 1, 0)
  
  model <- glm(Response ~ score + cohort,
               data=temp,
               family = "binomial")
  
  model_summmary <- summary(model)
  
  score_pval <- model_summmary$coefficients[2,4]
  
  cat(paste0(curr_signature," p-value: ", sprintf("%.6f",score_pval),"\n"))
  
}


summary(temp$score)

temp$predictions <- predict(model,temp,type="response")

ggplot(temp, aes(x = score, y = predictions)) +
  geom_smooth()

