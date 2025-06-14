args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

source("final_scripts/drug_treatment_functions.R")

set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"
# plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/"

#################################################################################

supercluster_signatures <- readRDS(paste0(dataDirectory,"genesets/supercluster_up_signatures.rds"))

all_signatures <- supercluster_signatures
geneset_title <- "Supercluster Signatures"

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

enlight_response_data <- enlight_response_data %>%
  filter(Dataset %in% drug_classes$small_molecules)

all_drugs <- unique(enlight_response_data$Dataset)

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
  ssgsea_res <- gsva(gsvapar, )
  
  ssgsea_res <- t(apply(ssgsea_res, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))
  # ssgsea_res <- t(scale(t(ssgsea_res)))
  
  colnames(ssgsea_res) <- colnames(curr_data)
  temp <- data.frame(t(ssgsea_res)) %>% 
    rownames_to_column()
  
  
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

saveRDS(df, paste0(dataDirectory, "dinstag_data_scores.rds"))


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


p <- ggplot(df)+
  geom_point(aes(x=mean_R, y=mean_NR,color=cohort,shape=geneset), size=2)+
  scale_shape_manual(values=c(15,16,17,18), name="Signature")+
  scale_x_continuous(expand=c(0,0),limits = c(0,1))+
  scale_y_continuous(expand=c(0,0),limits = c(0,1))+
  geom_abline(intercept = 0, slope = 1, size = .2)+
  xlab("Mean Responder Score")+
  ylab("Mean Non-Responder Score")+
  guides(color=guide_legend(title="Cohort"))+
  theme_classic()+
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.title = element_text(size=8),
        legend.text = element_text(size=4),
        legend.key.height = unit(3,"mm"),
        legend.key.width = unit(3,"mm"),
        axis.line = element_line(linewidth=.2),
        axis.ticks = element_line(linewidth = .2))



p


jpeg(paste0(plotDirectory,"figure_3e.jpg"), width=100, height = 80, units = "mm", res = 600)
print(p)
dev.off()


head(df)



t.test(df$mean_NR-df$mean_R, alternative = "greater")




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

  cat(model_summmary$coefficients[2,1],"\n")
  cat(paste0(curr_signature," p-value: ", sprintf("%.6f",score_pval),"\n"))

}


summary(temp$score)

temp$predictions <- predict(model,temp,type="response")

ggplot(temp, aes(x = score, y = predictions)) +
  geom_smooth()







