library(tidyverse)
library(ComplexHeatmap)
library(ggpubr)
source("../survival_analysis/score_gene_expression.R")


drug_classes <- list()
drug_classes <- append(drug_classes, list(c("MGH_Alpelisib","Sorafenib","Sorafenib_2","MK2206","Tipifarnib_1","Tipifarnib_2","Selinexor","Vismodegib")))
drug_classes <- append(drug_classes, list(c("Anti-PD1","Anti-PD1 +- Anti-CTLA4", "Anti-PD1_2", "Anti-PD1_3","Anti-PD1_4")))
drug_classes <- append(drug_classes, list(c("Bevacizumab","Bevacizumab_2","Bevacizumab_3","Bevacizumab_4","Trastuzumab","Trastuzumab_2","Trastuzumab_3","Trastuzumab_4","Trastuzumab_5","Cetuximab","Rituximab")))

names(drug_classes) <- c("small_molecules","immunotherapy","mAb")

# saveRDS(drug_classes, "C:/Users/ruoffcj/Documents/projects/drug_treatment/data/enlight_data/drug_classes.rds")


therapy_response_OR_by_drug <- function(data, drug_name, genesets_to_use, x){
  
  enlight_response_data <- read.csv("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/enlight_data/drug_response_classifications_with_type.csv")
  
  enlight_response_data %>% 
    filter
  
  #Score samples for given genesets
  cancer_gene_set_score_info <- compute_bulk_normalized_gene_set_scores(gene_exp_mat = data, gene_sets = genesets_to_use, num_controls = 5, num_bins = 10, q_thresh=0.95, gene_universe=NULL, use_median=F)
  
  #Subtract background geneset scores from foreground geneset scores
  normalized_gene_set_score_mat <- cancer_gene_set_score_info$fg - cancer_gene_set_score_info$bg
  
  
  df <- data.frame(matrix(NA, nrow=0, ncol=3))
  
  for(i in 1:length(genesets_to_use)){
    sorted_sample_scores <- sort(normalized_gene_set_score_mat[i,], decreasing = T)
    
    # x <- .1
    top_x_samples <- sorted_sample_scores[sorted_sample_scores > quantile(sorted_sample_scores, probs = (1-x))]
    bottom_x_samples <- sorted_sample_scores[sorted_sample_scores < quantile(sorted_sample_scores, probs = x)]
    
    names(top_x_samples) <- gsub("\\.","-",names(top_x_samples))
    names(bottom_x_samples) <- gsub("\\.","-",names(bottom_x_samples))
    
    #Number of Non-responder in top x%
    top_NR <- enlight_response_data %>% 
      filter(Sample.ID %in% names(top_x_samples)) %>% 
      pull(Response) 
    
    top_NR <- sum(top_NR=="Non-responder")
    
    #Number of responder in top x%
    top_R <- enlight_response_data %>% 
      filter(Sample.ID %in% names(top_x_samples)) %>% 
      pull(Response)
    
    top_R <- sum(top_R=="Responder")
    
    
    #Number of Non-responder in bottom x%
    bottom_NR <- enlight_response_data %>% 
      filter(Sample.ID %in% names(bottom_x_samples)) %>% 
      pull(Response)
    
    bottom_NR <- sum(bottom_NR=="Non-responder")
    
    #Number of responder in bottom x%
    bottom_R <- enlight_response_data %>% 
      filter(Sample.ID %in% names(bottom_x_samples)) %>% 
      pull(Response)
    
    bottom_R <- sum(bottom_R=="Responder")
    
    
    top_NR_percent <- top_NR/(top_NR+top_R)
    top_R_percent <- top_R/(top_NR+top_R)
    
    bottom_NR_percent <- bottom_NR/(bottom_NR+bottom_R)
    bottom_R_percent <- bottom_R/(bottom_NR+bottom_R)
    
    
    
    NR_OR <- (top_NR_percent+1)/(bottom_NR_percent+1)
    
    R_OR <- (bottom_R_percent+1)/(top_R_percent+1)
    
    curr_geneset_name <- names(genesets_to_use)[i]
    
    
    df <- rbind(df,c(curr_geneset_name,NR_OR, R_OR,drug_name))
  }
  
  colnames(df) <- c("geneset","NR_OR","R_OR","drug")
  
  df <- data.frame(df)
  
  df$NR_OR <- as.numeric(df$NR_OR)
  df$R_OR <- as.numeric(df$R_OR)
  df$geneset <- factor(df$geneset, levels = names(genesets_to_use))
  
  
  return(df)
}

therapy_response_OR <- function(data_list, genesets_to_use, x){
  
  enlight_response_data <- read.csv("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/enlight_data/drug_response_classifications_with_type.csv")
  or_df <- data.frame(matrix(NA, nrow=0, ncol=5))
  scores_df <- data.frame(matrix(NA, nrow=length(geneset_to_use), ncol=0))
  
  for(i in 1:length(data_list)){
    cat(names(data_list)[i], "\n")
    
    #Select current experiemnt data, drug, drug class, and cancer type
    curr_data <- data_list[[i]]
    curr_drug <- names(data_list)[i]
    curr_cancer_type <- enlight_response_data %>% 
      filter(Dataset == curr_drug) %>% 
      pull(cancer_type) %>% 
      unique()
    
    if(curr_drug %in% drug_classes[[1]]){
      curr_class <- names(drug_classes)[1]
    } else if(curr_drug %in% drug_classes[[2]]){
      curr_class <- names(drug_classes)[2]
    } else if(curr_drug %in% drug_classes[[3]]){
      curr_class <- names(drug_classes)[3]      
    } else{
      next
    }
    
    
    #Remove genesets that have no overlap with current data row names
    # temp <- lapply(genesets_to_use, FUN = function(x){ if(length(intersect(x, rownames(curr_data))) == 0){
    #   return(x)
    # }})
    
    #Score all samples for each geneset
    cancer_gene_set_score_info <- compute_bulk_normalized_gene_set_scores(gene_exp_mat = curr_data, gene_sets = genesets_to_use, num_controls = 5, num_bins = 10, q_thresh=0.95, gene_universe=NULL, use_median=F)
    
    #Subtract background geneset scores from foreground geneset scores
    normalized_gene_set_score_mat <- cancer_gene_set_score_info$fg - cancer_gene_set_score_info$bg
    
    scores_df <- cbind(scores_df, normalized_gene_set_score_mat)
    
    curr_drug_or_df <- data.frame(matrix(NA, nrow=0, ncol=6))
    
    #Calculate current drug OR for each geneset
    for(j in 1:length(genesets_to_use)){
      sorted_sample_scores <- sort(normalized_gene_set_score_mat[j,], decreasing = T)
      
      # x <- .3
      top_x_samples <- sorted_sample_scores[sorted_sample_scores > quantile(sorted_sample_scores, probs = (1-x))]
      bottom_x_samples <- sorted_sample_scores[sorted_sample_scores < quantile(sorted_sample_scores, probs = x)]
      
      names(top_x_samples) <- gsub("\\.","-",names(top_x_samples))
      names(bottom_x_samples) <- gsub("\\.","-",names(bottom_x_samples))
      
      if(curr_drug == "Selinexor"){
        names(top_x_samples) <- gsub("X","",names(top_x_samples))
        names(bottom_x_samples) <- gsub("X","",names(bottom_x_samples))
        
      }
      
      #Number of Non-responder in top x%
      top_NR <- enlight_response_data %>% 
        filter(Sample.ID %in% names(top_x_samples)) %>% 
        pull(Response) 
      
      top_NR <- sum(top_NR=="Non-responder")
      
      #Number of responder in top x%
      top_R <- enlight_response_data %>% 
        filter(Sample.ID %in% names(top_x_samples)) %>% 
        pull(Response)
      
      top_R <- sum(top_R=="Responder")
      
      
      #Number of Non-responder in bottom x%
      bottom_NR <- enlight_response_data %>% 
        filter(Sample.ID %in% names(bottom_x_samples)) %>% 
        pull(Response)
      
      bottom_NR <- sum(bottom_NR=="Non-responder")
      
      #Number of responder in bottom x%
      bottom_R <- enlight_response_data %>% 
        filter(Sample.ID %in% names(bottom_x_samples)) %>% 
        pull(Response)
      
      bottom_R <- sum(bottom_R=="Responder")
      
      
      # fisher_table <- matrix(c(1,2,3,4), nrow = 2,
      #                        dimnames =list(c("Top x%", "Bottom x%"),
      #                                       c("NR", "All")))
      # 
      # fisher_table <- matrix(c(top_NR,bottom_NR,(top_NR+top_R),(bottom_NR+bottom_R)), nrow = 2,
      #                        dimnames =
      #                          list(c("Top x%", "Bottom x%"),
      #                               c("NR", "All")))
      # 
      # #Run fisher test
      # fisher_table <- fisher_table+1
      # fisher_results <- fisher.test(fisher_table, alternative = "two.sided")
      # 
      # NR_OR <- fisher_results$estimate
      
      #Manual OR calculation
      top_NR_percent <- top_NR/(top_NR+top_R)
      top_R_percent <- top_R/(top_NR+top_R)
      
      bottom_NR_percent <- bottom_NR/(bottom_NR+bottom_R)
      bottom_R_percent <- bottom_R/(bottom_NR+bottom_R)
      
      
      NR_OR <- (top_NR_percent+.5)/(bottom_NR_percent+.5)
      
      R_OR <- (bottom_R_percent+.5)/(top_R_percent+.5)
      
      curr_geneset_name <- names(genesets_to_use)[j]
      
      
      curr_drug_or_df <- rbind(curr_drug_or_df,c(curr_geneset_name,NR_OR, R_OR, curr_drug, curr_cancer_type, curr_class))
    }
    
    colnames(curr_drug_or_df) <- c("geneset","NR_OR","R_OR","drug","cancer_type","class")
    
    curr_drug_or_df <- data.frame(curr_drug_or_df)
    
    curr_drug_or_df$NR_OR <- as.numeric(curr_drug_or_df$NR_OR)
    curr_drug_or_df$R_OR <- as.numeric(curr_drug_or_df$R_OR)
    curr_drug_or_df$geneset <- factor(curr_drug_or_df$geneset, levels = names(genesets_to_use))
    
    
    or_df <- rbind(or_df, curr_drug_or_df)
  }
  
  return(list(or_df,scores_df))
}
#################################################################################

enlight_response_data <- read.csv("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/enlight_data/drug_response_classifications_with_type.csv")

all_drugs <- unique(enlight_response_data$Dataset)

data_list <- list()
for(curr_drug in all_drugs){
  data_list[[curr_drug]] <- read.csv(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/enlight_data/", curr_drug,".csv"), row.names = 1)
}

data_list <- data_list[-13]
#################################################################################
curr_cell_line <- "K562"
signature_length <- 200

cluster_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_cluster_signatures_",signature_length,".rds"))

supercluster_signature <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/supercluster_signature.rds"))


raj_resistance_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_resistance_signatures.rds")

# A549_pkc_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/A549_PKC_signaling_cluster_signatures_200.rds")

raj_resistant_2017 <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_resistant_2017.rds")

raj_early_response_2017 <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_early_response_2017.rds")



#################################################################################

genesets_to_use <- active_cluster_signatures

# genesets_to_use <- append(genesets_to_use, raj_resistant_2017)
# genesets_to_use <- append(genesets_to_use, raj_early_response_2017)

response_result <- therapy_response_OR(data_list, genesets_to_use, .5)

df <- response_result[[1]]
original_df <- df


df <- original_df %>% 
  filter(original_df$class == "small_molecules")

ggplot(df, aes(x=drug, y=NR_OR, fill=cancer_type))+
  facet_wrap(~geneset)+
  geom_col()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_hline(yintercept=1)+
  ggtitle(paste0(curr_cell_line, " Cluster Signature Scores Predicting Non-Responder"))


ggplot(df, aes(x=drug, y=R_OR, fill=cancer_type))+
  facet_wrap(~geneset)+
  geom_col()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_hline(yintercept=1)

active <- list()
active <- append(active, list(c(4,9,13)))
active <- append(active, list(c(5,11)))
active <- append(active, list(c(5,8,13,16,17)))
names(active) <- c("A549","K562","MCF7")

active_clusters <- active[[curr_cell_line]]

active_clusters_names <- paste0(curr_cell_line,"_cluster",active_clusters,"_signature")

active_clusters_names <- paste0(active_clusters,"_active")

vec1 <- df %>% 
  filter(df$geneset %in% active_clusters_names & df$class == "small_molecules") %>% 
  pull(NR_OR)

vec2 <- df %>% 
  filter(!df$geneset %in% active_clusters_names & df$class == "small_molecules") %>% 
  pull(NR_OR)

boxplot_df <- data.frame("OR" = c(vec1,vec2),"group" = c(rep("Active", length(vec1)),rep("Inactive", length(vec2))))

wilcox_res <- wilcox.test(vec1,vec2)
pvalue <- wilcox_res$p.value
plot_title <- paste0(curr_cell_line, " Cluster Signatures (p value: ", sprintf("%.4f", pvalue), ")")

ggplot(boxplot_df)+
  geom_boxplot(aes(x=group,y=OR),fill=c("pink","lightblue"))+
  xlab("Cluster Resistance Group")+
  ylab("Non-reponder Odds Ratio")+
  ggtitle(plot_title)
#################################################################################

scores_df <- response_result[[2]]

sample_responses <- enlight_response_data %>% 
  select(Sample.ID,Response) %>% 
  rename("rowname"  = "Sample.ID")

temp <- data.frame(t(scores_df)) %>% 
  rownames_to_column() %>% 
  pivot_longer(!rowname, names_to = "signature",values_to = "score") %>% 
  merge(., sample_responses, by="rowname") 


ggplot(temp)+
  geom_boxplot(aes(x=Response,y=score))+
  facet_wrap(~signature)


summary(temp$score)

#################################################################################
emergent <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/trapnell_cluster_groups/emergent_clusters.rds")

emergent_clusters <- emergent[[curr_cell_line]]

emergent_clusters_names <- paste0(curr_cell_line,"_cluster",emergent_clusters,"_signature")
colnames(df)

df <- df %>% 
  select(geneset, NR_OR) %>% 
  mutate("emergent" = ifelse(geneset %in% emergent_clusters_names, geneset, "non-emergent")) %>% 
  group_by(emergent) %>% 
  mutate("emergent_median_OR" = median(NR_OR)) %>% 
  group_by(geneset) %>% 
  mutate("cluster_median_OR" = median(NR_OR))

p <- ggboxplot(df, x = "geneset", y = "NR_OR", fill="cluster_median_OR")+
  scale_fill_gradient2(low="blue",mid="white", high="red", midpoint = 1)+
  ggtitle(paste0(curr_cell_line," Small Molecule Non-responder Prediction"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#  Add p-value
p + stat_compare_means()




df$active <- ifelse(grepl("[0-9]a", df$geneset), "active","inactive")
p <- ggboxplot(df, x = "active", y = "NR_OR", fill="cluster_median_OR")+
  scale_fill_gradient2(low="blue",mid="white", high="red", midpoint = 1)+
  ggtitle(paste0(curr_cell_line," Small Molecule Non-responder Prediction"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#  Add p-value
p + stat_compare_means()






















cluster_median_ORs <- df %>%
  select(geneset, NR_OR) %>%
  group_by(geneset) %>%
  mutate("cluster_median_score" = median(NR_OR)) %>%
  arrange(geneset) %>%
  select(cluster_median_score) %>%
  unique()









p <- ggboxplot(df, x = "emergent", y = "NR_OR", fill="emergent_mean_OR")+
  scale_fill_gradientn(colours=c("blue","white","red"),na.value = "transparent",breaks=c(-1.5,0,1.5))+
  ggtitle(curr_cell_line)
#  Add p-value
p + stat_compare_means()




vec1 <- df %>% 
  filter(df$geneset %in% emergent_clusters_names & df$class == "small_molecules") %>% 
  pull(NR_OR)

vec2 <- df %>% 
  filter(!df$geneset %in% emergent_clusters_names & df$class == "small_molecules") %>% 
  pull(NR_OR)

boxplot_df <- data.frame("OR" = c(vec1,vec2),"group" = c(rep("emergent", length(vec1)),rep("Inemergent", length(vec2))))

wilcox_res <- wilcox.test(vec1,vec2)
pvalue <- wilcox_res$p.value
plot_title <- paste0(curr_cell_line, " Cluster Signatures (p value: ", sprintf("%.4f", pvalue), ")")

ggplot(boxplot_df)+
  geom_boxplot(aes(x=group,y=OR),fill=c("pink","lightblue"))+
  xlab("Cluster Group")+
  ylab("Non-reponder Odds Ratio")+
  ggtitle(plot_title)


#################################################################################
total_df <- matrix(NA,nrow=0,ncol=4)

for(i in 1:length(data_list)){
  df <- therapy_response_OR_by_drug(data_list[[i]], names(data_list)[i], cluster_signatures, .1)  
  
  
  total_df <- rbind(total_df,df)
}

ggplot(total_df, aes(x=drug,y=NR_OR, fill=NR_OR))+
  facet_wrap(~geneset)+
  geom_col()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))



active <- c(5,8,13,16,17)

active <- paste0(curr_cell_line,"_cluster",active,"_signature")

vec1 <- total_df %>% 
  filter(total_df$geneset %in% active) %>% 
  pull(NR_OR)

vec2 <- total_df %>% 
  filter(!total_df$geneset %in% active) %>% 
  pull(NR_OR)


boxplot(vec1,vec2)
wilcox.test(vec1,vec2)

