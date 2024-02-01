library(DESeq2)
library(tidyverse)
library(ggpubr)
library(GSVA)
# install.packages('corto')

#################################################################################
rac_type_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/global_rac_type_signatures.rds")
rac_type_signatures <- rac_type_signatures[grepl("type1",names(rac_type_signatures))]

rac_supercluster_consensus_signature <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_supercluster_consensus_signature.rds")
type1_supercluster_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/type1_supercluster_signatures.rds")
type1_supercluster_signatures <- type1_supercluster_signatures[1:2]

all_signatures <- type1_supercluster_signatures
geneset_title <- "RAC Type 1 Supercluster Signatures"
#################################################################################
drug_classes <- list()
drug_classes <- append(drug_classes, list(c("MGH_Alpelisib","Sorafenib","Sorafenib_2","MK2206","Tipifarnib_1","Tipifarnib_2","Selinexor","Vismodegib")))
drug_classes <- append(drug_classes, list(c("Anti-PD1","Anti-PD1 +- Anti-CTLA4", "Anti-PD1_2", "Anti-PD1_3","Anti-PD1_4")))
drug_classes <- append(drug_classes, list(c("Bevacizumab","Bevacizumab_2","Bevacizumab_3","Bevacizumab_4","Trastuzumab","Trastuzumab_2","Trastuzumab_3","Trastuzumab_4","Trastuzumab_5","Cetuximab","Rituximab")))
names(drug_classes) <- c("small_molecules","immunotherapy","mAb")

enlight_response_data <- read.csv("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/enlight_data/drug_response_classifications_with_type.csv")

rna_seq_datasets <- c("Trastuzumab_5","Bevacizumab_3", "Selinexor", "Vismodegib", "MGH_Ribociclib", "MGH_Alpelisib", "BRAFi_1", "Anti-PD1", "Anti-PD1_2", "Anti-PD1_5")

# enlight_response_data <- enlight_response_data %>%
#   filter(Dataset %in% rna_seq_datasets)#!Dataset %in% drug_classes$immunotherapy)

# enlight_response_data <- enlight_response_data %>%
#   filter(Dataset %in% drug_classes$small_molecules)

all_drugs <- unique(enlight_response_data$Dataset)

if("Anti-PD1 +- Anti-CTLA4" %in% all_drugs){
  all_drugs <- all_drugs[-which(all_drugs == "Anti-PD1 +- Anti-CTLA4")]  
}

all_sample_scores <- matrix(NA,nrow=0,ncol=(length(all_signatures)+1))

# curr_drug <- all_drugs[1]

for(curr_drug in all_drugs){
  cat(curr_drug, "\n")
  
  curr_data <- read.csv(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/enlight_data/", curr_drug,".csv"), row.names = 1)
  
  #Remove rows that are NAs
  curr_data <- curr_data[!is.na(curr_data[,1]),]
  
  #Normalize data if data are integers
  if(curr_data[which(curr_data != 0)[1],1]%%1 == 0){
    
    meta <- enlight_response_data %>% 
      filter(Dataset == curr_drug) %>% 
      select(Sample.ID,Response) %>% 
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
  # ssgsea_res <- ssgsea(as.matrix(curr_data), all_signatures, scale = TRUE)
  # ssgsea_res <- t(scale(t(ssgsea_res)))
  ssgsea_res <- t(scale(t(ssgsea_res)))
  
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

names1 <- sapply(as.character(df$geneset), FUN=function(x) gsub("_", " ",x))
names2 <- gsub("type1 ","",names1)
names3 <- gsub("supercluster","supercluster ",names2)
df$geneset <- str_to_title(names3)

df$drug <- sapply(df$Dataset, function(x) strsplit(x, "_")[[1]][1])

p <- ggboxplot(df, x = "drug", y = "score",
               color = "Response", palette = "jco",
               facet.by = c("cancer_type", "geneset"), short.panel.labs = FALSE,
               add="")


# Use only p.format as label. Remove method name.
p <- p + stat_compare_means(aes(group=Response),label = "p.format")+
  ggtitle(paste0(geneset_title," ssGSEA Scores"))+
  xlab("Drug")+
  ylab("Score")+
  coord_flip()+
  theme(legend.position="right",
        legend.title = element_blank(),
        strip.text = element_text(size=20),
        title = element_text(size=30),
        axis.text = element_text(size=18),
        legend.text = element_text(size=20),
        legend.key.size = unit(2, 'cm'))



p <- facet(p, facet.by = "geneset",
      short.panel.labs = T)

png("/data/ruoffcj/projects/drug_treatment/final_figures/figure_4e.png",
    width=22, height = 20, units="in",res=300)


print(p)

dev.off()


