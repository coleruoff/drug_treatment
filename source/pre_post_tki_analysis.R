library(readxl)
library(GSVA)
library(tidyverse)
library(ggpubr)

# Read in data and fix rownames
data <- read_xlsx("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/TKI_EGFR_rebiopsy_data/GSE165019_EarlyRebiopsy_genes.fpkm_table.xlsx")

data_rownames <- make.names(data$`Gene ID`, unique = T)

data <- data[,-1]

rownames(data) <- data_rownames

colnames(data) <- gsub(pattern = " ", "_" , colnames(data))

# Set up metadata dataset

data_metadata <- data.frame(matrix(NA, nrow=ncol(data),ncol=2))
data_metadata[,1] <- colnames(data)
data_metadata[,2] <- rep(c("pre","post"),4)
data_metadata[,3] <- c(rep("short", 8),rep("long",8))
colnames(data_metadata) <- c("sample", "treatment_stage","response")
rownames(data_metadata) <- colnames(data)

# Read in gene signatures

rac_type_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/global_rac_type_signatures.rds")
rac_type_signatures <- rac_type_signatures[grepl("type1",names(rac_type_signatures))]

rac_supercluster_consensus_signature <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_supercluster_consensus_signature.rds")
type1_supercluster_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/type1_supercluster_signatures.rds")
type1_supercluster_signatures <- type1_supercluster_signatures[1:2]

geneset_to_use <- c(rac_type_signatures[1],type1_supercluster_signatures)

geneset_title <- "Type 1 Signatures"

# Score samples

ssgsea_res <- gsva(as.matrix(data), geneset_to_use, method="ssgsea")

df <- data.frame(ssgsea_res) %>%
  rownames_to_column("geneset") %>% 
  pivot_longer(!geneset, names_to = "sample",values_to = "score")


df <- merge(df,data_metadata,by = "sample")


# Plot boxplots

p <- ggboxplot(df %>% 
                 filter(treatment_stage=="pre"), x = "response", y = "score",
               color = "response", palette = "jco",
               facet.by = "geneset", short.panel.labs = FALSE,
               add="")


# Use only p.format as label. Remove method name.
p + stat_compare_means(label = "p.format", method = "wilcox")+
  ggtitle(paste0(geneset_title," ssGSEA Scores (Pre-treatment)"))+
  xlab("")+
  ylab("Score")+
  theme(legend.position="right",
        strip.text = element_text(size=20),
        title = element_text(size=30))


p <- ggboxplot(df %>% 
                 filter(treatment_stage=="post"), x = "response", y = "score",
               color = "response", palette = "jco",
               facet.by = "geneset", short.panel.labs = FALSE,
               add="")


# Use only p.format as label. Remove method name.
p + stat_compare_means(label = "p.format", method = "wilcox")+
  ggtitle(paste0(geneset_title," ssGSEA Scores (Post-treatment)"))+
  xlab("")+
  ylab("Score")+
  theme(legend.position="right",
        strip.text = element_text(size=20),
        title = element_text(size=30))




p <- ggboxplot(df, x = "treatment_stage", y = "score",
               color = "response", palette = "jco",
               facet.by = "geneset", short.panel.labs = FALSE,
               add="")


# Use only p.format as label. Remove method name.
p + stat_compare_means(label = "p.format", method = "wilcox")+
  ggtitle(paste0(geneset_title," ssGSEA Scores (Pre-treatment)"))+
  xlab("")+
  ylab("Score")+
  theme(legend.position="right",
        strip.text = element_text(size=20),
        title = element_text(size=30))


