args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])
library(readxl)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GSVA")
library(GSVA)
library(tidyverse)
library(ggpubr)
library(matrixStats)
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"
# plotDirectory <- "/data/ruoffcj/projects/drug_treatment/"

################################################################################

# Read in data and fix rownames
data <- read_xlsx(paste0(dataDirectory,"raw_data/TKI_EGFR_rebiopsy_data/GSE165019_EarlyRebiopsy_genes.fpkm_table.xlsx"))

data_rownames <- make.names(data$`Gene ID`, unique = T)

data <- data[,-1]

rownames(data) <- data_rownames

colnames(data) <- gsub(pattern = " ", "_" , colnames(data))

# Set up metadata dataset

data_metadata <- data.frame(matrix(NA, nrow=ncol(data),ncol=2))
data_metadata[,1] <- colnames(data)
data_metadata[,2] <- rep(c("pre","post"),4)
data_metadata[,3] <- c(rep("Short", 8),rep("Long",8))
colnames(data_metadata) <- c("sample", "treatment_stage","response")
rownames(data_metadata) <- colnames(data)

# Read in gene signatures

rac_signatures <- readRDS(paste0(dataDirectory, "genesets/global_rac_signatures.rds"))

supercluster_signatures <- readRDS(paste0(dataDirectory,"genesets/rac_supercluster_signatures.rds"))

geneset_to_use <- c(rac_signatures[1],supercluster_signatures)

geneset_title <- "RAC Signatures"

#Filter out genes with low/zero variance
data_to_use <- as.matrix(data)
row_vars <- rowVars(data_to_use)
data_to_use <- data_to_use[which(row_vars > 1e-10),]

# Score samples
gsvapar <- ssgseaParam(data_to_use, geneset_to_use)
ssgsea_res <- gsva(gsvapar)

df <- data.frame(ssgsea_res) %>%
  rownames_to_column("geneset") %>% 
  pivot_longer(!geneset, names_to = "sample",values_to = "score")


df <- merge(df,data_metadata,by = "sample")


# Plot boxplots

df$Recurrence <- df$response

facet_names <- c("A549 Global RAC Signature","Supercluster 1 Signature","Supercluster 2 Signature")

p <- ggboxplot(df %>% 
                 filter(treatment_stage=="pre"), x = "Recurrence", y = "score",
               color = "Recurrence", palette = "jco",
               facet.by = "geneset", short.panel.labs = T,
               add="", panel.labs = list(geneset=facet_names), 
               panel.labs.font = list(size=20),ncol=4)



# Use only p.format as label. Remove method name.
pre_plot <- p + stat_compare_means(label = "p.format", method = "wilcox", size=3,
                                   label.x = 1.2, 
                                   label.y = 2.15)+
  xlab("")+
  ylab("Score")+
  ylim(1,2.2)+
  theme(legend.position="right",
        strip.text.x.top = element_text(size=6),
        title = element_text(size=10),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10),
        legend.key.height = unit(5,"mm"),
        legend.key.width = unit(5,"mm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=10))


pre_plot

tiff(paste0(plotDirectory,"figure_4c.tiff"), width=140, height = 80, units = "mm", res = 1000)

print(pre_plot)

dev.off()


# png(paste0(plotDirectory,"figure_4c.png"),
#     width=16, height = 6, units="in",res=300)
# 
# 
# print(pre_plot)
# 
# dev.off()


