args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

source("final_scripts/drug_treatment_functions.R")
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"
# plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_figures/"

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

supercluster_signatures <- readRDS(paste0(dataDirectory,"genesets/supercluster_up_signatures.rds"))

geneset_to_use <- c(rac_signatures[1],supercluster_signatures)

geneset_title <- "RAC Signatures"

#Filter out genes with low/zero variance
data_to_use <- as.matrix(data)
# row_vars <- rowVars(data_to_use)
# data_to_use <- data_to_use[which(row_vars > 1e-10),]

# Score samples
gsvapar <- ssgseaParam(data_to_use, geneset_to_use)
ssgsea_res <- gsva(gsvapar)

df <- data.frame(ssgsea_res) %>%
  rownames_to_column("geneset") %>% 
  pivot_longer(!geneset, names_to = "sample",values_to = "score")

df <- merge(df,data_metadata,by = "sample")

# Plot boxplots
df$Recurrence <- df$response

facet_names <- c("A549 Global RAC Signature","Supercluster 1 Signature","Supercluster 2 Signature","Supercluster 3 Signature")
# facet_names <- c("A549 Global RAC Signature","Supercluster 1 Signature","Supercluster 2 Signature")

p <- ggboxplot(df %>% 
                 filter(treatment_stage=="pre"), x = "Recurrence", y = "score",
               fill = "Recurrence", palette = "jco", size=.2, outlier.shape = NA,
               facet.by = "geneset", short.panel.labs = T,
               add="", panel.labs = list(geneset=facet_names), 
               panel.labs.font = list(size=20),ncol=4)


# Use only p.format as label. Remove method name.
pre_plot <- p + 
  stat_compare_means(label = "p.format", method = "wilcox", size=2, label.x = 1.2, label.y = 1.3,
                     method.args = list(alternative = "less"))+
  xlab("Recurrence")+
  ylab("Score")+
  ylim(0,1.35)+
  theme(legend.position="none",
        strip.text.x.top = element_text(size=4),
        title = element_text(),
        legend.text = element_text(size=4),
        legend.title = element_text(size=6),
        legend.key.height = unit(2,"mm"),
        legend.key.width = unit(2,"mm"),
        axis.text.x = element_text(size=6),
        axis.title = element_text(size=8),
        axis.text.y = element_text(size=6),
        axis.line = element_line(linewidth=.2),
        axis.ticks = element_line(linewidth = .2))


pre_plot



jpeg(paste0(plotDirectory,"figure_3c.jpg"), width=120, height = 50, units = "mm", res = 1000)
print(pre_plot)
dev.off()


