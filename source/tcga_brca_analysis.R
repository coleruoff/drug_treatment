if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
library(tidyverse)
library(SummarizedExperiment)
library(DESeq2)
library(survminer)
library(survival)


################################################################################
# Get clincal data from GDC

clincal_brca <- GDCquery_clinic('TCGA-BRCA')

any(colnames(clincal_brca) %in% c("vital_status","days_to_last_follow_up","days_to_death"))

filtered_clincal_brca <- clincal_brca[,which(colnames(clincal_brca) %in% c("vital_status","days_to_last_follow_up","days_to_death"))]
table(filtered_clincal_brca$vital_status)

clincal_brca$deceased <- ifelse(clincal_brca$vital_status == "Alive", F, T)

clincal_brca$OS <- ifelse(clincal_brca$vital_status == "Alive", 
                          clincal_brca$days_to_last_follow_up,
                          clincal_brca$days_to_death)

################################################################################
# Get gene expression data

# Get list of projects
gdcprojects <- getGDCprojects()

tcga_projects <- gdcprojects[grepl("TCGA", gdcprojects$id),]$id

all_samples <- c()
for(i in tcga_projects){
  query_tcga <- GDCquery(project = i,
                         data.category = 'Transcriptome Profiling',
                         data.type = 'Gene Expression Quantification',
                         experimental.strategy = "RNA-Seq",
                         workflow.type = 'STAR - Counts',
                         access = 'open')
  
  output_query_tcga <- getResults(query_tcga)
  
  tumor <- output_query_tcga[output_query_tcga$sample_type == "Primary Tumor", "cases"][1:30]
  
  all_samples <- append(all_samples, tumor)
}

query_tcga <- GDCquery(project = tcga_projects,
                       data.category = 'Transcriptome Profiling',
                       data.type = 'Gene Expression Quantification',
                       experimental.strategy = "RNA-Seq",
                       workflow.type = 'STAR - Counts',
                       access = 'open',
                       barcode = all_samples)

# Get project summary
getProjectSummary("TCGA-BRCA")

# Build a query
query_tcga <- GDCquery(project = 'TCGA-BRCA',
         data.category = 'Transcriptome Profiling',
         data.type = 'Gene Expression Quantification',
         experimental.strategy = "RNA-Seq",
         workflow.type = 'STAR - Counts',
         access = 'open')
         # barcode = c('TCGA-E9-A1RH-01A-21R-A169-07','TCGA-C8-A26W-01A-11R-A16F-07','TCGA-E9-A1RH-11A-34R-A169-07'))

output_query_tcga <- getResults(query_tcga)

tumor <- output_query_tcga[output_query_tcga$sample_type == "Primary Tumor", "cases"][1:30]


query_tcga <- GDCquery(project = 'TCGA-BRCA',
                       data.category = 'Transcriptome Profiling',
                       data.type = 'Gene Expression Quantification',
                       experimental.strategy = "RNA-Seq",
                       workflow.type = 'STAR - Counts',
                       access = 'open',
                       barcode = tumor)


# Download data - GDCdownload
GDCdownload(query_tcga)


# Prepare data
tcga_brca_data <- GDCprepare(query_tcga, summarizedExperiment = T)

brca_matrix <- assay(tcga_brca_data, 'unstranded')

# extract gene and sample metadata from summarizedExperiment object
gene_metadata <- as.data.frame(rowData(tcga_brca_data))
coldata <- as.data.frame(colData(tcga_brca_data))

# Normalize data with DESeq
dds <- DESeqDataSetFromMatrix(countData = brca_matrix,
                              colData = coldata,
                              design = ~ 1)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# variance stabalizing transformation
vsd <- vst(dds, blind = F)
brca_matrix_vst <- assay(vsd)


rownames(brca_matrix_vst) 

temp <- gene_metadata %>% 
  filter(gene_id %in% rownames(brca_matrix_vst)) %>% 
  select(gene_id,gene_name)


if(all(temp$gene_id == rownames(brca_matrix_vst))){
  rownames(brca_matrix_vst) <- temp$gene_name
}



ssgsea_result <- ssgsea(brca_matrix_vst,supercluster_signature)

ssgsea_result_norm <- t(scale(t(ssgsea_result)))


scores_mat <- matrix(NA, nrow=nrow(ssgsea_result_norm),ncol=ncol(ssgsea_result_norm))
colnames(scores_mat) <- colnames(ssgsea_result_norm)
rownames(scores_mat) <- rownames(ssgsea_result_norm)

for(i in 1:nrow(ssgsea_result_norm)){
  scores_mat[i,] <- ifelse(ssgsea_result_norm[i,] < quantile(ssgsea_result_norm[i,],prob=.25), 'low', ifelse(ssgsea_result_norm[i,] > quantile(ssgsea_result_norm[i,],prob=.75), 'high', "med"))
  
  # scores_mat[i,] <- ifelse(ssgsea_result_norm[i,] > quantile(ssgsea_result_norm[i,],prob=.5), 'high', 'low')
}

scores_mat <- t(scores_mat) %>% 
  as.data.frame() %>% 
  rownames_to_column("submitter_id")

scores_mat$submitter_id <- gsub("-01.*", "", scores_mat$submitter_id)

################################################################################

final_df <- merge(scores_mat,clincal_brca,by.x='submitter_id')

colnames(scores_mat)
# Fitting survival curve
fit <- survfit(Surv(OS, deceased) ~ supercluster_signature, data = final_df)


# fit

ggsurvplot(fit,
           data = final_df,
           pval = T,
           risk.table = F)
