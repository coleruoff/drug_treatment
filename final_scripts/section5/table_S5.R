args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
setwd(args[1])
library(xlsx)
#################################################################################
df <- list()

yeast_human_orthologs_up <- readRDS(paste0(dataDirectory, "genesets/yeast_human_orthologs_up.rds"))

ecoli_human_orthologs_up <- readRDS(paste0(dataDirectory, "genesets/ecoli_human_orthologs_up.rds"))


df[["yeast_upregulated_orthologs"]] <- yeast_human_orthologs_up$yeast_human_orthologs_up
df[["ecoli_upregulated_orthologs"]] <- ecoli_human_orthologs_up$ecoli_human_orthologs_up


max_len <- max(lengths(df))

df <- lapply(df, sort)

df <- lapply(df, FUN = function(x) append(x,rep("",max_len-length(x))))

df <- data.frame(df)

write.xlsx(df, file=paste0(dataDirectory, "supplementary_tables/table_S5.xlsx"), sheetName="ortholog_genesets", row.names=FALSE)

#################################################################################
df <- list()
#############################
# Yeast

df[["yeast_supercluster2"]] <- readRDS(paste0(dataDirectory, "genesets/yeast_sc2_overlap.rds"))

#############################
# E. coli

df[["ecoli_supercluster2"]] <- readRDS(paste0(dataDirectory, "genesets/ecoli_sc2_overlap.rds"))


max_len <- max(lengths(df))

df <- lapply(df, sort)

df <- lapply(df, FUN = function(x) append(x,rep("",max_len-length(x))))

df <- data.frame(df)

write.xlsx(df, file=paste0(dataDirectory, "supplementary_tables/table_S5.xlsx"), sheetName="supercluster_ortholog_overlap", append=T, row.names=FALSE)
