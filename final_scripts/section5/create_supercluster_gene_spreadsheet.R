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

write.xlsx(df, file=paste0(dataDirectory, "supplementary_tables/supplementary_table_5.xlsx"), sheetName="ortholog_genesets", row.names=FALSE)

#################################################################################
df <- list()
#############################
# Yeast

supercluster_signatures <- readRDS(paste0(dataDirectory, "genesets/rac_supercluster_signatures.rds"))

names1 <- gsub("type1_", "", names(supercluster_signatures))
names2 <- gsub("_", " ", names1)
names3 <- gsub("supercluster", "Supercluster ", names2)
names(supercluster_signatures) <- gsub("signature", "Signature", names3)

supercluster_signatures[[1]] <- intersect(supercluster_signatures[[1]], yeast_human_orthologs_up$yeast_human_orthologs_up)
supercluster_signatures[[2]] <- intersect(supercluster_signatures[[2]], yeast_human_orthologs_up$yeast_human_orthologs_up)

df[["yeast_supercluster1"]] <- supercluster_signatures[[1]]
df[["yeast_supercluster2"]] <- supercluster_signatures[[2]]

#############################
# E. coli

names1 <- gsub("type1_", "", names(supercluster_signatures))
names2 <- gsub("_", " ", names1)
names3 <- gsub("supercluster", "Supercluster ", names2)
names(supercluster_signatures) <- gsub("signature", "Signature", names3)

supercluster_signatures[[1]] <- intersect(supercluster_signatures[[1]], ecoli_human_orthologs_up$ecoli_human_orthologs_up)
supercluster_signatures[[2]] <- intersect(supercluster_signatures[[2]], ecoli_human_orthologs_up$ecoli_human_orthologs_up)

df[["ecoli_supercluster1"]] <- supercluster_signatures[[1]]
df[["ecoli_supercluster2"]] <- supercluster_signatures[[2]]


max_len <- max(lengths(df))

df <- lapply(df, sort)

df <- lapply(df, FUN = function(x) append(x,rep("",max_len-length(x))))

df <- data.frame(df)

write.xlsx(df, file=paste0(dataDirectory, "supplementary_tables/supplementary_table_5.xlsx"), sheetName="supercluster_ortholog_overlap", append=T, row.names=FALSE)
