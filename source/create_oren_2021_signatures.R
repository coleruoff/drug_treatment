library(readxl)

oren_signatures_df <- read_xlsx("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/oren_2021_gene_signatures.xlsx")

oren_signatures_df <- oren_signatures_df[,!grepl("HALLMARK",colnames(oren_signatures_df))]

oren_signatures <- list()
for(i in 1:ncol(oren_signatures_df)){
  
  temp <- as.vector(oren_signatures_df[,i])
  
  oren_signatures[[i]] <- temp[[1]][!is.na(temp[[1]])]
  
}
names(oren_signatures) <- colnames(oren_signatures_df)

to_remove <- c("S.GENES.ITAY","G2M.GENES.ITAY","CELL_CYCLE")

oren_signatures <- oren_signatures[!names(oren_signatures) %in% to_remove]

saveRDS(oren_signatures, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/oren_2021_gene_signatures.rds")


