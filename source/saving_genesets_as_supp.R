
length(ecoli_human_orthologs_up$ecoli_human_orthologs_up) <- length(yeast_human_orthologs_up$yeast_human_orthologs_up)
upregulated_orthologs <- list(ecoli_human_orthologs_up$ecoli_human_orthologs_up,yeast_human_orthologs_up$yeast_human_orthologs_up)
names(upregulated_orthologs) <- c("ecoli_orthologs_up","yeast_orthologs_up")


write.csv(upregulated_orthologs,"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/upregulated_orthologs.csv", 
          row.names = F)



yeast_up <- upregulated_genes


write.csv(yeast_up,"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/yeast_resistance_geneset.csv", 
          row.names = F)




ecoli_genesets <- append(ecoli_AMR_genesets_up, list("consensus" = ecoli_AMR_consensus_geneset_up))



max_len <- max(lengths(ecoli_genesets))

ecoli_genesets <- lapply(lapply(ecoli_genesets, unlist), "length<-", max_len)


write.csv(ecoli_genesets,"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/ecoli_resistance_genesets.csv", 
          row.names = F)



