


kegg <- gmtPathways("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt")

kegg_t2g <- matrix(NA, nrow=0,ncol=2)

for(i in 1:length(kegg)){
  
  curr_segment <- cbind(rep(names(kegg)[i], length(kegg[[i]])), kegg[[i]])
  
  kegg_t2g <- rbind(kegg_t2g, curr_segment) 
}

colnames(kegg_t2g) <- c("gs_name","human_gene_symbol")
kegg_t2g <- as.tibble(kegg_t2g)

saveRDS(kegg_t2g, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/kegg_t2g.rds")



kegg_tgf <- list("kegg_tgf_beta"= kegg$KEGG_TGF_BETA_SIGNALING_PATHWAY)

saveRDS(kegg_tgf,"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/kegg_tgf_beta.rds")
