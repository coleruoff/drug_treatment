

MPs <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ITH_meta_programs.rds")


MPs
mp_t2g <- matrix(NA, nrow=0,ncol=2)

for(i in 1:length(MPs)){

  curr_segment <- cbind(rep(names(MPs)[i], length(MPs[[i]])), MPs[[i]])
 
  mp_t2g <- rbind(mp_t2g, curr_segment) 
}

colnames(mp_t2g) <- c("gs_name","human_gene_symbol")
mp_t2g <- as.tibble(mp_t2g)

saveRDS(mp_t2g, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ith_meta_programs_t2g.rds")
