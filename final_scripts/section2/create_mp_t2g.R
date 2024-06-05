args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
setwd(args[1])
library(tidyverse)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

MPs <- readRDS(paste0(dataDirectory, "genesets/ITH_meta_programs.rds"))

mp_t2g <- matrix(NA, nrow=0,ncol=2)

for(i in 1:length(MPs)){

  curr_segment <- cbind(rep(names(MPs)[i], length(MPs[[i]])), MPs[[i]])
 
  mp_t2g <- rbind(mp_t2g, curr_segment) 
}

colnames(mp_t2g) <- c("gs_name","human_gene_symbol")
mp_t2g <- as_tibble(mp_t2g)

saveRDS(mp_t2g, paste0(dataDirectory,"genesets/ith_meta_programs_t2g.rds"))
