setwd("/data/ruoffcj/projects/drug_treatment/")
library(tidyverse)
library(decoupleR)
library(OmnipathR)

remotes::install_github('saezlab/decoupleR')

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

net <- get_collectri(organism='human', split_complexes=FALSE)

all_tfs <- unique(net$source)

regulons <- list()
for(curr_tf in all_tfs){
  
  curr_regulon <- net %>% 
    filter(source == curr_tf) %>% 
    pull(target)
  
  curr_regulon <- append(curr_regulon,curr_tf)
  
  regulons <- append(regulons, list(curr_regulon))
  
}



saveRDS(regulons, paste0(dataDirectory, "genesets/tf_regulons.rds"))


regulons <- readRDS(paste0(dataDirectory, "genesets/tf_regulons.rds"))
names(regulons) <- paste0(all_tfs,"_regulon")
