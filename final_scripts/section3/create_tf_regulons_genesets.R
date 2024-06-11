args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])
library(tidyverse)
library(decoupleR)
library(OmnipathR)

remotes::install_github('saezlab/decoupleR')

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

#################################################################################

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

names(regulons) <- paste0(all_tfs,"_regulon")

saveRDS(regulons, paste0(dataDirectory, "genesets/tf_regulons.rds"))



