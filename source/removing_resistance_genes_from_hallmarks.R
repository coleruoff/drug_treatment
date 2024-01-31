

hallmarks <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/hallmarks.rds")

resistance_signature <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_watermelon_resistance_signature.rds")

hallmarks_without_resistance <- list()
for(i in 1:length(hallmarks)){
  
  hallmarks_without_resistance <- append(hallmarks_without_resistance,list(hallmarks[[i]][!hallmarks[[i]] %in% resistance_signature[[1]]]))
  
}
names(hallmarks_without_resistance) <- names(hallmarks)
lengths(hallmarks)
lengths(hallmarks_without_resistance)

saveRDS(hallmarks_without_resistance,"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/hallmarks_without_resistance.rds")



MPs <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ITH_meta_programs.rds")


MPs_without_resistance <- list()
for(i in 1:length(MPs)){
  
  MPs_without_resistance <- append(MPs_without_resistance,list(MPs[[i]][!MPs[[i]] %in% resistance_signature[[1]]]))
  
}
names(MPs_without_resistance) <- names(MPs)
lengths(MPs)
lengths(MPs_without_resistance)

saveRDS(MPs_without_resistance,"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/MPs_without_resistance.rds")
