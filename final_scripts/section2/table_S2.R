args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

source("final_scripts/drug_treatment_functions.R")

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"



if (file.exists(paste0(dataDirectory, "supplementary_tables/table_S2.xlsx"))) {
  #Delete file if it exists
  file.remove(paste0(dataDirectory, "supplementary_tables/table_S2.xlsx"))
}

if(!file.exists(paste0(dataDirectory,"supplementary_tables/"))){
  dir.create(paste0(dataDirectory,"supplementary_tables/")) 
}

################################################################################
crispr_ko_data <- read_xlsx(paste0(dataDirectory,"gottesman_crispr_data/CRISPR_KO_summary_negative.xlsx"), sheet = 5)

crispr_ko_data[,5]

crispr_ko_ranks <- crispr_ko_data %>% 
  arrange(`neg|rank`) %>% 
  dplyr::select(id,`neg|rank`) %>% 
  mutate(ranks = 1:nrow(.))  %>% 
  select(id,ranks)

colnames(crispr_ko_ranks) <- c("gene","rank")

write.xlsx(as.data.frame(crispr_ko_ranks), file=paste0(dataDirectory, "supplementary_tables/table_S2.xlsx"), sheetName="CRISPR_gene_ranks", append=F, row.names=FALSE)

