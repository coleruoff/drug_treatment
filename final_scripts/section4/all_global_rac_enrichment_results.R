setwd("/data/ruoffcj/projects/drug_treatment/")
source("source/final_scripts/drug_treatment_functions.R")
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)
library(tidyverse)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"
plotDirectory <- "/data/ruoffcj/projects/drug_treatment/"

cell_lines <- c("A549","K562","MCF7")

# hallmark term2gene list
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, human_gene_symbol)

# ITH Meta-programs term2gene list
mp_t2g <- readRDS(paste0(dataDirectory, "genesets/ith_meta_programs_t2g.rds"))

#Remove 'specific' MPs
specifc_mps <- c("MP39 Metal-response","MP31 Alveolar","MP29 NPC/OPC","MP28 Oligo normal","MP27 Oligo Progenitor","MP38 Glutathione","MP41 Unassigned","MP35 Hemato-related-I","MP37 Hemato-related-II","MP32 Skin-pigmentation","MP36 IG","MP16 MES (glioma)","MP15 EMT IV")
mp_t2g <- mp_t2g %>% 
  filter(!gs_name %in% specifc_mps)
######################################################################################

hallmark_results <- list()
mp_results <- list()
go_results <- list()

curr_cell_line <- cell_lines[1]

all_ranks <- readRDS(paste0(dataDirectory, "genesets/global_rac_ranks.rds"))

for(curr_cell_line in cell_lines){
  
  cat(curr_cell_line, "\n")
  
  # curr_de <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/",curr_cell_line,"_cluster_all_markers_de.rds"))
  # curr_de <- readRDS(paste0(dataDirectory, "de_results/",curr_cell_line,"_cluster_de.rds"))
  
  
  
  curr_clusters <- 1:3
  
  hallmark_plots <- list()
  mp_plots <- list()
  go_plots <- list()
  
  for(i in curr_clusters){
    
    cat(i, "\n")
    
    ranks <- all_ranks[[curr_cell_line]][[i]]
    
    hallmark_gsea <- GSEA(ranks, TERM2GENE = m_t2g)
    mp_gsea <- GSEA(ranks, TERM2GENE = mp_t2g)
    go_gsea <- gseGO(geneList     = ranks,
                     OrgDb        = org.Hs.eg.db,
                     ont          = "BP",
                     minGSSize    = 100,
                     maxGSSize    = 500,
                     pvalueCutoff = 0.01,
                     verbose      = FALSE,
                     keyType = "SYMBOL")
    
    
    cluster_name <- paste0("Cluster ", i)
    hallmark_results[[curr_cell_line]][[cluster_name]] <- hallmark_gsea
    mp_results[[curr_cell_line]][[cluster_name]] <- mp_gsea
    go_results[[curr_cell_line]][[cluster_name]] <- go_gsea
    
  }
}



saveRDS(hallmark_results, paste0(dataDirectory, "enrichment_results/all_global_rac_hallmark_results.rds"))
saveRDS(mp_results, paste0(dataDirectory, "enrichment_results/all_global_rac_mp_results.rds"))
saveRDS(go_results, paste0(dataDirectory, "enrichment_results/all_global_rac_go_results.rds"))




