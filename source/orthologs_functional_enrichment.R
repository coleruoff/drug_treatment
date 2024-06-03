library(clusterProfiler)
library(tidyverse)
library(ggpubr)

# genesets_to_use <- ecoli_up_genesets
# genesets_title <- "E Coli Downregulated Genes"

GO_dotplots_for_genesets <- function(genesets_to_use, genesets_title){
  
  dotplots <- list()
  
  for(i in 1:length(genesets_to_use)){
    cat(i,"/",length(genesets_to_use),"\n")
    
    go_enrich <- enrichGO(gene = genesets_to_use[[i]],
                          OrgDb = "org.Hs.eg.db",
                          keyType = 'SYMBOL',
                          readable = T,
                          ont = "ALL",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.10)
    
    
    if(nrow(go_enrich)> 0){
      
      p <- dotplot(go_enrich)+ggtitle(names(genesets_to_use)[i])
      
      dotplots <- append(dotplots, list(p))
      
    }
  }
  
  
  plt <- egg::ggarrange(plots = dotplots, ncol = 4)
  
  main_title <- paste0("Functional Enrichment of ", genesets_title)
  plt <- annotate_figure(plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
  
  return(plt)
}

################################################################################

# E coli AMR orthologs
ecoli_AMR_orthologs <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ecoli_AMR_genesets_orthologs.rds")

ecoli_up_genesets <- ecoli_AMR_orthologs[grepl("up", names(ecoli_AMR_orthologs))]
plot(GO_dotplots_for_genesets(ecoli_up_genesets, "E Coli upregulated Orthologs"))

ecoli_down_genesets <- ecoli_AMR_orthologs[grepl("down", names(ecoli_AMR_orthologs))]
plot(GO_dotplots_for_genesets(ecoli_down_genesets, "E Coli Downregulated Orthologs"))


# E coli Stress response orthologs
ecoli_stress_orthologs <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ecoli_stress_orthologs.rds")

ecoli_stress_up_genesets <- ecoli_stress_orthologs[grepl("up", names(ecoli_stress_orthologs))]
plot(GO_dotplots_for_genesets(ecoli_stress_up_genesets, "E Coli Stress Up Orthologs"))

ecoli_stress_down_genesets <- ecoli_stress_orthologs[grepl("down", names(ecoli_stress_orthologs))]
plot(GO_dotplots_for_genesets(ecoli_stress_down_genesets, "E Coli Stress down Orthologs"))



# E coli AMR consensus orthologs
ecoli_AMR_consensus_orthologs <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ecoli_AMR_consensus_orthologs.rds")
plot(GO_dotplots_for_genesets(ecoli_AMR_consensus_orthologs, "E Coli AMR Consensus Orthologs"))

# E coli Stress response consensus orthologs
ecoli_stress_consensus_orthologs <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ecoli_stress_consensus_orthologs.rds")
plot(GO_dotplots_for_genesets(ecoli_stress_consensus_orthologs, "E Coli Stress Consensus Orthologs"))



# Yeast stress orthologs
yeast_stress_orthologs <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/yeast_upregulated_orthologs.rds")
plot(GO_dotplots_for_genesets(list(yeast_stress_orthologs), "Yeast Upregulated Orthologs"))

