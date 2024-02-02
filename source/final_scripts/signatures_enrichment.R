library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)
source("source/cole_functions.R")

genesets_characterization <- function(genesets_to_use, universe_to_use, num_pathways=10){
  m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::select(gs_name, human_gene_symbol)
  # hallmark_uni <- unique(c(unlist(all_signatures), m_t2g$human_gene_symbol))
  
  mp_t2g <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ith_meta_programs_t2g.rds")
  
  # mp_uni <- unique(c(unlist(all_signatures), mp_t2g$human_gene_symbol))
  
  hallmarks_dotplots <- list()
  mps_dotplots <- list()
  go_dotplots <- list()
  go_results <- list()
  
  for(i in 1:length(genesets_to_use)){

    curr_cluster <- names(genesets_to_use)[i]
    cat(curr_cluster, "\n")
    
    curr_geneset <- genesets_to_use[[i]]

    
    # Hallmarks Enrichment
    hallmark_enrichment_results <- enricher(curr_geneset, TERM2GENE=m_t2g,
                                            universe = universe_to_use)
    
    # MPs Enrichment
    mp_enrichment_results <- enricher(curr_geneset, TERM2GENE=mp_t2g,
                                      universe = universe_to_use)
    
    # GO Functional Enrichment
    ego <- enrichGO(gene          = curr_geneset,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = .01,
                    qvalueCutoff  = .05,
                    readable      = TRUE,
                    keyType = "SYMBOL",
                    universe = universe_to_use)
    
    go_results <- append(go_results,list(ego))
    
    if(nrow(hallmark_enrichment_results) > 0){
      p <- dotplot(hallmark_enrichment_results,
                   showCategory = num_pathways) + ggtitle(paste0("", curr_cluster))
      hallmarks_dotplots <- append(hallmarks_dotplots, list(p))
    }
    if(!is.null(mp_enrichment_results) && nrow(mp_enrichment_results) > 0){
      p <- dotplot(mp_enrichment_results,
                   showCategory = num_pathways)+ggtitle(paste0("", curr_cluster))
      mps_dotplots <- append(mps_dotplots, list(p))
    }
    if(nrow(ego) > 0){
      p <- dotplot(ego,
                   showCategory = num_pathways) + ggtitle(paste0("", curr_cluster))
      go_dotplots <- append(go_dotplots, list(p))
    }
  }
  
  
  ret_list <- list(hallmarks_dotplots, mps_dotplots, go_dotplots,go_results)
  names(ret_list) <- c("hallmarks_dotplots", "mps_dotplots", "go_dotplots","go_results")
 
  return(ret_list) 
}

cell_lines <- c("A549","K562","MCF7")

# cell_line_universes <- list()
# for(curr_cell_line in cell_lines){
#   
#   cat(curr_cell_line, "\n")
#   
#   #find genes expressed in > 1% of cells
#   data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/",curr_cell_line,"_processed_filtered.rds"))
#   gene_universe <- names(apply(data@assays$RNA@data, 1, FUN = function(x) sum(x>0) >.01*ncol(data)))
#   
#   
#   cell_line_universes <- append(cell_line_universes, list(gene_universe))
# }
# 
# names(cell_line_universes) <- cell_lines
# 
# gene_universe <- intersect(cell_line_universes[[1]], cell_line_universes[[2]])
# gene_universe <- intersect(gene_universe, cell_line_universes[[3]])
# 
# saveRDS(gene_universe, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/cell_line_gene_universe_intersection.rds")
# saveRDS(cell_line_universes, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/cell_line_universes.rds")

gene_universe_intersection <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/cell_line_gene_universe_intersection.rds")
cell_line_universes <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/cell_line_universes.rds")

################################################################################
# Plotting for RAC vs rest cells signature

rac_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_signatures.rds")

curr_cell_line <- cell_lines[1]

for(curr_cell_line in cell_lines){
  
  curr_rac_signatures <- rac_signatures[grepl(curr_cell_line, names(rac_signatures))]
  
  names1 <- gsub("_", " ", names(curr_rac_signatures))
  names2 <- gsub("cluster", "Cluster ", names1)
  names(curr_rac_signatures) <- gsub("signature", "Signature ", names2)
  
  all_dotplots <- genesets_characterization(curr_rac_signatures, universe_to_use = cell_line_universes[[curr_cell_line]])
  
  png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/rac_enrichment/",curr_cell_line,"_rac_hallmarks.png"), width = 1000,height = 1000)
  hallmarks_plt <- egg::ggarrange(plots = all_dotplots$hallmarks_dotplots, ncol = 2)
  main_title <- paste0("Cancer Hallmarks Enrichment of ", curr_cell_line ," RACs")
  hallmarks_plt <- annotate_figure(hallmarks_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
  plot(hallmarks_plt)
  dev.off()
  
  png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/rac_enrichment/",curr_cell_line,"_rac_mps.png"), width = 1000,height = 1000)
  mps_plt <- egg::ggarrange(plots = all_dotplots$mps_dotplots, ncol = 2)
  main_title <- paste0("ITH Meta-programs Enrichment of ", curr_cell_line ," RACs")
  mps_plt <- annotate_figure(mps_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
  plot(mps_plt)
  dev.off()
  
  png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/rac_enrichment/",curr_cell_line,"_rac_functional_enrichment.png"), width = 1000,height = 1000)
  go_plt <- egg::ggarrange(plots = all_dotplots$go_dotplots, ncol = 2)
  main_title <- paste0("GO Functional Enrichment of ", curr_cell_line ," RACs")
  go_plt <- annotate_figure(go_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
  plot(go_plt)
  dev.off()
}

################################################################################
# Plotting for RAC type 1 vs rest/inactive cells signature

global_rac_type_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/global_rac_type_signatures.rds")

global_rac_type1_signatures <- global_rac_type_signatures[grepl("type1",names(global_rac_type_signatures))]

names(global_rac_type1_signatures) <- cell_lines


all_dotplots <- list()
for(curr_cell_line in cell_lines){
  
  curr_signature <- global_rac_type1_signatures[curr_cell_line]
  curr_dotplots <- genesets_characterization(curr_signature, universe_to_use = cell_line_universes[[curr_cell_line]])  
  
  all_dotplots[["hallmarks_dotplots"]] <- append(all_dotplots[["hallmarks_dotplots"]],curr_dotplots[["hallmarks_dotplots"]])
  all_dotplots[["mps_dotplots"]] <- append(all_dotplots[["mps_dotplots"]],curr_dotplots[["mps_dotplots"]])
  all_dotplots[["go_dotplots"]] <- append(all_dotplots[["go_dotplots"]],curr_dotplots[["go_dotplots"]])
  all_dotplots[["go_results"]] <- append(all_dotplots[["go_results"]],curr_dotplots[["go_results"]])
  
}



A549_go <- all_dotplots$go_results[[1]]@result

strsplit(A549_go %>% 
  filter(Description == "cytoplasmic translation") %>% 
  pull(geneID), "/")[[1]]

K562_go <- all_dotplots$go_results[[2]]@result

strsplit(K562_go %>% 
           filter(Description == "cytoplasmic translation") %>% 
           pull(geneID), "/")[[1]]

MCF7_go <- all_dotplots$go_results[[3]]@result

strsplit(MCF7_go %>% 
           filter(Description == "cytoplasmic translation") %>% 
           pull(geneID), "/")[[1]]



png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/type1_vs_inactive_enrichment/type1_vs_inactive_hallmarks.png"), width = 1000,height = 1000)
hallmarks_plt <- egg::ggarrange(plots = all_dotplots$hallmarks_dotplots, ncol = 2)
main_title <- paste0("Cancer Hallmarks Enrichment of RAC Type 1 Cells")
hallmarks_plt <- annotate_figure(hallmarks_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(hallmarks_plt)
dev.off()

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/type1_vs_inactive_enrichment/type1_vs_inactive_mps.png"), width = 1000,height = 1000)
mps_plt <- egg::ggarrange(plots = all_dotplots$mps_dotplots, ncol = 2)
main_title <- paste0("ITH Meta-programs Enrichment of RAC Type 1 Cells")
mps_plt <- annotate_figure(mps_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(mps_plt)
dev.off()

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/type1_vs_inactive_enrichment/type1_vs_inactive_functional_enrichment.png"), width = 1000,height = 1000)
go_plt <- egg::ggarrange(plots = all_dotplots$go_dotplots, ncol = 2)
main_title <- paste0("GO Functional Enrichment of RAC Type 1 Cells")
go_plt <- annotate_figure(go_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(go_plt)
dev.off()

################################################################################
# Plotting for individual RAC type 1 vs rest/inactive cells signature

rac_type1_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_type1_signatures.rds")

all_dotplots <- genesets_characterization(rac_type1_signatures, universe_to_use = gene_universe_intersection)

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/type1_vs_inactive_enrichment/individual_type1_vs_inactive_hallmarks.png"), width = 2000,height = 2000)
hallmarks_plt <- egg::ggarrange(plots = all_dotplots$hallmarks_dotplots, ncol = 5)
main_title <- paste0("Cancer Hallmarks Enrichment of RAC Type 1 Clusters")
hallmarks_plt <- annotate_figure(hallmarks_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(hallmarks_plt)
dev.off()

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/type1_vs_inactive_enrichment/individual_type1_vs_inactive_mps.png"), width = 2000,height = 2000)
mps_plt <- egg::ggarrange(plots = all_dotplots$mps_dotplots, ncol = 5)
main_title <- paste0("ITH Meta-programs Enrichment of RAC Type 1 Clusters")
mps_plt <- annotate_figure(mps_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(mps_plt)
dev.off()

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/type1_vs_inactive_enrichment/individual_type1_vs_inactive_functional_enrichment.png"), width = 2000,height = 2000)
go_plt <- egg::ggarrange(plots = all_dotplots$go_dotplots, ncol = 5)
main_title <- paste0("GO Functional Enrichment of RAC Type 1 Clusters")
go_plt <- annotate_figure(go_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(go_plt)
dev.off()

################################################################################
# Plotting for DE genes between Type 1 and Type 2 cells within each RAC

for(curr_cell_line in cell_lines){
  
  curr_intra_rac_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_rac_within_cluster_de_signatures.rds"))
  
  names(curr_intra_rac_signatures) <- paste0(curr_cell_line, " Cluster ",names(curr_intra_rac_signatures))
  
  all_dotplots <- genesets_characterization(curr_intra_rac_signatures, universe_to_use = cell_line_universes[[curr_cell_line]])
  
  png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/intra_rac_enrichment/",curr_cell_line,"_intra_rac_hallmarks.png"), width = 1000,height = 1000)
  hallmarks_plt <- egg::ggarrange(plots = all_dotplots$hallmarks_dotplots, ncol = 2)
  main_title <- paste0("Cancer Hallmarks Enrichment in DE Genes Between Type 1 and Type 2 Cells Within ", curr_cell_line ," RACs")
  hallmarks_plt <- annotate_figure(hallmarks_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
  plot(hallmarks_plt)
  dev.off()

  png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/intra_rac_enrichment/",curr_cell_line,"_intra_rac_mps.png"), width = 1000,height = 1000)
  mps_plt <- egg::ggarrange(plots = all_dotplots$mps_dotplots, ncol = 2)
  main_title <- paste0("ITH Meta-programs Enrichment in DE Genes Between Type 1 and Type 2 Cells Within ", curr_cell_line ," RACs")
  mps_plt <- annotate_figure(mps_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
  plot(mps_plt)
  dev.off()
  if(curr_cell_line == "A549"){
    png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/intra_rac_enrichment/",curr_cell_line,"_intra_rac_functional_enrichment.png"), width = 1000,height = 1400)
  } else if (curr_cell_line == "MCF7"){
    png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/intra_rac_enrichment/",curr_cell_line,"_intra_rac_functional_enrichment.png"), width = 1000,height = 1400)
  } else {
    png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/intra_rac_enrichment/",curr_cell_line,"_intra_rac_functional_enrichment.png"), width = 1000,height = 1200)
  }
  
  go_plt <- egg::ggarrange(plots = all_dotplots$go_dotplots, ncol = 2)
  main_title <- paste0("GO Functional Enrichment in DE Genes Between Type 1 and Type 2 Cells Within ", curr_cell_line ," RACs")
  go_plt <- annotate_figure(go_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
  plot(go_plt)
  dev.off()
}
################################################################################
# Plotting for GLOBAL DE genes between Type 1 and Type 2 cells within each RAC

global_intra_rac_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/global_intra_rac_signatures")

names(global_intra_rac_signatures) <- cell_lines

all_dotplots <- list()
for(curr_cell_line in cell_lines){
  
  curr_signature <- global_intra_rac_signatures[curr_cell_line]
  curr_dotplots <- genesets_characterization(curr_signature, universe_to_use = cell_line_universes[[curr_cell_line]])  
  
  all_dotplots[["hallmarks_dotplots"]] <- append(all_dotplots[["hallmarks_dotplots"]],curr_dotplots[["hallmarks_dotplots"]])
  all_dotplots[["mps_dotplots"]] <- append(all_dotplots[["mps_dotplots"]],curr_dotplots[["mps_dotplots"]])
  all_dotplots[["go_dotplots"]] <- append(all_dotplots[["go_dotplots"]],curr_dotplots[["go_dotplots"]])
  all_dotplots[["go_results"]] <- append(all_dotplots[["go_results"]],curr_dotplots[["go_results"]])
  
}


png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/intra_rac_enrichment/global_intra_rac_hallmarks.png"), width = 1000,height = 1000)
hallmarks_plt <- egg::ggarrange(plots = all_dotplots$hallmarks_dotplots, ncol = 2)
main_title <- paste0("Cancer Hallmarks Enrichment of DE Genes Between RAC Type 1 and Type 2 Cells")
hallmarks_plt <- annotate_figure(hallmarks_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(hallmarks_plt)
dev.off()

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/intra_rac_enrichment/global_intra_rac_mps.png"), width = 1000,height = 1000)
mps_plt <- egg::ggarrange(plots = all_dotplots$mps_dotplots, ncol = 2)
main_title <- paste0("ITH Meta-programs Enrichment of of DE Genes Between RAC Type 1 and Type 2 Cells")
mps_plt <- annotate_figure(mps_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(mps_plt)
dev.off()

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/intra_rac_enrichment/global_intra_rac_functional_enrichment.png"), width = 1000,height = 1000)
go_plt <- egg::ggarrange(plots = all_dotplots$go_dotplots, ncol = 2)
main_title <- paste0("GO Functional Enrichment of of DE Genes Between RAC Type 1 and Type 2 Cells")
go_plt <- annotate_figure(go_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(go_plt)
dev.off()



################################################################################
# Plotting for RAC supercluster signature

rac_supercluster_consensus_signature <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_supercluster_consensus_signature.rds")

names(rac_supercluster_consensus_signature) <- NULL

all_dotplots <- genesets_characterization(rac_supercluster_consensus_signature,universe_to_use = gene_universe_intersection)

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/supercluster_enrichment/supercluster_hallmarks.png"), width = 1000,height = 1000)
hallmarks_plt <- egg::ggarrange(plots = all_dotplots$hallmarks_dotplots, ncol = 2)
main_title <- paste0("Cancer Hallmarks Enrichment of Supercluster Signature")
hallmarks_plt <- annotate_figure(hallmarks_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(hallmarks_plt)
dev.off()

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/supercluster_enrichment/supercluster_mps.png"), width = 1000,height = 1000)
mps_plt <- egg::ggarrange(plots = all_dotplots$mps_dotplots, ncol = 2)
main_title <- paste0("ITH Meta-programs Enrichment of Supercluster Signature")
mps_plt <- annotate_figure(mps_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(mps_plt)
dev.off()

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/supercluster_enrichment/supercluster_functional_enrichment.png"), width = 1000,height = 1000)
go_plt <- egg::ggarrange(plots = all_dotplots$go_dotplots, ncol = 3)
main_title <- paste0("GO Functional Enrichment of Supercluster Signature")
go_plt <- annotate_figure(go_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(go_plt)
dev.off()

################################################################################
# Plotting for RAC type 1 superclusters signatures
type1_supercluster_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/type1_supercluster_signatures.rds")

names1 <- gsub("type1_", "", names(type1_supercluster_signatures))
names2 <- gsub("_", " ", names1)
names3 <- gsub("supercluster", "Supercluster ", names2)
names(type1_supercluster_signatures) <- gsub("signature", "Signature ", names3)

all_dotplots <- genesets_characterization(type1_supercluster_signatures, universe_to_use = gene_universe_intersection, num_pathways=20)

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/supercluster_enrichment/type1_superclusters_hallmarks.png"), width = 1000,height = 800)
hallmarks_plt <- egg::ggarrange(plots = all_dotplots$hallmarks_dotplots, ncol = 2)
main_title <- paste0("Cancer Hallmarks Enrichment of RAC Type 1 Superclusters Signatures")
hallmarks_plt <- annotate_figure(hallmarks_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(hallmarks_plt)
dev.off()

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/supercluster_enrichment/type1_superclusters_mps.png"), width = 1000,height = 800)
mps_plt <- egg::ggarrange(plots = all_dotplots$mps_dotplots, ncol = 2)
main_title <- paste0("ITH Meta-programs Enrichment of RAC Type 1 Superclusters Signatures")
mps_plt <- annotate_figure(mps_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(mps_plt)
dev.off()

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/supercluster_enrichment/type1_superclusters_functional_enrichment.png"), width = 1000,height = 800)
go_plt <- egg::ggarrange(plots = all_dotplots$go_dotplots, ncol = 2)
main_title <- paste0("GO Functional Enrichment of RAC Type 1 Superclusters Signatures")
go_plt <- annotate_figure(go_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(go_plt)
dev.off()




temp <- all_dotplots$go_results[[1]]@result %>% 
  dplyr::select(Description,GeneRatio,p.adjust) 


temp$Count <- sapply(temp$GeneRatio, function(x) as.numeric(strsplit(x, "/")[[1]][1]))

temp <- temp %>% 
  filter(p.adjust<.05)


################################################################################
# Plotting for RAC type 1 vs rest/inactive cells DOWNREGULATED signature

global_rac_type_down_signatures<- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/global_rac_type_down_signatures.rds")

global_rac_type_down_signatures <- global_rac_type_down_signatures[grepl("type1",names(global_rac_type_down_signatures))]

names(global_rac_type_down_signatures) <- cell_lines
all_dotplots <- genesets_characterization(global_rac_type_down_signatures,universe_to_use = gene_universe_intersection)


png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/type1_vs_inactive_enrichment/type1_vs_inactive_hallmarks_DOWN.png"), width = 1000,height = 1000)
hallmarks_plt <- egg::ggarrange(plots = all_dotplots$hallmarks_dotplots, ncol = 2)
main_title <- paste0("Cancer Hallmarks Enrichment of RAC Type 1 Cells Downregulated Genes")
hallmarks_plt <- annotate_figure(hallmarks_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(hallmarks_plt)
dev.off()

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/type1_vs_inactive_enrichment/type1_vs_inactive_mps_DOWN.png"), width = 1000,height = 1000)
mps_plt <- egg::ggarrange(plots = all_dotplots$mps_dotplots, ncol = 2)
main_title <- paste0("ITH Meta-programs Enrichment of RAC Type 1 Cells Downregulated Genes")
mps_plt <- annotate_figure(mps_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(mps_plt)
dev.off()

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/type1_vs_inactive_enrichment/type1_vs_inactive_functional_enrichment_DOWN.png"), width = 1000,height = 1000)
go_plt <- egg::ggarrange(plots = all_dotplots$go_dotplots, ncol = 2)
main_title <- paste0("GO Functional Enrichment of RAC Type 1 Cells Downregulated Genes")
go_plt <- annotate_figure(go_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(go_plt)
dev.off()

################################################################################
# Plotting for RAC type 1 superclusters DOWNREGULATED signatures

type1_supercluster_down_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/type1_supercluster_down_signatures.rds")

names1 <- gsub("type1_", "", names(type1_supercluster_down_signatures))
names2 <- gsub("_", " ", names1)
names3 <- gsub("supercluster", "Supercluster ", names2)
names(type1_supercluster_down_signatures) <- gsub("signature", "Signature ", names3)

all_dotplots <- genesets_characterization(type1_supercluster_down_signatures, universe_to_use = gene_universe)

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/supercluster_enrichment/type1_superclusters_hallmarks_DOWN.png"), width = 1000,height = 1000)
hallmarks_plt <- egg::ggarrange(plots = all_dotplots$hallmarks_dotplots, ncol = 2)
main_title <- paste0("Cancer Hallmarks Enrichment of RAC Type 1 Superclusters Signatures Downregulated Genes")
hallmarks_plt <- annotate_figure(hallmarks_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(hallmarks_plt)
dev.off()

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/supercluster_enrichment/type1_superclusters_mps_DOWN.png"), width = 1000,height = 1000)
mps_plt <- egg::ggarrange(plots = all_dotplots$mps_dotplots, ncol = 2)
main_title <- paste0("ITH Meta-programs Enrichment of RAC Type 1 Superclusters Signatures Downregulated Genes")
mps_plt <- annotate_figure(mps_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(mps_plt)
dev.off()

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/supercluster_enrichment/type1_superclusters_functional_enrichment_DOWN.png"), width = 1000,height = 1000)
go_plt <- egg::ggarrange(plots = all_dotplots$go_dotplots, ncol = 2)
main_title <- paste0("GO Functional Enrichment of RAC Type 1 Superclusters Signatures Downregulated Genes")
go_plt <- annotate_figure(go_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(go_plt)
dev.off()

################################################################################
# Plotting for shared genes between RACs signatures and yeast orthologs

global_rac_type_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/global_rac_type_signatures.rds")

global_rac_type1_signatures <- global_rac_type_signatures[grepl("type1",names(global_rac_type_signatures))]

names(global_rac_type1_signatures) <- cell_lines

yeast_upregulated_orthologs <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/yeast_upregulated_orthologs.rds")

all_dotplots <- list()
for(curr_cell_line in cell_lines){
  
  curr_signature <- list(intersect(global_rac_type1_signatures[[curr_cell_line]],yeast_upregulated_orthologs))
  names(curr_signature) <- curr_cell_line
  curr_dotplots <- genesets_characterization(curr_signature, universe_to_use = cell_line_universes[[curr_cell_line]])  
  
  all_dotplots[["hallmarks_dotplots"]] <- append(all_dotplots[["hallmarks_dotplots"]],curr_dotplots[["hallmarks_dotplots"]])
  all_dotplots[["mps_dotplots"]] <- append(all_dotplots[["mps_dotplots"]],curr_dotplots[["mps_dotplots"]])
  all_dotplots[["go_dotplots"]] <- append(all_dotplots[["go_dotplots"]],curr_dotplots[["go_dotplots"]])
  all_dotplots[["go_results"]] <- append(all_dotplots[["go_results"]],curr_dotplots[["go_results"]])
  
}

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/yeast_figures/yeast_hallmarks.png"), width = 1000,height = 1000)
hallmarks_plt <- egg::ggarrange(plots = all_dotplots$hallmarks_dotplots, ncol = 2)
main_title <- paste0("Hallmarks Enrichment of Shared Genes Between RACs and Yeast Stress Orthologs")
hallmarks_plt <- annotate_figure(hallmarks_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(hallmarks_plt)
dev.off()

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/yeast_figures/yeast_mps.png"), width = 1000,height = 1000)
mps_plt <- egg::ggarrange(plots = all_dotplots$mps_dotplots, ncol = 2)
main_title <- paste0("ITH Meta-programs Enrichment of Shared Genes Between RACs and Yeast Stress Orthologs")
mps_plt <- annotate_figure(mps_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(mps_plt)
dev.off()

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/yeast_figures/yeast_functional_enrichment.png"), width = 1000,height = 1000)
go_plt <- egg::ggarrange(plots = all_dotplots$go_dotplots, ncol = 2)
main_title <- paste0("GO Functional Enrichment of Shared Genes Between RACs and Yeast Stress Orthologs")
go_plt <- annotate_figure(go_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(go_plt)
dev.off()

################################################################################
# Plotting for shared genes between RAC Type 1 signatures and e coli orthologs

global_rac_type_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/global_rac_type_signatures.rds")
global_rac_type1_signatures <- global_rac_type_signatures[grepl("type1",names(global_rac_type_signatures))]
names(global_rac_type1_signatures) <- cell_lines

ecoli_AMR_genesets_orthologs <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ecoli_AMR_genesets_orthologs.rds")
ecoli_AMR_consensus_orthologs <- find_consensus_geneset(ecoli_AMR_genesets_orthologs,2)


all_dotplots <- list()
for(curr_cell_line in cell_lines){
  
  curr_signature <- list(intersect(global_rac_type1_signatures[[curr_cell_line]],ecoli_AMR_consensus_orthologs))
  names(curr_signature) <- curr_cell_line
  curr_dotplots <- genesets_characterization(curr_signature, universe_to_use = cell_line_universes[[curr_cell_line]])  
  
  all_dotplots[["hallmarks_dotplots"]] <- append(all_dotplots[["hallmarks_dotplots"]],curr_dotplots[["hallmarks_dotplots"]])
  all_dotplots[["mps_dotplots"]] <- append(all_dotplots[["mps_dotplots"]],curr_dotplots[["mps_dotplots"]])
  all_dotplots[["go_dotplots"]] <- append(all_dotplots[["go_dotplots"]],curr_dotplots[["go_dotplots"]])
  all_dotplots[["go_results"]] <- append(all_dotplots[["go_results"]],curr_dotplots[["go_results"]])
  
}

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/ecoli_figures/ecoli_hallmarks.png"), width = 1000,height = 1000)
hallmarks_plt <- egg::ggarrange(plots = all_dotplots$hallmarks_dotplots, ncol = 2)
main_title <- paste0("Hallmarks Enrichment of Shared Genes Between RACs and E Coli Resistance Orthologs")
hallmarks_plt <- annotate_figure(hallmarks_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(hallmarks_plt)
dev.off()

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/ecoli_figures/ecoli_mps.png"), width = 1000,height = 1000)
mps_plt <- egg::ggarrange(plots = all_dotplots$mps_dotplots, ncol = 2)
main_title <- paste0("ITH Meta-programs Enrichment of Shared Genes Between RACs and E Coli Resistance Orthologs")
mps_plt <- annotate_figure(mps_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(mps_plt)
dev.off()

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/ecoli_figures/ecoli_functional_enrichment.png"), width = 1000,height = 1000)
go_plt <- egg::ggarrange(plots = all_dotplots$go_dotplots, ncol = 2)
main_title <- paste0("GO Functional Enrichment of Shared Genes Between RACs and E Coli Resistance Orthologs")
go_plt <- annotate_figure(go_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(go_plt)
dev.off()



##################################################################################
# RAC Type 1 Superclusters and Yeast
type1_supercluster_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/type1_supercluster_signatures.rds")

names1 <- gsub("type1_", "", names(type1_supercluster_signatures))
names2 <- gsub("_", " ", names1)
names3 <- gsub("supercluster", "Supercluster ", names2)
names(type1_supercluster_signatures) <- gsub("signature", "Signature ", names3)

yeast_upregulated_orthologs <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/yeast_upregulated_orthologs.rds")


supercluster_ovelaps <- list()
for(curr_supercluster in type1_supercluster_signatures){
  supercluster_ovelaps <- append(supercluster_ovelaps, list(intersect(curr_supercluster,yeast_upregulated_orthologs)))
  
}
  
names(supercluster_ovelaps) <- names(type1_supercluster_signatures)

all_dotplots <- genesets_characterization(supercluster_ovelaps, universe_to_use = gene_universe_intersection, num_pathways=20)

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/yeast_figures/type1_superclusters_yeast_hallmarks.png"), width = 1000,height = 800)
hallmarks_plt <- egg::ggarrange(plots = all_dotplots$hallmarks_dotplots, ncol = 2)
main_title <- paste0("Cancer Hallmarks Enrichment of Shared Genes Between RAC Type 1 Superclusters Signatures and Yeast Orthologs")
hallmarks_plt <- annotate_figure(hallmarks_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(hallmarks_plt)
dev.off()

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/yeast_figures/type1_superclusters_yeast_mps.png"), width = 1000,height = 800)
mps_plt <- egg::ggarrange(plots = all_dotplots$mps_dotplots, ncol = 2)
main_title <- paste0("ITH Meta-programs Enrichment of Shared Genes Between RAC Type 1 Superclusters Signatures and Yeast Orthologs")
mps_plt <- annotate_figure(mps_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(mps_plt)
dev.off()

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/yeast_figures/type1_superclusters_yeast_functional_enrichment.png"), width = 1000,height = 800)
go_plt <- egg::ggarrange(plots = all_dotplots$go_dotplots, ncol = 2)
main_title <- paste0("GO Functional Enrichment of Shared Genes Between RAC Type 1 Superclusters Signatures and Yeast Orthologs")
go_plt <- annotate_figure(go_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(go_plt)
dev.off()


##################################################################################
# RAC Type 1 Superclusters and E coli
type1_supercluster_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/type1_supercluster_signatures.rds")

names1 <- gsub("type1_", "", names(type1_supercluster_signatures))
names2 <- gsub("_", " ", names1)
names3 <- gsub("supercluster", "Supercluster ", names2)
names(type1_supercluster_signatures) <- gsub("signature", "Signature ", names3)

ecoli_AMR_genesets_orthologs <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ecoli_AMR_genesets_orthologs.rds")
ecoli_AMR_consensus_orthologs <- find_consensus_geneset(ecoli_AMR_genesets_orthologs,2)

supercluster_ovelaps <- list()
for(curr_supercluster in type1_supercluster_signatures){
  supercluster_ovelaps <- append(supercluster_ovelaps, list(intersect(curr_supercluster,ecoli_AMR_consensus_orthologs)))
  
}

names(supercluster_ovelaps) <- names(type1_supercluster_signatures)


all_dotplots <- genesets_characterization(supercluster_ovelaps, universe_to_use = gene_universe_intersection, num_pathways=20)

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/ecoli_figures/type1_superclusters_ecoli_hallmarks.png"), width = 1000,height = 800)
hallmarks_plt <- egg::ggarrange(plots = all_dotplots$hallmarks_dotplots, ncol = 2)
main_title <- paste0("Cancer Hallmarks Enrichment of Shared Genes Between RAC Type 1 Superclusters Signatures and E. coli Orthologs")
hallmarks_plt <- annotate_figure(hallmarks_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(hallmarks_plt)
dev.off()

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/ecoli_figures/type1_superclusters_ecoli_mps.png"), width = 1000,height = 800)
mps_plt <- egg::ggarrange(plots = all_dotplots$mps_dotplots, ncol = 2)
main_title <- paste0("ITH Meta-programs Enrichment of Shared Genes Between RAC Type 1 Superclusters Signatures and E. coli Orthologs")
mps_plt <- annotate_figure(mps_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(mps_plt)
dev.off()

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/ecoli_figures/type1_superclusters_ecoli_functional_enrichment.png"), width = 1000,height = 800)
go_plt <- egg::ggarrange(plots = all_dotplots$go_dotplots, ncol = 2)
main_title <- paste0("GO Functional Enrichment of Shared Genes Between RAC Type 1 Superclusters Signatures and E. coli Orthologs")
go_plt <- annotate_figure(go_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(go_plt)
dev.off()





