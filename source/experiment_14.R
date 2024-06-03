setwd("/data/ruoffcj/projects/drug_treatment/")
source("source/read_in_all_cell_lines.R")
library(tidyverse)
library(Seurat)



# Using FindALlMarkers
# 
# pre_rac_signatures <- list()
# for(curr_cell_line in cell_lines){
#   data <- all_data[[curr_cell_line]]
#   
#   pre_clusters <- data@meta.data %>%
#     count(treatment_stage,Cluster) %>%
#     filter(treatment_stage == "pre" & n > 10) %>%
#     pull(Cluster)
# 
#   pre_data <- data[,data$treatment_stage=='pre' & data$Cluster %in% pre_clusters]
#   
#   Idents(pre_data) <- pre_data$emergent_rac
#   
#   
#   de_res <- FindAllMarkers(pre_data)
#   
#   pre_rac_signature <- de_res %>% 
#     filter(cluster == "non_emergent_rac" & avg_log2FC > 0 & p_val_adj <0.05) %>% 
#     arrange(desc(avg_log2FC)) %>% 
#     pull(gene)
#   
#   
#   pre_rac_signature <- pre_rac_signature[!grepl("MT-",pre_rac_signature)]
#   
#   pre_rac_signatures[[curr_cell_line]] <- pre_rac_signature[1:200]
# }



# Using FindMarkers
# curr_cell_line <- "K562"
# 
# all_rac_type_signatures <- list()
# all_rac_type_ranks <- list()
# 
# for(curr_cell_line in cell_lines){
#   data <- all_data[[curr_cell_line]]
#   
#   post_data <- data[,data$treatment_stage=='post']
#   
#   
#   Idents(post_data) <- post_data$emergent_rac
#   
#   #####
#   # Pre Existing RACs vs Non-RACs
#   cat("1\n")
#   pre_vs_non_de_res <- FindMarkers(post_data, ident.1="non_emergent_rac",ident.2 = "non_rac")
#   
#   curr_signature <- pre_vs_non_de_res %>% 
#     filter(avg_log2FC > 0 & p_val_adj <0.05) %>% 
#     arrange(desc(avg_log2FC)) %>% 
#     rownames_to_column() %>% 
#     pull(rowname)
#   
#   curr_ranks <- pre_vs_non_de_res %>% 
#     arrange(desc(avg_log2FC)) %>% 
#     rownames_to_column() %>% 
#     dplyr::select(rowname,avg_log2FC) %>% 
#     deframe()
#   
#   
#   curr_signature <- curr_signature[!grepl("MT-",curr_signature)]
#   curr_ranks <- curr_ranks[!grepl("MT-",names(curr_ranks))]
#   
#   all_rac_type_signatures[[curr_cell_line]][["pre_vs_non"]] <- curr_signature
#   all_rac_type_ranks[[curr_cell_line]][["pre_vs_non"]] <- curr_ranks
#  
#   #####
#   # Emergent RACs vs Non-RACs
#   cat("2\n")
#   em_vs_non_de_res <- FindMarkers(post_data, ident.1="emergent_rac",ident.2 = "non_rac")
#   
#   curr_signature <- em_vs_non_de_res %>% 
#     filter(avg_log2FC > 0 & p_val_adj <0.05) %>% 
#     arrange(desc(avg_log2FC)) %>% 
#     rownames_to_column() %>% 
#     pull(rowname)
#   
#   curr_ranks <- em_vs_non_de_res %>% 
#     arrange(desc(avg_log2FC)) %>% 
#     rownames_to_column() %>% 
#     dplyr::select(rowname,avg_log2FC) %>% 
#     deframe()
#   
#   
#   curr_signature <- curr_signature[!grepl("MT-",curr_signature)]
#   curr_ranks <- curr_ranks[!grepl("MT-",names(curr_ranks))]
#   
#   all_rac_type_signatures[[curr_cell_line]][["em_vs_non"]] <- curr_signature
#   all_rac_type_ranks[[curr_cell_line]][["em_vs_non"]] <- curr_ranks
#   
#   #####
#   # Emergent RACs vs Pre Existing 
#   cat("3\n")
#   em_vs_pre_de_res <- FindMarkers(post_data, ident.1="emergent_rac",ident.2 = "non_emergent_rac")
#   
#   curr_signature <- em_vs_pre_de_res %>% 
#     filter(avg_log2FC > 0 & p_val_adj <0.05) %>% 
#     arrange(desc(avg_log2FC)) %>% 
#     rownames_to_column() %>% 
#     pull(rowname)
#   
#   curr_ranks <- em_vs_pre_de_res %>% 
#     arrange(desc(avg_log2FC)) %>% 
#     rownames_to_column() %>% 
#     dplyr::select(rowname,avg_log2FC) %>% 
#     deframe()
#   
#   
#   curr_signature <- curr_signature[!grepl("MT-",curr_signature)]
#   curr_ranks <- curr_ranks[!grepl("MT-",names(curr_ranks))]
#   
#   all_rac_type_signatures[[curr_cell_line]][["em_vs_pre"]] <- curr_signature
#   all_rac_type_ranks[[curr_cell_line]][["em_vs_pre"]] <- curr_ranks
#   
#   #####
#   # Pre Existing vs Emergent RACs
#   cat("4\n")
#   pre_vs_em_de_res <- FindMarkers(post_data, ident.1="non_emergent_rac",ident.2 = "emergent_rac")
#   
#   curr_signature <- pre_vs_em_de_res %>% 
#     filter(avg_log2FC > 0 & p_val_adj <0.05) %>% 
#     arrange(desc(avg_log2FC)) %>% 
#     rownames_to_column() %>% 
#     pull(rowname)
#   
#   curr_ranks <- pre_vs_em_de_res %>% 
#     arrange(desc(avg_log2FC)) %>% 
#     rownames_to_column() %>% 
#     dplyr::select(rowname,avg_log2FC) %>% 
#     deframe()
#   
#   
#   curr_signature <- curr_signature[!grepl("MT-",curr_signature)]
#   curr_ranks <- curr_ranks[!grepl("MT-",names(curr_ranks))]
#   
#   all_rac_type_signatures[[curr_cell_line]][["pre_vs_em"]] <- curr_signature
#   all_rac_type_ranks[[curr_cell_line]][["pre_vs_em"]] <- curr_ranks
#   
# }
# 
# 
# saveRDS(all_rac_type_signatures, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/experiment_data/experiment_14_all_rac_type_signatures.rds")
# saveRDS(all_rac_type_ranks, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/experiment_data/experiment_14_all_rac_type_ranks.rds")


#################################################################################
# Plot enrichment results


all_rac_type_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/experiment_data/experiment_14_all_rac_type_signatures.rds")
all_rac_type_ranks <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/experiment_data/experiment_14_all_rac_type_ranks.rds")



plot_gsea_enrichment <- function(ranks,plot_title){
  # Meta-programs
  cat("Meta-programs\n")
  mp_gsea <- GSEA(ranks, TERM2GENE = mp_t2g)
  
  df <- mp_gsea@result %>% 
    filter(p.adjust < 0.05) %>% 
    arrange(desc(NES)) %>% 
    dplyr::select(Description,NES,p.adjust)
  
  
  if(nrow(df)> 10){
    df <- df[1:10,]
  }
  
  mp_plot <- ggplot(df)+
    geom_col(aes(x=(fct_reorder(Description,NES)),y=NES, fill=NES))+
    scale_fill_continuous(high="red",low="pink")+
    xlab("Meta-Program")+
    ggtitle(paste0(curr_cell_line, " ", plot_title))+
    theme(plot.title = element_text(size=30),
          axis.text = element_text(size=15),
          axis.title = element_text(size=20,face="bold"))+
    coord_flip()
  
  
  # Hallmarks
  cat("Hallmarks\n")
  hallmark_gsea <- GSEA(ranks, TERM2GENE = m_t2g)
  
  df <- hallmark_gsea@result %>%
    filter(p.adjust < 0.05) %>%
    arrange(desc(NES)) %>%
    dplyr::select(Description,NES,p.adjust)
  
  
  hallmark_plot <- ggplot(df[1:10,])+
    geom_col(aes(x=(fct_reorder(Description,NES)),y=NES, fill=NES))+
    scale_fill_continuous(high="red",low="pink")+
    xlab("Cancer Hallmark")+
    ggtitle(paste0(curr_cell_line, " ", plot_title))+
    theme(plot.title = element_text(size=30),
          axis.text = element_text(size=15),
          axis.title = element_text(size=20,face="bold"))+
    coord_flip()
  
  
  # GO Pathways
  cat("GO Pathways\n")
  go_gsea <- gseGO(geneList     = ranks,
                   OrgDb        = org.Hs.eg.db,
                   ont          = "BP",
                   minGSSize    = 100,
                   maxGSSize    = 500,
                   pvalueCutoff = 0.01,
                   verbose      = FALSE,
                   keyType = "SYMBOL")
  
  
  df <- go_gsea@result %>% 
    filter(p.adjust < 0.05) %>% 
    arrange(desc(NES)) %>% 
    dplyr::select(Description,NES,p.adjust)
  
  if(nrow(df)> 10){
    df <- df[1:10,]
  }
  
  go_plot <- ggplot(df)+
    geom_col(aes(x=(fct_reorder(Description,NES)),y=NES, fill=NES))+
    scale_fill_continuous(high="red",low="pink")+
    xlab("GO Pathways")+
    ggtitle(paste0(curr_cell_line, " ", plot_title))+
    theme(plot.title = element_text(size=30),
          axis.text = element_text(size=15),
          axis.title = element_text(size=20,face="bold"))+
    coord_flip()
  
  
  
  plot_list <- list(hallmark_plot,mp_plot,go_plot)
  
  
  p <- ggarrange(plotlist = plot_list,ncol=3,nrow=1)
  
  
  return(p)
}

plot_ora_enrichment <- function(signature,plot_title){

  all_dotplots <- genesets_characterization(signature, cell_line_universes[[curr_cell_line]])
  
  hallmarks_plt <- ggarrange(plotlist = all_dotplots$hallmarks_dotplots, ncol = 1, common.legend = T, legend=c("right"))
  main_title <- paste0("\nCancer Hallmarks")
  hallmarks_plt <- annotate_figure(hallmarks_plt, top = text_grob(main_title, color = "black", face = "bold", size = 30))
  
  
  mps_plt <- ggarrange(plotlist = all_dotplots$mps_dotplots, ncol = 1, common.legend = T,legend=c("right"))
  main_title <- paste0("\nITH Meta-programs")
  mps_plt <- annotate_figure(mps_plt, top = text_grob(main_title, color = "black", face = "bold", size = 30))
  
  
  go_plt <- ggarrange(plotlist = all_dotplots$go_dotplots, ncol = 1, common.legend = T,legend=c("right"))
  main_title <- paste0("\nGO Pathways")
  go_plt <- annotate_figure(go_plt, top = text_grob(main_title, color = "black", face = "bold", size = 30))
  
  
  plots <- list(hallmarks_plt,mps_plt, go_plt)
  
  figure <- ggarrange(plotlist = plots, ncol=3, nrow=1, common.legend = T,legend=c("right"))
  
  p <- annotate_figure(figure, left = text_grob("", rot = 90, vjust = 1, size=35, face="bold"),
                        bottom = text_grob("", size=35, face="bold"),
                        top=text_grob(paste0(plot_title, " Enrichment"), size=40, face="bold"))
  
  return(p)
}

genesets_characterization <- function(genesets_to_use, universe_to_use, num_pathways=10){
  m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::select(gs_name, human_gene_symbol)
  # hallmark_uni <- unique(c(unlist(all_signatures), m_t2g$human_gene_symbol))
  
  mp_t2g <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ith_meta_programs_t2g.rds")
  
  specifc_mps <- c("MP39 Metal-response","MP31 Alveolar","MP29 NPC/OPC","MP28 Oligo normal","MP27 Oligo Progenitor","MP38 Glutathione","MP41 Unassigned","MP35 Hemato-related-I","MP37 Hemato-related-II","MP32 Skin-pigmentation","MP36 IG","MP16 MES (glioma)","MP15 EMT IV")
  
  mp_t2g <- mp_t2g %>% 
    filter(!gs_name %in% specifc_mps)
  
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
      p <- barplot(hallmark_enrichment_results,
                   showCategory = num_pathways, font.size=20) + 
        ggtitle(paste0("", curr_cluster))+
        theme(plot.title = element_text(size=20))
      
      hallmarks_dotplots <- append(hallmarks_dotplots, list(p))
    }
    if(!is.null(mp_enrichment_results) && nrow(mp_enrichment_results) > 0){
      p <- barplot(mp_enrichment_results,
                   showCategory = num_pathways, font.size=20)+
        ggtitle(paste0("", curr_cluster))+
        theme(plot.title = element_text(size=20))
      
      mps_dotplots <- append(mps_dotplots, list(p))
    }
    if(nrow(ego) > 0){
      p <- barplot(ego,
                   showCategory = num_pathways, font.size=20) + 
        ggtitle(paste0("", curr_cluster))+
        theme(plot.title = element_text(size=20))
      
      go_dotplots <- append(go_dotplots, list(p))
    }
  }
  
  
  ret_list <- list(hallmarks_dotplots, mps_dotplots, go_dotplots,go_results)
  names(ret_list) <- c("hallmarks_dotplots", "mps_dotplots", "go_dotplots","go_results")
  
  return(ret_list) 
}

# GSEA
for(curr_cell_line in cell_lines){
  
  ranks <- all_rac_type_ranks[[curr_cell_line]]$pre_vs_non
  plot_title <- "pre-rac vs non-rac"
  
  p1 <- plot_gsea_enrichment(ranks, plot_title)
  
  
  ranks <- all_rac_type_ranks[[curr_cell_line]]$em_vs_non
  plot_title <- "em-rac vs non-rac"
  
  p2 <- plot_gsea_enrichment(ranks, plot_title)
  
  
  ranks <- all_rac_type_ranks[[curr_cell_line]]$em_vs_pre
  plot_title <- "em-rac vs pre-rac"
  
  p3 <- plot_gsea_enrichment(ranks, plot_title)
  
  
  ranks <- all_rac_type_ranks[[curr_cell_line]]$pre_vs_em
  plot_title <- "pre-rac vs em-rac"
  
  p4 <- plot_gsea_enrichment(ranks, plot_title)
  
  
  
  plot_list <- list(p1,p2,p3,p4)
  
  
  p_final <- ggarrange(plotlist = plot_list,ncol=1,nrow=4)
  

  p_final
}






# ORA
cell_line_universes <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/cell_line_universes.rds")
for(curr_cell_line in cell_lines){
  
  signature <- list(all_rac_type_signatures[[curr_cell_line]]$pre_vs_non)
  plot_title <- "pre-rac vs non-rac"
  
  p1 <- plot_ora_enrichment(signature,plot_title)
  
  
  signature <- list(all_rac_type_signatures[[curr_cell_line]]$em_vs_non)
  plot_title <- "em-rac vs non-rac"
  
  
  p2 <- plot_ora_enrichment(signature,plot_title)
  
  
  signature <- list(all_rac_type_signatures[[curr_cell_line]]$em_vs_pre)
  plot_title <- "em-rac vs pre-rac"
  
  p3 <- plot_ora_enrichment(signature,plot_title)
  
  
  signature <- list(all_rac_type_signatures[[curr_cell_line]]$pre_vs_em)
  plot_title <- "pre-rac vs em-rac"
  
  p4 <- plot_ora_enrichment(signature,plot_title)
  
  
  
  plot_list <- list(p1,p2,p3,p4)
  
  
  p_final <- ggarrange(plotlist = plot_list,ncol=1,nrow=4)
  
  
  p_final
}







































































