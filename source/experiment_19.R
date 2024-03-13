library(tidyverse)
library(Seurat)
library(msigdbr)
library(org.Hs.eg.db)
library(ggpubr)
source("source/read_in_all_cell_lines.R")
source("source/cole_functions.R")

cell_lines <- c("A549","K562","MCF7")
RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- cell_lines

create_barplots <- function(ranks,signature_name){
  m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::select(gs_name, human_gene_symbol)
  # hallmark_uni <- unique(c(unlist(all_signatures), m_t2g$human_gene_symbol))
  
  mp_t2g <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ith_meta_programs_t2g.rds")
  
  specifc_mps <- c("MP39 Metal-response","MP31 Alveolar","MP29 NPC/OPC","MP28 Oligo normal","MP27 Oligo Progenitor","MP38 Glutathione","MP41 Unassigned","MP35 Hemato-related-I","MP37 Hemato-related-II","MP32 Skin-pigmentation","MP36 IG","MP16 MES (glioma)","MP15 EMT IV")
  
  mp_t2g <- mp_t2g %>% 
    filter(!gs_name %in% specifc_mps)
  # Meta-programs
  # mp_gsea <- GSEA(ranks, TERM2GENE = mp_t2g)
  # 
  # df <- mp_gsea@result %>% 
  #   filter(p.adjust < 0.05) %>% 
  #   arrange(desc(NES)) %>% 
  #   dplyr::select(Description,NES,p.adjust)
  # 
  # 
  # if(nrow(df)> 10){
  #   df <- df[1:10,]
  # }
  # mp_plot <- ggplot(df)+
  #   geom_col(aes(x=(fct_reorder(Description,NES)),y=NES, fill=NES))+
  #   scale_fill_continuous(high="red",low="pink")+
  #   xlab("Meta-Program")+
  #   ggtitle(paste0(curr_cell_line, " ", signature_name))+
  #   theme(plot.title = element_text(size=30),
  #         axis.text = element_text(size=15),
  #         axis.title = element_text(size=20,face="bold"))+
  #   coord_flip()
  # 
  # 
  # # Hallmarks
  # hallmark_gsea <- GSEA(ranks, TERM2GENE = m_t2g)
  # 
  # df <- hallmark_gsea@result %>%
  #   filter(p.adjust < 0.05) %>%
  #   arrange(desc(NES)) %>%
  #   dplyr::select(Description,NES,p.adjust)
  # 
  # 
  # hallmark_plot <- ggplot(df[1:10,])+
  #   geom_col(aes(x=(fct_reorder(Description,NES)),y=NES, fill=NES))+
  #   scale_fill_continuous(high="red",low="pink")+
  #   xlab("Cancer Hallmark")+
  #   ggtitle(paste0(curr_cell_line, " ", signature_name))+
  #   theme(plot.title = element_text(size=30),
  #         axis.text = element_text(size=15),
  #         axis.title = element_text(size=20,face="bold"))+
  #   coord_flip()
  
  
  # GO Pathways
  go_gsea <- gseGO(geneList     = ranks,
                   OrgDb        = org.Hs.eg.db,
                   ont          = "ALL",
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
    ggtitle(paste0(curr_cell_line, " ", signature_name))+
    theme(plot.title = element_text(size=30),
          axis.text = element_text(size=15),
          axis.title = element_text(size=20,face="bold"))+
    coord_flip()
  
  
  kegg_t2g <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/kegg_t2g.rds")
  kegg_gsea <- GSEA(ranks, TERM2GENE = kegg_t2g)

  df <- kegg_gsea@result %>%
    filter(p.adjust < 0.05) %>%
    arrange(desc(NES)) %>%
    dplyr::select(Description,NES,p.adjust)


  if(nrow(df)> 10){
    df <- df[1:10,]
  }
  kegg_plot <- ggplot(df)+
    geom_col(aes(x=(fct_reorder(Description,NES)),y=NES, fill=NES))+
    scale_fill_continuous(high="red",low="pink")+
    xlab("KEGG Pathways")+
    ggtitle(paste0(curr_cell_line, " ", signature_name))+
    theme(plot.title = element_text(size=30),
          axis.text = element_text(size=15),
          axis.title = element_text(size=20,face="bold"))+
    coord_flip()

  
  
  # plot_list <- list(hallmark_plot,mp_plot,go_plot)
  plot_list <- list(go_plot,kegg_plot)
  
  p <- ggarrange(plotlist = plot_list,ncol=2,nrow=1)
  
  
  return(p)
}

####################################################################################

all_signatures <- list()
curr_cell_line <- cell_lines[2]
for(curr_cell_line in cell_lines){
  # curr_de <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/",curr_cell_line,"_cluster_all_markers_de.rds"))
  curr_de <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/",curr_cell_line, "_cluster_all_markers_de.rds"))
  
  curr_clusters <- unique(curr_de$cluster)
  
  all_plots <- list()
  signatures <- list()
  for(i in RACs[[curr_cell_line]]){
    
    cat(i,"\n")
    
    temp <- curr_de %>% 
      filter(cluster == i & avg_log2FC > 0 & p_val_adj < 0.05) %>% 
      arrange(desc(avg_log2FC)) %>% 
      pull(gene) 
    
    signatures <- append(signatures, list(temp))
    
    
    temp <- curr_de %>% 
      filter(cluster == i) %>% 
      arrange(desc(avg_log2FC)) %>% 
      dplyr::select(gene,avg_log2FC) %>% 
      deframe()
    
    signature_name <- paste0("Cluster ", i)
    all_plots <- append(all_plots,list(create_barplots(temp,signature_name)))
    
  }
  
  
  p <- ggarrange(plotlist = all_plots,ncol=2,nrow=4)
  p
  
  
  
  all_signatures[[curr_cell_line]] <- curr_signatures
}


library(tidyverse)

data <- all_data[["MCF7"]]

data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")



ggplot(data@meta.data) +
  geom_boxplot(aes(x=Cluster, y =percent.mt))
  

