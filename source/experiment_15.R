library(tidyverse)
library(Seurat)
source("source/cole_functions.R")

cell_lines <- c("A549","K562","MCF7")


all_signatures <- list()
for(curr_cell_line in cell_lines){
  curr_de <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/",curr_cell_line,"_cluster_all_markers_de.rds"))
  
  curr_clusters <- unique(curr_de$cluster)
  
  curr_signatures <- list()
  for(i in curr_clusters){
    
    curr_signatures[[i]] <- curr_de %>% 
      filter(cluster == i & avg_log2FC > 0 & p_val_adj < 0.05) %>% 
      pull(gene)
    
  }
  
  all_signatures[[curr_cell_line]] <- curr_signatures
}


supercluster1_signatures <- c(list(all_signatures[["A549"]][[9]]),list(all_signatures[["K562"]][[5]]),list(all_signatures[["MCF7"]][[8]]))

supercluster1_signature <- list("supercluster1_signature" = find_consensus_geneset(supercluster1_signatures,2))

##############

supercluster2_signatures <- c(list(all_signatures[["A549"]][[19]]),list(all_signatures[["K562"]][[11]]),list(all_signatures[["MCF7"]][[5]]))

supercluster2_signature <- list("supercluster2_signature" = find_consensus_geneset(supercluster2_signatures,2))

##############

supercluster3_signatures <- c(list(all_signatures[["A549"]][[14]]),list(all_signatures[["K562"]][[9]]),list(all_signatures[["MCF7"]][[13]]))

supercluster3_signature <- list("supercluster3_signature" = find_consensus_geneset(supercluster3_signatures,2))

supercluster_signatures <- c(supercluster1_signature,supercluster2_signature,supercluster3_signature)

saveRDS(supercluster_signatures,"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_supercluster_signatures.rds")

######################################################################################
all_ranks <- list()
for(curr_cell_line in cell_lines){
  curr_de <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/",curr_cell_line,"_cluster_all_markers_de.rds"))
  
  curr_clusters <- unique(curr_de$cluster)
  
  curr_signatures <- list()
  for(i in curr_clusters){
    
    curr_signatures[[i]] <- curr_de %>% 
      filter(cluster==i) %>% 
      arrange(desc(avg_log2FC)) %>% 
      dplyr::select(gene,avg_log2FC) %>% 
      deframe()
    
  }
  
  all_ranks[[curr_cell_line]] <- curr_signatures
}


curr_cell_line <- "K562"
cluster_num <- 5
ranks <- all_ranks[[curr_cell_line]][[cluster_num]]

# Meta-programs
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
  ggtitle(paste0(curr_cell_line, " Cluster ",cluster_num," Enrichment"))+
  theme(plot.title = element_text(size=30),
        axis.text = element_text(size=15),
        axis.title = element_text(size=20,face="bold"))+
  coord_flip()


# Hallmarks
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
#   ggtitle(paste0(curr_cell_line, " emergent RAC Enrichment"))+
#   theme(plot.title = element_text(size=30),
#         axis.text = element_text(size=15),
#         axis.title = element_text(size=20,face="bold"))+
#   coord_flip()


# GO Pathways
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
  ggtitle("")+
  theme(plot.title = element_text(size=30),
        axis.text = element_text(size=15),
        axis.title = element_text(size=20,face="bold"))+
  coord_flip()




plot_list <- list(mp_plot,go_plot)


p <- ggarrange(plotlist = plot_list,ncol=2,nrow=1)


p






