library(tidyverse)
library(Seurat)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)

cell_lines <- c("A549","K562","MCF7")


curr_cell_line <- "MCF7"

all_ranks <- list()
rac_signatures <- list()
for(curr_cell_line in cell_lines){
  data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))
  # RACs <- list(c(4,9,12,13,14,16,18),c(4,5,11),c(5,8,12,13,17))
  RACs <- list(c(12),c(11),c(8))
  names(RACs) <- c("A549","K562","MCF7")
  clusters_of_interest <- RACs[[curr_cell_line]]
  
  
  scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
  threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
  active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
  clusters_of_interest <- RACs[[curr_cell_line]]
  
  
  data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, data$Cluster, "Non-RAC"), col.name = "RAC")
  #Add metadata for RAC and Cell Group
  data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
  data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, "1", ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), 2, 0)), col.name = "cell_group")
  data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, paste0(data$Cluster, "_1"), ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), paste0(data$Cluster, "_2"), paste0(data$Cluster, "_0"))), col.name = "cell_cluster_group")
  
  
  # DimPlot(data, group.by='Cluster',split.by = "rac")
  
  
  Idents(data) <- data$rac
  de_res <- FindAllMarkers(data)
  
  
  ranks <- de_res %>% 
    filter(cluster == "rac") %>% 
    arrange(desc(avg_log2FC)) %>% 
    dplyr::select(gene,avg_log2FC) %>% 
    deframe()
  
  
  all_ranks[[curr_cell_line]] <- ranks
  
  curr_signature <- de_res %>% 
    filter(p_val_adj<0.05 & avg_log2FC > 0 & cluster=="rac") %>% 
    arrange(desc(avg_log2FC)) %>% 
    pull(gene)
  
  
  rac_signatures[[curr_cell_line]] <- curr_signature[1:200]
  
  
}


shared_genes <- intersect(intersect(names(all_ranks[[1]]),names(all_ranks[[2]])),names(all_ranks[[3]]))

all_ranks[[1]] <- all_ranks[[1]][names(all_ranks[[1]]) %in% shared_genes]
all_ranks[[2]] <- all_ranks[[2]][names(all_ranks[[2]]) %in% shared_genes]
all_ranks[[3]] <- all_ranks[[3]][names(all_ranks[[3]]) %in% shared_genes]


A549_ranks <- data.frame(names(all_ranks[[1]]),all_ranks[[1]])
K562_ranks <- data.frame(names(all_ranks[[2]]),all_ranks[[2]])
MCF7_ranks <- data.frame(names(all_ranks[[3]]),all_ranks[[3]])

colnames(A549_ranks) <- c("gene","value")
colnames(K562_ranks) <- c("gene","value")
colnames(MCF7_ranks) <- c("gene","value")



colnames(temp)

temp <- merge(A549_ranks,K562_ranks,by="gene")

df <- merge(temp,MCF7_ranks,by="gene")


ranks <- df %>% 
  column_to_rownames("gene") %>% 
  rowMeans()

ranks <- sort(ranks,decreasing = T)

# Meta-programs
mp_gsea <- GSEA(ranks, TERM2GENE = mp_t2g)

df <- mp_gsea@result %>% 
  filter(p.adjust < 0.05) %>% 
  arrange(desc(NES)) %>% 
  dplyr::select(Description,NES,p.adjust)



mp_plot <- ggplot(df[1:10,])+
  geom_col(aes(x=(fct_reorder(Description,NES)),y=NES, fill=NES))+
  scale_fill_continuous(high="red",low="pink")+
  xlab("Meta-Program")+
  ggtitle(paste0(curr_cell_line, " RAC Enrichment"))+
  theme(plot.title = element_text(size=30),
        axis.text = element_text(size=15),
        axis.title = element_text(size=20,face="bold"))+
  coord_flip()


# Hallmarks
hallmark_gsea <- GSEA(ranks, TERM2GENE = m_t2g)

df <- hallmark_gsea@result %>% 
  filter(p.adjust < 0.05) %>% 
  arrange(desc(NES)) %>% 
  dplyr::select(Description,NES,p.adjust)


hallmark_plot <- ggplot(df[1:10,])+
  geom_col(aes(x=(fct_reorder(Description,NES)),y=NES, fill=NES))+
  scale_fill_continuous(high="red",low="pink")+
  xlab("Cancer Hallmark")+
  ggtitle(paste0(curr_cell_line, " RAC Enrichment"))+
  theme(plot.title = element_text(size=30),
        axis.text = element_text(size=15),
        axis.title = element_text(size=20,face="bold"))+
  coord_flip()


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

go_plot <- ggplot(df[1:10,])+
  geom_col(aes(x=(fct_reorder(Description,NES)),y=NES, fill=NES))+
  scale_fill_continuous(high="red",low="pink")+
  xlab("GO Pathways")+
  ggtitle(paste0(curr_cell_line, " RAC Enrichment"))+
  theme(plot.title = element_text(size=30),
        axis.text = element_text(size=15),
        axis.title = element_text(size=20,face="bold"))+
  coord_flip()




plot_list <- list(hallmark_plot,mp_plot,go_plot)


p <- ggarrange(plotlist = plot_list,ncol=3,nrow=1)


p




supercluster_signature <- find_consensus_geneset(rac_signatures,1)




curr_geneset <- supercluster_signature
universe_to_use <- gene_universe_intersection

# Hallmarks Enrichment
hallmark_enrichment_results <- enricher(curr_geneset, TERM2GENE=m_t2g,
                                        universe = universe_to_use)


hallmark_plot <- barplot(hallmark_enrichment_results)


# MPs Enrichment
mp_enrichment_results <- enricher(curr_geneset, TERM2GENE=mp_t2g,
                                  universe = universe_to_use)



mp_plot <- barplot(mp_enrichment_results)

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


go_plot <- barplot(ego)


plot_list <- list(hallmark_plot,mp_plot,go_plot)


p <- ggarrange(plotlist=plot_list,ncol=3)

annotate_figure(p)




