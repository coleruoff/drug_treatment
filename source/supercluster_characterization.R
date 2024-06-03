library(clusterProfiler)
library(msigdbr)
source("source/cole_functions.R")

cell_lines <- c("A549","K562","MCF7")

A549_active <- c(4,9,13)

K562_active <- c(5)

MCF7_active <-  c(5,8,17)

supercluster_components <- list(A549_active, K562_active, MCF7_active)
names(supercluster_components) <- cell_lines

signature_length <- 200

curr_cell_line <- "A549"
A549_cluster_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_cluster_signatures_",signature_length,".rds"))


curr_cell_line <- "K562"
K562_cluster_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_cluster_signatures_",signature_length,".rds"))


curr_cell_line <- "MCF7"
MCF7_cluster_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_cluster_signatures_",signature_length,".rds"))

cell_line_cluster_signatures <- c(list(A549_cluster_signatures), list(K562_cluster_signatures), list(MCF7_cluster_signatures))
names(cell_line_cluster_signatures) <- cell_lines

A549_temp <- unlist(A549_cluster_signatures[A549_active])

K562_temp <- unlist(K562_cluster_signatures[K562_active])

MCF7_temp <- unlist(MCF7_cluster_signatures[MCF7_active])

supercluster_markers <- unique(intersect(intersect(A549_temp,K562_temp),MCF7_temp))

temp <- unique(c(A549_temp,K562_temp,MCF7_temp))

supercluster_markers <- supercluster_markers[!grepl("^MT",supercluster_markers)]

supercluster_signature <- list("supercluster_signature" = supercluster_markers)

supercluster_consensus_signature <- find_consensus_geneset(c(A549_cluster_signatures[A549_active],K562_cluster_signatures[K562_active],MCF7_cluster_signatures[MCF7_active]),2)

supercluster_consensus_signature <- supercluster_consensus_signature[!grepl("^MT",supercluster_consensus_signature)]

supercluster_consensus_signature <- supercluster_consensus_signature[!grepl("^RP",supercluster_consensus_signature)]

saveRDS(supercluster_signature, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/supercluster_signature.rds")


#################################################################################
# GO Functional Enrichment of shared markers
#################################################################################

go_enrich <- enrichGO(gene = supercluster_markers,
                      OrgDb = "org.Hs.eg.db",
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10)

dotplot(go_enrich,showCategory=15)+
  ggtitle("Supercluster Functional Enrichment")+
  theme(plot.title = element_text(size=30))

#################################################################################
# Hallmark Enrichment of shared markers
#################################################################################

hallmark_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, human_gene_symbol)
head(hallmark_t2g)


em <- enricher(supercluster_markers, TERM2GENE=hallmark_t2g)

dotplot(em)

#################################################################################
# Common GO pathways in each component cluster of supercluster
#################################################################################

all_pathways <- list()

for(curr_cell_line in cell_lines){
  
  cat(curr_cell_line, "\n")
  clusters <- supercluster_components[[curr_cell_line]]
  
  signatures <- cell_line_cluster_signatures[[curr_cell_line]]
  
  
  for(curr_cluster in clusters){
    
    cat(curr_cluster, "\n")
    
    go_enrich <- enrichGO(gene = signatures[[curr_cluster]],
                          OrgDb = "org.Hs.eg.db",
                          keyType = 'SYMBOL',
                          readable = T,
                          ont = "ALL",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.10)
    
    
    curr_pathways <- go_enrich@result %>% 
      filter(p.adjust < 0.05) %>% 
      arrange(desc(Count)) %>% 
      pull(Description)
    
    curr_pathways <- curr_pathways[1:10]
    
    all_pathways <- append(all_pathways, list(curr_pathways))
    
  }
}

pathway_freq_table <- sort(table(unlist(all_pathways)), decreasing = T)

pathway_freq_table <- pathway_freq_table

pathway_freq_table <- data.frame(cbind(names(pathway_freq_table),pathway_freq_table))

rownames(pathway_freq_table) <- NULL
colnames(pathway_freq_table) <- c("pathway","frequency")


ggplot(pathway_freq_table,aes(x = fct_rev(fct_reorder(pathway, frequency)), y = frequency)) +
  geom_col() +
  geom_col(aes(x=pathway,y=frequency, fill=frequency))+
  xlab("GO Pathways")+
  ylab("Frequency")+
  ggtitle("Frequency of Pathways in Top 10 Pathways of Component Clusters")+
  theme(axis.text.x = element_text(size=15,angle=45,hjust=1,vjust=1),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        title = element_text(size=25))+
  NoLegend()
  



supercluster_markers <- unique(intersect(intersect(A549_temp,K562_temp),MCF7_temp))




