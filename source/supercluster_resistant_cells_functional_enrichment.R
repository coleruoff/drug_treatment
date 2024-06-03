library(clusterProfiler)
source("source/cole_functions.R")

cell_lines <- c("A549","K562","MCF7")

all_active_signatures <- list()
for(curr_cell_line in cell_lines){
  de_result <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/",curr_cell_line,"_resistant_clusters_de.rds"))
  
  # Intracluster DE Upregulated genes
  intracluster_active_de_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/",curr_cell_line, "intracluster_active_de_signatures.rds"))
  
  RACs <- list(c(4,9,12,13,14,16,18),c(4,5,9,11),c(5,8,12,13,17))
  names(RACs) <- c("A549","K562","MCF7")
  clusters_of_interest <- RACs[[curr_cell_line]]
  
  active_signatures <- list()
  
  for(curr_cluster in clusters_of_interest){
    
    cat(curr_cluster, "\n")
    
    inactive_signature <- de_result %>% 
      filter(cluster == curr_cluster & p_val_adj < .05 & avg_log2FC > .5) %>% 
      arrange(desc(avg_log2FC)) %>% 
      pull(gene)
    
    active_signature <- de_result %>% 
      filter(cluster == paste0(curr_cluster,"_resistant") & p_val_adj < .05 & avg_log2FC > .5) %>% 
      arrange(desc(avg_log2FC)) %>% 
      pull(gene)
    
    #Find shared genes between active and inactive componenents
    shared_genes <- intersect(active_signature, inactive_signature)
    
    #Remove shared genes
    active_signature <- active_signature[!active_signature %in% shared_genes]
    
    #Add genes from DE between active and inactive within this cluster
    active_signature <- unique(append(active_signature, intracluster_active_de_signatures[[paste0(curr_cluster,"_active_signature")]]))
    
    active_signatures <- append(active_signatures, list(active_signature))
    
  }
  
  names(active_signatures) <- paste0(curr_cell_line, "_", clusters_of_interest)
  
  all_active_signatures <- append(all_active_signatures, active_signatures)
}



################################################################################

supercluster1_signatures <- all_active_signatures[c(paste0("A549_",c(4,9)),paste0("K562_",c(4,9,5)),paste0("MCF7_",c(8,12)))]

supercluster1_signature <- find_consensus_geneset(supercluster1_signatures,2)

go_enrich <- enrichGO(gene = supercluster1_signature,
                      OrgDb = "org.Hs.eg.db",
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10)


dotplot(go_enrich) + ggtitle("supercluster 1")



supercluster1_signature[grepl("RPS3",supercluster1_signature)]

##############

supercluster2_signatures <- all_active_signatures[c(paste0("A549_",c(19)),paste0("K562_",c(11)),paste0("MCF7_",c(5)))]

supercluster2_signature <- find_consensus_geneset(supercluster2_signatures, 2)

go_enrich <- enrichGO(gene = supercluster2_signature,
                      OrgDb = "org.Hs.eg.db",
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10)


dotplot(go_enrich) + ggtitle("supercluster 2")

##############

supercluster3_signatures <- all_active_signatures[c(paste0("A549_",c(14)),paste0("MCF7_",c(13)))]

supercluster3_signature <- find_consensus_geneset(supercluster3_signatures, 2)

go_enrich <- enrichGO(gene = supercluster3_signature,
                      OrgDb = "org.Hs.eg.db",
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10)


dotplot(go_enrich) + ggtitle("supercluster 3")

##############

supercluster4_signatures <- all_active_signatures[c(paste0("A549_",c(13,17)),paste0("MCF7_",c(17)))]

supercluster4_signature <- find_consensus_geneset(supercluster4_signatures, 2)

go_enrich <- enrichGO(gene = supercluster4_signature,
                     OrgDb = "org.Hs.eg.db",
                     keyType = 'SYMBOL',
                     readable = T,
                     ont = "ALL",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.10)


dotplot(go_enrich) + ggtitle("supercluster 4")

##############

supercluster5_signatures <- all_active_signatures[c(paste0("A549_",c(12,16)))]

supercluster5_signature <- find_consensus_geneset(supercluster5_signatures, 2)

go_enrich <- enrichGO(gene = supercluster5_signature,
                      OrgDb = "org.Hs.eg.db",
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10)


dotplot(go_enrich) + ggtitle("supercluster 5")
################################################################################

# Supercluster 1
supercluster1_list <- list()
supercluster1_list <- append(supercluster1_list,A549_active_signatures[paste0(c(9),"_active_signature")])
supercluster1_list <- append(supercluster1_list,K562_active_signatures[paste0(c(5),"_active_signature")])
supercluster1_list <- append(supercluster1_list,MCF7_active_signatures[paste0(c(8),"_active_signature")])

supercluster1_signature <- find_consensus_geneset(supercluster1_list,2)
go_enrich <- enrichGO(gene = supercluster1_signature,
                      OrgDb = "org.Hs.eg.db",
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10)


dotplot(go_enrich) + ggtitle("supercluster 1")

# Supercluster 2
supercluster2_list <- list()
supercluster2_list <- append(supercluster2_list,A549_active_signatures[paste0(c(4,9),"_active_signature")])
supercluster2_list <- append(supercluster2_list,K562_active_signatures[paste0(c(5),"_active_signature")])
supercluster2_list <- append(supercluster2_list,MCF7_active_signatures[paste0(c(8,12),"_active_signature")])

supercluster2_signature <- find_consensus_geneset(supercluster2_list,2)
go_enrich <- enrichGO(gene = supercluster2_signature,
                      OrgDb = "org.Hs.eg.db",
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10)


dotplot(go_enrich) + ggtitle("supercluster 2")

# Supercluster 3
supercluster3_list <- list()
supercluster3_list <- append(supercluster3_list,A549_active_signatures[paste0(c(4,9),"_active_signature")])
supercluster3_list <- append(supercluster3_list,K562_active_signatures[paste0(c(4,5,9),"_active_signature")])
supercluster3_list <- append(supercluster3_list,MCF7_active_signatures[paste0(c(8,12),"_active_signature")])

supercluster3_signature <- find_consensus_geneset(supercluster3_list,2)
go_enrich <- enrichGO(gene = supercluster3_signature,
                      OrgDb = "org.Hs.eg.db",
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10)


dotplot(go_enrich) + ggtitle("supercluster 3")



supercluster_signatures <- list(supercluster1_signature,supercluster2_signature,supercluster3_signature)
names(supercluster_signatures) <- c("supercluster1","supercluster2","supercluster3")


# Resistant subpopulation supercluster 1
genes_to_use <- intersect(intersect(A549_active_signatures[[paste0(9,"_active_signature")]],K562_active_signatures[[paste0(5,"_active_signature")]]),MCF7_active_signatures[[paste0(8,"_active_signature")]])


go_enrich <- enrichGO(gene = genes_to_use,
                      OrgDb = "org.Hs.eg.db",
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10)


dotplot(go_enrich) + ggtitle("supercluster 1")


# Resistant subpopulation supercluster 2

genes_to_use <- Reduce(intersect, list(A549_active_signatures[[paste0(9,"_active_signature")]],A549_active_signatures[[paste0(4,"_active_signature")]],K562_active_signatures[[paste0(5,"_active_signature")]],MCF7_active_signatures[[paste0(8,"_active_signature")]],MCF7_active_signatures[[paste0(12,"_active_signature")]]))


go_enrich <- enrichGO(gene = genes_to_use,
                      OrgDb = "org.Hs.eg.db",
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10)


dotplot(go_enrich) + ggtitle("supercluster 2")

# Resistant subpopulation supercluster 3

genes_to_use <- Reduce(intersect, list(A549_active_signatures[[paste0(9,"_active_signature")]],A549_active_signatures[[paste0(4,"_active_signature")]],K562_active_signatures[[paste0(4,"_active_signature")]],K562_active_signatures[[paste0(5,"_active_signature")]],K562_active_signatures[[paste0(9,"_active_signature")]],MCF7_active_signatures[[paste0(8,"_active_signature")]],MCF7_active_signatures[[paste0(12,"_active_signature")]]))


go_enrich <- enrichGO(gene = genes_to_use,
                      OrgDb = "org.Hs.eg.db",
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10)


dotplot(go_enrich) + ggtitle("supercluster 3")


