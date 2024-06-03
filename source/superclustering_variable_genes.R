library(Seurat)
library(ComplexHeatmap)
library(tidyverse)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

#Set cell line
curr_cell_line <- "A549"
#Read in cell line data
A549.data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

#Set cell line
curr_cell_line <- "K562"
#Read in cell line data
K562.data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

#Set cell line
curr_cell_line <- "MCF7"
#Read in cell line data
MCF7.data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))


#################################################################################
# Cluster based on mean expression of variable genes
#################################################################################

A549.data <- FindVariableFeatures(A549.data)
A549_variable_genes <- VariableFeatures(A549.data)

K562.data <- FindVariableFeatures(K562.data)
K562_variable_genes <- VariableFeatures(K562.data)

MCF7.data <- FindVariableFeatures(MCF7.data)
MCF7_variable_genes <- VariableFeatures(MCF7.data)

all_variable_genes <- unique(c(A549_variable_genes, K562_variable_genes, MCF7_variable_genes))

all_variable_genes <- all_variable_genes[all_variable_genes %in% rownames(A549.data)]
all_variable_genes <- all_variable_genes[all_variable_genes %in% rownames(K562.data)]
all_variable_genes <- all_variable_genes[all_variable_genes %in% rownames(MCF7.data)]

# Remove cell line specific markers
# cell_line_de <- readRDS("processed_data/cell_line_de.rds")
# 
# A549_cell_line_markers <- cell_line_de %>% 
#   filter(cluster == "A549" & p_val_adj < 0.05) %>% 
#   arrange(desc(avg_log2FC)) %>% 
#   pull(gene) %>% 
#   head(200)
# 
# K562_cell_line_markers <- cell_line_de %>% 
#   filter(cluster == "K562" & p_val_adj < 0.05) %>% 
#   arrange(desc(avg_log2FC)) %>% 
#   pull(gene) %>% 
#   head(200)
# 
# 
# MCF7_cell_line_markers <- cell_line_de %>% 
#   filter(cluster == "MCF7" & p_val_adj < 0.05) %>% 
#   arrange(desc(avg_log2FC)) %>% 
#   pull(gene) %>% 
#   head(200)
# 
# all_variable_genes <- all_variable_genes[!all_variable_genes %in% A549_cell_line_markers]
# 
# all_variable_genes <- all_variable_genes[!all_variable_genes %in% K562_cell_line_markers]
# 
# all_variable_genes <- all_variable_genes[!all_variable_genes %in% MCF7_cell_line_markers]

cell_lines <- c("A549", "K562", "MCF7")

heatmap <- matrix(NA, ncol=49, nrow=length(all_variable_genes))
colnames(heatmap) <- rep("", 49)

j <- 1
for(cell_line in cell_lines){
  
  cat(cell_line,"\n")
  data <- eval(parse(text=paste0(cell_line, ".data")))
  
  #total_mean <- rowMeans(data@assays$originalexp@data[all_variable_genes,])
  
  clusters <- 1:nlevels(data$Cluster)
  
  for(curr_cluster in clusters){
    cat(curr_cluster,"\n")
    
    cluster_mean <- rowMeans(data@assays$RNA@data[all_variable_genes, data$Cluster == curr_cluster])
    
    total_mean <- rowMeans(data@assays$RNA@data[all_variable_genes, data$Cluster != curr_cluster])
    
    #heatmap[,j] <- cluster_mean-total_mean
    heatmap[,j] <- total_mean-cluster_mean
    colnames(heatmap)[j] <- paste0(cell_line,"_",curr_cluster)
    
    j <- j + 1
  }
}

active_clusters <- c("A549_4","A549_9","A549_13", "K562_5","K562_11","MCF7_5","MCF7_8","MCF7_13","MCF7_16","MCF7_17")

emergent_cluster <- c("A549_14","A549_15","A549_16", "A549_17","A549_18","A549_19","K562_9","MCF7_13","MCF7_15","MCF7_18")

clusters_of_interest <- emergent_cluster

cor_heatmap <- cor(heatmap, method = "spearman")

row_idx <- which(rownames(cor_heatmap) %in% clusters_of_interest)
fontsizes <- rep(10, nrow(cor_heatmap))
fontsizes[row_idx] <- 18
fontcolors <- rep('black', nrow(cor_heatmap))
fontcolors[row_idx] <- 'red'
fontfaces <- rep('plain',nrow(cor_heatmap))
fontfaces[row_idx] <- 'bold'

rowAnno <- rowAnnotation(rows = anno_text(rownames(cor_heatmap), gp = gpar(fontsize = fontsizes, fontface = fontfaces, col = fontcolors)))


RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")
clusters_of_interest <- paste0("A549_",RACs[["A549"]])
clusters_of_interest <- append(clusters_of_interest, paste0("K562_",RACs[["K562"]]))
clusters_of_interest <- append(clusters_of_interest, paste0("MCF7_",RACs[["MCF7"]]))

rac_ha <- rowAnnotation(RAC = c(ifelse(colnames(cor_heatmap) %in% clusters_of_interest,"RAC","Non-RAC")),
                            col = list(RAC = c("RAC" = "darkgreen", "Non-RAC" = "lightblue")))


ht <- Heatmap(cor_heatmap, name="Spearman Correlation", cluster_rows = T, cluster_columns = T,
              right_annotation = rac_ha)


# file_name <- "documents/figures/supercluster_heatmap.png"
# png(file_name, width= 1200, height= 1200, units="px")

draw(ht, column_title="Correlations of Cluster Differential Mean Expression of Most Variable Genes", column_title_gp = gpar(fontsize = 26))

dev.off()











#################################################################################
# Cluster based on mean expression of union of emergent cluster signatures
#################################################################################

A549_variable_genes <- VariableFeatures(A549.data)

K562_variable_genes <- VariableFeatures(K562.data)

MCF7_variable_genes <- VariableFeatures(MCF7.data)

all_variable_genes <- unique(c(A549_variable_genes, K562_variable_genes, MCF7_variable_genes))

# Remove cell line specific markers
cell_line_de <- readRDS("processed_data/cell_line_de.rds")

A549_cell_line_markers <- cell_line_de %>% 
  filter(cluster == "A549" & p_val_adj < 0.05) %>% 
  arrange(desc(avg_log2FC)) %>% 
  pull(gene) %>% 
  head(200)

K562_cell_line_markers <- cell_line_de %>% 
  filter(cluster == "K562" & p_val_adj < 0.05) %>% 
  arrange(desc(avg_log2FC)) %>% 
  pull(gene) %>% 
  head(200)


MCF7_cell_line_markers <- cell_line_de %>% 
  filter(cluster == "MCF7" & p_val_adj < 0.05) %>% 
  arrange(desc(avg_log2FC)) %>% 
  pull(gene) %>% 
  head(200)

all_variable_genes <- all_variable_genes[!all_variable_genes %in% A549_cell_line_markers]

all_variable_genes <- all_variable_genes[!all_variable_genes %in% K562_cell_line_markers]

all_variable_genes <- all_variable_genes[!all_variable_genes %in% MCF7_cell_line_markers]

cell_lines <- c("A549", "K562", "MCF7")

heatmap <- matrix(NA, ncol=49, nrow=length(all_variable_genes))
colnames(heatmap) <- rep("", 49)

j <- 1
for(cell_line in cell_lines){
  
  cat(cell_line,"\n")
  data <- eval(parse(text=paste0(cell_line, ".data")))
  
  clusters <- 1:nlevels(data$Cluster)
  
  for(curr_cluster in clusters){
    cat(curr_cluster,"\n")
    
    cluster_mean <- rowMeans(data@assays$originalexp@data[all_variable_genes, data$Cluster == curr_cluster])
    
    total_mean <- rowMeans(data@assays$originalexp@data[all_variable_genes, data$Cluster != curr_cluster])
    
    heatmap[,j] <- cluster_mean-total_mean
    # heatmap[,j] <- total_mean-cluster_mean
    colnames(heatmap)[j] <- paste0(cell_line,"_",curr_cluster)
    
    j <- j + 1
  }
}


emergent_clusters <- c("A549_14","A549_15","A549_16","A549_17","A549_18","A549_19", "K562_9","K562_11","MCF7_13","MCF7_15","MCF7_18")

cor_heatmap <- cor(heatmap, method = "spearman")

row_idx <- which(rownames(cor_heatmap) %in% emergent_clusters)
fontsizes <- rep(10, nrow(cor_heatmap))
fontsizes[row_idx] <- 18
fontcolors <- rep('black', nrow(cor_heatmap))
fontcolors[row_idx] <- 'red'
fontfaces <- rep('plain',nrow(cor_heatmap))
fontfaces[row_idx] <- 'bold'

rowAnno <- rowAnnotation(rows = anno_text(rownames(cor_heatmap), gp = gpar(fontsize = fontsizes, fontface = fontfaces, col = fontcolors)))

ht <- Heatmap(cor_heatmap, name="Spearman\nCorrelation", cluster_rows = T, cluster_columns = T,
              right_annotation = rowAnno, show_row_names = F, heatmap_legend_param = list(legend_gp = gpar(fontsize = 12), labels_gp = gpar(fontsize=15),grid_width = unit(1, "cm")))


# file_name <- "documents/figures/supercluster_heatmap.png"
# png(file_name, width= 1200, height= 1200, units="px")

draw(ht, column_title="Correlations of Cluster Differential Mean Expression of Most Variable Genes", column_title_gp = gpar(fontsize = 40))

dev.off()

clustered_order <- colnames(cor_heatmap)[row_order(ht)]

#Select range for supercluster of interest
temp <- sort(clustered_order[which(clustered_order == "MCF7_16"):which(clustered_order == "A549_14")])


A549_clusters <- as.numeric(gsub("A549_","", temp[grepl("A549",temp)]))
K562_clusters <- as.numeric(gsub("K562_","", temp[grepl("K562",temp)]))
MCF7_clusters <- as.numeric(gsub("MCF7_","", temp[grepl("MCF7",temp)]))



A549_temp <- unlist(A549_cluster_signatures[A549_clusters])

K562_temp <- unlist(K562_cluster_signatures[K562_clusters])

MCF7_temp <- unlist(MCF7_cluster_signatures[MCF7_clusters])



supercluster_markers <- unique(intersect(intersect(A549_temp,K562_temp),MCF7_temp))



supercluster_markers <- supercluster_markers[!grepl("^MT",supercluster_markers)]




go_enrich <- enrichGO(gene = supercluster_markers,
                      OrgDb = "org.Hs.eg.db",
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10,
)

dotplot(go_enrich,showCategory=15)+
  ggtitle("Supercluster Functional Enrichment")+
  theme(plot.title = element_text(size=30))




