library(ComplexHeatmap)
library(Seurat)
library(circlize)
library(ggpubr)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/"

geneset_group_matrix <- function(data, ident_to_use, genesets_to_use){
  
  Idents(data) <- ident_to_use
  
  num_groups <- nlevels(data)
  
  groups <- levels(data)
  
  cat(genesets_to_use, "\n")
  
  scores <- readRDS(paste0(dataDirectory, "data/aucell_score_objects/",curr_cell_line, "_processed_filtered_", genesets_to_use, "_aucell_scores.rds"))
  threshold <- readRDS(paste0(dataDirectory, "data/aucell_score_objects/",curr_cell_line, "_processed_filtered_", genesets_to_use, "_aucell_thresholds.rds"))
  
  # scores <- scale(scores)
  
  geneset_names <- colnames(scores)
  
  score_matrix <- matrix(NA, ncol=length(groups), nrow=length(geneset_names))
  pvalue_matrix <- matrix(NA, ncol=length(groups), nrow=length(geneset_names))
  
  for(curr_group in groups){
    
    cat(curr_group, "\n")
    
    curr_cell_names <- colnames(data)[data[[ident_to_use]] == curr_group]
    rest_names <- colnames(data)[data[[ident_to_use]] != curr_group]
    
    for(curr_geneset in geneset_names){
      
      curr_scores <- scores[curr_cell_names, curr_geneset]
      
      rest_scores <- scores[rest_names, curr_geneset]
      
      wilcox_res <- wilcox.test(curr_scores,rest_scores)
      
      i <- which(curr_geneset == geneset_names)
      j <- which(curr_group == groups)
      
      score_matrix[i,j] <- mean(curr_scores)
      pvalue_matrix[i,j] <- wilcox_res$p.value
      
    }
  }
  
  scaled_matrix <- t(scale(t(score_matrix)))
  
  colnames(scaled_matrix) <- groups
  rownames(scaled_matrix) <- geneset_names
  
  colnames(score_matrix) <- groups
  rownames(score_matrix) <- geneset_names
  
  return(list(scaled_matrix, pvalue_matrix,score_matrix))
}


################################################################################

curr_cell_line <- "MCF7"

data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")
clusters_of_interest <- RACs[[curr_cell_line]]

data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, data$Cluster, "Non-RAC"), col.name = "RAC")
#Add metadata for RAC and Cell Group
data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, "1", ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), 2, 0)), col.name = "cell_group")
data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, paste0(data$Cluster, "_1"), ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), paste0(data$Cluster, "_2"), paste0(data$Cluster, "_0"))), col.name = "cell_cluster_group")

genesets_name <- "ecoli_AMR_genesets_orthologs"

ret_mat <- geneset_group_matrix(data, "cell_cluster_group", genesets_name)

#Set all non-significant heatmap cells to 0
heatmap_matrix <- ifelse(ret_mat[[2]]<.05,ret_mat[[1]], 0)

colnames(heatmap_matrix) <- paste0("Cluster ",colnames(ret_mat[[1]]))
rownames(heatmap_matrix) <- rownames(ret_mat[[1]])

heatmap_matrix_up <- heatmap_matrix[!grepl("down",rownames(heatmap_matrix)),]

rownames(heatmap_matrix_up) <- sapply(X = rownames(heatmap_matrix_up), FUN = function(x) {strsplit(x,split = "_")[[1]][1]})


df <- data.frame(cbind(colnames(heatmap_matrix),heatmap_matrix[1,]))

colnames(df) <- c("cluster","value")
df$value <- as.numeric(df$value)

rac_ha <- HeatmapAnnotation(cell_group = c(ifelse(grepl("_0",colnames(heatmap_matrix)),"Non-RAC", ifelse(grepl("_1",colnames(heatmap_matrix)), "RAC Type 1","RAC Type 2"))),
                            col = list(cell_group = c("RAC Type 1" = "red","RAC Type 2" = "orange", "Non-RAC" = "lightblue")), show_annotation_name = F,
                            annotation_legend_param = list(legend_gp=gpar(fontsize=20), grid_height=unit(1,"cm"),grid_width=unit(1,"cm"),
                                                           title="Cluster Type", labels_gp = gpar(fontsize = 14)))

colnames(heatmap_matrix) <- gsub("_0", "", colnames(heatmap_matrix))
colnames(heatmap_matrix) <- gsub("_1", " (Type 1)", colnames(heatmap_matrix))
colnames(heatmap_matrix) <- gsub("_2", " (Type 2)", colnames(heatmap_matrix))


ht <- Heatmap(heatmap_matrix_up, name="Z-Score", cluster_rows = F, cluster_columns = T,
              bottom_annotation = rac_ha, column_title = "", column_title_side = "bottom",
              row_title = "Antimicrobial Drugs\n", row_title_side = "left", row_title_gp = gpar(fontsize=25),
              row_names_side = "left", column_names_rot = 45, 
              row_names_gp = gpar(fontsize=20),
              column_names_gp = gpar(fontsize=20),
              heatmap_legend_param = list(legend_gp = gpar(fontsize = 12),legend_height = unit(3, "cm"), grid_width=unit(1,"cm"),
                                          labels_gp = gpar(fontsize = 14)))

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/ecoli_figures/",curr_cell_line, "_ecoli_heatmap.png"), width = 1200,height=800)

draw(ht, column_title = paste0("                  E. Coli Antimicrobial Resistance Orthologs Mean Cluster AUCell Score (", curr_cell_line, ")\n"), 
     column_title_gp = gpar(fontsize = 25, fontface = "bold"),  padding = unit(c(2, 2, 10, 2), "mm"),
     heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend=T)

dev.off()

#################################################################################

# RACs BOXPLOTS

type1_cell_names <- colnames(data)[data$cell_group == "1"]
type2_cell_names <- colnames(data)[data$cell_group == "2"]
non_rac_cell_names <- colnames(data)[data$cell_group == "0"]

type1_cell_values <- scores[type1_cell_names,]
type2_cell_values <- scores[type2_cell_names,]
non_rac_values <- scores[non_rac_cell_names,]

boxplot_df <- as.data.frame(cbind(scores,ifelse(rownames(scores) %in% non_rac_cell_names,"Non-RAC", ifelse(rownames(scores) %in% type1_cell_names, "RAC Type 1","RAC Type 2"))))

colnames(boxplot_df) <- c("value","group")
boxplot_df$value <- as.numeric(boxplot_df$value)
boxplot_df$group <- factor(boxplot_df$group, levels = c("RAC Type 1","RAC Type 2","Non-RAC"))

my_comparisons <- list( c("Non-RAC", "RAC Type 2"), c("Non-RAC", "RAC Type 1"),c("RAC Type 1","RAC Type 2"))

p <- ggboxplot(boxplot_df, x = "group", y = "value", fill="group",short.panel.labs = FALSE)

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/ecoli_figures/",curr_cell_line, "_ecoli_boxplots_racs.png"), width = 1200,height=800)
# Use only p.format as label. Remove method name.
p + stat_compare_means(comparisons = my_comparisons,label = "p.format", method = "wilcox", size=8)+
  ggtitle(paste0("E. Coli Antimicrobial Resistance Orthologs AUCell Scores (",curr_cell_line,")"))+
  xlab("")+
  ylab("AUCell Score")+
  scale_fill_manual(values=c("red", "orange","lightblue"),name = "Cell Groups")+
  theme(legend.position="right",
        title = element_text(size=28, face="bold"),
        axis.text = element_text(size=30),
        legend.text = element_text(size=20))


dev.off()

#################################################################################


#################################################################################

# Functional Enrichment of Overlapping Genes

RAC_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/active_vs_inactive_signatures.rds")

ecoli_orthologs <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ecoli_AMR_genesets_orthologs.rds")

RAC_overlap_with_ecoli_signatures <- list()
for(i in 1:length(RAC_signatures)){
  
  temp <- c()
  for(j in 1:length(ecoli_orthologs)){
    
    temp <- append(temp, intersect(RAC_signatures[[i]],ecoli_orthologs[[j]]))
  }
  RAC_overlap_with_ecoli_signatures <- append(RAC_overlap_with_ecoli_signatures, list(unique(temp)))
  
}

names(RAC_overlap_with_ecoli_signatures) <- c("A549","K562","MCF7")

all_dotplots <- genesets_characterization(RAC_overlap_with_ecoli_signatures)


curr_width <- 1000
curr_height <- 1000

# Plotting for shared genes between RACs signatures and e coli orthologs
png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/ecoli_figures/ecoli_hallmarks.png"), width = curr_width,height = curr_height)
hallmarks_plt <- egg::ggarrange(plots = all_dotplots$hallmarks_dotplots, ncol = 2)
main_title <- paste0("Cancer Hallmarks Enrichment of Shared Genes Between RAC Active Cells Signature and ecoli Stress Orthologs")
hallmarks_plt <- annotate_figure(hallmarks_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(hallmarks_plt)
dev.off()

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/ecoli_figures/ecoli_mps.png"), width = curr_width,height = curr_height-500)
mps_plt <- egg::ggarrange(plots = all_dotplots$mps_dotplots, ncol = 2)
main_title <- paste0("ITH Meta-programs Enrichment of Shared Genes Between RAC Active Cells Signature and ecoli Stress Orthologs")
mps_plt <- annotate_figure(mps_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(mps_plt)
dev.off()

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/ecoli_figures/ecoli_functional_enrichment.png"), width = curr_width,height = curr_height)
go_plt <- egg::ggarrange(plots = all_dotplots$go_dotplots, ncol = 2)
main_title <- paste0("GO Functional Enrichment of Shared Genes Between RAC Active Cells Signature and ecoli Stress Orthologs")
go_plt <- annotate_figure(go_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(go_plt)
dev.off()
