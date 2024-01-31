library(Seurat)
source("source/cole_functions.R")

watermelon_data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/watermelon_data/pc9_processed.rds"))

################################################################################
# Calculate FC for each time point
################################################################################

Idents(watermelon_data) <- watermelon_data$time_point

time_points <- levels(Idents(watermelon_data))

time_points_fc <- list()
for(curr_time_point in time_points){
  
  curr_fc_result <- FoldChange(watermelon_data, ident.1 = curr_time_point)
  
  time_points_fc <- append(time_points_fc, list(curr_fc_result))
  
}


names(time_points_fc) <- paste0("pc9_day",time_points,"_fc")

################################################################################
# Run GSEA of raj resistance signature for each time point
################################################################################

#Create ranked genelists for each time point
genesets_with_ranks <- list()
for(curr_fc_res in time_points_fc){
  
  ranks <- curr_fc_res$avg_log2FC
  
  names(ranks) <- rownames(curr_fc_res)
  
  ranks <- sort(ranks, decreasing = T)
  
  genesets_with_ranks <- append(genesets_with_ranks, list(ranks))
  
}

names(genesets_with_ranks) <- paste0("pc9_day",time_points,"_fc")

#Read in resistance signature
curr_cell_line <- "MCF7"
signature_length <- 200

trapnell_cluster_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_cluster_signatures_",signature_length,".rds"))


genesets2 <- trapnell_cluster_signatures

# Run GSEA

result <- create_GSEA_matrix(genesets_with_ranks, genesets2)


A549_active <- c(4,9,13)

K562_active <- c(5,11)

MCF7_active <-  c(5,8,13,16,17)

active <- eval(parse(text=paste0(curr_cell_line,"_active")))


gsea_matrix <- result[[1]]

inactive <- (1:nrow(gsea_matrix))[-active]

colnames(gsea_matrix) <- paste0("Day ", time_points)
rownames(gsea_matrix) <- names(genesets2)
plot_title <- paste0(curr_cell_line," Cluster Markers Enrichment Along Resistance Development in PC9")


col_fun <-colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
Heatmap(gsea_matrix, cluster_columns = F, cluster_rows=F, column_title = plot_title, name="NES", 
        column_names_rot = 45, col = col_fun)



wilcox.test(gsea_matrix[,4][active],gsea_matrix[,4][inactive])
boxplot(gsea_matrix[,4][active],gsea_matrix[,4][inactive])

        



