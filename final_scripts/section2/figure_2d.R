args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

source("final_scripts/drug_treatment_functions.R")

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"
# plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_figures/"

################################################################################
# plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/"

################################################################################

cell_lines <- c("A549","K562","MCF7")

supercluster_components <- readRDS(paste0(dataDirectory, "processed_data/supercluster_components.rds"))

all_data <- readRDS(paste0(dataDirectory, "processed_data/sciplex_data/all_cell_lines_data.rds"))


drug_classes <- c("Antioxidant","Apoptotic regulation","Cell cycle regulation",
                  "DNA damage & DNA repair","Epigenetic regulation","Focal adhesion signaling",
                  "HIF signaling","JAK/STAT signaling","Metabolic regulation",
                  "Neuronal signaling","Nuclear receptor signaling","PKC signaling",
                  "Protein folding & Protein degradation","TGF/BMP signaling",
                  "Tyrosine kinase signaling")



post_data <- list()
for(curr_cell_line in cell_lines){
  data <- all_data[[curr_cell_line]]
  
  post_data[[curr_cell_line]] <- data[,data$treatment_stage == "post"]
  
}

heatmap <- matrix(NA, ncol=3,nrow=length(drug_classes))
colnames(heatmap) <- paste0("Supercluster ",1:3)
rownames(heatmap) <- drug_classes

for(curr_sc in 1:length(supercluster_components)){
  
  cat(curr_sc, "\n")
  
  for(curr_drug_class in drug_classes){
  
    a <- 0
    b <- 0
    c <- 0
    d <- 0
    
    for(curr_cell_line in cell_lines){
      
      data <- post_data[[curr_cell_line]]
      
      clusters <- supercluster_components[[curr_sc]][[curr_cell_line]]
      
      a <- a + sum(data$Cluster %in% clusters & data$pathway_level_1 == curr_drug_class)
      b <- b + sum(!(data$Cluster %in% clusters) & data$pathway_level_1 == curr_drug_class)
      c <- c + sum(data$Cluster %in% clusters & !(data$pathway_level_1 == curr_drug_class))
      d <- d + sum(!(data$Cluster %in% clusters & data$pathway_level_1 == curr_drug_class))
      
    }
    
    mat <- matrix(c(a,c,b,d),nrow=2)
    
    res <- fisher.test(mat)
    
    # heatmap[which(curr_drug_class == drug_classes),curr_sc] <- res$estimate
    heatmap[which(curr_drug_class == drug_classes),curr_sc] <- (a*d)/(b*c)
    
  }
}


col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

ht <- Heatmap(log(heatmap), cluster_rows = F,cluster_columns = F, name="log(OR)",
        rect_gp = gpar(col = "white", lwd = 1), column_names_gp = gpar(fontsize=4), column_names_rot = 45,
        row_names_gp = gpar(fontsize=4), col = col_fun,
        heatmap_legend_param = list(title_gp = gpar(fontsize = 4),legend_height = unit(2, "mm"), grid_width=unit(2,"mm"),
                            labels_gp = gpar(fontsize = 4),legend_gp = gpar(lwd = .5)))





jpeg(paste0(plotDirectory,"figure_2d.jpg"), width=60, height = 60, units = "mm", res = 1000)
draw(ht, heatmap_legend_side = "left")
dev.off()






