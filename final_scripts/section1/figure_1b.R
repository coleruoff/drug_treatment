args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

source("final_scripts/drug_treatment_functions.R")

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"
# plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_figures/"

###################################################################################

cell_lines <- c("A549","K562","MCF7")

all_data <- list()
for(curr_cell_line in cell_lines){
  # Read in current cell line data
  data <- readRDS(paste0(dataDirectory, "processed_data/sciplex_data/", curr_cell_line, "_processed_filtered.rds"))
  
  # Read in resistance signature AUCell scores
  scores <- readRDS(paste0(dataDirectory, "aucell_score_objects/",curr_cell_line, "_processed_filtered_resistance_signature_aucell_scores.rds"))
  
  # Add AUCell scores to metadata then subset Cluster and scores for boxplot plotting
  data <- AddMetaData(data, metadata = scores, col.name = "curr_geneset")
  
  all_data[[curr_cell_line]] <- data
}


plots <- list()

for(curr_cell_line in cell_lines){
  
  data <- all_data[[curr_cell_line]]
  
  df <- data@meta.data %>% 
    filter(treatment_stage == "post") %>% 
    dplyr::select(Cluster, curr_geneset) %>% 
    group_by(Cluster) %>% 
    mutate("group_median_score" = median(curr_geneset))
  
  # Plot current cell line's boxplots
  plot_title <- curr_cell_line
  
  if(curr_cell_line == "K562"){
    p_label_x <- 1.5
    p_label_y <- .15
  } else {
    p_label_x <- 2
    p_label_y <- .15
  }
  
  p <- ggboxplot(df, x = "Cluster", y = "curr_geneset", fill="group_median_score", outlier.shape = NA,size=.1)+
    stat_compare_means(label.x = p_label_x, label.y=p_label_y, size=2)+
    scale_fill_gradient(low="white", high="orange", name="      Median\nCluster Score")+
    ggtitle(plot_title)+
    xlab("")+
    ylab("")+
    theme(legend.position="right",
          axis.text = element_text(size=6),
          legend.text = element_text(size=6),
          legend.title = element_text(size=6),
          legend.key.height = unit(5,"mm"),
          legend.key.width = unit(3,"mm"),
          title = element_text(size=8),
          axis.line = element_line(linewidth=.2),
          axis.ticks = element_line(linewidth = .2))
  
  # p
  # Add current cell line boxplots to final list of plots
  plots <- append(plots,list(p))
  

}


# Arrange each cell line plot, annotate figure, and plot
figure <- ggarrange(plotlist = plots, nrow=3, ncol=1, common.legend = T, legend=c("right"))

p <- annotate_figure(figure, left = textGrob("Score", rot = 90, vjust = 1, gp = gpar(fontsize=8)),
                     bottom = textGrob("Clusters", gp = gpar(fontsize=8)))



jpeg(paste0(plotDirectory,"figure_1b.jpg"), width=100, height=120, units = "mm", res = 1000)
print(p)
dev.off()





