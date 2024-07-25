args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

library(ggpubr)
library(Seurat)
library(tidyverse)
library(grid)
library(lmtest) 

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"
plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/"
###################################################################################

cell_lines <- c("A549","K562","MCF7")

plots <- list()

for(curr_cell_line in cell_lines){
  
  # Read in current cell line data
  data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))
  
  # Read in resistance signature AUCell scores
  scores <- readRDS(paste0(dataDirectory, "aucell_score_objects/",curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
  
  # Add AUCell scores to metadata then subset Cluster and scores for boxplot plotting
  data <- AddMetaData(data, metadata = scores, col.name = "curr_geneset")
  
  df <- data@meta.data %>% 
    filter(treatment_stage == "post") %>% 
    dplyr::select(Cluster, curr_geneset) %>% 
    group_by(Cluster) %>% 
    mutate("group_median_score" = median(curr_geneset))
  
  # Plot current cell line's boxplots
  plot_title <- curr_cell_line
  p <- ggboxplot(df, x = "Cluster", y = "curr_geneset", fill="group_median_score",outlier.shape = 20,size=.2)+
    scale_fill_gradient(low="white", high="orange", name="      Median\nCluster Score")+
    ggtitle(plot_title)+
    xlab("")+
    ylab("")+
    theme(legend.position="right",
          axis.text = element_text(size=8),
          legend.text = element_text(size=10),
          legend.title = element_text(size=10),
          legend.key.height = unit(5,"mm"),
          legend.key.width = unit(5,"mm"))
  
  # Add current cell line boxplots to final list of plots
  plots <- append(plots,list(p))
  
  # Comparing real data to null model with likelihood-ratio test
  model <- lm(curr_geneset ~ Cluster, data = df)
  
  summary(model)
  
  
  null_model <- lm(curr_geneset ~ 1, data = df)
  
  summary(null_model)
  
  
  likelihood_ratio_test <- lrtest(null_model, model) 
  cat(likelihood_ratio_test$`Pr(>Chisq)`[2], "\n")
}

# Arrange each cell line plot, annotate figure, and plot
figure <- ggarrange(plotlist = plots, nrow=1,ncol=3, common.legend = T,legend=c("right"))

p <- annotate_figure(figure, left = textGrob("Score", rot = 90, vjust = 1, gp = gpar(fontsize=20, fontface="bold")),
                bottom = textGrob("Clusters", gp = gpar(fontsize=20, fontface="bold")))


tiff(paste0(plotDirectory,"figure_1b.tiff"),width=400, height = 100, units = "mm", res = 1000)

print(p)

dev.off()

# png(paste0(plotDirectory,"figure_1b.png"),
#     width=20, height=20, units = "in", res = 300)
# 
# print(p)
# 
# dev.off()



# Comparing real data to null model with likelihood-ratio test
model <- lm(curr_geneset ~ Cluster, data = df)

summary(model)


null_model <- lm(curr_geneset ~ 1, data = df)

summary(null_model)


likelihood_ratio_test <- lrtest(null_model, model) 
likelihood_ratio_test$`Pr(>Chisq)`[2]


