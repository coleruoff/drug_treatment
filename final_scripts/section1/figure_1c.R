args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

source("final_scripts/drug_treatment_functions.R")

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"
# plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/"

################################################################################

cell_lines <- c("A549","K562","MCF7")

all_racs <- list()
plot_data <- list()

for(curr_cell_line in cell_lines){
  cat(curr_cell_line,"\n")
  
  data <- readRDS(paste0(dataDirectory,"processed_data/sciplex_data/", curr_cell_line, "_processed_filtered.rds"))
  
  scores <- readRDS(paste0(dataDirectory, "aucell_score_objects/", curr_cell_line, "_processed_filtered_resistance_signature_aucell_scores.rds"))
  threshold <- readRDS(paste0(dataDirectory, "aucell_score_objects/", curr_cell_line, "_processed_filtered_resistance_signature_aucell_thresholds.rds"))
  
  data <- AddMetaData(object = data, 
                      metadata = ifelse(scores[,1] > threshold$threshold, 1, 0), 
                      col.name = 'active')
  
  data <- data[,data$treatment_stage=="post"]
  
  all_clusters <- as.numeric(levels(data))
  
  ORs <- c()
  pvals <- c()
  final_df <- list()
  
  for(i in all_clusters){
    
    in_cluster <- ifelse(data$Cluster == i, 1, 0)
    
    df <- as.data.frame(cbind(data$active,in_cluster))
    
    df$in_cluster <- factor(df$in_cluster, levels=c(1,0))
    df$V1 <- factor(df$V1, levels=c(1,0))
    
    colnames(df) <- c("active", "in_cluster")
    
    contin_table <- table(df)
    
    res <- fisher.test(contin_table)
    
    ORs <- append(ORs, res$estimate)
    pvals <- append(pvals, res$p.value)
    
  }
  
  # Plot OR bar plots
  plot_df <- as.data.frame(cbind(all_clusters,ORs,p.adjust(pvals)))
  colnames(plot_df) <- c("cluster","or","pval")
  
  plot_df <- plot_df %>% 
    mutate(signif = ifelse(pval < 0.05 & or > 1, "*","")) %>% 
    mutate(rac = ifelse(or > 1 & pval < 0.05, "RAC","Non-RAC"))
  
  plot_df$cluster <- factor(as.character(plot_df$cluster), levels=all_clusters)
  
  
  plot_data[[curr_cell_line]] <- plot_df
  
  
  curr_racs <- plot_df %>% 
    dplyr::filter(or > 1 & pval < 0.05) %>% 
    pull(cluster) %>% 
    as.numeric()
  
  all_racs <- append(all_racs, list(curr_racs))
}

names(all_racs) <- cell_lines

saveRDS(all_racs, paste0(dataDirectory, "processed_data/all_racs.rds"))

plots <- list()
for(curr_cell_line in cell_lines){
  plot_df <- plot_data[[curr_cell_line]]
  
  plot_title <- curr_cell_line
  
  p <- ggplot(plot_df)+
    geom_col(aes(x=cluster,y=or, fill=rac))+
    geom_text(aes(x=cluster, y=or,label=signif),size=3)+
    scale_fill_manual(name="Cluster Type",values=c("lightblue","orange"))+
    geom_hline(yintercept=1, linetype="dashed", color = "red", linewidth=.1)+
    ggtitle(plot_title)+
    xlab("")+
    ylab("")+
    theme_classic()+
    theme(legend.position="right",
          axis.text = element_text(size=6, color="black"),
          legend.text = element_text(size=6),
          legend.title = element_text(size=6),
          legend.key.height = unit(3,"mm"),
          legend.key.width = unit(3,"mm"),
          title = element_text(size=8),
          axis.line = element_line(linewidth=.2),
          axis.ticks = element_line(linewidth = .2))
  
  
  plots <- append(plots,list(p))
}

figure <- ggarrange(plotlist = plots, nrow=3, common.legend = T, legend=c("right"))

p <- annotate_figure(figure, left = text_grob("Odds Ratio", rot = 90, vjust = 1, size=8),
                     bottom = text_grob("Clusters", size=8))

p
jpeg(paste0(plotDirectory,"figure_1c.jpg"), width=100, height=120, units = "mm", res = 600)
print(p)
dev.off()






