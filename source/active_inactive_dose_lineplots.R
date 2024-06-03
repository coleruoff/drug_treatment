library(tidyverse)

curr_cell_line <- "A549"

cell_lines <- c("A549","K562","MCF7")

plots <- list()
for(curr_cell_line in cell_lines){
  
  cat(curr_cell_line,"\n")
  
  data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))
  
  scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
  threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
  
  active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
  active_cell_names <- active_cell_names[active_cell_names %in% colnames(data)]
  
  all_doses <- sort(unique(data$dose))
  all_clusters <- sort(unique(data$Cluster))
  
  # total_num_active <- length(active_cell_names)
  # total_num_inactive <- ncol(data)-total_num_active
  
  ORs <- c()
  pvals <- c()
  
  total_df <- data.frame(matrix(NA,nrow=0,ncol=7))
  
  for(i in all_clusters){
    
    for(j in all_doses){
      curr_dose_names <- colnames(data)[data$dose == j]
      
      total_num_active <- sum(active_cell_names %in% curr_dose_names)
      total_num_inactive <- length(curr_dose_names)-total_num_active
      
      curr_cluster_names <- colnames(data)[data$Cluster == i & data$dose == j]
      
      curr_num_active <- sum(curr_cluster_names %in% active_cell_names)
      curr_num_inactive <- length(curr_cluster_names) - curr_num_active
      
      
      contin_table <- matrix(c(curr_num_active,total_num_active,curr_num_inactive,total_num_inactive), ncol=2)
      
      res <- fisher.test(contin_table)
      
      ORs <- append(ORs, res$estimate)
      pvals <- append(pvals, res$p.value)
      
      active_percent <- curr_num_active/(curr_num_inactive+curr_num_active)
      inactive_percent <- 1-active_percent
      
      cluster_size_percent <- length(curr_cluster_names)/length(curr_dose_names)
      
      prolif_index <- log(mean(data$proliferation_index[colnames(data) %in% curr_cluster_names])/mean(data$proliferation_index))
      
      
      total_df <- rbind(total_df,c(i,j,res$estimate,active_percent,inactive_percent,cluster_size_percent,prolif_index))
    }
  }
  
  colnames(total_df) <- c("cluster","dose","OR","percent","inactive_percent","cluster_size_percent","prolif_index")
  total_df$OR <- as.numeric(total_df$OR)
  total_df$percent <- as.numeric(total_df$percent)
  total_df$inactive_percent <- as.numeric(total_df$inactive_percent)
  total_df$cluster_size_percent <- as.numeric(total_df$cluster_size_percent)
  total_df$prolif_index <- as.numeric(total_df$prolif_index)
  
  total_df$cluster <- paste0("Cluster ", total_df$cluster)
  all_clusters <- paste0("Cluster ", all_clusters)
  
  
  
  # png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/active_inactive_or_figures/",curr_cell_line,"_dose_lineplots.png"),
  #     width=1000, height = 1000)
  
  # plot_title <- paste0("Resistant Active/Inactive Odds Ratio (", curr_cell_line, ")")
  plot_title <- curr_cell_line
  p <- ggplot(total_df,aes(x=dose, y=OR, group=1))+
    geom_line()+
    geom_hline(yintercept=1, linetype="dashed", color = "red", size=1)+
    geom_point()+
    xlab("")+
    ylab("")+
    ggtitle(plot_title)+
    facet_wrap(~factor(cluster, levels=paste0("Cluster ", 1:19)))+
    theme(axis.text.x = element_text(angle=45, hjust=1),
          legend.position="right",
          title = element_text(size=28, face="bold"),
          axis.text = element_text(size=18),
          legend.text = element_text(size=14),
          legend.title = element_text(size=12),
          strip.text = element_text(size=20))
  
  plots <- append(plots,list(p))
}

jpeg(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/active_inactive_or_figures/figure_1d.jpeg"),
     width=3000, height=1000)

figure <- ggarrange(plotlist = plots, ncol=3, nrow=1,common.legend = F)

p <- annotate_figure(figure, left = text_grob("Odds Ratio", rot = 90, vjust = 1, size=35, face="bold"),
                bottom = text_grob("Dose", size=35, face="bold"),
                top=text_grob("Resistant Active/Inactive Odds Ratio", size=40, face="bold"))



print(p)


dev.off()

























# 
# 
# plot_title <- paste0(curr_cell_line, " Resistant Active Percent")
# ggplot(total_df,aes(x=dose, y=percent,group=1))+
#   geom_line()+
#   geom_point()+
#   xlab("Dose")+
#   ylab("Percent")+
#   ggtitle(plot_title)+
#   theme(axis.text.x = element_text(angle=45, hjust=1))+
#   facet_wrap(~factor(cluster, levels=all_clusters))+
#   geom_hline(yintercept=mean(total_df$percent[!is.na(total_df$percent)]), linetype="dashed", color = "red", size=1)
# 
# 
# 
# plot_title <- paste0(curr_cell_line, " Resistant Inactive Percent")
# ggplot(total_df,aes(x=dose, y=inactive_percent,group=1))+
#   geom_line()+
#   geom_point()+
#   xlab("Dose")+
#   ylab("Percent")+
#   ggtitle(plot_title)+
#   theme(axis.text.x = element_text(angle=45, hjust=1))+
#   facet_wrap(~factor(cluster, levels=all_clusters))+
#   geom_hline(yintercept=mean(total_df$inactive_percent[!is.na(total_df$inactive_percent)]), linetype="dashed", color = "red", size=1)
# 
# 
# plot_title <- paste0(curr_cell_line, " Cluster Size Percent")
# ggplot(total_df,aes(x=dose, y=cluster_size_percent,group=1))+
#   geom_line()+
#   geom_point()+
#   xlab("Dose")+
#   ylab("Percent")+
#   ggtitle(plot_title)+
#   theme(axis.text.x = element_text(angle=45, hjust=1))+
#   facet_wrap(~factor(cluster, levels=all_clusters))
# 
# 
# plot_title <- paste0(curr_cell_line, " Cluster Proliferation Index Relative to Total")
# ggplot(total_df,aes(x=dose, y=prolif_index,group=1))+
#   geom_line()+
#   geom_point()+
#   xlab("Dose")+
#   ylab("log(Cluster Proliferation Index/Total Proliferation Index)")+
#   ggtitle(plot_title)+
#   theme(axis.text.x = element_text(angle=45, hjust=1))+
#   facet_wrap(~factor(cluster, levels=all_clusters))+
#   geom_hline(yintercept=0, linetype="dashed", color = "red", size=1)
# 
# 
# 
# 
# 
# 
# 
# df <- data.frame(cbind(paste0("Dose ", all_clusters),ORs))
# colnames(df) <- c("cluster","or")
# df$cluster <- factor(df$cluster, levels = paste0("Dose ", all_clusters))
# df$or <- as.numeric(df$or)
# 
# df$color <- ifelse(df$or > 1, "RAC","Non-RAC")
# 
# 
# cor_val <- cor(1:5,df$or, method = "spearman")
# plot_title <- paste0(curr_cell_line, " Resistant Active/Inactive Odds Ratio (Spearman: ", cor_val, ")")
# 
# 
# ggplot(df,aes(x=cluster, y=or,group=1))+
#   geom_line()+
#   geom_hline(yintercept=1, linetype="dashed", color = "red", size=1)+
#   geom_point()+
#   xlab("Dose")+
#   ylab("Odds Ratio")+
#   ggtitle(plot_title)+
#   theme(axis.text.x = element_text(angle=45, hjust=1))
# 
# all_clusters[ORs > 1.5]
# 
# 
# 
# 
