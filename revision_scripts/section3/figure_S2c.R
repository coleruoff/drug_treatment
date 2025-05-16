args = commandArgs(trailingOnly=TRUE)
# dataDirectory <- paste0(args[1],"final_data/")
setwd(args[1])

source("revision_scripts/drug_treatment_functions.R")
set.seed(42)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"
plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_figures/"

################################################################################
# data <- readRDS(paste0(dataDirectory, "breast_premalignancy_data/datasets_for_boxplot_fisheroverlap.Rds"))
# data <- readRDS(paste0(dataDirectory, "breast_premalignancy_data/datasets_for_boxplot_newSignatures.Rds"))
data <- readRDS(paste0(dataDirectory, "breast_premalignancy_data/datasets_for_boxplot_fisheroverlap_newSignatures.Rds"))

plot_data <- data[[3]]

plot_data$log_odds[which(plot_data$log_odds == -Inf)] <- -5

plot_data$variable <- factor(plot_data$variable, levels = c("MCF7","supercluster1_signature", "supercluster2_signature", "supercluster3_signature"))
# 
# facet.labs <- c("MCF7 Global RAC Signature","Supercluster 1 Signature", "Supercluster 2 Signature", "Supercluster 3 Signature")
# 
# p <- ggbarplot(plot_data, x="variable",y="log_odds",fill="Gene_class", position = position_dodge(.7),
#           palette = "Dark2")+
#   coord_flip()+
#   geom_text(aes(x=variable, y=log_odds,label=sprintf(fmt = "%.2f", FDR)),size=3)+
#   scale_x_discrete(labels= facet.labs)+
#   labs(x = "",
#        y = "Log Odds",
#        fill="Gene Class")+
#   theme(legend.position="top",
#         strip.text.x.top = element_text(size=4),
#         title = element_text(),
#         legend.text = element_text(size=4),
#         legend.title = element_text(size=6),
#         legend.key.height = unit(2,"mm"),
#         legend.key.width = unit(2,"mm"),
#         axis.text.x = element_text(size=6),
#         axis.title = element_text(size=8),
#         axis.text.y = element_text(size=8),
#         axis.line = element_line(linewidth=.2),
#         axis.ticks = element_line(linewidth = .2))
# 
# p
# 
# 
# jpeg(paste0(plotDirectory,"figure_S2a.jpg"), width=120, height = 80, units = "mm", res = 1000)
# print(p)
# dev.off()



log_odd_ht <- plot_data %>%
  select(Gene_class,variable,log_odds) %>%
  pivot_wider(names_from = variable,values_from = log_odds) %>%
  column_to_rownames("Gene_class") %>%
  t()

fdr_ht <- plot_data %>%
  select(Gene_class,variable,FDR) %>%
  pivot_wider(names_from = variable,values_from = FDR) %>%
  column_to_rownames("Gene_class") %>%
  t()


colnames(log_odd_ht) <- gsub("_"," ",colnames(log_odd_ht))
colnames(log_odd_ht)[1] <- gsub("pr","Pr",colnames(log_odd_ht)[1])

rownames(log_odd_ht) <- gsub("supercluster","Supercluster ",rownames(log_odd_ht))
rownames(log_odd_ht) <- gsub("_signature"," Signature",rownames(log_odd_ht))
rownames(log_odd_ht)[1] <- "MCF7 Global RAC Signature"

ht <- Heatmap(log_odd_ht, name="log(Odds)", cluster_rows = F,cluster_columns = F,
              rect_gp = gpar(col = "white", lwd = 2), column_names_rot=45,
              row_names_gp = gpar(fontsize=6),column_names_gp = gpar(fontsize=6),
              heatmap_legend_param = list(title_gp = gpar(fontsize = 6),legend_height = unit(2, "mm"), grid_width=unit(2,"mm"),
                                          labels_gp = gpar(fontsize = 4),legend_gp = gpar(lwd = .5)),
              
              cell_fun = function(j, i, x, y, w, h, fill) {
                if(fdr_ht[i, j] < 0.001) {
                  grid.text("***", x, y)
                } else if(fdr_ht[i, j] < 0.01) {
                  grid.text("**", x, y)
                }})


jpeg(paste0(plotDirectory,"figure_S2c.jpg"), width=120, height = 80, units = "mm", res = 1000)
draw(ht,heatmap_legend_side = "left")
dev.off()

