args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

source("final_scripts/drug_treatment_functions.R")
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"
# plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_figures/"

################################################################################
data <- readRDS(paste0(dataDirectory, "breast_premalignancy_data/datasets_for_boxplot_fisheroverlap_newSignatures.Rds"))

plot_data <- data[[3]]

plot_data$log_odds[which(plot_data$log_odds == -Inf)] <- -5

plot_data$variable <- factor(plot_data$variable, levels = c("MCF7","supercluster1_signature", "supercluster2_signature", "supercluster3_signature"))



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

