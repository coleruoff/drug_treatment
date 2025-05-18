args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

source("final_scripts/drug_treatment_functions.R")
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"
# plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_figures/"

################################################################################
# data <- readRDS(paste0(dataDirectory, "breast_premalignancy_data/datasets_for_boxplot_fisheroverlap.Rds"))
# data <- readRDS(paste0(dataDirectory, "breast_premalignancy_data/datasets_for_boxplot_newSignatures.Rds"))
data <- readRDS(paste0(dataDirectory, "breast_premalignancy_data/datasets_for_boxplot_fisheroverlap_newSignatures.Rds"))
                
plot_data <- data[[1]]

plot_data$variable <- factor(plot_data$variable, levels = c("MCF7","supercluster1_signature", "supercluster2_signature", "supercluster3_signature"))

facet.labs <- c("MCF7 Global RAC Signature","Supercluster 1 Signature", "Supercluster 2 Signature", "Supercluster 3 Signature")

colnames(plot_data)[2] <- "Progression"

p <- ggboxplot(plot_data,x = "Progression", y="value", facet.by = "variable", scales="free", nrow=1, panel.labs = list(variable = facet.labs),
          fill = "Progression", palette = "Dark2", size=.2, outlier.shape = NA)+
  stat_compare_means(label = "p.format", method = "wilcox", size=2, method.args = list(alternative = "less"), label.y.npc = "top",vjust = 1)+
  labs(x = "Progression",
       y = "Score")+
  theme(legend.position="none",
        strip.text.x.top = element_text(size=4),
        title = element_text(),
        legend.text = element_text(size=4),
        legend.title = element_text(size=6),
        legend.key.height = unit(2,"mm"),
        legend.key.width = unit(2,"mm"),
        axis.text.x = element_text(size=6),
        axis.title = element_text(size=8),
        axis.text.y = element_text(size=4),
        axis.line = element_line(linewidth=.2),
        axis.ticks = element_line(linewidth = .2))

p


jpeg(paste0(plotDirectory,"figure_3d.jpg"), width=120, height = 50, units = "mm", res = 1000)
print(p)
dev.off()


###############################################################################

