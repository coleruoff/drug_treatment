args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])
source("final_scripts/drug_treatment_functions.R")
library(Seurat)
library(tidyverse)

if (file.exists(paste0(plotDirectory,"figure_1a.png"))) {
  #Delete file if it exists
  file.remove(paste0(plotDirectory,"figure_1a.png"))
}

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"
source("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_scripts/drug_treatment_functions.R")
plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/"

################################################################################

# Read in resistance signature 
initial_resistance_signature <- readRDS(paste0(dataDirectory, "genesets/initial_resistance_signature.rds"))

# Calculate GSEA enrichment at each time point for resistance signature
gsea_results <- watermelon_validation(initial_resistance_signature)

# Plot NES of resistance signature at each time point as line plot
df <- as.data.frame(gsea_results)
colnames(df) <- c("day","value")



# df$value <- as.numeric(df$value)
df$day <- factor(df$day, levels = df$day)
df$dot_color <- ifelse(df$value > 0, "pos","neg")

p <- ggplot(df, aes(x=day,y=NES))+
  geom_line(aes(group=1), linewidth=2)+
  geom_point(aes(color=dot_color),size=5)+
  scale_color_manual(values=c("pos" = "red","neg" ="blue"))+
  xlab("Time Point")+
  ylab("NES")+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=15),
        plot.title = element_text(size=30,face="bold"))+
  NoLegend()+
  geom_errorbar(aes(ymin=NES-log2err, ymax=NES+log2err), width=.2,
                position=position_dodge(.9)) 


p


png(paste0(plotDirectory,"figure_1a.png"),width=12,height=5, units = "in", res = 300)

print(p)

dev.off()


