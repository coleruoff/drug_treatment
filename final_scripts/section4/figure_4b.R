args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])
library(tidyverse)
library(ggpubr)
library(patchwork)
library(grid)
set.seed(42)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"
# plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/"

################################################################################

hazard_ratio_df <- readRDS(paste0(dataDirectory, "survival_analysis_hazard_ratios/pan_cancer_hazard_ratios_OS.rds"))
hazard_ratio_df <- as.data.frame(hazard_ratio_df)

df <- hazard_ratio_df
df$hazard_ratio <- log(df$hazard_ratio)
df$hazard_ratio.high <- log(df$hazard_ratio.high)
df$hazard_ratio.low <- log(df$hazard_ratio.low)
df$project_signature <- paste0(df$project,df$signature)
df$project_signature <- factor(df$project_signature)


plots <- list()
for(i in 1:length(unique(df$signature))){
  
  curr_signature <- unique(df$signature)[i]
  
  curr_df <- df %>% 
    filter(signature == curr_signature)
  
  color_to_use <- c('slateblue2','firebrick1')[i]
  
  curr_df$project <- ifelse(curr_df$p_value < 0.05, paste0("* ", curr_df$project), curr_df$project)
  
  
  p <- ggplot(curr_df, aes(x=hazard_ratio,y=reorder(project,hazard_ratio,mean),color=signature))+
    geom_point(position = position_dodge(width = 1),size=3)+
    geom_errorbar(aes(xmin=hazard_ratio.low, xmax=hazard_ratio.high), width=1, size=2,position = "dodge") +
    scale_color_manual(name="",values=c(color_to_use))+
    geom_vline(xintercept = 0, linetype="dashed")+
    labs(x="",y="")+
    facet_wrap(~signature)+
    theme_classic()+
    theme(legend.position="none",
          strip.text = element_text(size=30),
          axis.text.x = element_text(size=28),
          axis.text.y = element_text(size=18),
          plot.margin = unit(c(1,1,1,1), "cm"))
  
  plots <- append(plots, list(p))
  
}

final_plot <- ggarrange(plotlist = plots)
  

p <- annotate_figure(final_plot, 
                left = textGrob("TCGA Project", rot = 90, vjust = 1, gp = gpar(fontsize = 32)),
                bottom = textGrob("Log(Hazard Ratio)\n", gp = gpar(fontsize = 32)))




png(paste0(plotDirectory, "figure_4b.png"),
    width=26, height = 16, units="in",res=300)

print(p)

dev.off()


# count projects with at least one significant and positive log(HR)s 
for(curr_signature in unique(df$signature)){
  cat(curr_signature, "\n")
  num_gr_1 <- df %>% 
    filter(signature == curr_signature & hazard_ratio > 0) %>% 
    nrow()
  
  num_gr_1_sig <- df %>% 
    filter(signature == curr_signature & hazard_ratio > 0 & p_value < 0.05) %>% 
    nrow()
  
  cat(paste0("HR > 1:                 ", num_gr_1, "/33 projects\n"))
  cat(paste0("HR > 1 and pval < 0.05: ", num_gr_1_sig, "/33 projects\n\n"))
}



