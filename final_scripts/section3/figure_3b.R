args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

source("final_scripts/drug_treatment_functions.R")
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"
# plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_figures/"

################################################################################

hazard_ratio_df <- readRDS(paste0(dataDirectory, "survival_analysis_hazard_ratios/pan_cancer_hazard_ratios_OS.rds"))
hazard_ratio_df <- as.data.frame(hazard_ratio_df)


df <- hazard_ratio_df
df$hazard_ratio <- log(df$hazard_ratio)
df$hazard_ratio.high <- log(df$hazard_ratio.high)
df$hazard_ratio.low <- log(df$hazard_ratio.low)
df$project_signature <- paste0(df$project,df$signature)
df$project_signature <- factor(df$project_signature)
i <- 2

plots <- list()

for(i in 1:length(unique(df$signature))){
  
  curr_signature <- unique(df$signature)[i]
  
  curr_df <- df %>% 
    filter(signature == curr_signature)
  
  res <- t.test(curr_df$hazard_ratio, alternative = "greater")
  cat(res$p.value, "\n")
  
  color_to_use <- c('slateblue2','firebrick1')[i]
  
  # curr_df$project <- ifelse(curr_df$p_value < 0.05, paste0("* ", curr_df$project), curr_df$project)
  
  p <- ggplot(curr_df, aes(x=hazard_ratio,y=project,color=signature))+
    geom_point(position = position_dodge(width = .5),size=.2)+
    geom_errorbar(aes(xmin=hazard_ratio.low, xmax=hazard_ratio.high), width=.5, size=.2, position = "dodge") +
    scale_color_manual(name="",values=c(color_to_use))+
    geom_vline(xintercept = 0, linetype="dashed", size = .2)+
    labs(x="",y="")+
    facet_wrap(~signature)+
    theme_classic()+
    theme(legend.position="none",
          strip.text = element_text(size=6),
          axis.text = element_text(size=6),
          axis.title = element_text(size=6),
          axis.line = element_line(linewidth=.2),
          axis.ticks = element_line(linewidth = .2))
  
  plots <- append(plots, list(p))
  
}

final_plot <- ggarrange(plotlist = plots, nrow=1)


p <- annotate_figure(final_plot, 
                     left = textGrob("TCGA Project", rot = 90, vjust = 1, gp = gpar(fontsize = 6)),
                     bottom = textGrob("Log(Hazard Ratio)", gp = gpar(fontsize = 6)))



p


jpeg(paste0(plotDirectory,"figure_3b.jpg"), width=180, height = 100, units = "mm", res = 1000)
print(p)
dev.off()



# length(tcga_projects)

# count projects with at least one significant and positive log(HR)s 
for(curr_signature in unique(df$signature)){
  cat(curr_signature, "\n")
  num_gr_1 <- df %>% 
    filter(signature == curr_signature & hazard_ratio > 0) %>% 
    nrow()
  
  num_gr_1_sig <- df %>% 
    filter(signature == curr_signature & hazard_ratio > 0 & p_value < 0.05) %>% 
    nrow()
  
  cat(paste0("HR > 1:                 ", num_gr_1, "/18 projects\n"))
  cat(paste0("HR > 1 and pval < 0.05: ", num_gr_1_sig, "/18 projects\n\n"))
}

# Count projects with at least one HR > 1
count <- 0
for(curr_project in unique(hazard_ratio_df$project)){
  
  
  # hr_values <- hazard_ratio_df %>% 
  #   filter(project == curr_project) %>% 
  #   pull(hazard_ratio)
  # 
  # count <- ifelse(sum(hr_values > 1) > 0, count+1,count)
  
  temp <- hazard_ratio_df %>% 
    filter(project == curr_project & hazard_ratio > 0 & p_value < 0.05) %>% 
    nrow()
  
  count <- ifelse(temp > 0, count+1,count)
  
}


count/length(unique(hazard_ratio_df$project))


