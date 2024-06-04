args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])
library(tidyverse)
library(ggpubr)
library(patchwork)
library(grid)
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"
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
          strip.text = element_text(size=20),
          axis.text.x = element_text(size=18),
          axis.text.y = element_text(size=18))
  
  plots <- append(plots, list(p))
  
}

final_plot <- ggarrange(plotlist = plots)+
  theme(plot.margin = margin(.5,0.1,.1,.1, "cm")) 
  

p <- annotate_figure(final_plot, 
                left = textGrob("\nTCGA Project", rot = 90, vjust = 1, gp = gpar(fontsize = 26)),
                bottom = textGrob("Log(Hazard Ratio)\n", gp = gpar(fontsize = 26)),
                top = textGrob("Hazard Ratios for Supercluster Signatures Across TCGA Cancer Types", gp=gpar(fontsize = 46)))



png(paste0(plotDirectory, "figure_4c.png"),
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





# 
# # wrangle results into pre-plotting table form
# res_plot <- df  |>
#   # round estimates and 95% CIs to 2 decimal places for journal specifications
#   mutate(across(
#     c(hazard_ratio, hazard_ratio.low, hazard_ratio.high),
#     ~ str_pad(
#       round(.x, 2),
#       width = 4,
#       pad = "0",
#       side = "right"
#     )
#   ),
#   
#   # add an "-" between HR estimate confidence intervals
#   estimate_lab = paste0(hazard_ratio, " (", hazard_ratio.low, " - ", hazard_ratio.high, ")")) |>
#   # round p-values to two decimal places, except in cases where p < .001
#   mutate(p_value = case_when(
#     p_value < .001 ~ "<0.001",
#     round(p_value, 2) == .05 ~ as.character(round(p_value,3)),
#     p_value < .01 ~ str_pad( # if less than .01, go one more decimal place
#       as.character(round(p_value, 3)),
#       width = 4,
#       pad = "0",
#       side = "right"
#     ),
#     TRUE ~ str_pad( # otherwise just round to 2 decimal places and pad string so that .2 reads as 0.20
#       as.character(round(p_value, 2)),
#       width = 4,
#       pad = "0",
#       side = "right")
#   )) |>
#   # add a row of data that are actually column names which will be shown on the plot in the next step
#   bind_rows(
#     data.frame(
#       project = "Project",
#       project_signature="Project",
#       estimate_lab = "Hazard Ratio (95% CI)",
#       conf.low = "",
#       conf.high = "",
#       p_value = "p-value"
#     )
#   ) |>  mutate(project_signature = fct_rev(fct_relevel(project_signature, "Project")))
# 
# glimpse(res_plot)
# 
# #Remove every other project name
# res_plot$project
# new_project <- c()
# for(i in 1:length(res_plot$project)){
#   if(grepl("1",res_plot$signature[i]) | grepl("Project",res_plot$project[i])){
#     new_project <- append(new_project, res_plot$project[i])
#   } else{
#     new_project <- append(new_project,"")
#   }
# }
# 
# res_plot$project <- new_project
# 
# 
# mid_df <- res_plot %>% filter(!project_signature=="Project")
# mid_df$hazard_ratio <- as.numeric(mid_df$hazard_ratio)
# mid_df$hazard_ratio.high <- as.numeric(mid_df$hazard_ratio.high)
# mid_df$hazard_ratio.low <- as.numeric(mid_df$hazard_ratio.low)
# 
# p <- ggplot(mid_df, aes(x=hazard_ratio,y=project_signature,color=signature))+
#   geom_point(position = position_dodge(width = 1),size=3)+
#   geom_errorbar(aes(xmin=hazard_ratio.low, xmax=hazard_ratio.high), width=1, size=2,position = "dodge") +
#   scale_color_manual(name="",values=c('slateblue2','firebrick1'))+
#   geom_vline(xintercept = 0, linetype="dashed")+
#   theme_classic()+
#   coord_cartesian(ylim=c(1,67),xlim=c(-20,50))+
#   labs(x="Log(Hazard Ratio)",y="")
# 
# p_mid <- p +
#   theme(axis.line.y = element_blank(),
#         axis.ticks.y= element_blank(),
#         axis.text.y= element_blank(),
#         axis.title.y= element_blank(),
#         axis.text.x= element_text(size=20), 
#         axis.title.x= element_text(size=20), 
#         legend.text = element_text(size=16),
#         legend.key.height = unit(1.5,"cm"),
#         legend.key.width = unit(1.5,"cm"),
#         legend.position = c(.85,.7))
# 
# 
# p_left <-
#   res_plot  |>
#   ggplot(aes(y = project_signature))
# 
# p_left <-
#   p_left +
#   geom_text(aes(x = 0, label = project), hjust = 0, fontface = "bold")
# 
# p_left <-
#   p_left +
#   geom_text(
#     aes(x = 1, label = estimate_lab),
#     hjust = 0,
#     fontface = ifelse(res_plot$estimate_lab == "Hazard Ratio (95% CI)", "bold", "plain")
#   )
# 
# p_left <-
#   p_left +
#   theme_void() +
#   coord_cartesian(xlim = c(0, 4))
# 
# p_left
# 
# 
# res_plot$p_value <- ifelse(res_plot$p_value < 0.05, paste0(res_plot$p_value, " *"), res_plot$p_value)
# 
# # right side of plot - pvalues
# p_right <-
#   res_plot  |>
#   ggplot() +
#   geom_text(
#     aes(x = 0, y = project_signature, label = p_value),
#     hjust = 0,
#     fontface = ifelse(res_plot$p_value == "p-value", "bold", "plain")) +
#   theme_void()
# 
# p_right
# 
# 
# 
# layout <- c(
#   area(t = 0, l = 0, b = 30, r = 3), # left plot, starts at the top of the page (0) and goes 30 units down and 3 units to the right
#   area(t = 1, l = 2, b = 30, r = 4), # middle plot starts a little lower (t=1) because there's no title. starts 1 unit right of the left plot (l=4, whereas left plot is r=3), goes to the bottom of the page (30 units), and 6 units further over from the left plot (r=9 whereas left plot is r=3)
#   area(t = 0, l = 5, b = 30, r = 11) # right most plot starts at top of page, begins where middle plot ends (l=9, and middle plot is r=9), goes to bottom of page (b=30), and extends two units wide (r=11)
# )
# 
# # final plot arrangement
# p <- p_left + p_right + p_mid  + plot_layout(design = layout) + plot_annotation("Hazard Ratios for Supercluster Signatures Across TCGA Cancer Types\n",theme=theme(plot.title=element_text(hjust=0.5,face='bold',size=30)))


# png(paste0(plotDirectory, "final_figures/figure_4c.png"),
#     width=16, height = 18, units="in",res=300)
# 
# print(p)
# 
# dev.off()





