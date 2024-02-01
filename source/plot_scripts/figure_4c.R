setwd("/data/ruoffcj/projects/drug_treatment/")
library(tidyverse)
library(ggpubr)
library(patchwork)

hazard_ratio_df <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/survival_analysis_hazard_ratios/pan_cancer_hazard_ratios.rds")

df <- hazard_ratio_df


df$hazard_ratio <- log(df$hazard_ratio)
df$hazard_ratio.high <- log(df$hazard_ratio.high)
df$hazard_ratio.low <- log(df$hazard_ratio.low)


df$project_signature <- paste0(df$project,df$signature)



p <- ggplot(df, aes(x=hazard_ratio,color=signature,y=project))+
  geom_point(position = position_dodge(width = 1))+
  geom_errorbar(aes(xmin=hazard_ratio.low, xmax=hazard_ratio.high), width=1, position = "dodge") +
  geom_vline(xintercept = 0, linetype="dashed")+
  theme_classic()+
  coord_cartesian(ylim=c(1,34),xlim=c(-20,25))+
  labs(x="Log(Hazard Ratio)",y="")



p

p_mid <- p +
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        legend.text = element_text(size=16),
        legend.key.height = unit(1.5,"cm"),
        legend.key.width = unit(1.5,"cm"))+
  scale_color_discrete(name = "")



# wrangle results into pre-plotting table form
res_plot <- df  |>
  # round estimates and 95% CIs to 2 decimal places for journal specifications
  mutate(across(
    c(hazard_ratio, hazard_ratio.low, hazard_ratio.high),
    ~ str_pad(
      round(.x, 2),
      width = 4,
      pad = "0",
      side = "right"
    )
  ),
  # add an "-" between HR estimate confidence intervals
  estimate_lab = paste0(hazard_ratio, " (", hazard_ratio.low, "-", hazard_ratio.high, ")")) |>
  # round p-values to two decimal places, except in cases where p < .001
  mutate(p_value = case_when(
    p_value < .001 ~ "<0.001",
    round(p_value, 2) == .05 ~ as.character(round(p_value,3)),
    p_value < .01 ~ str_pad( # if less than .01, go one more decimal place
      as.character(round(p_value, 3)),
      width = 4,
      pad = "0",
      side = "right"
    ),
    TRUE ~ str_pad( # otherwise just round to 2 decimal places and pad string so that .2 reads as 0.20
      as.character(round(p_value, 2)),
      width = 4,
      pad = "0",
      side = "right")
  )) |>
  # add a row of data that are actually column names which will be shown on the plot in the next step
  bind_rows(
    data.frame(
      project = "Project",
      estimate_lab = "Hazard Ratio (95% CI)",
      conf.low = "",
      conf.high = "",
      p_value = "p-value"
    )
  ) |>
  mutate(project_signature = fct_rev(fct_relevel(project_signature, "Model")))

glimpse(res_plot)

new_project <- c()
for(i in 1:length(res_plot$project)){
  if(i%%2==1){
    new_project <- append(new_project, res_plot$project[i])
  } else{
    new_project <- append(new_project,"")
  }
}

res_plot$project <- new_project

p_left <-
  res_plot  |>
  ggplot(aes(y = project_signature))

p_left <-
  p_left +
  geom_text(aes(x = 0, label = project), hjust = 0, fontface = "bold")

p_left <-
  p_left +
  geom_text(
    aes(x = 1, label = estimate_lab),
    hjust = 0,
    fontface = ifelse(res_plot$estimate_lab == "Hazard Ratio (95% CI)", "bold", "plain")
  )

p_left <-
  p_left +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))

p_left

# right side of plot - pvalues
p_right <-
  res_plot  |>
  ggplot() +
  geom_text(
    aes(x = 0, y = project_signature, label = p_value),
    hjust = 0,
    fontface = ifelse(res_plot$p_value == "p-value", "bold", "plain")) +
  theme_void()

p_right



layout <- c(
  area(t = 0, l = 0, b = 30, r = 3), # left plot, starts at the top of the page (0) and goes 30 units down and 3 units to the right
  area(t = 1, l = 4, b = 30, r = 9), # middle plot starts a little lower (t=1) because there's no title. starts 1 unit right of the left plot (l=4, whereas left plot is r=3), goes to the bottom of the page (30 units), and 6 units further over from the left plot (r=9 whereas left plot is r=3)
  area(t = 0, l = 9, b = 30, r = 11) # right most plot starts at top of page, begins where middle plot ends (l=9, and middle plot is r=9), goes to bottom of page (b=30), and extends two units wide (r=11)
)
# final plot arrangement



# figure1 <- p_left + p_mid + p_right + plot_layout(design = layout)

plots <- list(p_left,p_mid,p_right)

figure <- ggarrange(plotlist = plots,ncol=3,widths = c(0.5, 1.4,0.5))

p <- annotate_figure(figure, top = text_grob("Hazard Ratios for Supercluster Signatures Across TCGA Cancer Types\n", face = "bold", size = 36))

png("/data/ruoffcj/projects/drug_treatment/final_figures/figure_4c.png",
    width=20, height = 16, units="in",res=300)

print(p)

dev.off()

