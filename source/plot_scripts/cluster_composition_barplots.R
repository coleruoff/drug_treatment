

curr_cell_line <- "A549"


data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))








df <- data@meta.data %>% 
  count(treatment_stage,Cluster) %>% 
  filter(n > 10)




df$n <- ifelse(df$n < 10, 0, df$n) 


plot_title <- paste0("Treatment Stage Cell Counts Across Clusters (",curr_cell_line,")")

ggplot(df)+
  geom_col(aes(x=Cluster, y=n , fill=treatment_stage), position = 'dodge')+
  ggtitle(plot_title)+
  xlab("Cluster")+
  ylab("Cell Counts")+
  scale_fill_discrete(name = "Treatment Stage")+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=15),
        plot.title = element_text(size=30,face="bold"),
        legend.title = element_text(size=20),
        legend.text = element_text(size=15))




df$n <- (df$n/(total_pre+total_post))*100

total_pre <- sum(df %>% 
  filter(treatment_stage=="pre") %>% 
  pull(n))


total_post <- sum(df %>% 
                   filter(treatment_stage=="post") %>% 
                   pull(n))


pre_df <- df %>% 
  filter(treatment_stage=="pre")

pre_df$n <- (pre_df$n/total_pre)*100


post_df <- df %>% 
  filter(treatment_stage=="post")

post_df$n <- (post_df$n/total_post)*100

df_2 <- rbind(pre_df,post_df)

plot_title <- paste0("Treatment Stage Cluster Percentages (",curr_cell_line,")")

ggplot(df_2)+
  geom_col(aes(x=Cluster, y=n , fill=treatment_stage), position = 'dodge')+
  ggtitle(plot_title)+
  xlab("Cluster")+
  ylab("Percentage of Treatment Stage")+
  scale_fill_discrete(name = "Treatment Stage")+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=15),
        plot.title = element_text(size=30,face="bold"),
        legend.title = element_text(size=20),
        legend.text = element_text(size=15))
