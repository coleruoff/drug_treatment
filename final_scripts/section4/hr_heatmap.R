dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

hazard_ratio_df <- readRDS(paste0(dataDirectory, "survival_analysis_hazard_ratios/pan_cancer_hazard_ratios_OS.rds"))
hazard_ratio_df <- as.data.frame(hazard_ratio_df)

df <- hazard_ratio_df

df$hazard_ratio <- log(df$hazard_ratio)
df$hazard_ratio.high <- log(df$hazard_ratio.high)
df$hazard_ratio.low <- log(df$hazard_ratio.low)
df$project_signature <- paste0(df$project,df$signature)
df$project_signature <- factor(df$project_signature)

random_hazard_ratio_df <- readRDS(paste0(dataDirectory, "survival_analysis_hazard_ratios/pan_cancer_hazard_ratios_OS_random.rds"))
random_hazard_ratio_df <- as.data.frame(random_hazard_ratio_df)

control_df <- random_hazard_ratio_df

control_df$hazard_ratio <- log(control_df$hazard_ratio)
control_df$hazard_ratio.high <- log(control_df$hazard_ratio.high)
control_df$hazard_ratio.low <- log(control_df$hazard_ratio.low)
control_df$project_signature <- paste0(control_df$project,control_df$signature)

control_df$signature <- paste0(control_df$signature, " control")

colnames(df)
colnames(control_df)

temp <- rbind(df, control_df)

temp <- temp %>% 
  dplyr::select("project","hazard_ratio","p_value","signature")


temp$color <- ifelse(temp$hazard_ratio > 0, "positive","negative")

temp$color <- ifelse(temp$hazard_ratio > 0 & temp$p_value < 0.05, "positive + significant", temp$color)
temp$color <- ifelse(temp$hazard_ratio < 0 & temp$p_value < 0.05, "negative + significant", temp$color)

ggplot(temp, aes(x=signature,y=project,fill=color))+
  geom_tile(color="black")+
  scale_fill_manual(values=c("royalblue1","royalblue4","red1","red4"))




real_sc1_hrs <- df %>% 
  filter(signature == "Supercluster 1") %>% 
  arrange(project) %>% 
  pull(hazard_ratio)

control_sc1_hrs <- control_df %>% 
  filter(signature == "Supercluster 1 control") %>% 
  arrange(project) %>% 
  pull(hazard_ratio)


wilcox.test(real_sc1_hrs,control_sc1_hrs, alternative = "greater")

boxplot(real_sc1_hrs,control_sc1_hrs)


real_sc2_hrs <- df %>% 
  filter(signature == "Supercluster 2") %>% 
  arrange(project) %>% 
  pull(hazard_ratio)

control_sc2_hrs <- control_df %>% 
  filter(signature == "Supercluster 2 control") %>% 
  arrange(project) %>% 
  pull(hazard_ratio)


wilcox.test(real_sc2_hrs,control_sc2_hrs, alternative = "greater")

boxplot(real_sc2_hrs,control_sc2_hrs)


real_sc3_hrs <- df %>% 
  filter(signature == "Supercluster 3") %>% 
  arrange(project) %>% 
  pull(hazard_ratio)

control_sc3_hrs <- control_df %>% 
  filter(signature == "Supercluster 3 control") %>% 
  arrange(project) %>% 
  pull(hazard_ratio)


wilcox.test(real_sc3_hrs,control_sc3_hrs, alternative = "greater")

boxplot(real_sc3_hrs,control_sc3_hrs)





sum(real_sc1_hrs>0)


