library(readxl)
library(tidyverse)
library(decoupleR)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"

# net <- get_collectri(organism='human', split_complexes=FALSE)
# full_tf_list <- unique(net$source)

supercluster_top_tfs <- readRDS(paste0(dataDirectory, "genesets/rac_supercluster_top_tf_list.rds"))
supercluster_bottom_tfs <- readRDS(paste0(dataDirectory, "genesets/rac_supercluster_bottom_tf_list.rds"))

crispr_ko_data <- read_xlsx("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/gottesman_crispr_data/CRISPR_KO_summary_negative.xlsx", sheet = 5)
plot(density(crispr_ko_data$`neg|score`))


quantile(crispr_ko_data$`neg|score`, probs = .39)

final_df <- list()

for(j in 1:2){
  
  crispr_ko_top <- crispr_ko_data %>% 
    filter(`neg|score` < quantile(crispr_ko_data$`neg|score`, probs = .25)) %>% 
    pull(id)
  
  a <- sum(supercluster_top_tfs[[j]] %in% crispr_ko_top)
  b <- sum(!supercluster_top_tfs[[j]] %in% crispr_ko_top)
  
  #all other TFs
  # c <- sum(full_tf_list %in% crispr_ko_top)
  # d <- sum(!full_tf_list %in% crispr_ko_top)
  
  # all other genes
  # c <- sum(crispr_ko_data$id %in% crispr_ko_top)
  # d <- sum(!crispr_ko_data$id %in% crispr_ko_top)
  
  #bottom TFs
  c <- sum(supercluster_bottom_tfs[[j]] %in% crispr_ko_top)
  d <- sum(!supercluster_bottom_tfs[[j]] %in% crispr_ko_top)
  
  contin_table <- matrix(c(a,c,b,d),ncol=2)
  
  fisher_res <- fisher.test(contin_table)
  
  final_df[["OR"]] <- append(final_df[["OR"]], fisher_res$estimate)
  final_df[["pval"]] <- append(final_df[["pval"]], fisher_res$p.value)
  final_df[["sc"]] <- append(final_df[["sc"]], j)
  
}


df <- data.frame(final_df)

ggplot(df, aes(x=sc,y=OR))+
  geom_bar(stat = "identity", color="black",linewidth=.5)+
  geom_hline(yintercept = 1,color="red",linewidth=3)



