library(readxl)
library(tidyverse)
library(decoupleR)

head(crispr_ko_data)
dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"

net <- get_collectri(organism='human', split_complexes=FALSE)
full_tf_list <- unique(net$source)

supercluster_top_tfs <- readRDS(paste0(dataDirectory, "genesets/rac_supercluster_top_tf_list.rds"))
supercluster_bottom_tfs <- readRDS(paste0(dataDirectory, "genesets/rac_supercluster_bottom_tf_list.rds"))

crispr_ko_data <- read_xlsx("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/gottesman_crispr_data/CRISPR_KO_summary_negative.xlsx", sheet = 5)
plot(density(crispr_ko_data$`neg|score`))

final_df <- list()
cutoff <- .05
j <- 2
for(cutoff in seq(.05,1,by=.05)){

  for(j in 1:2){
    
    crispr_ko_top <- crispr_ko_data %>% 
      filter(`neg|score` <= quantile(crispr_ko_data$`neg|score`, probs = cutoff)) %>% 
      pull(id)
    
    supercluster_top_tfs[[j]] <- supercluster_top_tfs[[j]][supercluster_top_tfs[[j]] %in% crispr_ko_data$id]
    supercluster_bottom_tfs[[j]] <- supercluster_bottom_tfs[[j]][supercluster_bottom_tfs[[j]] %in% crispr_ko_data$id]
    
    a <- sum(supercluster_top_tfs[[j]] %in% crispr_ko_top)
    b <- length(supercluster_top_tfs[[j]])
    
    #all other TFs
    c <- sum(full_tf_list %in% crispr_ko_top)
    d <- length(full_tf_list)
    
    # all other genes
    # c <- sum(crispr_ko_data$id %in% crispr_ko_top)
    # d <- length(crispr_ko_data$id)
    # 
    #bottom TFs
    # c <- sum(supercluster_bottom_tfs[[j]] %in% crispr_ko_top)
    # d <- length(supercluster_bottom_tfs[[j]])
    
    contin_table <- matrix(c(a,c,b,d),ncol=2)
    
    fisher_res <- fisher.test(contin_table)
    
    final_df[["cutoff"]] <- append(final_df[["cutoff"]], cutoff)
    final_df[["OR"]] <- append(final_df[["OR"]], fisher_res$estimate)
    final_df[["pval"]] <- append(final_df[["pval"]], fisher_res$p.value)
    final_df[["sc"]] <- append(final_df[["sc"]], j)
    
  }
}

df <- data.frame(final_df)

ggplot(df, aes(x=cutoff,y=OR))+
  facet_wrap(~sc,nrow=2)+
  geom_bar(stat = "identity", color="black",linewidth=.5)+
  geom_hline(yintercept = 1,color="red",linewidth=3)


ggplot(df)
























