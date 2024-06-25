args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])
library(readxl)
library(tidyverse)
library(ggpubr)
library(decoupleR)
library(poolr)
library(rstatix) 
set.seed(42)
# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

#################################################################################

# net <- get_collectri(organism='human', split_complexes=FALSE)
# full_tf_list <- unique(net$source)

supercluster_top_tfs <- readRDS(paste0(dataDirectory, "genesets/rac_supercluster_top_tf_list.rds"))
supercluster_bottom_tfs <- readRDS(paste0(dataDirectory, "genesets/rac_supercluster_bottom_tf_list.rds"))

# supercluster_top_tfs <- readRDS(paste0(dataDirectory, "genesets/supercluster_ranks_top_tfs.rds"))
# supercluster_bottom_tfs <- readRDS(paste0(dataDirectory, "genesets/supercluster_ranks_bottom_tfs.rds"))

final_df <- list()

crispr_ko_data <- read_xlsx("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/gottesman_crispr_data/CRISPR_KO_summary_negative.xlsx", sheet = 5)
j <- 1
for(j in 1:2){
  crispr_ko_ranks <- crispr_ko_data %>% 
    arrange(`neg|fdr`) %>% 
    dplyr::select(id,`neg|rank`) 
  
  crispr_ko_ranks$ranks <- 1:nrow(crispr_ko_ranks)
  
  top_ranks <- crispr_ko_ranks %>% 
    filter(id %in% (supercluster_top_tfs[[j]])) %>% 
    pull(ranks)
  
  
  # down_ranks <- crispr_ko_ranks %>%
  #   filter(id %in% (supercluster_bottom_tfs[[j]])) %>%
  #   pull(ranks)
  
  rest_ranks <- crispr_ko_ranks %>%
    filter(!id %in% (supercluster_top_tfs[[j]])) %>%
    pull(ranks)
  

  
  final_df[["rank"]] <- append(final_df[["rank"]], c(top_ranks,rest_ranks))
  final_df[["top"]] <- append(final_df[["top"]], c(rep("Top",length(top_ranks)),rep("Bottom",length(rest_ranks))))
  final_df[["geneset"]] <- append(final_df[["geneset"]], rep(paste0("Supercluster ", j),length(c(top_ranks,rest_ranks))))
  
}


crispr_ko_ranks$id[crispr_ko_ranks$id %in% (supercluster_top_tfs[[j]])]
crispr_ko_ranks$id[!crispr_ko_ranks$id %in% (supercluster_top_tfs[[j]])]

sort(top_ranks)[1:10]
sort(rest_ranks)[1:10]


final_df <- data.frame(final_df)

final_df$top <- factor(final_df$top,levels = c("Top","Bottom"))

saveRDS(final_df, paste0(dataDirectory, "supercluster_tf_top_bottom_crispr_ranks_data.rds"))

p <- ggboxplot(final_df, x = "top", y = "rank",
               fill = "top", short.panel.labs = T, ncol=1)+
  facet_wrap(~geneset, strip.position = "bottom")


stat.test <- final_df %>%
  group_by(geneset) %>%
  wilcox_test(rank ~ top) %>%
  adjust_pvalue() %>%
  add_significance()

# stat.test$p.adj[1] <- "ns"

p <- p + stat_pvalue_manual(stat.test, label = "p.adj", y.position = 1000,size=8)+
  scale_y_reverse()+
  scale_fill_discrete(name = "TFs")+
  ylab("Gene Rank")+
  xlab("")+
  theme(strip.text = element_text(size=20),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid"),
        legend.position = "right",
        legend.title = element_text(size=26),
        legend.text = element_text(size=20),
        axis.title = element_text(size=30),
        axis.text = element_text(size=20),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

p


print(p)




