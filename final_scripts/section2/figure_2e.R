args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

source("final_scripts/drug_treatment_functions.R")

set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"
# plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_figures/"

#################################################################################

supercluster_up_genes <- readRDS(paste0(dataDirectory, "genesets/supercluster_up_signatures.rds"))
supercluster_down_genes <- readRDS(paste0(dataDirectory, "genesets/supercluster_down_signatures.rds"))

final_df <- list()

crispr_ko_data <- read_xlsx(paste0(dataDirectory,"gottesman_crispr_data/CRISPR_KO_summary_negative.xlsx"), sheet = 5)

crispr_ko_data[,5]
j <- 1
for(j in 1:length(supercluster_up_genes)){
  crispr_ko_ranks <- crispr_ko_data %>% 
    arrange(`neg|rank`) %>% 
    dplyr::select(id,`neg|rank`) %>% 
    mutate(ranks = 1:nrow(.))  
  
  
  top_ranks <- crispr_ko_ranks %>%
    filter(id %in% supercluster_up_genes[[j]]) %>%
    pull(ranks)

  rest_ranks <- crispr_ko_ranks %>%
    filter(id %in% supercluster_down_genes[[j]]) %>%
    pull(ranks)
  
  # rest_ranks <- crispr_ko_ranks %>%
  #   filter(id %in% full_tf_list & !id %in% supercluster_top_tfs[[j]]) %>%
  #   pull(`neg|score`)

  # rest_ranks <- crispr_ko_ranks %>%
  #   filter(!id %in% supercluster_top_tfs[[j]]) %>%
  #   pull(`neg|score`)
  
  
  final_df[["rank"]] <- append(final_df[["rank"]], c(top_ranks,rest_ranks))
  final_df[["top"]] <- append(final_df[["top"]], c(rep("Upregulated",length(top_ranks)),rep("Downregulated",length(rest_ranks))))
  final_df[["geneset"]] <- append(final_df[["geneset"]], rep(paste0("Supercluster ", j),length(c(top_ranks,rest_ranks))))
  
}


final_df <- data.frame(final_df)

final_df$top <- factor(final_df$top,levels = c("Upregulated","Downregulated"))

# saveRDS(final_df, paste0(dataDirectory, "supercluster_tf_top_bottom_crispr_ranks_data.rds"))
# 


p <- ggviolin(final_df, x = "top", y = "rank",
              fill = "top", short.panel.labs = T, ncol=1,
              add="boxplot", size=.2, add.params = list(size=.2))+
  facet_wrap(~geneset, strip.position = "top",nrow=1)

p
stat.test <- final_df %>%
  group_by(geneset) %>%
  wilcox_test(rank ~ top) %>%
  adjust_pvalue() %>%
  add_significance()

stat.test$p.adj <- paste0("p = ", sprintf("%.5f", stat.test$p.adj))
# stat.test$p.adj[c(2,3)] <- "ns"


p <- p + stat_pvalue_manual(stat.test, label = "p.adj", y.position = 0, size=1, vjust = -2, label.y = .5)+
  scale_y_reverse(breaks=c(1, 10000,20000))+
  scale_fill_discrete(name = "Genes")+
  ylab("CRISPR Rank")+
  xlab("")+
  theme(strip.text = element_text(size=6),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid"),
        legend.position = "right",
        legend.title = element_text(size=5),
        legend.text = element_text(size=4),
        legend.key.height = unit(3,"mm"),
        legend.key.width = unit(3,"mm"),
        axis.title = element_text(size=6),
        axis.text = element_text(size=6),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        title = element_text(size=8),
        axis.line = element_line(linewidth=.2),
        axis.ticks = element_line(linewidth = .2))


p

jpeg(paste0(plotDirectory,"figure_2e.jpg"), width=100, height = 60, units = "mm", res = 1000)
print(p)
dev.off()





j <- 1

top_ranks <- crispr_ko_ranks %>%
  filter(id %in% supercluster_up_genes[[j]]) %>%
  arrange(ranks) %>% 
  pull(id)

bottom_ranks <- crispr_ko_ranks %>%
  filter(id %in% supercluster_down_genes[[j]]) %>%
  arrange(ranks) %>% 
  pull(id)

top_ranks[1:25]
which("ABCB1" == top_ranks)



