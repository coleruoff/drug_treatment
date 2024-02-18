ecoli_AMR_genesets <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/ecoli_data/ecoli_AMR_genesets.rds")



all_genes <- unique(unlist(ecoli_AMR_genesets))

count_dist <- c()
for(curr_gene in all_genes){
  count_dist <- append(count_dist, sum(sapply(ecoli_AMR_genesets, FUN = function(x) curr_gene %in% x)))
}


df <- as.data.frame(cbind(all_genes,count_dist))

sum(df$count_dist >= 2)

df$count_dist <- as.numeric(df$count_dist)

df <- df %>% 
  filter(count_dist > 1)

ggplot(df)+
  geom_bar(aes(x=count_dist))+
  scale_x_continuous(breaks=c(2:9))+
  ylab("Count")+
  xlab("# of Experiments")+
  ggtitle("Frequency of E. coli genes across 9 experiments")



