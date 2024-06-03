#!/bin/bash

module load R/4.3

working_directory="/data/CDSL_hannenhalli/Cole/projects/drug_treatment/"
source_directory="/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_scripts/section5/"

echo "create_candida_auris_orthologs"
Rscript --vanilla ${source_directory}create_candida_auris_orthologs.R $working_directory

echo "create_ecoli_orthologs"
Rscript --vanilla ${source_directory}create_ecoli_orthologs.R $working_directory

#AUCell scoring
echo "A549 scoring"
Rscript --vanilla ${source_directory}../AUCell_precompute_thresholds.R $working_directory "sciPlex_data/A549_processed_filtered" "RNA" "yeast_human_orthologs_up" $SLURM_CPUS_PER_TASK
Rscript --vanilla ${source_directory}../AUCell_precompute_thresholds.R $working_directory "sciPlex_data/A549_processed_filtered" "RNA" "ecoli_human_orthologs_up" $SLURM_CPUS_PER_TASK

echo "K562 scoring"
Rscript --vanilla ${source_directory}../AUCell_precompute_thresholds.R $working_directory "sciPlex_data/K562_processed_filtered" "RNA" "yeast_human_orthologs_up" $SLURM_CPUS_PER_TASK
Rscript --vanilla ${source_directory}../AUCell_precompute_thresholds.R $working_directory "sciPlex_data/K562_processed_filtered" "RNA" "ecoli_human_orthologs_up" $SLURM_CPUS_PER_TASK

echo "MCF7 scoring"
Rscript --vanilla ${source_directory}../AUCell_precompute_thresholds.R $working_directory "sciPlex_data/MCF7_processed_filtered" "RNA" "yeast_human_orthologs_up" $SLURM_CPUS_PER_TASK
Rscript --vanilla ${source_directory}../AUCell_precompute_thresholds.R $working_directory "sciPlex_data/MCF7_processed_filtered" "RNA" "ecoli_human_orthologs_up" $SLURM_CPUS_PER_TASK

echo "figure 5a"
Rscript --vanilla ${source_directory}figure_5a.R $working_directory

echo "figure 5b"
Rscript --vanilla ${source_directory}figure_5b.R $working_directory

echo "figure 5c"
Rscript --vanilla ${source_directory}figure_5c.R $working_directory

echo "figure 5d"
Rscript --vanilla ${source_directory}figure_5c.R $working_directory

echo "create_supercluster_gene_spreadsheet"
Rscript --vanilla ${source_directory}create_supercluster_gene_spreadsheet.R $working_directory
